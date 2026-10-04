#!/usr/bin/env python3
"""Build release kits and exercise them from isolated consumer directories."""
import argparse
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
import re
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'release/manifest.json'
sys.path.insert(0, str(ROOT / 'benchmarks'))
from publication import PublicSanitizer
from wasm_publication import sanitize_javascript, sanitize_wasm


def run(args, cwd, env=None):
    # Logs are local/CI diagnostics, never copied into release assets.
    subprocess.run([str(a) for a in args], cwd=cwd, env=env, check=True)


def tracked_sources(root, paths):
    names = subprocess.check_output(['git', 'ls-files', '-z', '--', *paths], cwd=root).decode().split('\0')
    result = {}
    for name in filter(None, names):
        path = root / name
        if path.is_symlink() or not path.is_file():
            raise ValueError('source kits require regular files: ' + name)
        result[name] = path.read_bytes()
    return result


def write_archive(destination, files):
    """Deterministic gzip/tar metadata without local users or filesystem paths."""
    with destination.open('wb') as stream, gzip.GzipFile(filename='', mode='wb', fileobj=stream, mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode='w') as archive:
            for name, data in sorted(files.items()):
                if Path(name).is_absolute() or '..' in Path(name).parts:
                    raise ValueError('archive paths must be relative')
                info = tarfile.TarInfo(name)
                info.size = len(data)
                info.mode = 0o644
                archive.addfile(info, io.BytesIO(data))


def unpack(archive, destination):
    """Extract only the regular, relative members produced by this pipeline."""
    with tarfile.open(archive) as package:
        for entry in package:
            if not entry.isfile() or Path(entry.name).is_absolute() or '..' in Path(entry.name).parts:
                raise ValueError('unexpected archive member')
            path = destination / entry.name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(package.extractfile(entry).read())


def table(manifest):
    lines = ['# API compatibility: ' + manifest['version'], '',
             'This is the shared float64 contract for these release kits. A check mark is',
             'a declared capability validated by the shared conformance suite in CI.',
             'Package installation tests separately exercise solving, multiplication, SVD and its inverse APIs, CSR multiplication, conjugate gradient, GMRES, ILU(0) reuse, RCM/AMD permutations and sparse Cholesky and IC(0)-CG reuse.', '',
             '| Port | Required runtime/toolchain | ' + ' | '.join(manifest['capabilities']) + ' |',
             '| --- | --- | ' + ' | '.join('---' for _ in manifest['capabilities']) + ' |']
    for name, port in manifest['ports'].items():
        lines.append('| ' + name + ' | ' + port['requires'] + ' | ' + ' | '.join('✓' if c in port['capabilities'] else '—' for c in manifest['capabilities']) + ' |')
    lines += ['', 'C# entries apply to `Matrix<double>` on .NET 10. The netstandard2.1',
              'target contains the older compatibility classes and does not promise this',
              'complete contract. ARM64 and WebAssembly use the C API/kernels; C++ shares',
              'the C general-eigenvalue, SVD and factorization kernels.', '']
    return '\n'.join(lines)


def check_manifest(manifest):
    version = manifest['version']
    if not re.fullmatch(r'\d+\.\d+\.\d+', version):
        raise ValueError('expected a numeric release version')
    api_readme = ROOT / 'docs/api' / version / 'README.md'
    if not api_readme.is_file():
        raise ValueError('missing API contract for release version ' + version)
    api_lines = api_readme.read_text().splitlines()
    if not api_lines or api_lines[0] != '# linear-A ' + version + ' API contract':
        raise ValueError('API contract title differs from release manifest')
    api_link = re.search(r'\[Versioned API contract\]\(([^)]+)\)', (ROOT / 'README.md').read_text())
    if not api_link or api_link.group(1) != 'docs/api/' + version + '/README.md':
        raise ValueError('root README API contract link differs from release manifest')
    for name in ('python', 'rust', 'julia'):
        path = ROOT / 'ports' / name / {'python':'pyproject.toml','rust':'Cargo.toml','julia':'Project.toml'}[name]
        if not re.search(r'^version\s*=\s*"' + re.escape(version) + '"', path.read_text(), re.M):
            raise ValueError(name + ' package version differs from release manifest')
    if json.loads((ROOT / 'ports/typescript/package.json').read_text())['version'] != version:
        raise ValueError('TypeScript version differs from release manifest')
    for port in manifest['ports'].values():
        if set(port['capabilities']) != set(manifest['capabilities']):
            raise ValueError('shared capability declarations differ')
        if not (ROOT / 'release/examples' / port['example']).is_file():
            raise ValueError('missing consumer example')
    if (ROOT / 'docs/api' / version / 'compatibility.md').read_text() != table(manifest):
        raise ValueError('regenerate compatibility.md with --write-table')


def kit_files(name, port, version, provenance=None):
    files = tracked_sources(ROOT, port['source_roots'])
    for path in ['LICENSE', 'Directory.Build.props', 'global.json', 'ports/README.md', 'ports/SOLVING.md', 'ports/SVD.md', 'ports/SPARSE.md', 'ports/CHOLESKY.md', 'ports/AMD.md', 'ports/IC0.md']:
        files[path] = (ROOT / path).read_bytes()
    for path in (ROOT / 'docs/api' / version).glob('*.md'):
        files['docs/api/' + version + '/' + path.name] = path.read_bytes()
    files['example/' + port['example']] = (ROOT / 'release/examples' / port['example']).read_bytes()
    files['release.json'] = json.dumps({'version': version, 'port': name, **(provenance or {}), **port}, indent=2).encode()
    return files


def smoke(name, port, kit, output, version):
    with tempfile.TemporaryDirectory(prefix='linear-a-consumer-') as temporary:
        base = Path(temporary)
        source, consumer = base / 'source', base / 'consumer'
        unpack(kit, source)
        consumer.mkdir()
        example = consumer / port['example']
        shutil.copyfile(source / 'example' / port['example'], example)
        env = os.environ.copy()
        # Prevent repository imports from disguising incomplete packages.
        env.pop('PYTHONPATH', None)
        env.pop('NODE_PATH', None)
        env['PYTHONNOUSERSITE'] = '1'
        if name in ('c', 'cpp', 'assembly'):
            if name == 'assembly' and platform.machine().lower() not in ('arm64','aarch64'):
                raise RuntimeError('assembly consumer requires an ARM64 host')
            if name == 'cpp':
                command = [os.environ.get('CXX','c++'), '-std=c++17', '-I', source/'ports/cpp', example]
            else:
                command = [os.environ.get('CC','cc'), '-std=c11', '-I', source/'ports/c', example, source/'ports/c/matrix.c']
                if name == 'assembly': command += ['-DMATRIX_USE_ASM', source/'ports/assembly/kernels.S']
            run([*command, '-O2', '-ffp-contract=off', '-lm', '-o', consumer/'example'], consumer, env)
            run([consumer/'example'], consumer, env)
        elif name == 'python':
            run([sys.executable, '-m', 'venv', base/'venv'], consumer, env)
            python = base/'venv/bin/python'
            package = source/'ports/python'
            shutil.copyfile(source/'LICENSE', package/'LICENSE')
            metadata = package/'pyproject.toml'
            metadata.write_text(metadata.read_text().replace('[project]', '[project]\nlicense = {file = \"LICENSE\"}'))
            run([python, '-m', 'pip', 'wheel', '--no-deps', '--wheel-dir', output, source/'ports/python'], consumer, env)
            wheel = next(output.glob('linear_a_python-*.whl'))
            run([python, '-m', 'pip', 'install', '--no-index', '--no-deps', wheel], consumer, env)
            run([python, '-I', example], consumer, env)
        elif name in ('typescript','wasm'):
            package = source/'ports'/name
            if name == 'typescript':
                run(['npm','ci','--ignore-scripts'], package, env)
                run(['npm','run','build'], package, env)
                meta = json.loads((package/'package.json').read_text())
                meta['files'] = ['dist/*.js','dist/*.d.ts','README.md','LICENSE']
                # Ship only library modules; tests and protocol executables are not API.
                for path in (package/'dist').glob('*'):
                    if '.test.' in path.name or path.name.startswith('runner.'): path.unlink()
                (package/'package.json').write_text(json.dumps(meta,indent=2)+'\n')
            else:
                runtime = package/'runtime'
                run([sys.executable, package/'build.py', '--output', runtime], source, env)
                sanitizer = PublicSanitizer(root=source)
                (runtime/'matrix.mjs').write_text(sanitize_javascript((runtime/'matrix.mjs').read_text(),sanitizer))
                (runtime/'matrix.wasm').write_bytes(sanitize_wasm((runtime/'matrix.wasm').read_bytes(),sanitizer))
                # Keep the development wrapper's override while supplying a relocatable default.
                (package/'index.mjs').write_text('import { createMatrixAPI as create } from "./matrix.mjs";\nexport const createMatrixAPI = (options = {}) => create({moduleUrl: new URL("./runtime/matrix.mjs", import.meta.url), ...options});\n')
                (package/'package.json').write_text(json.dumps({'name':'linear-a-wasm','version':version,'type':'module','exports':'./index.mjs','files':['index.mjs','matrix.mjs','runtime','README.md','LICENSE'],'license':'MIT'}))
            shutil.copyfile(source/'LICENSE',package/'LICENSE')
            run(['npm','pack','--ignore-scripts','--pack-destination',output], package, env)
            package_file = next(output.glob('linear-a-'+name+'-*.tgz'))
            (consumer/'package.json').write_text('{"private":true,"type":"module"}')
            run(['npm','install','--ignore-scripts','--no-audit','--no-fund',package_file], consumer, env)
            if name == 'typescript':
                run([package/'node_modules/.bin/tsc', '--target', 'ES2022', '--module', 'NodeNext', '--strict', '--outDir', consumer/'dist', example], consumer, env)
                example = consumer/'dist/typescript.mjs'
            run(['node',example], consumer, env)
        elif name == 'rust':
            package = source/'ports/rust'
            shutil.copyfile(source/'LICENSE',package/'LICENSE')
            run(['cargo','package','--offline','--allow-dirty','--manifest-path',package/'Cargo.toml'], consumer, env)
            crate = next((package/'target/package').glob('*.crate'))
            shutil.copyfile(crate,output/crate.name)
            installed = base/'installed'
            unpack(crate,installed)
            (consumer/'Cargo.toml').write_text('[package]\nname="release-consumer"\nversion="0.0.0"\nedition="2021"\n[dependencies]\nlinear-a-rust={path='+json.dumps(str(installed/('linear-a-rust-'+version)))+'}\n[[bin]]\nname="consumer"\npath="rust.rs"\n')
            run(['cargo','run','--offline','--quiet'], consumer, env)
        elif name == 'go':
            (consumer/'go.mod').write_text('module release-consumer\n\ngo 1.22\nrequire github.com/jmcummings77/linear-A/ports/go v'+version+'\nreplace github.com/jmcummings77/linear-A/ports/go => '+str(source/'ports/go')+'\n')
            run(['go','run','.'], consumer, env)
        elif name == 'julia':
            (consumer/'Project.toml').write_text('[deps]\nLinearAMatrices = "a8a34c5f-5193-4f59-8cc0-d599ed85a1f5"\n')
            # A path manifest installs the extracted package without a registry lookup.
            (consumer/'Manifest.toml').write_text('manifest_format = "2.0"\n[[deps.LinearAMatrices]]\npath = '+json.dumps(str(source/'ports/julia'))+'\nuuid = "a8a34c5f-5193-4f59-8cc0-d599ed85a1f5"\nversion = "'+version+'"\n')
            env['JULIA_LOAD_PATH']='@:@stdlib'
            run(['julia','--startup-file=no','--project='+str(consumer),example], consumer, env)
        elif name in ('csharp','fsharp'):
            extension = 'cs' if name == 'csharp' else 'fs'
            project = source/('linear-A/linear-A/linear-A.csproj' if name=='csharp' else 'ports/fsharp/LinearA.FSharp.fsproj')
            readme = project.parent/'PACKAGE-README.md'
            shutil.copyfile(source/'docs/api'/version/'README.md',readme)
            project.write_text(project.read_text().replace('</Project>', '<PropertyGroup><PackageReadmeFile>PACKAGE-README.md</PackageReadmeFile></PropertyGroup><ItemGroup><None Include="PACKAGE-README.md" Pack="true" PackagePath="/" /></ItemGroup></Project>'))
            package_id = 'LinearA.' + ('CSharp' if name=='csharp' else 'FSharp')
            run(['dotnet','pack',project,'-c','Release','-o',output,'-p:PackageId='+package_id,'-p:Version='+version,'-p:Authors=J.M. Cummings','-p:PackageLicenseExpression=MIT','-p:IncludeSymbols=false','-p:DebugType=none','-p:NuGetAudit=false'], consumer, env)
            (consumer/('Consumer.'+extension+'proj')).write_text('<Project Sdk="Microsoft.NET.Sdk"><PropertyGroup><OutputType>Exe</OutputType><TargetFramework>net10.0</TargetFramework><NuGetAudit>false</NuGetAudit></PropertyGroup><ItemGroup>'+ ('<Compile Include="fsharp.fs" />' if name=='fsharp' else '') + '<PackageReference Include="'+package_id+'" Version="'+version+'" /></ItemGroup></Project>')
            env['NUGET_PACKAGES']=str(base/'nuget')
            run(['dotnet','restore','--source',output,'--source','https://api.nuget.org/v3/index.json'],consumer,env)
            run(['dotnet','run','-c','Release','--no-restore'],consumer,env)


def audit_assets(paths, private_paths=None):
    """Reject build-location leaks in source members, wheels, and compiled assets."""
    private_paths = private_paths or [str(ROOT), str(Path.home()), 'linear-a-consumer-']
    needles = [value.encode(encoding) for value in private_paths for encoding in ('utf-8','utf-16-le')]
    def check(data):
        if any(needle in data for needle in needles):
            raise ValueError('release asset contains an identifying local path')
    for path in paths:
        if path.name.endswith(('.tar.gz','.tgz','.crate')):
            with tarfile.open(path) as archive:
                for entry in archive:
                    check(entry.name.encode())
                    if entry.isfile(): check(archive.extractfile(entry).read())
        elif path.suffix in ('.whl','.nupkg'):
            with zipfile.ZipFile(path) as archive:
                for name in archive.namelist():
                    check(name.encode())
                    check(archive.read(name))
        elif path.is_file():
            check(path.read_bytes())


def main():
    manifest=json.loads(MANIFEST.read_text())
    version=manifest['version']
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=version,
                        help='print the release version from manifest.json and exit')
    parser.add_argument('--port', action='append', choices=list(manifest['ports']))
    parser.add_argument('--output', type=Path, default=None)
    parser.add_argument('--tag', help='require this release tag to match the version and a clean checkout')
    parser.add_argument('--source-only', action='store_true', help='prepare source kits without installation verification')
    parser.add_argument('--write-table', action='store_true')
    args=parser.parse_args()
    if args.write_table:
        (ROOT/'docs/api'/version/'compatibility.md').write_text(table(manifest)); return
    check_manifest(manifest)
    output=(args.output or ROOT/'.build/releases'/version).resolve()
    output.mkdir(parents=True,exist_ok=True)
    if any(output.iterdir()): parser.error('output directory must be empty (prevents mixing releases)')
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True))
    if args.tag and (args.tag != 'v'+version or dirty or args.source_only or args.port):
        parser.error('tagged releases require a matching version, clean checkout, all ports and installation checks')
    ports=args.port or list(manifest['ports'])
    for name in ports:
        print('Packaging '+name,flush=True)
        port=manifest['ports'][name]
        archive=output/('linear-a-'+name+'-'+version+'-source.tar.gz')
        write_archive(archive,kit_files(name,port,version,{'revision':revision,'dirty':dirty}))
        if not args.source_only: smoke(name,port,archive,output,version)
    audit_assets(output.iterdir())
    assets=[{'file':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(output.iterdir()) if p.is_file()]
    (output/'release.json').write_text(json.dumps({'version':version,'revision':revision,'dirty':dirty,'created_utc':datetime.now(timezone.utc).isoformat(),'ports':ports,'consumer_checks':'not-run' if args.source_only else 'passed','assets':assets},indent=2)+'\n')
    index = output/'release.json'
    checksums = assets + [{'file':index.name, 'sha256':hashlib.sha256(index.read_bytes()).hexdigest()}]
    (output/'SHA256SUMS').write_text(''.join(a['sha256']+'  '+a['file']+'\n' for a in checksums))
    print('Release assets ready: '+str(output))

if __name__ == '__main__': main()
