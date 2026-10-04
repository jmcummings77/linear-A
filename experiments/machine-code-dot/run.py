#!/usr/bin/env python3
"""Isolated ARM64 scalar dot-product experiment. Python standard library only."""
import argparse
import html
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import shutil
import statistics
import struct
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent
BUILD = ROOT / '.build'
SRC = ROOT / 'src'
RESULTS = ROOT / 'results'
sys.path.insert(0, str(ROOT.parents[1] / 'benchmarks'))
LABELS = {'machine-code': 'Hand-encoded machine code', 'assembly': 'Assembly (identical instructions)',
          'c': 'C', 'cpp': 'C++', 'rust': 'Rust', 'go': 'Go', 'csharp': 'C#', 'fsharp': 'F#',
          'typescript': 'TypeScript / Node', 'python': 'Python', 'julia': 'Julia'}
ENV = dict(os.environ, DOTNET_CLI_HOME=str(BUILD / 'cli'), NUGET_PACKAGES=str(BUILD / 'nuget'),
           DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_GENERATE_ASPNET_CERTIFICATE='false',
           DOTNET_TieredCompilation='0', DOTNET_NOLOGO='1', GOCACHE=str(BUILD / 'go-cache'),
           GOPATH=str(BUILD / 'go'), GOTOOLCHAIN='local', JULIA_NUM_THREADS='1',
           JULIA_DEPOT_PATH=str(BUILD / 'julia-depot'))

def command(args, timeout=120):
    result = subprocess.run([str(x) for x in args], cwd=ROOT, env=ENV,
                            text=True, capture_output=True, timeout=timeout)
    if result.returncode:
        raise RuntimeError('Command failed: ' + ' '.join(map(str, args)) + '\n' + result.stdout + result.stderr)
    return result.stdout.strip()

def executable(name, override=None):
    found = override or shutil.which(name)
    if not found:
        raise RuntimeError('Missing tool: ' + name)
    return str(Path(found).resolve())

def build(args, output=RESULTS):
    if platform.system() != 'Darwin' or platform.machine() != 'arm64':
        raise RuntimeError('This experiment targets ARM64 macOS.')
    BUILD.mkdir(exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    for directory in ('cli', 'nuget', 'go-cache', 'go', 'julia-depot'):
        (BUILD / directory).mkdir(exist_ok=True)
    clang, cpp, rust, go, node = [executable(n) for n in ('clang', 'clang++', 'rustc', 'go', 'node')]
    dotnet = executable('dotnet', args.dotnet)
    julia = executable('julia', args.julia)
    versions = {key: command(cmd).splitlines()[0] for key, cmd in {
        'clang': [clang, '--version'], 'rust': [rust, '--version'], 'go': [go, 'version'],
        'node': [node, '--version'], 'python': [sys.executable, '--version'],
        'dotnet': [dotnet, '--version'], 'julia': [julia, '--version']}.items()}
    flags = ['-O3', '-ffp-contract=off', '-Wall', '-Wextra']
    if not args.no_build:
        command([clang, *flags, '-DRAW_MACHINE_CODE', SRC / 'native_runner.c', '-o', BUILD / 'machine-code'])
        command([clang, *flags, '-c', SRC / 'native_runner.c', '-o', BUILD / 'harness.o'])
        command([clang, *flags, '-c', SRC / 'dot.c', '-o', BUILD / 'c.o'])
        command([cpp, *flags, '-std=c++17', '-c', SRC / 'dot.cpp', '-o', BUILD / 'cpp.o'])
        command([clang, '-c', SRC / 'assembly.S', '-o', BUILD / 'assembly.o'])
        for key in ('c', 'cpp', 'assembly'):
            command([clang, BUILD / 'harness.o', BUILD / (key + '.o'), '-o', BUILD / key])
        command([rust, '-C', 'opt-level=3', '-C', 'target-cpu=native', SRC / 'dot.rs', '-o', BUILD / 'rust'])
        command([go, 'build', '-o', BUILD / 'go-runner', SRC / 'dot.go'])
        sdk_version = versions['dotnet']
        sdk_source = Path(dotnet).parent / 'sdk' / sdk_version / 'FSharp' / 'library-packs'
        for language, extension in (('csharp', 'csproj'), ('fsharp', 'fsproj')):
            cmd = [dotnet, 'build', ROOT / language / ('Dot.' + extension), '-c', 'Release',
                   '--disable-build-servers', '-m:1', '-p:NuGetAudit=false']
            # FSharp.Core ships with the SDK; this avoids an unnecessary network dependency.
            if sdk_source.is_dir():
                cmd.extend(['--source', sdk_source])
            command(cmd)
    repository_prefix = str(ROOT.parents[1]) + os.sep
    disassembly = command(['xcrun', 'llvm-objdump', '--disassemble', BUILD / 'assembly.o']).replace(repository_prefix, '')
    assembled_words = [int(x, 16) for x in re.findall(r'^\s*[0-9a-f]+:\s+([0-9a-f]{8})\s', disassembly, re.M)]
    manual_words = [int(x, 16) for x in re.findall(r'0x([0-9a-f]{8})', (SRC / 'raw_dot.h').read_text())]
    if assembled_words != manual_words or len(manual_words) != 9:
        raise RuntimeError('Assembly control does not match the hand-encoded 36-byte function.')
    (output / 'raw_dot.bin').write_bytes(struct.pack('<9I', *manual_words))
    (output / 'disassembly.txt').write_text(disassembly)
    (output / 'c-disassembly.txt').write_text(command(['xcrun', 'llvm-objdump', '--disassemble', BUILD / 'c.o']).replace(repository_prefix, ''))
    runners = {key: [str(BUILD / key)] for key in ('machine-code', 'assembly', 'c', 'cpp', 'rust')}
    runners.update({'go': [str(BUILD / 'go-runner')], 'typescript': [node, str(SRC / 'dot.ts')],
                    'python': [sys.executable, str(SRC / 'dot.py')],
                    'julia': [julia, '--startup-file=no', '--history-file=no', str(SRC / 'dot.jl')]})
    for key in ('csharp', 'fsharp'):
        runners[key] = [dotnet, str(ROOT / key / 'bin' / 'Release' / 'net10.0' / 'Dot.dll')]
    return runners, versions, manual_words

_REFERENCE = {}
def numerators(n, seed):
    key = n, seed
    if key not in _REFERENCE:
        # Exact integer oracle: every input is an integer divided by 16.
        _REFERENCE[key] = [sum(((i * 17 + seed * 13) % 101 - 50) *
                                  ((i * 29 + (seed + 1) * 7) % 103 - 51)
                                  for i in range(batch*n, (batch+1)*n)) for batch in range(32)]
    return _REFERENCE[key]

def invoke(runner, n, iterations, seed):
    data = json.loads(command(runner + [str(n), str(iterations), str(seed)]))
    parts = numerators(n, seed)
    expected = ((iterations // 32) * sum(parts) + sum(parts[:iterations % 32])) / 256
    if data['iterations'] != iterations or not math.isfinite(data['checksum']):
        raise RuntimeError('Invalid benchmark output: ' + repr(data))
    if data['checksum'] != expected:
        raise RuntimeError('Exact checksum mismatch: expected %r, got %r' % (expected, data['checksum']))
    if not math.isfinite(data['elapsed_ns']) or data['elapsed_ns'] < 0:
        raise RuntimeError('Invalid elapsed time')
    return data

def render(result, destination):
    from report_design import apply_report_design
    data = json.dumps(result, allow_nan=False).replace('<', '\\u003c').replace('&', '\\u0026')
    template = (ROOT / 'report.html').read_text()
    replacements = {'__DATA__': data, '__NAMES__': json.dumps(LABELS), '__SCRIPT__': (ROOT / 'report.mjs').read_text()}
    page = re.sub('|'.join(replacements), lambda match: replacements[match[0]], template)
    Path(destination).write_text(apply_report_design(page))


def report(result, output=RESULTS):
    render(result, output / 'index.html')
    (output / 'results.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dotnet')
    parser.add_argument('--julia')
    parser.add_argument('--no-build', action='store_true')
    parser.add_argument('--samples', type=int, default=5)
    parser.add_argument('--sizes', type=int, nargs='+', default=[16, 256, 4096, 65536])
    parser.add_argument('--output', type=Path, default=RESULTS)
    args = parser.parse_args()
    if args.samples < 3 or any(n < 1 or n > 1048576 for n in args.sizes): parser.error('At least 3 samples and lengths 1..1048576 required')
    revision = command(['git', 'rev-parse', 'HEAD'])
    dirty = bool(command(['git', 'status', '--porcelain', '--untracked-files=all']))
    runners, versions, words = build(args, args.output)
    checks = 0
    for key, runner in runners.items():
        for n, count, seed in ((0,33,7),(1,35,0),(3,67,7),(16,97,13),(255,65,991),(4096,33,7)):
            invoke(runner,n,count,seed); checks += 1
        print('Verified ' + LABELS[key], flush=True)
    timings = []
    rng = random.Random(7312)
    for n in args.sizes:
        iterations, samples = {}, {key: [] for key in runners}
        order = list(runners); rng.shuffle(order)
        for key in order:
            initial = invoke(runners[key],n,32,7)
            estimate = max(1, initial['elapsed_ns']/32)
            iterations[key] = max(32, min(10000000, int(30000000/estimate)))
        for sample in range(args.samples):
            rng.shuffle(order)
            for key in order:
                value = invoke(runners[key],n,iterations[key],7)
                checks += 1
                samples[key].append(value['elapsed_ns']/iterations[key])
            print('Length %d: sample %d/%d complete' % (n,sample+1,args.samples),flush=True)
        for key in runners:
            timings.append(dict(implementation=key,n=n,iterations=iterations[key],
                                ns_per_call=samples[key],median_ns=statistics.median(samples[key])))
        (args.output / 'timing-checkpoint.json').write_text(json.dumps(timings, indent=2))
    result = dict(revision=revision, dirty=dirty,
                  metadata=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),
                               os=platform.platform(),architecture=platform.machine(),versions=versions,
                               cpu=platform.processor() or platform.machine(),
                               compiler_flags='Native: -O3 -ffp-contract=off; Rust: opt-level=3 target-cpu=native',
                               dotnet_tiered_compilation=False,seed=7,order_seed=7312,
                               scope='Separate scalar dot kernels, not full matrix-library comparisons'),
                  correctness_checks=checks,samples=args.samples,
                  words=['0x%08x'%word for word in words],timings=timings)
    report(result, args.output)
    print('Report: '+str(args.output/'index.html'),flush=True)
    for n in args.sizes:
        print('\nLength %d'%n)
        for row in sorted((r for r in timings if r['n']==n),key=lambda r:r['median_ns']):
            print('%-32s %12.2f ns/call' % (LABELS[row['implementation']], row['median_ns']))

if __name__ == '__main__':
    main()
