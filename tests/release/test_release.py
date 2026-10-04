import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
import zipfile
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('release_build',ROOT/'release/build.py')
build=importlib.util.module_from_spec(spec)
spec.loader.exec_module(build)

class ReleaseTests(unittest.TestCase):
    def test_manifest_versions_examples_and_checked_table(self):
        manifest=json.loads(build.MANIFEST.read_text())
        self.assertEqual(len(manifest['ports']),11)
        build.check_manifest(manifest)

    def test_archive_is_reproducible_and_has_no_host_metadata(self):
        with tempfile.TemporaryDirectory() as d:
            a,b=Path(d)/'a.tar.gz',Path(d)/'b.tar.gz'
            files={'z/module.py':b'pass\n','LICENSE':b'license'}
            build.write_archive(a,files);build.write_archive(b,files)
            self.assertEqual(a.read_bytes(),b.read_bytes())
            with tarfile.open(a) as archive:
                for member in archive:
                    self.assertEqual((member.uid,member.gid,member.uname,member.gname,member.mtime),(0,0,'','',0))
            build.unpack(a,Path(d)/'installed')
            self.assertEqual((Path(d)/'installed/z/module.py').read_bytes(),b'pass\n')

    def test_unsafe_archive_members_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            for name in ['/escape','../escape']:
                with self.assertRaises(ValueError):
                    build.write_archive(Path(d)/'bad.tar.gz',{name:b'x'})
            archive=Path(d)/'symlink.tar'
            with tarfile.open(archive,'w') as output:
                member=tarfile.TarInfo('link');member.type=tarfile.SYMTYPE;member.linkname='/tmp'
                output.addfile(member)
            with self.assertRaises(ValueError): build.unpack(archive,Path(d)/'out')

    def test_python_package_includes_lazy_solver_dependency(self):
        # The solver is imported lazily; basic import tests missed its omission.
        metadata=(ROOT/'ports/python/pyproject.toml').read_text()
        self.assertIn('"solve"',metadata)
        self.assertIn('"ilu"',metadata)

    def test_audit_rejects_paths_inside_compiled_packages(self):
        with tempfile.TemporaryDirectory() as d:
            package=Path(d)/'package.nupkg'
            with zipfile.ZipFile(package,'w') as archive:
                archive.writestr('lib/library.dll','/home/private-user/source'.encode('utf-16-le'))
            with self.assertRaises(ValueError):
                build.audit_assets([package], ['/home/private-user'])
            with zipfile.ZipFile(package,'w') as archive:
                archive.writestr('lib/library.dll',b'portable library')
            build.audit_assets([package], ['/home/private-user'])

    def test_tagged_release_rejects_unverified_source_only_build(self):
        with tempfile.TemporaryDirectory() as d:
            result=subprocess.run([sys.executable,str(ROOT/'release/build.py'),
                '--tag','v0.10.0','--source-only','--output',d],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('tagged releases require',result.stderr)
            self.assertFalse((Path(d)/'release.json').exists())

    def test_native_kits_include_shared_headers_and_license(self):
        manifest=json.loads(build.MANIFEST.read_text())
        for port in ['c','cpp','assembly','wasm']:
            files=build.kit_files(port,manifest['ports'][port],manifest['version'])
            self.assertIn('ports/c/solve_core.h',files)
            self.assertIn('ports/c/sparse_core.h',files)
            self.assertIn('ports/c/ilu_core.h',files)
            self.assertIn('ports/SPARSE.md',files)
            self.assertIn('ports/c/general_eigen.h',files)
            self.assertIn('LICENSE',files)
            self.assertIn('docs/api/0.10.0/README.md',files)
            self.assertFalse(any('.build/' in name or 'node_modules/' in name or '.git/' in name for name in files))

if __name__=='__main__': unittest.main()
