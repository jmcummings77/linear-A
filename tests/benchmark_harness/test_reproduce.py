"""Publication captures measure a fresh, stable committed source tree."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'benchmarks'))
import reproduce


RUNNER = '''import argparse, json, subprocess
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument('--output', type=Path, required=True)
p.add_argument('--require-all', action='store_true')
p.add_argument('--change-source', action='store_true')
p.add_argument('--fail', action='store_true')
args = p.parse_args()
assert args.require_all
assert not Path('.build/stale-binary').exists(), 'Reused the caller build'
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
data = dict(schema_version=1, revision=revision, dirty=False, source_sha256='a'*64,
            implementations=[dict(id='python', status='passed')], results=[], profiles=[])
args.output.mkdir(parents=True)
(args.output / 'results.json').write_text(json.dumps(data))
(args.output / 'index.html').write_text('<script type="application/json" id="data">' + json.dumps(data) + '</script>')
if args.change_source: Path('tracked.txt').write_text('modified while measuring')
if args.fail: raise SystemExit(1)
'''


class CaptureTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.root = self.directory / 'repository'
        self.root.mkdir()
        self.output = self.directory / 'output'
        self.git('init', '-q')
        self.git('config', 'user.name', 'Benchmark test')
        self.git('config', 'user.email', 'benchmark@example.invalid')
        (self.root / 'benchmarks').mkdir()
        (self.root / 'benchmarks/run.py').write_text(RUNNER)
        (self.root / '.gitignore').write_text('.build/\n')
        (self.root / 'tracked.txt').write_text('committed source')
        self.git('add', '.')
        self.git('commit', '-qm', 'Fixture')

    def git(self, *arguments):
        return reproduce.git(self.root, *arguments)

    def test_capture_uses_fresh_build_tree_and_records_exact_source_and_command(self):
        (self.root / '.build').mkdir()
        (self.root / '.build/stale-binary').write_text('stale')
        before = self.git('status', '--porcelain')
        manifest = reproduce.capture(self.output, root=self.root)
        self.assertEqual(manifest['revision'], self.git('rev-parse', 'HEAD'))
        self.assertEqual(manifest['source_tree'], self.git('rev-parse', 'HEAD^{tree}'))
        self.assertEqual(manifest['command'],
                         ['python3', 'benchmarks/run.py', '--require-all', '--output', '<output>'])
        self.assertEqual(manifest, json.loads((self.output / 'provenance.json').read_text()))
        self.assertEqual(before, self.git('status', '--porcelain'))
        self.assertEqual((self.root / '.build/stale-binary').read_text(), 'stale')

    def test_staged_unstaged_and_untracked_source_changes_are_rejected(self):
        for state in ('unstaged', 'staged', 'untracked'):
            with self.subTest(state=state):
                if state == 'untracked':
                    (self.root / 'new-source.py').write_text('new source')
                else:
                    (self.root / 'tracked.txt').write_text('changed')
                    if state == 'staged':
                        self.git('add', 'tracked.txt')
                with self.assertRaisesRegex(ValueError, 'clean checkout'):
                    reproduce.capture(self.output, root=self.root)
                self.assertFalse(self.output.exists())
                self.git('reset', '--hard', 'HEAD')
        self.assertTrue((self.root / 'new-source.py').exists())

    def test_source_mutation_or_failed_run_never_promotes_partial_output(self):
        for argument, error in (('--change-source', ValueError), ('--fail', subprocess.CalledProcessError)):
            with self.subTest(argument=argument), self.assertRaises(error):
                reproduce.capture(self.output, arguments=[argument], root=self.root)
            self.assertFalse(self.output.exists())
            self.assertEqual((self.root / 'tracked.txt').read_text(), 'committed source')

    def test_existing_output_is_not_overwritten(self):
        self.output.mkdir()
        sentinel = self.output / 'results.json'
        sentinel.write_text('keep existing results')
        with self.assertRaisesRegex(ValueError, 'never overwritten'):
            reproduce.capture(self.output, root=self.root)
        self.assertEqual(sentinel.read_text(), 'keep existing results')

    def test_harness_option_abbreviations_cannot_bypass_fresh_capture(self):
        for option in ('--no-build', '--no-b', '--no', '--out', '--output=elsewhere',
                       '--render-o', '--render-only=old.json', '--help', '-h', '--amd-lib=x', '--'):
            with self.subTest(option=option), self.assertRaises(ValueError):
                reproduce.checked_arguments('run.py', [option])


if __name__ == '__main__':
    unittest.main()
