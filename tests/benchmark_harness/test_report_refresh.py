"""Presentation refresh preserves public metadata regardless of the current host."""
import contextlib
import io
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'benchmarks'))
import amd_bench
import cholesky_bench
import ic0_bench
import multigrid_bench
import ordering_bench
import refresh_reports
import report
import report_design
import sparse
from publication import PublicSanitizer

ROOT = Path(__file__).resolve().parents[2]


class PublishedReportRefreshTests(unittest.TestCase):
    def test_all_saved_reports_refresh_unchanged_for_ci_runner_account(self):
        # Exercise the real dispatch, including experiments and the comparison
        # catalog, but put every output in a disposable copy of the saved pages.
        with tempfile.TemporaryDirectory() as directory, \
             patch('publication.Path.home', return_value=Path('/home/runner')), \
             patch.object(report, 'live_bundle', side_effect=AssertionError('Must retain saved runtime')), \
             patch.object(sparse, 'live_bundle', side_effect=AssertionError('Must retain saved runtime')):
            destination = Path(directory)
            for relative in ('benchmarks/reports', 'experiments/matmul-locality/results',
                             'experiments/machine-code-dot/results'):
                shutil.copytree(ROOT / relative, destination / relative)
            for relative in ('applications/index.html', 'benchmarks/comparison/index.html'):
                target = destination / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(ROOT / relative, target)
            (destination / 'benchmarks/pages').mkdir()
            modules = {relative: refresh_reports.module(name, relative) for name, relative in (
                ('refresh_applications', 'applications/build.py'),
                ('refresh_locality', 'experiments/matmul-locality/run.py'),
                ('refresh_dot', 'experiments/machine-code-dot/run.py'))}
            with patch.object(refresh_reports, 'ROOT', destination), \
                 patch.object(refresh_reports, 'module', side_effect=lambda name, relative: modules[relative]), \
                 patch.object(refresh_reports, 'render_directory', side_effect=lambda:
                     report_design.render_directory(destination / 'benchmarks/pages/index.html')), \
                 contextlib.redirect_stdout(io.StringIO()):
                # refresh() asserts exact equality of all embedded data and
                # executable bundles and unchanged measurement file hashes.
                refresh_reports.refresh()

    def test_fresh_inputs_keep_host_sanitization_even_with_a_saved_runtime(self):
        data = {'schema_version': 1, 'machine': {}, 'methodology': {},
                'implementations': [{'id': 'c', 'name': 'C', 'status': 'passed',
                    'toolchain': 'private compiler for runner',
                    'build_commands': [['/home/runner/bin/clang', '/home/runner/private.c']]}],
                'results': [], 'profiles': []}
        renderers = (report, sparse, ordering_bench, cholesky_bench, amd_bench,
                     ic0_bench, multigrid_bench)
        with tempfile.TemporaryDirectory() as directory, \
             patch('publication.Path.home', return_value=Path('/home/runner')):
            for renderer in renderers:
                with self.subTest(renderer=renderer.__name__):
                    target = Path(directory) / (renderer.__name__ + '.html')
                    renderer.render(data, target, live_override={'available': False})
                    saved = refresh_reports.embedded(target, 'data')
                    implementation = saved['implementations'][0]
                    self.assertEqual(implementation['toolchain'], 'private compiler for [user]')
                    self.assertEqual(implementation['build_commands'], [['clang', 'private.c']])

    def test_published_mode_keeps_generic_path_and_credential_redaction(self):
        with patch('publication.Path.home', return_value=Path('/home/runner')):
            sanitizer = PublicSanitizer.for_published()
            self.assertEqual(sanitizer.text('benchmarks/dotnet/Runner.csproj'),
                             'benchmarks/dotnet/Runner.csproj')
            self.assertEqual(sanitizer.text('/home/private/secret.c'), 'secret.c')
            self.assertEqual(sanitizer.text('password=secret-value'), '[redacted credential]')


if __name__ == '__main__':
    unittest.main()
