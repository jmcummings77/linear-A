"""End-to-end checks of the independently checkable nonnormal experiment."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'experiments/nonnormal-gmres/run.py'


class NonnormalGmresTests(unittest.TestCase):
    def test_worked_example_restart_loss_and_preconditioning(self):
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run([sys.executable, str(SCRIPT), '--output', directory],
                           check=True, capture_output=True, text=True)
            data = json.loads((Path(directory)/'results.json').read_text())
            report = (Path(directory)/'RESULTS.md').read_text()
        self.assertEqual(len(data['runs']), 20)
        rows = {(r['s'], r['restart'], r['right_jacobi']): r for r in data['runs']}
        plain = rows[.25, 1, False]
        self.assertAlmostEqual(plain['history'][1]['x'][1], 4/5, places=13)
        self.assertAlmostEqual(plain['history'][2]['x'][0], -32/65, places=13)
        self.assertAlmostEqual(plain['history'][2]['x'][1], 68/65, places=13)
        self.assertAlmostEqual(plain['history'][2]['true_residual'], 13**.5/65, places=13)
        # At step two GMRES(1) is in a new cycle, not a degree-two Krylov space.
        self.assertEqual(plain['history'][2]['cycle_start'], 1)
        self.assertEqual(plain['history'][2]['local_step'], 1)
        self.assertAlmostEqual(plain['history'][2]['cycle_envelope'], .5/(5**.5))
        for s in (0, .25, .75, 2, 4):
            for jacobi in (False, True):
                self.assertEqual(rows[s, 2, jacobi]['reason'], 'converged')
                self.assertLessEqual(rows[s, 2, jacobi]['iterations'], 2)
            a, b = rows[s, 1, False], rows[s, 1, True]
            self.assertEqual(a['iterations'], b['iterations'])
            for x, y in zip(a['history'], b['history']):
                self.assertAlmostEqual(x['true_residual'], y['true_residual'], places=13)
        self.assertNotEqual(rows[4, 1, False]['reason'], 'converged')
        self.assertGreater(rows[4, 1, False]['history'][-1]['true_residual'], .99)
        self.assertIn('not a certified error bound or stopping criterion', report)
        self.assertNotIn(str(ROOT), json.dumps(data))

    def test_verification_does_not_write_reports_and_bounded_run_keeps_limit(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory)/'absent'
            subprocess.run([sys.executable, str(SCRIPT), '--verify-only', '--limit', '2',
                            '--output', str(destination)], check=True, capture_output=True)
            self.assertFalse(destination.exists())
            subprocess.run([sys.executable, str(SCRIPT), '--limit', '2',
                            '--output', str(destination)], check=True, capture_output=True)
            data = json.loads((destination/'results.json').read_text())
        for row in data['runs']:
            if row['restart'] == 1 and row['s']:
                self.assertEqual(row['reason'], 'iteration_limit')
                self.assertEqual(row['iterations'], 2)

    def test_independent_projection_check_rejects_self_consistent_wrong_step(self):
        # Isolate port imports from identically named benchmark modules.
        code = '''import importlib.util
from unittest.mock import patch
spec = importlib.util.spec_from_file_location('experiment', %r)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
original = module.gmres
def broken(*args, **kwargs):
    result = original(*args, **kwargs)
    result.iterates[1][1] = 0.0
    result.residuals[1] = 1.0
    result.estimated_residuals[1] = 1.0
    return result
with patch.object(module, 'gmres', broken):
    try:
        module.trajectory(.25, 2, False)
    except AssertionError as error:
        assert 'exact minimal-residual projection' in str(error)
    else:
        raise AssertionError('corrupted projection was accepted')
''' % str(SCRIPT)
        subprocess.run([sys.executable, '-c', code], check=True, capture_output=True, text=True)

    def test_invalid_limit_is_rejected_before_output(self):
        result = subprocess.run([sys.executable, str(SCRIPT), '--limit', '101'],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)


if __name__ == '__main__':
    unittest.main()
