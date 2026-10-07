"""Analytic graph metadata and independent checks in the untimed report."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    'nonnormal_graph_report', ROOT / 'experiments/nonnormal-gmres/graph_report.py')
report = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(report)


class NonnormalGraphReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = report.graph_diagnostics()
        cls.cases = {case['id']: case for case in cls.data['cases']}

    def test_known_spectra_full_conditioning_and_ramanujan_scope(self):
        self.assertEqual(set(self.cases),
                         {'petersen', 'triangular-prism', 'prism-8', 'prism-16', 'prism-32'})
        petersen = self.cases['petersen']
        self.assertEqual((petersen['n'], petersen['nnz'], petersen['edge_count']), (10, 40, 15))
        self.assertEqual(petersen['adjacency_spectrum'], [-2.0]*4 + [1.0]*5 + [3.0])
        self.assertEqual(petersen['spectral_gap'], 2.0)
        self.assertEqual(petersen['kappa2'], 21.0)
        triangular = self.cases['triangular-prism']
        for actual, expected in zip(triangular['adjacency_spectrum'], [-2, -2, 0, 0, 1, 3]):
            self.assertAlmostEqual(actual, expected, places=13)
        self.assertAlmostEqual(triangular['spectral_gap'], 2.0)
        self.assertAlmostEqual(triangular['kappa2'], 21.0)
        gaps = []
        for m in (8, 16, 32):
            case = self.cases['prism-%d' % m]
            self.assertEqual((case['n'], case['nnz']), (2*m, 8*m))
            self.assertTrue(case['bipartite'])
            self.assertAlmostEqual(case['adjacency_min'], -3.0)
            self.assertEqual(case['shifted_lambda_min'], .25)
            self.assertAlmostEqual(case['shifted_lambda_max'], 6.25)
            self.assertAlmostEqual(case['kappa2'], 25.0)
            gaps.append(case['spectral_gap'])
        self.assertGreater(gaps[0], gaps[1])
        self.assertGreater(gaps[1], gaps[2])
        self.assertAlmostEqual(gaps[2], .03842943919353914, places=13)
        self.assertEqual([self.cases[key]['is_ramanujan'] for key in
                          ('petersen', 'triangular-prism', 'prism-8', 'prism-16', 'prism-32')],
                         [True, True, True, False, False])
        self.assertIn('constant eigenvector', self.data['scope']['conditioning'])
        self.assertIn('scalar scaling', self.data['scope']['jacobi'])
        self.assertIn('not a Ramanujan family', self.data['scope']['graphs'])

    def test_actual_solutions_and_residuals_are_consistent(self):
        for case in self.data['cases']:
            with self.subTest(graph=case['id']):
                self.assertGreater(case['energy'], 0)
                self.assertAlmostEqual(case['energy_from_rhs'], case['energy'], places=12)
                self.assertEqual([solve['jacobi'] for solve in case['solves']], [False, True])
                truth = case['manufactured_solution']
                self.assertGreater(len(set(truth)), 1)
                rhs = case['rhs']
                csr = case['csr']
                for solve in case['solves']:
                    self.assertEqual(solve['reason'], 'converged')
                    self.assertTrue(solve['converged'])
                    self.assertLessEqual(solve['iterations'], 2000)
                    self.assertEqual(len(solve['history']), solve['iterations'] + 1)
                    # Recompute from the serialized matrix to cross-check the
                    # report's independent edge-based physical residual.
                    x = solve['x']
                    product = [sum(csr['values'][p] * x[csr['indices'][p]]
                                   for p in range(csr['offsets'][i], csr['offsets'][i+1]))
                               for i in range(case['n'])]
                    residual = math.hypot(*(b - ax for b, ax in zip(rhs, product)))
                    self.assertLessEqual(abs(residual - solve['true_residual']),
                                         solve['rounding_allowance'])
                    self.assertLessEqual(residual, solve['stopping_threshold'] +
                                         solve['rounding_allowance'])
                    error = math.hypot(*(a - b for a, b in zip(x, truth)))
                    self.assertAlmostEqual(error, solve['solution_error'], places=15)
                    self.assertLess(solve['relative_solution_error'], 1e-8)
        encoded = json.dumps(self.data, allow_nan=False)
        self.assertEqual(json.loads(encoded), self.data)
        self.assertNotIn(str(ROOT), encoded)
        self.assertNotIn('elapsed', encoded)

    def test_corrupt_csr_cannot_pass_using_its_own_product(self):
        original = report.graph_reference.shifted_laplacian
        def broken(*args, **kwargs):
            csr = copy.deepcopy(original(*args, **kwargs))
            csr['values'][0] += 1.0
            return csr
        with patch.object(report.graph_reference, 'shifted_laplacian', broken):
            with self.assertRaisesRegex(AssertionError, 'independent graph edge action'):
                report.graph_diagnostics()

    def test_unexpected_nonconvergence_fails_the_report(self):
        original = report.CSRMatrix.conjugate_gradient
        def broken(*args, **kwargs):
            return original(*args, **kwargs)._replace(converged=False, reason='iteration_limit')
        with patch.object(report.CSRMatrix, 'conjugate_gradient', broken):
            with self.assertRaisesRegex(AssertionError, 'unexpectedly failed to converge'):
                report.graph_diagnostics()

    def test_self_consistent_zero_claim_fails_independent_residual_check(self):
        original = report.CSRMatrix.conjugate_gradient
        def broken(*args, **kwargs):
            result = original(*args, **kwargs)
            result.x[:] = [0.0] * len(result.x)
            result.iterates[-1][:] = result.x
            result.residuals[-1] = 0.0
            return result
        with patch.object(report.CSRMatrix, 'conjugate_gradient', broken):
            with self.assertRaisesRegex(AssertionError, 'independent graph true residual'):
                report.graph_diagnostics()


if __name__ == '__main__':
    unittest.main()
