import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'benchmarks'))
from solve_reference import exact_solve, exact_least_squares, solve_fixtures
from reference import matrix

class SolveReferenceTests(unittest.TestCase):
    def test_gauss_jordan_hand_solution_and_pivot(self):
        self.assertEqual(exact_solve(matrix(2,2,[0,2,1,3]),matrix(2,2,[4,2,7,4]))['values'],[1,1,2,1])
        with self.assertRaises(ValueError):
            exact_solve(matrix(2,2,[1,2,2,4]),matrix(2,1,[1,2]))

    def test_exact_fit_with_noisy_rhs(self):
        # Best constant is the arithmetic mean; QR production never forms AᵀA.
        self.assertEqual(exact_least_squares(matrix(3,1,[1,1,1]),matrix(3,1,[1,2,6]))['values'],[3])

    def test_fixture_names_and_reference_shapes(self):
        cases=solve_fixtures()
        self.assertEqual(len({c['name'] for c in cases}),len(cases))
        for case in cases:
            if not case['invalid'] and case['op']!='rcond':
                self.assertEqual(case['expected']['rows'],case['a']['cols'])
                self.assertEqual(case['expected']['cols'],case['b']['cols'])
