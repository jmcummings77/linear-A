import unittest
from matrix import Matrix

class PseudoinverseTests(unittest.TestCase):
    def test_missing_rhs_cannot_request_an_inverse_accidentally(self):
        a=Matrix(1,2,[1,1])
        for rhs in [None, [], 0]:
            with self.assertRaises(TypeError):a.solve_minimum_norm(rhs)
        self.assertEqual(a.solve_minimum_norm(Matrix(1,1,[0])).values,[0,0])

    def test_cutoff_and_scaling_failures_are_explicit(self):
        a=Matrix(2,2,[4,0,0,1])
        self.assertEqual(a.pseudoinverse(.25).values,[.25,0,0,0])
        self.assertEqual(a.spectral_diagnostics(.25).rank,1)
        with self.assertRaises(ArithmeticError):a.solve_minimum_norm(Matrix(2,1,[1e300,1e-300]))
