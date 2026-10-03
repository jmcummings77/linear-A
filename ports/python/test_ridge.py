import math
import unittest
from matrix import Matrix

class RidgeTests(unittest.TestCase):
    def test_options_and_owned_result(self):
        a,b=Matrix(1,1,[2]),Matrix(1,1,[3])
        for value in [-1,math.nan,math.inf]:
            with self.assertRaises(ValueError):a.solve_ridge(b,value)
        with self.assertRaises(TypeError):a.solve_ridge(None,1)
        self.assertAlmostEqual(a.solve_ridge(b,1)[0,0],1.2)
        self.assertEqual(a.values,[2]);self.assertEqual(b.values,[3])
        self.assertEqual(a.solve_ridge(b,0).values,a.solve_minimum_norm(b).values)

    def test_scaled_filter_preserves_tiny_representable_result(self):
        self.assertAlmostEqual(Matrix(1,1,[1e-300]).solve_ridge(Matrix(1,1,[1e300]),1e300)[0,0]/1e-300,1)
