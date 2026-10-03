import math
import unittest
from matrix import Matrix

class SvdTests(unittest.TestCase):
    def test_options_nonconvergence_and_nonfinite(self):
        a=Matrix(3,3,[1,2,3,4,5,7,6,8,9])
        for tolerance,sweeps in [(0,100),(math.nan,100),(1,100),(1e-12,0),(1e-12,10001),(1e-12,True)]:
            with self.assertRaises(ValueError): a.svd(tolerance,sweeps)
        with self.assertRaises(ArithmeticError): a.svd(max_sweeps=1)
        a._values[0]=math.nan
        with self.assertRaises(ValueError):a.svd()

    def test_results_are_owned_and_input_is_unchanged(self):
        a=Matrix(2,3,[3,0,0,0,4,0]);original=a.values
        r=a.svd()
        self.assertEqual(a.values,original)
        r.u._values[0]=99;r.values[0]=99;r.vt._values[0]=99
        self.assertEqual(a.values,original)
        self.assertEqual(a.svd().values,[4,3])

    def test_rank_deficient_null_basis(self):
        r=Matrix(8,8,[1.]*64).svd()
        self.assertAlmostEqual(r.values[0],8)
        for p in range(8):
            for q in range(8):
                self.assertAlmostEqual(sum(r.u.values[i*8+p]*r.u.values[i*8+q] for i in range(8)),float(p==q),places=10)

if __name__=='__main__':unittest.main()
