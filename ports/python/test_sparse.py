import unittest
from matrix import Matrix, CSRMatrix

class SparseTests(unittest.TestCase):
    def test_owned_canonical_arrays(self):
        rp,ci,v=[0,1],[0],[2]
        a=CSRMatrix(1,1,rp,ci,v)
        rp[1]=0;ci[0]=3;v[0]=9
        a.values[0]=99
        self.assertEqual(a.matvec([3]),[6])
        self.assertEqual(CSRMatrix.from_dense(Matrix(2,2,[0,2,0,0])).row_offsets,[0,1,1])
    def test_frames_do_not_alias_solution_or_each_other(self):
        a=CSRMatrix(2,2,[0,2,4],[0,1,0,1],[4,1,1,3])
        r=a.conjugate_gradient([6,7],capture=True)
        self.assertTrue(r.converged)
        r.x[0]=99
        self.assertAlmostEqual(r.iterates[-1][0],1)
        r.iterates[0][0]=99
        self.assertAlmostEqual(r.iterates[-1][0],1)
    def test_nonfinite_and_invalid_inputs(self):
        for v in [float('nan'),float('inf')]:
            with self.assertRaises(ValueError):CSRMatrix(1,1,[0,1],[0],[v])
        a=CSRMatrix(1,1,[0,1],[0],[1])
        for options in [dict(rtol=float('nan')),dict(atol=float('inf')),dict(max_iterations=-1)]:
            with self.assertRaises(ValueError):a.conjugate_gradient([1],**options)
        with self.assertRaises(ValueError):a.matvec([])
