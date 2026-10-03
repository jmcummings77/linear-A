import unittest
from sparse import CSRMatrix
from cholesky import SparseCholeskySymbolic

class CholeskyOwnershipTests(unittest.TestCase):
    def test_reuse_and_pattern_contract(self):
        a=CSRMatrix(2,2,[0,2,4],[0,1,0,1],[4,2,2,3])
        s=SparseCholeskySymbolic(a);f=s.factorize(a)
        s.row_offsets[1]=99;s.column_indices[0]=99;s.fill_steps[0]=99
        scaled=CSRMatrix(2,2,[0,2,4],[0,1,0,1],[8,4,4,6])
        g=s.factorize(scaled)
        for actual in (f.solve([8,8]),g.solve([16,16])):
            self.assertAlmostEqual(actual[0],1);self.assertAlmostEqual(actual[1],2)
        with self.assertRaises(ValueError):s.factorize(CSRMatrix(2,2,[0,1,2],[0,1],[4,3]))
        with self.assertRaises(ValueError):f.solve([1])
        with self.assertRaises(ValueError):f.solve([float('inf'),1])
        with self.assertRaises(ArithmeticError):f.solve([1.7976931348623157e308,-1.7976931348623157e308])
        del a,s
        self.assertEqual(f.solve([0,0]),[0,0])
        lower=f.lower;lower.values[0]=99
        self.assertAlmostEqual(f.solve([8,8])[0],1)
