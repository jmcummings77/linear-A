import unittest,sys,random,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT/'benchmarks'))
from ic0_reference import reference,factor,apply,fixtures,BREAKDOWN
from sparse_reference import dense_csr,multiply
from sparse import CSRMatrix
from cholesky import IC0,SparseCholeskySymbolic

def matrix(a):return CSRMatrix(a['rows'],a['cols'],a['offsets'],a['indices'],a['values'])
class IC0Tests(unittest.TestCase):
 def test_random_rational_oracle_and_masked_reconstruction(self):
  rng=random.Random(2026)
  for n in range(1,13):
   a=[[0.]*n for _ in range(n)]
   for i in range(n):
    for j in range(i):
     if rng.random()<.3:a[i][j]=a[j][i]=rng.choice([-2.,-1.,1.,2.])
   for i in range(n):a[i][i]=1+sum(abs(x) for x in a[i])
   d=dense_csr(a);f=IC0(matrix(d));l=f.lower;rp,ci,v=factor(d)
   self.assertEqual(l.row_offsets,rp);self.assertEqual(l.column_indices,ci)
   for x,y in zip(l.values,v):self.assertAlmostEqual(x,y,places=12)
   lower=[[0.]*n for _ in range(n)]
   for i in range(n):
    for p in range(rp[i],rp[i+1]):lower[i][ci[p]]=l.values[p]
   for i in range(n):
    for p in range(rp[i],rp[i+1]):
     j=ci[p];self.assertAlmostEqual(sum(lower[i][k]*lower[j][k] for k in range(n)),a[i][j],places=11)
   for b in ([1.]*n,[float(i-3) for i in range(n)]):
    for x,y in zip(f.apply(b),apply(d,b)):self.assertAlmostEqual(x,y,places=11)
 def test_spd_breakdown_is_not_indefiniteness(self):
  # B=A-1.5 I has B²=2I, hence A's minimum eigenvalue is 1.5-sqrt(2)>0.
  b=[[BREAKDOWN[i][j]-(1.5 if i==j else 0) for j in range(4)] for i in range(4)]
  for i in range(4):
   for j in range(4):self.assertEqual(sum(b[i][k]*b[k][j] for k in range(4)),2*(i==j))
  a=matrix(dense_csr(BREAKDOWN));SparseCholeskySymbolic(a).factorize(a)
  with self.assertRaises(ValueError):reference(dense_csr(BREAKDOWN))
  with self.assertRaises(ArithmeticError):IC0(a)
 def test_ownership_conflicts_and_true_residuals(self):
  a=matrix(dense_csr([[4,1],[1,3]]));f=IC0(a);f.lower.values[0]=99;a.values[0]=99
  for rhs,truth in [([6,7],[1,2]),([3,-2],[1,-1])]:
   r=a.conjugate_gradient(rhs,preconditioner=f,capture=True)
   self.assertTrue(r.converged);self.assertEqual(r.iterations,1)
   for x,y in zip(r.x,truth):self.assertAlmostEqual(x,y)
  with self.assertRaises(ValueError):a.conjugate_gradient([0,0],jacobi=True,preconditioner=f)
  with self.assertRaises(ValueError):a.conjugate_gradient([0,0],preconditioner=IC0(matrix(dense_csr([[1]]))))
