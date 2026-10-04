import math
from pathlib import Path
import sys
import unittest
from multigrid import GeometricMultigrid
sys.path.append(str(Path(__file__).resolve().parents[2]/'benchmarks'))
from multigrid_reference import reference,transfer,truth

class MultigridTests(unittest.TestCase):
    def test_dense_oracle_and_owned_results(self):
        for w in (1,3,7):
            m=GeometricMultigrid(w);b=[float(i%5-2) for i in range(m.size)];saved=b[:]
            expected=list(map(float,reference(w,b)))
            for _ in range(2):
                x=m.apply(b)
                self.assertLess(max(abs(a-v) for a,v in zip(x,expected)),1e-13)
                x[0]=999
            self.assertEqual(b,saved)
            a=m.matrix;a.values[0]=999
            self.assertEqual(m.matrix.values[0],4)

    def test_preconditioner_is_symmetric_positive_definite(self):
        # Form the complete linear map from basis applications; independent LDL^T.
        m=GeometricMultigrid(3);n=m.size
        columns=[m.apply([float(i==j) for i in range(n)]) for j in range(n)]
        a=list(map(list,zip(*columns)));l=[[float(i==j) for j in range(n)] for i in range(n)];d=[]
        for i in range(n):
            for j in range(i):
                self.assertAlmostEqual(a[i][j],a[j][i],places=14)
                l[i][j]=(a[i][j]-sum(l[i][k]*d[k]*l[j][k] for k in range(j)))/d[j]
            d.append(a[i][i]-sum(l[i][k]**2*d[k] for k in range(i)))
            self.assertGreater(d[-1],0)
        b=list(range(n));x=m.apply(b)
        for i in range(n):self.assertAlmostEqual(x[i],sum(a[i][j]*b[j] for j in range(n)),places=13)

    def test_trace_is_actual_level_residual(self):
        m=GeometricMultigrid(7);b=[float(i%7-2) for i in range(m.size)];t=m.trace(b)
        self.assertEqual(t['x'],m.apply(b))
        self.assertEqual([(f['width'],f['phase']) for f in t['frames']],[(7,'enter'),(7,'pre_smooth'),(3,'enter'),(3,'pre_smooth'),(1,'enter'),(1,'coarse_solve'),(3,'correct'),(3,'post_smooth'),(7,'correct'),(7,'post_smooth')])
        for f in t['frames']:
            ax=GeometricMultigrid(f['width']).matrix.matvec(f['x'])
            for r,v,q in zip(f['residual'],f['b'],ax):self.assertAlmostEqual(r,v-q,places=13)
        t['frames'][0]['x'][0]=9
        self.assertEqual(t['frames'][1]['x'][0],m.trace(b)['frames'][1]['x'][0])

    def test_scaling_and_invalid_inputs(self):
        steps=[]
        for w in (7,15,31):
            m=GeometricMultigrid(w);a=m.matrix;x=truth(w);b=a.matvec(x)
            r=a.conjugate_gradient(b,preconditioner=m)
            self.assertTrue(r.converged);steps.append(r.iterations)
            self.assertLess(math.dist(a.matvec(r.x),b),math.hypot(*b)*1e-10)
        self.assertLessEqual(max(steps)-min(steps),5)
        for w in (0,2,4,256,-1,3.0,True):
            with self.assertRaises(ValueError):GeometricMultigrid(w)
        m=GeometricMultigrid(3)
        for b in ([1], [float('nan')]*9, [float('inf')]*9):
            with self.assertRaises(ValueError):m.apply(b)
        with self.assertRaises(ArithmeticError):m.apply([1.7e308]*9)
        with self.assertRaises(ValueError):m.matrix.conjugate_gradient([0]*9,jacobi=True,preconditioner=m)

if __name__=='__main__':unittest.main()
