import json
import math
from pathlib import Path
import subprocess
import sys
import unittest

from matrix import Matrix, SymmetricEigenDecomposition


class SymmetricEigenTests(unittest.TestCase):
    def assert_eigenpairs(self, a, expected):
        before = a.values
        result = a.eigen_symmetric()
        self.assertIsInstance(result, SymmetricEigenDecomposition)
        n, q = a.rows, result.vectors
        self.assertEqual((q.rows, q.cols), (n, n))
        self.assertEqual(result.values, sorted(result.values))
        scale = max([abs(value) for value in a.values] + [1e-300])
        for col, value in enumerate(expected):
            self.assertLessEqual(abs(result.values[col]/scale-value/scale), 1e-9)
            for row in range(n):
                transformed = math.fsum(a[row,k]/scale*q[k,col] for k in range(n))
                self.assertLessEqual(abs(transformed-result.values[col]/scale*q[row,col]), 1e-9)
                inner = math.fsum(q[k,row]*q[k,col] for k in range(n))
                self.assertLessEqual(abs(inner-int(row == col)), 1e-9)
        self.assertEqual(a.values, before)
        if n:
            result.values[0] = 99
            q[0,0] = 99
            self.assertEqual(a.values, before)

    def test_analytic_spectra_scaling_and_repeated_eigenspaces(self):
        for scale in (1, 1e-200, 1e200):
            self.assert_eigenpairs(Matrix(2,2,[2*scale,scale,scale,2*scale]), [scale,3*scale])
        self.assert_eigenpairs(Matrix(2,2,[0,1,1,0]), [-1,1])
        self.assert_eigenpairs(Matrix(2,2,[1,1,1,1]), [0,2])
        self.assert_eigenpairs(Matrix(3,3,[6,-2,-2,-2,3,4,-2,4,3]), [-1,(13-math.sqrt(33))/2,(13+math.sqrt(33))/2])
        self.assert_eigenpairs(Matrix(3,3,[0]*9), [0,0,0])
        self.assert_eigenpairs(Matrix.identity(4), [1,1,1,1])
        self.assert_eigenpairs(Matrix(), [])
        for n in (3,8,16):
            a = Matrix(n,n,[2 if i == j else -1 if abs(i-j) == 1 else 0 for i in range(n) for j in range(n)])
            self.assert_eigenpairs(a, [2-2*math.cos(k*math.pi/(n+1)) for k in range(1,n+1)])

    def test_diagonal_magnitudes_are_preserved_exactly(self):
        self.assertEqual(Matrix(2,2,[1e300,0,0,1e-300]).eigen_symmetric().values, [1e-300,1e300])
        self.assertEqual(Matrix(1,1,[5e-324]).eigen_symmetric().values, [5e-324])
        # Squaring the off-diagonal entries would underflow and falsely report convergence.
        result = Matrix(2,2,[1,1e-200,1e-200,2]).eigen_symmetric(tolerance=1e-250)
        self.assertNotEqual(result.vectors[1,0], 0)

    def test_invalid_domain_options_nonconvergence_and_range(self):
        for a in (Matrix(2,3), Matrix(2,2,[1,2,0,1]), Matrix(2,2,[1,5e-324,0,1])):
            with self.assertRaises(ValueError):
                a.eigen_symmetric()
        a = Matrix(3,3,[4,1,2,1,3,.5,2,.5,5])
        for tolerance in (0,-1,1,float('nan'),float('inf'),True):
            with self.assertRaises(ValueError):
                a.eigen_symmetric(tolerance=tolerance)
        for sweeps in (0,-1,1.5,True):
            with self.assertRaises(ValueError):
                a.eigen_symmetric(max_sweeps=sweeps)
        before = a.values
        with self.assertRaisesRegex(ValueError, 'converge'):
            a.eigen_symmetric(max_sweeps=1)
        self.assertEqual(a.values, before)
        with self.assertRaises(ValueError):
            Matrix(2,2,[1e308]*4).eigen_symmetric()

    def test_runner_eigenpair_shape_and_benchmark_consumption(self):
        runner = str(Path(__file__).with_name('runner.py'))
        result = subprocess.run([sys.executable,runner,'check','eigen_symmetric','2','2'],
                                input='2 1 1 2',text=True,capture_output=True,check=True)
        payload = json.loads(result.stdout)
        self.assertEqual(set(payload), {'eigenvalues','eigenvectors'})
        self.assertEqual(payload['eigenvalues'], [1,3])
        self.assertEqual(payload['eigenvectors']['rows'], 2)
        result = subprocess.run([sys.executable,runner,'bench','eigen_symmetric','3','2','17'],
                                text=True,capture_output=True,check=True)
        payload = json.loads(result.stdout)
        expected = 2*(sum(k*(2-2*math.cos(k*math.pi/4)) for k in range(1,4))+3)
        self.assertAlmostEqual(payload['checksum'], expected, places=9)
        self.assertGreater(payload['elapsed_ns'], 0)
