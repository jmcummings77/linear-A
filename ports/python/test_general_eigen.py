import math
import random
import sys
import unittest
from pathlib import Path
from matrix import Matrix
from runner import result_json, benchmark
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'benchmarks'))
from general_eigen_reference import general_eigen_fixtures, assert_general_eigen_result, general_eigen_checksum


class GeneralEigenTests(unittest.TestCase):
    def test_independent_spectra_and_complex_residuals(self):
        for case in general_eigen_fixtures():
            with self.subTest(case=case['name']):
                a=case['a']; m=Matrix(a['rows'],a['cols'],a['values'])
                if case['invalid']:
                    with self.assertRaises(ValueError): m.eigen_general()
                else:
                    result=m.eigen_general()
                    assert_general_eigen_result(result_json(result),a,case['expected_complex_eigenvalues'],
                                                spectrum_scale=case.get('spectrum_scale'),componentwise=case.get('componentwise',False))
                    self.assertEqual(m.values,a['values'])

    def test_random_dense_inputs_and_storage(self):
        rng=random.Random(413)
        for n in range(1,10):
            for _ in range(8):
                data=[rng.uniform(-10,10) for _ in range(n*n)]
                a=Matrix(n,n,data); result=a.eigen_general()
                assert_general_eigen_result(result_json(result),dict(rows=n,cols=n,values=data))
                roots=[complex(r,i) for r,i in zip(result.values_real,result.values_imag)]
                self.assertLess(abs(sum(roots)-a.trace()),1e-9)
                trace_square=sum(data[i*n+j]*data[j*n+i] for i in range(n) for j in range(n))
                self.assertLess(abs(sum(x*x for x in roots)-trace_square),1e-8)
                determinant=a.determinant()
                self.assertLess(abs(math.prod(roots)-determinant),1e-8*max(1,abs(determinant)))
                result.vectors_real[0,0]=900
                self.assertEqual(a.values,data)
                self.assertLessEqual(abs(a.eigen_general().vectors_real[0,0]),1)

    def test_options_and_nonconvergence(self):
        a=Matrix(4,4,[4,1,2,3,0,5,1,2,2,1,6,1,3,0,1,7])
        for bound in (0,-1,1.5,True,100001,float('inf')):
            with self.assertRaises(ValueError): a.eigen_general(bound)
        with self.assertRaisesRegex(ValueError,'converge'): a.eigen_general(1)
        for invalid in (math.nan,math.inf,-math.inf):
            with self.assertRaises(ValueError): Matrix(1,1,[invalid]).eigen_general()
        with self.assertRaises(ValueError): Matrix(2,2,[1e308]*4).eigen_general()

    def test_runner_consumes_complete_eigenpairs(self):
        for n in (1,2,3,8):
            result=benchmark(['eigen_general',str(n),'2','17'])
            self.assertAlmostEqual(result['checksum'],2*general_eigen_checksum(n,17),places=10)
