"""Independent complex eigenpair validation must catch false spectra and vectors."""
import copy
import math
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'benchmarks'))
from general_eigen_reference import *


def rotation_pairs():
    s=math.sqrt(.5)
    return dict(eigenvalues_real=[0,0],eigenvalues_imag=[-1,1],
                eigenvectors_real=matrix(2,2,[s,s,0,0]),eigenvectors_imag=matrix(2,2,[0,0,s,-s]))


class GeneralEigenReferenceTests(unittest.TestCase):
    def test_complex_rotation_and_arbitrary_phase(self):
        a=matrix(2,2,[0,-1,1,0]); r=rotation_pairs()
        assert_general_eigen_result(r,a,[[0,-1],[0,1]])
        r['eigenvectors_real'],r['eigenvectors_imag']=r['eigenvectors_imag'],r['eigenvectors_real']
        r['eigenvectors_imag']['values']=[-x for x in r['eigenvectors_imag']['values']]
        assert_general_eigen_result(r,a,[[0,-1],[0,1]])

    def test_rejects_wrong_spectrum_and_zero_or_wrong_vectors(self):
        a=matrix(2,2,[0,-1,1,0])
        for key,value in [('eigenvalues_imag',[0,0]),('eigenvalues_real',[1,1]),
                           ('eigenvectors_real',matrix(2,2,[0]*4)),
                           ('eigenvectors_imag',matrix(2,2,[0]*4)),
                           ('eigenvectors_imag',matrix(1,4,[0]*4)),('eigenvalues_imag',[True,1])]:
            bad=rotation_pairs();bad[key]=value
            with self.assertRaises(AssertionError): assert_general_eigen_result(bad,a,[[0,-1],[0,1]])

    def test_defective_columns_are_allowed_but_multiplicity_is_checked(self):
        a=matrix(2,2,[2,1,0,2])
        r=dict(eigenvalues_real=[2,2],eigenvalues_imag=[0,0],
               eigenvectors_real=matrix(2,2,[1,1,0,0]),eigenvectors_imag=matrix(2,2,[0]*4))
        assert_general_eigen_result(r,a,[[2,0],[2,0]])
        with self.assertRaises(AssertionError): assert_general_eigen_result(r,a,[[2,0],[3,0]])

    def test_balanced_spectrum_cannot_hide_behind_huge_input_norm(self):
        a=matrix(2,2,[0,1e300,-1e-300,0])
        r=dict(eigenvalues_real=[0,0],eigenvalues_imag=[0,0],
               eigenvectors_real=matrix(2,2,[1,1,0,0]),eigenvectors_imag=matrix(2,2,[0]*4))
        with self.assertRaisesRegex(AssertionError,'spectrum'):
            assert_general_eigen_result(r,a,[[0,-1],[0,1]],spectrum_scale=1)

    def test_componentwise_case_rejects_lost_tiny_vector_component(self):
        a=matrix(2,2,[0,1e300,-1e-300,0])
        r=dict(eigenvalues_real=[0,0],eigenvalues_imag=[-1,1],
               eigenvectors_real=matrix(2,2,[1,1,0,0]),eigenvectors_imag=matrix(2,2,[0]*4))
        with self.assertRaisesRegex(AssertionError,'Componentwise'):
            assert_general_eigen_result(r,a,[[0,-1],[0,1]],spectrum_scale=1,componentwise=True)

    def test_fixture_contract_and_analytic_workload(self):
        cases=general_eigen_fixtures()
        self.assertEqual(len(cases),20)
        self.assertEqual(len({c['name'] for c in cases}),20)
        for n in (1,2,3,8,16):
            a=generated_general_eigen(n,17); spectrum=general_eigen_spectrum(n,17)
            self.assertAlmostEqual(sum(re for re,im in spectrum),sum(a['values'][i*n+i] for i in range(n)))
            self.assertAlmostEqual(general_eigen_checksum(n,17),sum((i+1)*(re+abs(im)) for i,(re,im) in enumerate(spectrum))+n)
