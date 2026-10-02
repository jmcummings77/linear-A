"""Verify analytic spectra and reject malformed or numerically false eigensystems."""
import copy
import math
from pathlib import Path
import sys
from unittest import TestCase

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'benchmarks'))
from eigen_reference import assert_eigen_result, eigen_checksum, eigen_fixtures, eigen_spectrum, generated_eigen
from reference import matrix


def sine_eigensystem(n, seed):
    factor = math.sqrt(2 / (n+1))
    return {"eigenvalues": eigen_spectrum(n, seed), "eigenvectors": matrix(n, n,
        [factor*math.sin((i+1)*(j+1)*math.pi/(n+1)) for i in range(n) for j in range(n)])}


class EigenReferenceTests(TestCase):
    def test_toeplitz_generator_has_analytic_sine_eigenbasis(self):
        for n in (0, 1, 2, 3, 8, 24, 48):
            for seed in (0, 9, 17, 2147483646):
                a = generated_eigen(n, seed)
                result = sine_eigensystem(n, seed)
                assert_eigen_result(result, a, eigen_spectrum(n, seed))
                self.assertAlmostEqual(sum(result['eigenvalues']), sum(a['values'][i*n+i] for i in range(n)), places=10)
                self.assertAlmostEqual(eigen_checksum(n, seed),
                    math.fsum((i+1)*value for i,value in enumerate(result['eigenvalues'])) +
                    math.fsum(value*value for value in result['eigenvectors']['values']), places=10)

    def test_fixture_shapes_symmetry_and_invalid_domain(self):
        cases = eigen_fixtures()
        self.assertGreaterEqual(len(cases), 20)
        self.assertEqual(len({case['name'] for case in cases}), len(cases))
        for case in cases:
            self.assertEqual(case['op'], 'eigen_symmetric')
            a, n = case['a'], case['a']['rows']
            self.assertEqual(len(a['values']), n*a['cols'])
            if not case['invalid']:
                self.assertEqual(a['cols'], n)
                self.assertEqual(len(case['expected_eigenvalues']), n)
                for i in range(n):
                    for j in range(n):
                        self.assertEqual(a['values'][i*n+j], a['values'][j*n+i])
        self.assertEqual({case['name'] for case in cases if case['invalid']},
                         {'nonsquare eigensystem', 'nonsymmetric eigensystem', 'tiny asymmetry is rejected'})

    def test_signs_and_repeated_eigenspace_basis_are_not_fixed(self):
        a = matrix(3, 3, [2,0,0,0,2,0,0,0,2])
        rotated = {'eigenvalues':[2,2,2], 'eigenvectors':matrix(3,3,[.6,-.8,0,.8,.6,0,0,0,-1])}
        assert_eigen_result(rotated, a, [2,2,2])
        for scale in (1e-200, 1e200):
            scaled = matrix(3,3,[value*scale for value in a['values']])
            result = copy.deepcopy(rotated)
            result['eigenvalues'] = [2*scale]*3
            assert_eigen_result(result, scaled, [2*scale]*3)

    def test_rejects_correct_spectrum_with_wrong_eigenvectors(self):
        a = generated_eigen(3, 0)
        result = sine_eigensystem(3, 0)
        result['eigenvectors'] = matrix(3,3,[1,0,0,0,1,0,0,0,1])
        with self.assertRaisesRegex(AssertionError, 'residual'):
            assert_eigen_result(result, a, eigen_spectrum(3, 0))

    def test_rejects_dependent_vectors_even_when_residual_is_zero(self):
        result = {'eigenvalues':[2,2], 'eigenvectors':matrix(2,2,[1,1,0,0])}
        with self.assertRaisesRegex(AssertionError, 'orthonormal'):
            assert_eigen_result(result, matrix(2,2,[2,0,0,2]), [2,2])

    def test_tiny_inputs_do_not_get_an_absolute_error_allowance(self):
        a = matrix(2,2,[2e-200,1e-200,1e-200,2e-200])
        bad = {'eigenvalues':[0,0], 'eigenvectors':matrix(2,2,[1,0,0,1])}
        with self.assertRaisesRegex(AssertionError, 'residual'):
            assert_eigen_result(bad, a, [1e-200,3e-200])
        bad['eigenvalues'] = [0,5e-324]
        with self.assertRaisesRegex(AssertionError, 'residual'):
            assert_eigen_result(bad, matrix(2,2,[0,0,0,0]), [0,0])

    def test_valid_empty_and_subnormal_scalar(self):
        assert_eigen_result({'eigenvalues':[], 'eigenvectors':matrix(0,0,[])}, matrix(0,0,[]), [])
        assert_eigen_result({'eigenvalues':[5e-324], 'eigenvectors':matrix(1,1,[1])},
                            matrix(1,1,[5e-324]), [5e-324])

    def test_rejects_malformed_nonfinite_and_unsorted_results(self):
        a = generated_eigen(2, 0)
        good = sine_eigensystem(2, 0)
        bad_results = []
        for replacement in ([1], [1, True], [1, float('nan')], [1, float('inf')], [3,1]):
            bad = copy.deepcopy(good)
            bad['eigenvalues'] = replacement
            bad_results.append(bad)
        for replacement in ([1,0,0], [1,0,0,True], [1,0,0,float('nan')], [1e308,0,0,1e308]):
            bad = copy.deepcopy(good)
            bad['eigenvectors']['values'] = replacement
            bad_results.append(bad)
        bad = copy.deepcopy(good)
        bad['eigenvectors']['rows'] = True
        bad_results.append(bad)
        bad_results += [None, {}, {**good, 'extra':0}]
        for bad in bad_results:
            with self.subTest(result=bad), self.assertRaises(AssertionError):
                assert_eigen_result(bad, a, eigen_spectrum(2, 0))

    def test_invalid_generator_sizes(self):
        for bad in (-1, 1.5, True):
            with self.assertRaises(ValueError):
                generated_eigen(bad, 0)
            with self.assertRaises(ValueError):
                eigen_spectrum(bad, 0)
