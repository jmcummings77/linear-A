"""Independent math and routing checks for the algorithm comparison suite."""
import json
from pathlib import Path
import sys
from unittest import TestCase
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'benchmarks'))
from reference import generated_spd, benchmark_determinant, benchmark_checksum, determinant
from determinants import algorithm_cases
from run import measure


class DeterminantComparisonTests(TestCase):
    def test_spd_generator_is_exact_symmetric_and_strictly_dominant(self):
        for n in (1, 2, 3, 5, 16, 48):
            for seed in (0, 17, 2147483646):
                a = generated_spd(n, seed)['values']
                for i in range(n):
                    self.assertGreater(a[i*n+i], sum(abs(a[i*n+j]) for j in range(n) if j != i))
                    for j in range(n):
                        self.assertEqual(a[i*n+j], a[j*n+i])
                        self.assertEqual(a[i*n+j]*32, int(a[i*n+j]*32))

    def test_exact_spd_reference_agrees_with_permutation_oracle(self):
        for n in range(6):
            a = generated_spd(n, 17)
            self.assertEqual(benchmark_determinant(a, 32), determinant(a))
        self.assertEqual(benchmark_checksum('determinant_cholesky', 2, 17), 57303/1024)
        for n in (1, 4, 16):
            expected = benchmark_checksum('determinant_spd_lu', n, 17)
            self.assertEqual(expected, benchmark_checksum('determinant_cholesky', n, 17))
            self.assertEqual(expected, benchmark_checksum('determinant_spd', n, 17))

    def test_general_algorithms_receive_the_same_checksum(self):
        for op in ('determinant_lu', 'determinant_cofactor', 'determinant_small'):
            self.assertEqual(benchmark_checksum(op, 4, 17), benchmark_checksum('determinant', 4, 17))

    def test_cholesky_checks_include_invalid_domains(self):
        cases = algorithm_cases('cholesky')
        self.assertEqual({case['name'] for case in cases if case['invalid']},
                         {'nonsymmetric', 'indefinite', 'semidefinite', 'nonsquare'})
        self.assertTrue(all(case['op'] == 'determinant_cholesky' for case in cases))

    def test_timing_routes_algorithm_without_changing_expected_checksum(self):
        implementation = {'runner':['runner'], 'env':{}, 'operation_map':{'determinant_spd':'determinant_cholesky'}}
        with patch('run.require', return_value=json.dumps({'elapsed_ns':100, 'iterations':2, 'checksum':16})) as command:
            result = measure(implementation, 'determinant_spd', 2, 2, 17, 8)
        self.assertEqual(command.call_args[0][0], ['runner','bench','determinant_cholesky','2','2','17'])
        self.assertEqual(result['ns_per_op'], 50)
