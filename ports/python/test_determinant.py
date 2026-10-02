from fractions import Fraction
from itertools import permutations
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import unittest

from matrix import Matrix


def exact_determinant(values, size):
    total = Fraction(0)
    for order in permutations(range(size)):
        term = Fraction((-1) ** sum(order[i] > order[j] for i in range(size) for j in range(i + 1, size)))
        for row, col in enumerate(order):
            term *= Fraction(values[row * size + col])
        total += term
    return float(total)


class DeterminantTests(unittest.TestCase):
    def test_auto_and_lu_against_independent_exact_permutation_formula(self):
        random_source = random.Random(741)
        for size in range(7):
            for _ in range(5):
                values = [random_source.randrange(-8, 9) / 4 for _ in range(size * size)]
                matrix = Matrix(size, size, values)
                expected = exact_determinant(values, size)
                for algorithm in ("auto", "lu"):
                    with self.subTest(size=size, algorithm=algorithm):
                        self.assertTrue(math.isclose(matrix.determinant(algorithm), expected, rel_tol=1e-11, abs_tol=1e-9))
                        self.assertEqual(matrix.values, values)

    def test_mixed_exponents_subnormal_products_and_underflowing_elimination_factor(self):
        cases = [
            (1, [math.ulp(0.0)], math.ulp(0.0)),
            (2, [1e200] * 4, 0.0),
            (2, [1e-200, 1e-200, 3e-124, 6e-124], math.ulp(0.0)),
            (2, [1e308, 1e308, 1e-308, 2e-308], 1.0),
            (2, [1e308, 1e-308, 1e308, 2e-308], 1.0),
            (3, [1e308, 1e308, 1e-308, 1e-308, 2e-308, 1e308, 0, 0, 1e-308], 1e-308),
        ]
        for size, values, expected in cases:
            matrix = Matrix(size, size, values)
            for algorithm in ("auto", "lu"):
                with self.subTest(size=size, algorithm=algorithm, expected=expected):
                    actual = matrix.determinant(algorithm)
                    if expected in (0.0, math.ulp(0.0)):
                        self.assertEqual(actual, expected)
                    else:
                        self.assertTrue(math.isclose(actual, expected, rel_tol=1e-12, abs_tol=0.0))
                    self.assertEqual(matrix.values, values)
        for diagonal in ([1e-200, 1e-200, 1e200, 1e200], [1e200, 1e200, 1e-200, 1e-200]):
            matrix = Matrix(4, 4, [diagonal[row] if row == col else 0 for row in range(4) for col in range(4)])
            for algorithm in ("auto", "lu", "cholesky"):
                self.assertTrue(math.isclose(matrix.determinant(algorithm), 1.0, rel_tol=1e-12))

    def test_cholesky_spd_contract_and_nonmutation(self):
        # L=[[2,0,0],[1,3,0],[-1,2,4]] gives A=L*L^T and det(A)=24^2.
        values = [4, 2, -2, 2, 10, 5, -2, 5, 21]
        matrix = Matrix(3, 3, values)
        self.assertTrue(math.isclose(matrix.determinant("cholesky"), 576.0, rel_tol=1e-12))
        self.assertEqual(matrix.values, values)
        self.assertEqual(Matrix().determinant("cholesky"), 1.0)
        self.assertEqual(Matrix(1, 1, [math.ulp(0.0)]).determinant("cholesky"), math.ulp(0.0))
        for values in ([2, 100, 1, 2], [1, 2, 2, 1], [1, 1, 1, 1], [-1, 0, 0, -1]):
            matrix = Matrix(2, 2, values)
            with self.assertRaises(ValueError):
                matrix.determinant("cholesky")
            self.assertEqual(matrix.values, values)
        for algorithm in ("auto", "lu", "cholesky"):
            with self.assertRaises(ValueError):
                Matrix(2, 3).determinant(algorithm)
            with self.assertRaises(ValueError):
                Matrix(2, 2, [1e308, 0, 0, 1e308]).determinant(algorithm)
            self.assertEqual(Matrix(2, 2, [1e-200, 0, 0, 1e-200]).determinant(algorithm), 0.0)
        with self.assertRaises(ValueError):
            Matrix().determinant("cofactor")

    def test_algorithm_runner_names_and_shared_spd_benchmark_inputs(self):
        runner = str(Path(__file__).with_name("runner.py"))
        for operation in ("determinant_lu", "determinant_cholesky", "determinant_spd_lu"):
            result = subprocess.run([sys.executable, runner, "check", operation, "2", "2"],
                                    input="4 2 2 3", text=True, capture_output=True, check=True)
            self.assertTrue(math.isclose(json.loads(result.stdout)["value"], 8.0, rel_tol=1e-12))
        for operation in ("determinant_cholesky", "determinant_spd_lu"):
            result = subprocess.run([sys.executable, runner, "bench", operation, "3", "2", "17"],
                                    text=True, capture_output=True, check=True)
            data = json.loads(result.stdout)
            self.assertEqual(data["iterations"], 2)
            self.assertGreater(data["elapsed_ns"], 0)
            self.assertTrue(math.isclose(data["checksum"], 2 * 28161419 / 16384, rel_tol=1e-12))
        general = [161/16, -14/16, 3/16, 20/16, 229/16, -47/16, -30/16, -13/16, 196/16]
        for operation in ("determinant", "determinant_lu"):
            result = subprocess.run([sys.executable, runner, "bench", operation, "3", "2", "17"],
                                    text=True, capture_output=True, check=True)
            self.assertTrue(math.isclose(json.loads(result.stdout)["checksum"], 2 * exact_determinant(general, 3), rel_tol=1e-12))
        result = subprocess.run([sys.executable, runner, "check", "determinant_cholesky", "2", "2"],
                                input="1 2 2 1", text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        result = subprocess.run([sys.executable, runner, "bench", "determinant_cofactor", "2", "1", "17"],
                                text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")


if __name__ == "__main__":
    unittest.main()
