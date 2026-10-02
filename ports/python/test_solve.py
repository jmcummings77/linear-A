"""Independent exact references and numerical solver regressions."""
import math
import random
import unittest
from fractions import Fraction
from matrix import Matrix


def exact_solve(a, b):
    n, p = len(a), len(b[0]) if b else 0
    work = [[Fraction(v) for v in ar + br] for ar, br in zip(a, b)]
    for k in range(n):
        pivot = next(i for i in range(k, n) if work[i][k])
        work[k], work[pivot] = work[pivot], work[k]
        divisor = work[k][k]
        work[k] = [v / divisor for v in work[k]]
        for i in range(n):
            if i != k:
                coefficient = work[i][k]
                work[i] = [v - coefficient * w for v, w in zip(work[i], work[k])]
    return [float(v) for row in work for v in row[n:n+p]]


def matrix(rows):
    return Matrix(len(rows), len(rows[0]) if rows else 0, [v for row in rows for v in row])


class SolverTests(unittest.TestCase):
    def close(self, actual, expected, tolerance=2e-12):
        self.assertEqual(len(actual.values), len(expected))
        for a, e in zip(actual.values, expected):
            self.assertTrue(math.isclose(a, e, abs_tol=tolerance, rel_tol=tolerance), (a, e))

    def test_exact_generated_systems_and_reuse(self):
        rng = random.Random(741)
        for n in range(1, 7):
            for _ in range(10):
                a = [[rng.randrange(-5, 6) + (20 if i == j else 0) for j in range(n)] for i in range(n)]
                b = [[rng.randrange(-10, 11) for _ in range(3)] for _ in range(n)]
                expected = exact_solve(a, b)
                source, rhs = matrix(a), matrix(b)
                for factor in (source.factor_lu(), source.factor_qr()):
                    self.close(factor.solve(rhs), expected)
                    self.close(factor.solve(rhs), expected)
                self.assertEqual(source.values, [v for row in a for v in row])
                self.assertEqual(rhs.values, [v for row in b for v in row])

    def test_pivoting_and_snapshot(self):
        a, b = matrix([[0, 2], [1, 3]]), matrix([[4], [7]])
        factor = a.factor_lu()
        a[0, 0] = 100
        self.close(factor.solve(b), [1, 2])
        self.close(factor.solve(matrix([[2], [4]])), [1, 1])

    def test_spd_and_uniform_extreme_scaling(self):
        for scale in (1, 1e-300, 1e300):
            a = matrix([[4*scale, scale], [scale, 3*scale]])
            b = matrix([[6*scale], [7*scale]])
            for factor in (a.factor_lu(), a.factor_cholesky(), a.factor_qr()):
                self.close(factor.solve(b), [1, 2])
                self.assertAlmostEqual(factor.reciprocal_condition(), 11/25)

    def test_least_squares_exact_normal_equation_reference(self):
        # Normal equations are used only as an exact rational test oracle.
        a = [[1, 10], [1, 20], [1, 30], [1, 40]]
        b = [[2, 1], [4, 5], [5, 2], [8, 9]]
        ata = [[sum(row[i]*row[j] for row in a) for j in range(2)] for i in range(2)]
        atb = [[sum(a[k][i]*b[k][j] for k in range(4)) for j in range(2)] for i in range(2)]
        self.close(matrix(a).least_squares(matrix(b)), exact_solve(ata, atb))

    def test_shapes_empty_and_condition(self):
        for factor in (Matrix(0, 0, []).factor_lu(), Matrix(0, 0, []).factor_cholesky(), Matrix(0, 0, []).factor_qr()):
            self.assertEqual(factor.solve(Matrix(0, 3, [])).values, [])
            self.assertEqual(factor.reciprocal_condition(), 1)
        result = Matrix(3, 0, []).least_squares(Matrix(3, 2, [1]*6))
        self.assertEqual((result.rows, result.cols), (0, 2))
        self.assertEqual(matrix([[1, 0], [0, 4]]).factor_lu().reciprocal_condition(), 0.25)
        self.assertEqual(matrix([[1, 0], [0, 4]]).solve(Matrix(2, 0, [])).values, [])
        with self.assertRaises(ValueError):
            matrix([[1], [2]]).factor_qr().reciprocal_condition()
        with self.assertRaises(ValueError):
            matrix([[1]]).solve(matrix([[1], [2]]))

    def test_rejections(self):
        for values, method in [([[1, 2], [2, 4]], 'factor_lu'), ([[1, 2], [2, 4]], 'factor_qr'),
                               ([[1, 2], [0, 1]], 'factor_cholesky'), ([[1, 2], [2, 1]], 'factor_cholesky'),
                               ([[1, 2]], 'factor_qr'), ([[1], [2]], 'factor_lu'),
                               ([[1, 0], [0, 1e-18]], 'factor_qr'),
                               ([[1e300, 0], [0, 1e-300]], 'factor_lu')]:
            with self.subTest(method=method, values=values), self.assertRaises(ValueError):
                getattr(matrix(values), method)()
        with self.assertRaises(ValueError):
            matrix([[1e-300]]).solve(matrix([[1e300]]))
        self.close(matrix([[1, 0], [0, 1e-18]]).solve(matrix([[1], [1e-18]])), [1, 1])

if __name__ == '__main__':
    unittest.main()
