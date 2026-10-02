import json
from pathlib import Path
import subprocess
import sys
import unittest
from matrix import Matrix


class MatrixTests(unittest.TestCase):
    def test_arithmetic_and_storage(self):
        values = [1, 2, 3, 4, 5, 6]
        a = Matrix(2, 3, values)
        values[0] = 99
        self.assertEqual(a.transpose().values, [1, 4, 2, 5, 3, 6])
        self.assertEqual((a @ Matrix(3, 2, [7, 8, 9, 10, 11, 12])).values, [58, 64, 139, 154])
        self.assertEqual(a.scale(0.5).values, [0.5, 1, 1.5, 2, 2.5, 3])
        self.assertEqual((a + a - a).values, a.values)
        self.assertEqual((Matrix.identity(2) @ a).values, a.values)
        copy = a.copy()
        copy[0, 0] = 9
        row, col, exported = a.row(0), a.column(1), a.values
        row[0] = col[0] = exported[0] = 9
        self.assertEqual(a.values, [1, 2, 3, 4, 5, 6])

    def test_known_determinants_and_structure(self):
        cases = [(1, [-7], -7), (2, [1, 2, 3, 4], -2), (3, [6, 1, 1, 4, -2, 5, 2, 8, 7], -306),
                 (4, [3, 2, 0, 1, 4, 0, 1, 2, 3, 0, 2, 1, 9, 2, 3, 1], 24),
                 (3, [0, 1, 0, 0, 0, 1, 1, 0, 0], 1), (2, [1, 2, 2, 4], 0)]
        for size, values, expected in cases:
            a = Matrix(size, size, values)
            self.assertAlmostEqual(a.determinant(), expected)
            self.assertEqual(a.values, values)
        a = Matrix(3, 3, [0.5, 7, -1, 0, -1.5, 4, 0, 0, 2.5])
        self.assertEqual(a.triangular(), (True, False))
        self.assertEqual(a.transpose().triangular(), (False, True))
        self.assertEqual(a.trace(), 1.5)
        self.assertAlmostEqual(a.determinant(), -1.875)
        self.assertEqual(Matrix(2, 3).triangular(), (False, False))

    def test_empty_shapes(self):
        self.assertEqual(Matrix().determinant(), 1)
        self.assertEqual(Matrix().trace(), 0)
        self.assertEqual(Matrix().triangular(), (True, True))
        self.assertEqual((Matrix(2, 0) @ Matrix(0, 3)).values, [0] * 6)
        self.assertEqual(Matrix(2, 0).row(1), [])

    def test_errors_leave_inputs_unchanged(self):
        a = Matrix(1, 1, [1e308])
        with self.assertRaises(ValueError):
            a.scale(2)
        self.assertEqual(a[0, 0], 1e308)
        for bad in [float("nan"), float("inf")]:
            with self.assertRaises(ValueError):
                a[0, 0] = bad
        with self.assertRaises(IndexError):
            _ = a[-1, 0]
        with self.assertRaises(ValueError):
            Matrix(-1, 2)
        with self.assertRaises(ValueError):
            Matrix(2, 2, [1])
        for action in [lambda: a + Matrix(2, 2), lambda: a @ Matrix(2, 1), lambda: Matrix(2, 3).trace(), lambda: Matrix(2, 3).determinant()]:
            with self.assertRaises(ValueError):
                action()

    def test_runner_and_deterministic_benchmark(self):
        runner = str(Path(__file__).with_name("runner.py"))
        result = subprocess.run([sys.executable, runner, "check", "multiply", "1", "2", "2", "1"], input="2 3 4 5", text=True, capture_output=True, check=True)
        self.assertEqual(json.loads(result.stdout), {"rows": 1, "cols": 1, "values": [23]})
        result = subprocess.run([sys.executable, runner, "bench", "scale", "2", "7", "2147483646"], text=True, capture_output=True, check=True)
        payload = json.loads(result.stdout)
        values = [((i * 17 + 2147483646 * 13) % 101 - 50) / 16 for i in range(4)]
        self.assertEqual(payload["checksum"], (values[0] + values[2] + values[3]) * 1.25 * 7)
        self.assertGreater(payload["elapsed_ns"], 0)
        result = subprocess.run([sys.executable, runner, "check", "trace", "1", "1"], input="1 2", text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")


if __name__ == "__main__":
    unittest.main()
