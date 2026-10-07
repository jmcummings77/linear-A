import cmath
import copy
import importlib.util
import json
import math
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'benchmarks'))
from graph_reference import (SHIFT, fixtures, graph_energy, graph_product,
                             manufactured_solution, petersen_graph, prism_graph,
                             prism_spectrum, shifted_laplacian)
from sparse_reference import diffusion, multiply

spec = importlib.util.spec_from_file_location('graph_sparse_harness', ROOT / 'benchmarks/sparse.py')
harness = importlib.util.module_from_spec(spec)
spec.loader.exec_module(harness)


class GraphReferenceTests(unittest.TestCase):
    def test_petersen_exact_spectral_certificate(self):
        n, edges = petersen_graph()
        neighbors = [set() for _ in range(n)]
        for u, v in edges:
            neighbors[u].add(v)
            neighbors[v].add(u)
        self.assertEqual(len(edges), 15)
        self.assertEqual([len(row) for row in neighbors], [3] * 10)
        # H^2 + H = 2I + J exactly: the constant eigenvalue is 3;
        # on its orthogonal complement eigenvalues are 1 or -2.
        # Trace(H)=0 determines multiplicities 5 and 4 respectively.
        for i in range(n):
            for j in range(n):
                self.assertEqual(len(neighbors[i] & neighbors[j]) + int(j in neighbors[i]),
                                 2 * int(i == j) + 1)
        spectrum = [3] + [1] * 5 + [-2] * 4
        self.assertEqual(sum(spectrum), 0)
        self.assertLess(max(abs(x) for x in spectrum[1:]), 2 * math.sqrt(2))
        # The explicitly present pentagon independently certifies nonbipartiteness.
        self.assertTrue(all(tuple(sorted((i, (i + 1) % 5))) in edges for i in range(5)))

    def test_prism_fourier_spectrum_and_family_limit(self):
        for m in (3, 4, 5, 17):
            n, edges = prism_graph(m)
            neighbors = [[] for _ in range(n)]
            for u, v in edges:
                neighbors[u].append(v)
                neighbors[v].append(u)
            self.assertEqual([len(row) for row in neighbors], [3] * n)
            for k in range(m):
                for sign in (-1, 1):
                    eigenvalue = 2 * math.cos(2 * math.pi * k / m) + sign
                    vector = [sign**layer * cmath.exp(2j * math.pi * k * i / m)
                              for layer in (0, 1) for i in range(m)]
                    # Unit-magnitude Fourier entries; this allows trig roundoff,
                    # not an eigensolver or a numerical rank decision.
                    residual = max(abs(sum(vector[j] for j in neighbors[i]) - eigenvalue * vector[i])
                                   for i in range(n))
                    self.assertLess(residual, 1e-12)
        expected = [-2, -2, 0, 0, 1, 3]
        for actual, wanted in zip(prism_spectrum(3), expected):
            self.assertAlmostEqual(actual, wanted)
        self.assertLess(max(abs(x) for x in prism_spectrum(3)[:-1]), 2 * math.sqrt(2))
        self.assertGreater(prism_spectrum(17)[-2], 2 * math.sqrt(2))

    def test_edge_energy_and_csr_action_agree_independently(self):
        for n, edges in (petersen_graph(), prism_graph(3), prism_graph(17)):
            a = shifted_laplacian(n, edges)
            self.assertEqual(len(a['values']), 4 * n)
            for i in range(n):
                row = a['indices'][a['offsets'][i]:a['offsets'][i + 1]]
                self.assertEqual(row, sorted(set(row)))
            for vector in ([1.] * n, manufactured_solution(n), [(-1.)**i for i in range(n)]):
                result = graph_product(n, edges, vector)
                self.assertEqual(multiply(a, vector), result)
                self.assertEqual(math.fsum(x*y for x, y in zip(vector, result)),
                                 graph_energy(n, edges, vector))
                self.assertGreater(graph_energy(n, edges, vector), 0)
            self.assertEqual(graph_product(n, edges, [1.] * n), [SHIFT] * n)

    def test_oracle_detects_wrong_edge_or_shift(self):
        n, edges = petersen_graph()
        truth = manufactured_solution(n)
        expected = graph_product(n, edges, truth)
        wrong = copy.deepcopy(shifted_laplacian(n, edges))
        wrong['values'][0] += 1
        self.assertNotEqual(multiply(wrong, truth), expected)
        wrong = shifted_laplacian(n, edges[:-1])
        self.assertNotEqual(multiply(wrong, truth), expected)

    def test_invalid_graphs_cannot_create_psd_or_malformed_fixtures(self):
        for size in (True, 2, 3.5):
            with self.assertRaises(ValueError):
                prism_graph(size)
        for shift in (0, -1, math.inf, math.nan):
            with self.assertRaises(ValueError):
                shifted_laplacian(2, [(0, 1)], shift)
        for edges in ([(0, 0)], [(1, 0)], [(0, 2)], [(0, 1), (0, 1)]):
            with self.assertRaises(ValueError):
                shifted_laplacian(2, edges)

    def test_shared_graph_cases_have_manufactured_answers(self):
        cases = fixtures()
        self.assertEqual(len(cases), 6)
        self.assertEqual({case['op'] for case in cases}, {'spmv', 'cg'})
        self.assertTrue(all(case in harness.fixtures() for case in cases))
        for case in cases:
            a = case['a']
            if case['op'] == 'cg':
                self.assertEqual(multiply(a, case['expected']), case['b'])
            else:
                self.assertEqual(multiply(a, case['b']), case['expected'])


class GraphBenchmarkTests(unittest.TestCase):
    def test_default_benchmark_inputs_preserve_diffusion(self):
        for operation, case, metadata in harness.benchmark_cases(4):
            self.assertEqual(case['a'], diffusion(4, 2))
            self.assertEqual(metadata, {'contrast': 2})
            if case['op'] == 'cg':
                self.assertEqual(case['expected'], [1.] * 16)

    def test_prism_benchmark_has_nongrid_size_and_nonconstant_rhs(self):
        for operation, case, metadata in harness.benchmark_cases(17, 'prism'):
            self.assertEqual(case['a']['rows'], 34)
            self.assertEqual(len(case['a']['values']), 136)
            self.assertEqual(metadata, {'problem': 'shifted-prism'})
            self.assertGreater(len(set(case['b'])), 1)
            if case['op'] == 'cg':
                self.assertEqual(multiply(case['a'], case['expected']), case['b'])
            else:
                self.assertEqual(multiply(case['a'], case['b']), case['expected'])

    def test_prism_render_preserves_measurements_and_explains_separate_live_grid(self):
        data = dict(suite='sparse-prism-v1', results=[dict(implementation='python', operation='cg',
                    size=17, unknowns=34, nnz=136, problem='shifted-prism', median_ns=123)],
                    implementations=[], machine={})
        original = copy.deepcopy(data)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'index.html'
            harness.render(data, path, live_override={'available': False})
            html = path.read_text()
        self.assertIn('Cycle length m means 2m', html)
        self.assertIn('Cycle length / unknowns', html)
        self.assertIn('Number of unknowns · twice the cycle length', html)
        self.assertIn('separate grid diffusion and flow example', html)
        self.assertNotIn('Number of unknowns · grid width²', html)
        self.assertNotIn('data-report="sparse" aria-current="page"', html)
        embedded = re.search(r'<script id="data" type="application/json">(.*?)</script>', html, re.S).group(1)
        self.assertEqual(json.loads(embedded)['results'], data['results'])
        self.assertEqual(data, original)

    def test_cli_graph_metadata_describes_the_selected_workload(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / 'graph-run'
            with patch.object(sys, 'argv', ['sparse.py', '--workload', 'prism', '--sizes', '3',
                                          '--languages', 'python', '--samples', '1', '--output', str(destination)]), \
                 patch.object(harness.run, 'build_one', return_value={'id': 'python', 'status': 'passed'}), \
                 patch.object(harness, 'verify', return_value=True), \
                 patch.object(harness, 'benchmark', return_value=[]) as benchmark, \
                 patch.object(harness.run, 'require', return_value='test-revision'), \
                 patch.object(harness.run, 'fingerprint', return_value='test-source'), \
                 patch.object(harness, 'render'):
                self.assertEqual(harness.main(), 0)
            data = json.loads((destination / 'results.json').read_text())
        self.assertEqual(data['suite'], 'sparse-prism-v1')
        self.assertIn('not a scalable Ramanujan family', data['methodology']['matrix_layout'])
        self.assertIn('construction is excluded', data['methodology']['limitations'])
        self.assertEqual(benchmark.call_args.args[-1], 'prism')


if __name__ == '__main__':
    unittest.main()
