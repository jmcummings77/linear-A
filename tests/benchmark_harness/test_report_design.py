"""Report organization follows the recorded evidence without altering that evidence."""
import copy
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'benchmarks'))
import report


class Elements(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = {}
        self.links = []

    def handle_starttag(self, tag, attributes):
        attrs = dict(attributes)
        if 'id' in attrs:
            if attrs['id'] in self.ids:
                raise AssertionError('Duplicate ID: ' + attrs['id'])
            self.ids[attrs['id']] = attrs
        if tag == 'a':
            self.links.append(attrs)


class ReportDesignTests(unittest.TestCase):
    def fixture(self, operation='multiply'):
        return {'schema_version': 1, 'created_at': '2026-10-02T12:00:00+00:00',
                'revision': 'a' * 40, 'source_sha256': 'b' * 64, 'dirty': False,
                'seed': 17, 'suite': 'quick', 'machine': {'os': 'Darwin', 'architecture': 'arm64'},
                'methodology': {}, 'implementations': [{'id': 'c', 'name': 'C', 'status': 'passed'}],
                'results': [{'implementation': 'c', 'operation': operation, 'size': 8,
                             'status': 'passed', 'median_ns': 10, 'min_ns': 9, 'max_ns': 11,
                             'mad_ns': 1, 'iterations': 1, 'samples': [{'ns_per_op': 10}]}], 'profiles': []}

    def render(self, data, **kwargs):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / 'index.html'
            report.render(data, destination, **kwargs)
            html = destination.read_text()
        parser = Elements()
        parser.feed(html)
        return html, parser

    def test_titles_follow_recorded_report_scope(self):
        cases = [('multiply', 'Matrix performance'), ('determinant_spd', 'Determinant algorithms'),
                 ('eigen_symmetric', 'Eigenvalues &amp; eigenvectors'), ('cross', 'Vectors &amp; rotations')]
        for operation, title in cases:
            with self.subTest(operation=operation):
                html, parsed = self.render(self.fixture(operation), include_live=False)
                self.assertIn('<h1>' + title + '</h1>', html)
                self.assertIn('<title>linear-A · ' + title + '</title>', html)
                self.assertNotIn('hidden', parsed.ids['timings-panel'])

    def test_saved_evidence_precedes_live_tools_and_navigation_targets_exist(self):
        data = self.fixture()
        data['profiles'] = [{'implementation': 'c', 'status': 'available', 'stacks': [['multiply']], 'weights': [1]}]
        live = {'available': True, 'geometry_worker_source': 'geometry', 'accuracy_worker_source': 'accuracy'}
        html, parsed = self.render(data, live_override=live)
        ids = ['timings-panel', 'profiles-panel', 'correctness-panel', 'method-panel', 'explore', 'live-panel', 'accuracy-panel', 'geometry-panel']
        positions = [html.index('id="' + identifier + '"') for identifier in ids]
        self.assertEqual(positions, sorted(positions))
        self.assertTrue(all('hidden' not in parsed.ids[identifier] for identifier in ids))
        for link in parsed.links:
            target = link.get('href', '')
            if target.startswith('#'):
                self.assertIn(target[1:], parsed.ids)
        self.assertIn('id="profile-search"', html)
        self.assertIn('type="search"', html)

    def test_verification_hides_unavailable_sections_and_their_links(self):
        data = self.fixture()
        data['results'] = []
        html, parsed = self.render(data, include_live=False)
        self.assertIn('<h1>Correctness verification</h1>', html)
        for identifier in ['timings-panel', 'profiles-panel', 'explore', 'live-panel', 'accuracy-panel', 'geometry-panel', 'csv']:
            self.assertIn('hidden', parsed.ids[identifier])
        for link in parsed.links:
            if link.get('href') in {'#timings-panel', '#profiles-panel', '#explore'}:
                self.assertIn('hidden', link)
        self.assertIn('correctness checks only', html)

    def test_failed_timings_do_not_masquerade_as_verification_only(self):
        data = self.fixture()
        data['results'][0]['status'] = 'failed'
        html, parsed = self.render(data, include_live=False)
        self.assertIn('<h1>Matrix performance</h1>', html)
        self.assertIn('No successful timing samples', html)
        self.assertIn('hidden', parsed.ids['timings-panel'])

    def test_single_wasm_report_identity_follows_checks_timings_and_profiles(self):
        for kind, title in [('verification', 'WebAssembly verification'),
                            ('timings', 'WebAssembly performance'),
                            ('profiles', 'WebAssembly profiling check')]:
            with self.subTest(kind=kind):
                data = self.fixture()
                data['implementations'] = [{'id': 'wasm', 'name': 'WebAssembly', 'status': 'passed'}]
                data['results'][0]['implementation'] = 'wasm'
                if kind == 'verification':
                    data['results'] = []
                elif kind == 'profiles':
                    data['profiles'] = [{'implementation': 'wasm', 'status': 'available',
                                         'stacks': [['multiply']], 'weights': [1]}]
                html, _ = self.render(data, include_live=False)
                self.assertIn('<h1>' + title + '</h1>', html)
                self.assertIn('<title>linear-A · ' + title + '</title>', html)
                self.assertIn('for the WebAssembly implementation.', html)

    def test_regeneration_keeps_measurements_and_existing_bundle_without_building(self):
        data = self.fixture()
        live = {'available': True, 'sha256': 'c' * 64, 'wasm_base64': 'AGFzbQEAAAA=',
                'worker_source': '/* keep the published worker unchanged */', 'fixtures': []}
        original_data, original_live = copy.deepcopy(data), copy.deepcopy(live)
        with patch.object(report, 'live_bundle', side_effect=AssertionError('Must retain the saved bundle')):
            html, _ = self.render(data, live_override=live)
        for identifier, expected in [('data', original_data), ('live-data', original_live)]:
            embedded = re.search(r'<script id="' + identifier + r'" type="application/json">(.*?)</script>', html, re.S)
            self.assertEqual(json.loads(embedded.group(1)), expected)
        self.assertEqual(data, original_data)
        self.assertEqual(live, original_live)
        self.assertIn('<style id="report-design">', html)
        self.assertNotRegex(html, r'__(?:REPORT|TIMING|PROFILE|LIVE|ACCURACY|GEOMETRY)[A-Z_]*__')


if __name__ == '__main__':
    unittest.main()
