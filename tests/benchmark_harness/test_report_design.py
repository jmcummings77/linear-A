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
import report_design


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


class DirectoryElements(Elements):
    def __init__(self):
        super().__init__()
        self.cards = []
        self.current_card = None
        self.in_metadata = False
        self.div_families = []

    def handle_starttag(self, tag, attributes):
        super().handle_starttag(tag, attributes)
        attrs = dict(attributes)
        classes = attrs.get('class', '').split()
        if tag == 'div':
            family = attrs.get('id') if 'report-family' in classes else (self.div_families[-1] if self.div_families else None)
            self.div_families.append(family)
        if tag == 'a' and 'report' in classes:
            self.current_card = {'href': attrs.get('href'), 'family': self.div_families[-1] if self.div_families else None, 'metadata': ''}
            self.cards.append(self.current_card)
        if tag == 'span' and 'run-meta' in classes:
            self.in_metadata = True

    def handle_endtag(self, tag):
        if tag == 'div' and self.div_families:
            self.div_families.pop()
        if tag == 'span':
            self.in_metadata = False
        if tag == 'a':
            self.current_card = None

    def handle_data(self, text):
        if self.current_card is not None and self.in_metadata:
            self.current_card['metadata'] += text


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


class ReportDirectoryTests(unittest.TestCase):
    def render(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / 'directory.html'
            report_design.render_directory(destination)
            html = destination.read_text()
        parsed = DirectoryElements()
        parsed.feed(html)
        return html, parsed

    def test_every_current_benchmark_card_has_metadata_from_its_own_measurements(self):
        expected_names = {'latest', 'determinants', 'eigen', 'vectors', 'sparse',
                          'gmres', 'ilu', 'ordering', 'cholesky', 'amd', 'ic0'}
        html, parsed = self.render()
        cards = {card['href'].rstrip('/'): card for card in parsed.cards if card['metadata']}
        self.assertEqual(set(cards), expected_names)
        self.assertNotRegex(html, r'__[A-Z0-9_]+_META__')
        for name in expected_names:
            with self.subTest(report=name):
                data = json.loads((report_design.HERE / 'reports' / name / 'results.json').read_text())
                metadata = cards[name]['metadata']
                self.assertTrue(metadata.startswith(data['created_at'].split('T')[0] + ' · '))
                passed_implementations = [item for item in data['implementations'] if item['status'] == 'passed']
                passed_results = [row for row in data['results'] if row['status'] == 'passed']
                self.assertIn(f'{len(passed_implementations)} verified implementations', metadata)
                self.assertTrue(metadata.endswith(f'{len(passed_results)} timing results'))

    def test_directory_navigation_targets_exist_and_solver_cards_follow_the_learning_sequence(self):
        _, parsed = self.render()  # Parsing also rejects duplicate IDs.
        for link in parsed.links:
            href = link.get('href', '')
            if href.startswith('#'):
                with self.subTest(anchor=href):
                    self.assertIn(href[1:], parsed.ids)
        solver_cards = [card['href'] for card in parsed.cards
                        if card['family'] == 'sparse-solvers']
        self.assertEqual(solver_cards, ['sparse/', 'gmres/', 'ilu/'])
        self.assertIn('ilu/#explore', [link.get('href') for link in parsed.links])

    def test_new_template_marker_discovers_hyphenated_report_and_excludes_unsuccessful_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'pages').mkdir()
            (root / 'reports/future-solver2').mkdir(parents=True)
            (root / 'reports/empty').mkdir()
            (root / 'pages/template.html').write_text(
                '<html><head><style></style></head><body>'
                '<a class="report" href="future-solver2/"><span class="run-meta">__FUTURE_SOLVER2_META__</span></a>'
                '<a class="report" href="empty/"><span class="run-meta">__EMPTY_META__</span></a>'
                '</body></html>')
            data = {'created_at': '2031-04-05T23:59:00+00:00',
                    'implementations': [{'status': status} for status in ['passed', 'failed', 'passed', 'unavailable']],
                    'results': [{'status': status} for status in ['failed', 'passed', 'unavailable', 'passed', 'failed']]}
            source = root / 'reports/future-solver2/results.json'
            source.write_text(json.dumps(data))
            (root / 'reports/empty/results.json').write_text('{}')
            # Styling is independently tested; this fixture exercises discovery and metadata.
            with patch.object(report_design, 'HERE', root), patch.object(report_design, 'apply_report_design', side_effect=lambda html: html):
                html, parsed = self.render()
            self.assertNotRegex(html, r'__[A-Z0-9_]+_META__')
            self.assertEqual(parsed.cards[0]['metadata'], '2031-04-05 · 2 verified implementations · 2 timing results')
            self.assertEqual(parsed.cards[1]['metadata'], 'Undated run · 0 verified implementations · 0 timing results')
            self.assertEqual(json.loads(source.read_text()), data)


if __name__ == '__main__':
    unittest.main()
