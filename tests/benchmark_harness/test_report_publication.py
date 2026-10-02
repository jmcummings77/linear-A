"""Public artifacts keep recorded measurements while omitting local metadata."""
import copy
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'benchmarks'))
import report
from publication import PublicSanitizer


class PublicationPipelineTests(unittest.TestCase):
    def test_publish_copies_public_json_html_and_profiles_without_measuring(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'input'
            output = root / 'output'
            profile = source / 'profiles/python/samples.json'
            profile.parent.mkdir(parents=True)
            raw = {'samples': [{'at_ns': 123456789, 'stack': [
                'multiply (/Users/private-person/project/ports/python/matrix.py:52)']}],
                'end_ns': 123458901, 'environment': {'API_TOKEN': 'secret-value'},
                'hostname': 'private-machine'}
            profile.write_text(json.dumps(raw))
            data = {'schema_version': 1, 'created_at': '2026-10-02T14:02:03+00:00',
                'revision': 'ab' * 20, 'dirty': True, 'source_sha256': 'cd' * 32,
                'machine': {'os': 'Darwin', 'architecture': 'arm64', 'logical_cpus': 8,
                            'hostname': 'private-machine'},
                'environment': {'TOKEN': 'secret-value'}, 'seed': 17, 'suite': 'quick',
                'implementations': [{'id': 'python', 'name': 'Python', 'status': 'passed',
                    'toolchain': 'Python 3.13.7', 'build_commands': [[
                        '/Users/private-person/bin/python3.13', '-m', 'compileall',
                        '/Users/private-person/project/ports/python']]}],
                'results': [{'implementation': 'python', 'operation': 'multiply', 'size': 48,
                    'status': 'passed', 'iterations': 4, 'median_ns': 90.125,
                    'samples': [{'elapsed_ns': 360.5, 'iterations': 4,
                                 'checksum': -11.0625, 'ns_per_op': 90.125}]}],
                'profiles': [{'implementation': 'python', 'status': 'available',
                    'chronological': True, 'stacks': [raw['samples'][0]['stack']], 'weights': [0.002112],
                    'raw_file': 'profiles/python/samples.json'}]}
            original = copy.deepcopy(data)
            with patch.object(report, 'PublicSanitizer', return_value=PublicSanitizer(
                    root='/Users/private-person/project', home='/Users/private-person')), \
                 patch('subprocess.run', side_effect=AssertionError('Publication must not execute commands')):
                published = report.publish(data, output / 'index.html', source, include_live=False)
            self.assertEqual(data, original)
            self.assertEqual(json.loads(profile.read_text()), raw)
            public_json = json.loads((output / 'results.json').read_text())
            self.assertEqual(public_json, published)
            html = (output / 'index.html').read_text()
            embedded = json.loads(re.search(r'<script id="data" type="application/json">(.*?)</script>', html).group(1))
            self.assertEqual(embedded, published)
            self.assertEqual(published['results'], data['results'])
            for key in ('created_at', 'revision', 'dirty', 'source_sha256', 'seed'):
                self.assertEqual(published[key], data[key])
            self.assertEqual(published['implementations'][0]['toolchain'], 'Python 3.13.7')
            clean_raw = json.loads((output / 'profiles/python/samples.json').read_text())
            self.assertEqual(clean_raw['end_ns'], raw['end_ns'])
            self.assertEqual(clean_raw['samples'][0]['at_ns'], raw['samples'][0]['at_ns'])
            self.assertIn('multiply (ports/python/matrix.py:52)', clean_raw['samples'][0]['stack'])
            for path in output.rglob('*'):
                if path.is_file():
                    text = path.read_text()
                    for private in ('/Users/', 'private-person', 'private-machine', 'secret-value', 'API_TOKEN'):
                        self.assertNotIn(private, text)
            before = {p.relative_to(output): p.read_bytes() for p in output.rglob('*') if p.is_file()}
            report.publish(published, output / 'index.html', include_live=False)
            self.assertEqual(before, {p.relative_to(output): p.read_bytes() for p in output.rglob('*') if p.is_file()})

    def test_direct_html_render_also_sanitizes_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / 'index.html'
            data = {'implementations': [], 'results': [], 'profiles': [],
                    'env': {'TOKEN': 'secret-value'},
                    'machine': {'hostname': 'private-machine', 'architecture': 'arm64'}}
            report.render(data, destination, include_live=False)
            self.assertNotIn('private-machine', destination.read_text())
            self.assertNotIn('secret-value', destination.read_text())


if __name__ == '__main__':
    unittest.main()
