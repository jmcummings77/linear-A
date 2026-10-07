"""Saved report integrity and presentation-only regeneration."""
import copy
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'experiments/nonnormal-gmres/run.py'
SPEC = importlib.util.spec_from_file_location('nonnormal_html_report', SCRIPT.with_name('report.py'))
report = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(report)


class NonnormalReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run([sys.executable, str(SCRIPT), '--output', directory],
                           check=True, capture_output=True, text=True)
            cls.data = json.loads((Path(directory) / 'results.json').read_text())

    def test_all_recorded_results_survive_rendering_and_escaping(self):
        data = copy.deepcopy(self.data)
        data['test_note'] = '</script><script>alert("bad")</script>&@@RUN_ROWS@@'
        original = copy.deepcopy(data)
        rendered = report.render_html(data)
        payload = re.findall(r'<script id="data" type="application/json">(.*?)</script>', rendered, re.S)
        self.assertEqual(len(payload), 1)
        self.assertEqual(json.loads(payload[0]), data)
        self.assertEqual(data, original)
        self.assertNotIn('<script>alert("bad")', rendered)
        for identifier in ('results', 'explore', 'polynomials', 'graphs', 'method'):
            self.assertIn('id="' + identifier + '"', rendered)
        self.assertIn('Prism C_32 x K_2', rendered)
        self.assertIn('iteration limit', rendered)
        self.assertIn('stagnation', rendered)
        self.assertIn('exact zeros', rendered.lower())
        self.assertNotRegex(rendered, r'<(?:script|link)[^>]+(?:src|href)="https?://')

    def test_restyle_preserves_saved_data_and_uses_recorded_limits(self):
        data = copy.deepcopy(self.data)
        data['max_iterations'] = 2
        rendered = report.render_html(data)
        self.assertIn('-- --limit 2', rendered)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'results.json'
            raw = json.dumps(self.data, indent=2).encode()
            path.write_bytes(raw)
            subprocess.run([sys.executable, str(SCRIPT), '--render-only', str(path)],
                           check=True, capture_output=True, text=True)
            self.assertEqual(path.read_bytes(), raw)
            self.assertEqual((path.parent / 'index.html').read_text(), report.render_html(self.data))
            self.assertFalse((path.parent / 'RESULTS.md').exists())


if __name__ == '__main__':
    unittest.main()
