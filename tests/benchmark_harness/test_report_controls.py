"""Shared control assets stay local and leave embedded benchmark evidence intact."""
import json
from pathlib import Path
import re
import sys
import unittest
from urllib.parse import unquote
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'benchmarks'))
from report_design import apply_report_design, HERE


class ReportControlsTests(unittest.TestCase):
    def test_controls_are_inlined_once_and_preserve_recorded_data_and_worker_code(self):
        data = {'results': [{'median_ns': 1.25, 'status': 'passed'}]}
        worker = {'worker_source': 'self.onmessage = () => postMessage(42)', 'wasm_base64': 'AGFzbQ=='}
        source = ('<html><head><style>button { color: red }</style></head><body>'
                  '<button id="play" disabled>Play</button>'
                  f'<script id="data" type="application/json">{json.dumps(data)}</script>'
                  f'<script id="live" type="application/json">{json.dumps(worker)}</script>'
                  '</body></html>')
        once = apply_report_design(source)
        twice = apply_report_design(once)
        for html in [once, twice]:
            self.assertEqual(html.count('<script id="report-controls">'), 1)
            self.assertEqual(html.count('<style id="report-design">'), 1)
            self.assertEqual(html.count('<script id="report-appearance">'), 1)
            self.assertIn((HERE / 'report-ui.js').read_text(), html)
            self.assertIn('<button id="play" disabled>Play</button>', html)
            for identifier, expected in [('data', data), ('live', worker)]:
                embedded = re.search(r'<script id="' + identifier + '" type="application/json">(.*?)</script>', html, re.S)
                self.assertEqual(json.loads(embedded[1]), expected)
            self.assertNotRegex(html, r'<(?:script|link)[^>]+(?:src|href)="https?://')

    def test_implicit_body_reports_receive_the_same_inline_controls(self):
        html = apply_report_design('<!doctype html><title>Report</title><style></style><main>Saved evidence</main>')
        self.assertIn('<script id="report-controls">', html)
        self.assertTrue(html.rstrip().endswith('</script>'))

    def test_svg_masks_are_self_contained_and_have_no_accessible_text(self):
        css = (HERE / 'report-ui.css').read_text()
        masks = re.findall(r'\[data-report-icon="([a-z-]+)"\] \{ --report-icon: url\("data:image/svg\+xml,([^\"]+)"\); \}', css)
        self.assertGreaterEqual(len(masks), 15)
        self.assertEqual(len({name for name, _ in masks}), len(masks))
        for name, encoded in masks:
            with self.subTest(icon=name):
                svg = ET.fromstring(unquote(encoded))
                self.assertEqual(svg.attrib['viewBox'], '0 0 24 24')
                self.assertFalse(''.join(svg.itertext()).strip())
                self.assertFalse(any('href' in attr for node in svg.iter() for attr in node.attrib))
        self.assertIn('content: "";', css)
        self.assertIn('button:disabled', css)
        self.assertIn('button[aria-pressed="true"]', css)


if __name__ == '__main__':
    unittest.main()
