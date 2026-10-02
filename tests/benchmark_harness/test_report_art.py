"""Report artwork stays passive, self-contained, and separate from measurements."""
import copy
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "benchmarks"))
import report
from report_art import load_report_art


SVG_NAMESPACE = "http://www.w3.org/2000/svg"
SYMBOL_IDS = {"report-art-" + name for name in
              ("cube", "woven", "flow", "spiral", "mesh", "sphere")}


def sprite(contents, attributes=""):
    return (f'<svg xmlns="{SVG_NAMESPACE}" class="report-art-sprite" width="0" height="0" '
            f'aria-hidden="true" focusable="false" {attributes}>'
            f"<defs>{contents}</defs></svg>")


class ReportElements(HTMLParser):
    def __init__(self):
        super().__init__()
        self.elements = []

    def handle_starttag(self, tag, attrs):
        self.elements.append((tag, dict(attrs)))


class ReportArtTests(unittest.TestCase):
    def load_temporary(self, source):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "art.svg"
            path.write_text(source, encoding="utf-8")
            return load_report_art(path)

    def test_bundled_art_has_six_symbols_and_only_scoped_local_references(self):
        artwork = load_report_art()
        root = ET.fromstring(artwork)
        self.assertEqual(root.tag, f"{{{SVG_NAMESPACE}}}svg")
        for name, value in (("width", "0"), ("height", "0"),
                            ("aria-hidden", "true"), ("focusable", "false")):
            self.assertEqual(root.get(name), value)
        symbols = root.findall(f".//{{{SVG_NAMESPACE}}}symbol")
        self.assertEqual({symbol.get("id") for symbol in symbols}, SYMBOL_IDS)
        self.assertEqual(len(symbols), 6)
        ids = [element.get("id") for element in root.iter() if element.get("id")]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertTrue(all(identifier.startswith("report-art-") for identifier in ids))
        references = [value for element in root.iter() for name, value in element.attrib.items()
                      if name.rsplit("}", 1)[-1] == "href"]
        self.assertTrue(references)
        self.assertTrue(all(reference.startswith("#report-art-") for reference in references))
        self.assertTrue(all(reference[1:] in ids for reference in references))
        self.assertNotIn("url(", artwork.lower())

    def test_passive_geometry_and_approved_theme_paints_are_accepted(self):
        source = sprite('''
          <symbol id="report-art-example" viewBox="0 0 10 10">
            <g fill="none" stroke="currentColor" transform="translate(1, 2)">
              <path d="M0 0 L1 1"/>
              <line x1="0" y1="0" x2="1" y2="1"/>
              <polygon points="0,0 1,0 0,1" fill="var(--art-face-1)"/>
              <polyline points="0,0 1,0 0,1"/>
              <rect x="0" y="0" width="1" height="2" fill="var(--art-face-2)"/>
              <circle cx="1" cy="1" r="1"/>
              <ellipse cx="1" cy="1" rx="1" ry="2"/>
            </g>
          </symbol>
          <use href="#report-art-example"/>
        ''')
        root = ET.fromstring(self.load_temporary(source))
        self.assertIsNotNone(root.find(f".//{{{SVG_NAMESPACE}}}ellipse"))
        self.assertEqual(root.find(f".//{{{SVG_NAMESPACE}}}use").get("href"),
                         "#report-art-example")

    def test_active_elements_attributes_and_metadata_are_rejected(self):
        variants = {
            "script": "<script>alert(1)</script>",
            "foreign_object": "<foreignObject><div>text</div></foreignObject>",
            "style_element": "<style>path { fill: red; }</style>",
            "style_attribute": '<path d="M0 0" style="fill:currentColor"/>',
            "event_attribute": '<path d="M0 0" onload="alert(1)"/>',
            "mixed_case_event": '<path d="M0 0" onClick="alert(1)"/>',
            "animation": '<animate attributeName="href" to="https://example.invalid/a"/>',
            "external_image": '<image href="https://example.invalid/image.svg"/>',
            "metadata": "<metadata>private creator information</metadata>",
            "text_element": "<text>private creator information</text>",
            "text_node": "<g>private creator information</g>",
        }
        for name, contents in variants.items():
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.load_temporary(sprite(contents))

    def test_external_references_and_unapproved_paints_are_rejected(self):
        for reference in ("https://example.invalid/art.svg#shape", "//example.invalid/shape",
                          "data:image/svg+xml,svg", "javascript:alert(1)",
                          "file:///private/art.svg", "other.svg#report-art-example"):
            with self.subTest(reference=reference), self.assertRaises(ValueError):
                self.load_temporary(sprite(f'<use href="{reference}"/>'))
        with self.assertRaises(ValueError):
            self.load_temporary(sprite('<use xlink:href="https://example.invalid/art.svg"/>',
                                       'xmlns:xlink="http://www.w3.org/1999/xlink"'))
        for paint in ("url(https://example.invalid/paint)", "url(#report-art-example)",
                      "var(--unapproved)"):
            for attribute in ("fill", "stroke"):
                with self.subTest(paint=paint, attribute=attribute), self.assertRaises(ValueError):
                    self.load_temporary(sprite(f'<path d="M0 0" {attribute}="{paint}"/>'))

    def test_doctype_and_entity_declarations_are_rejected(self):
        variants = (
            '<!DOCTYPE svg SYSTEM "https://example.invalid/art.dtd">' + sprite(""),
            '<!DOCTYPE svg [<!ENTITY creator "private creator information">]>'
            + sprite("<g>&creator;</g>"),
            '<!DOCTYPE svg [<!ENTITY source SYSTEM "file:///private/art.svg">]>'
            + sprite("<g>&source;</g>"),
        )
        for source in variants:
            with self.subTest(source=source[:80]), self.assertRaises(ValueError):
                self.load_temporary(source)

    def test_duplicate_unscoped_and_dangling_identifiers_are_rejected(self):
        variants = {
            "duplicate": '<g id="report-art-same"/><path id="report-art-same" d="M0 0"/>',
            "unscoped_id": '<g id="unscoped"/>',
            "unscoped_reference": '<use href="#unscoped"/>',
            "dangling_reference": '<use href="#report-art-missing"/>',
        }
        for name, contents in variants.items():
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.load_temporary(sprite(contents))

    def test_render_embeds_artwork_without_changing_recorded_data(self):
        data = {"schema_version": 1, "created_at": "2026-10-02T12:34:56+00:00",
                "revision": "a" * 40, "source_sha256": "b" * 64, "dirty": True,
                "seed": 17, "suite": "quick", "machine": {}, "methodology": {},
                "implementations": [{"id": "python", "name": "Python", "status": "passed",
                                     "toolchain": "Python 3.13.7"}],
                "results": [{"implementation": "python", "operation": "multiply", "size": 8,
                             "status": "passed", "median_ns": 1234.5,
                             "samples": [{"elapsed_ns": 3703.5, "ns_per_op": 1234.5,
                                          "iterations": 3, "checksum": -12.75}]}],
                "profiles": [{"implementation": "python", "status": "available",
                              "chronological": True, "weights": [1.25],
                              "stacks": [["multiply (ports/python/matrix.py:52)"]]}]}
        original = copy.deepcopy(data)
        artwork = load_report_art()
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "index.html"
            with patch("subprocess.run", side_effect=AssertionError("Rendering must not rerun benchmarks")):
                report.render(data, destination, include_live=False)
            html = destination.read_text(encoding="utf-8")
        self.assertEqual(data, original)
        self.assertNotIn("__REPORT_ART__", html)
        self.assertIn(artwork, html)
        self.assertLess(html.index("<body>"), html.index(artwork))
        self.assertLess(html.index(artwork), html.index("</body>"))
        embedded = re.search(r'<script id="data" type="application/json">(.*?)</script>', html, re.S)
        self.assertIsNotNone(embedded)
        self.assertEqual(json.loads(embedded.group(1)), original)
        parsed = ReportElements()
        parsed.feed(html)
        self.assertFalse(any(tag == "script" and "src" in attrs for tag, attrs in parsed.elements))
        self.assertFalse(any(tag == "link" and attrs.get("rel") == "stylesheet"
                             for tag, attrs in parsed.elements))
        self.assertFalse(any(tag in ("img", "image") for tag, _ in parsed.elements))
        used_symbols = {attrs.get("href", "")[1:] for tag, attrs in parsed.elements if tag == "use"}
        self.assertTrue(SYMBOL_IDS.issubset(used_symbols))


if __name__ == "__main__":
    unittest.main()
