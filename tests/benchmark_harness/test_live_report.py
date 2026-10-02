"""The offline report embeds a validated, independent browser benchmark bundle."""
import base64
import copy
from html.parser import HTMLParser
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "benchmarks"))
import report


class LiveDataParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.in_live_data = False
        self.parts = []
        self.matches = 0

    def handle_starttag(self, tag, attrs):
        if tag == "script" and dict(attrs).get("id") == "live-data":
            self.in_live_data = True
            self.matches += 1

    def handle_endtag(self, tag):
        if tag == "script":
            self.in_live_data = False

    def handle_data(self, value):
        if self.in_live_data:
            self.parts.append(value)


class LiveReportTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.build = self.root / ".build/wasm"
        self.build.mkdir(parents=True)
        self.wrapper = self.root / "ports/wasm/matrix.mjs"
        self.worker = self.root / "benchmarks/live-worker.mjs"
        self.ui = self.root / "benchmarks/live-report.mjs"
        self.wrapper.parent.mkdir(parents=True)
        self.worker.parent.mkdir(parents=True)
        self.wrapper.write_text("export async function createMatrixAPI() {}\n")
        self.worker.write_text("self.onmessage = () => {};\n")
        self.ui.write_text("void 0;\n")
        (self.root / "benchmarks/accuracy-worker.mjs").write_text("export function computeAccuracy() {}\n")
        (self.root / "benchmarks/accuracy-report.mjs").write_text("void 0;\n")
        self.module = self.build / "matrix.mjs"
        self.module.write_text("export default async function createModule() { return {}; }\n")
        self.binary = self.build / "matrix.wasm"
        self.binary.write_bytes(b"\0asm\x01\0\0\0")
        self.data = {"implementations": [{"id": "wasm", "status": "passed"}],
                     "results": [], "profiles": []}
        self.root_patch = patch.object(report, "ROOT", self.root)
        self.root_patch.start()
        self.addCleanup(self.root_patch.stop)

    def assert_unavailable(self, bundle):
        self.assertIs(bundle["available"], False)
        self.assertIsInstance(bundle.get("reason"), str)
        self.assertTrue(bundle["reason"].strip())

    def embedded_live_data(self, path):
        parser = LiveDataParser()
        parser.feed(path.read_text())
        self.assertEqual(parser.matches, 1)
        return json.loads("".join(parser.parts))

    def test_valid_default_bundle_contains_real_fixtures_and_independent_sources(self):
        original = copy.deepcopy(self.data)
        bundle = report.live_bundle(self.data)
        self.assertIs(bundle["available"], True, bundle.get("reason"))
        self.assertEqual(bundle["wrapper_source"], self.wrapper.read_text())
        self.assertEqual(bundle["module_source"], self.module.read_text())
        self.assertEqual(bundle["worker_source"], self.worker.read_text())
        self.assertEqual(base64.b64decode(bundle["wasm_base64"], validate=True), self.binary.read_bytes())
        self.assertEqual(bundle["byte_length"], len(self.binary.read_bytes()))
        self.assertRegex(bundle["sha256"], r"^[0-9a-f]{64}$")
        self.assertEqual(len(bundle["fixtures"]), 39 + len(report.eigen_fixtures()) + len(report.general_eigen_fixtures()) + len(report.vector_fixtures()) + len(report.solve_fixtures()))
        self.assertTrue(any(case["op"] == "eigen_symmetric" for case in bundle["fixtures"]))
        # The independent oracle returns Fraction values for multiplication.
        # Packaging must normalize them while retaining JSON booleans as bools.
        normalized = json.loads(json.dumps(bundle["fixtures"], allow_nan=False))
        self.assertEqual(normalized, bundle["fixtures"])
        multiplication = next(case for case in normalized if case["op"] == "multiply" and not case["invalid"])
        self.assertTrue(all(type(value) in (int, float) for value in multiplication["expected"]["values"]))
        triangular = next(case for case in normalized if case["op"] == "triangular")
        self.assertIs(type(triangular["expected"]["upper"]), bool)
        self.assertEqual(self.data, original)
        bundle["fixtures"][0]["a"]["values"][0] = 99999
        again = report.live_bundle(self.data)
        self.assertNotEqual(again["fixtures"][0]["a"]["values"][0], 99999)
        self.assertEqual(self.data, original)

    def test_hash_tracks_embedded_module_and_binary_contents(self):
        first = report.live_bundle(self.data, wasm_directory=self.build)
        self.module.write_text(self.module.read_text() + "// changed generated module\n")
        second = report.live_bundle(self.data, wasm_directory=self.build)
        self.assertNotEqual(first["sha256"], second["sha256"])
        # A valid, empty type section changes retained bytes without compiling.
        # Optional build/debug metadata is removed during publication.
        self.binary.write_bytes(self.binary.read_bytes() + b"\x01\x01\0")
        third = report.live_bundle(self.data, wasm_directory=self.build)
        self.assertNotEqual(second["sha256"], third["sha256"])
        self.assertEqual(base64.b64decode(third["wasm_base64"]), self.binary.read_bytes())

    def test_explicit_build_directory_works_without_the_default_build(self):
        custom = self.root / "custom-build"
        custom.mkdir()
        (custom / "matrix.mjs").write_text(self.module.read_text())
        (custom / "matrix.wasm").write_bytes(self.binary.read_bytes())
        self.binary.unlink()
        self.assert_unavailable(report.live_bundle(self.data))
        bundle = report.live_bundle(self.data, wasm_directory=custom)
        self.assertIs(bundle["available"], True, bundle.get("reason"))

    def test_live_execution_requires_a_passed_wasm_implementation_and_opt_in(self):
        for data, include_live in ((self.data, False),
                                   ({"implementations": []}, True),
                                   ({"implementations": [{"id": "wasm", "status": "failed"}]}, True),
                                   ({"implementations": [{"id": "wasm", "status": "unavailable"}]}, True)):
            with self.subTest(data=data, include_live=include_live):
                self.assert_unavailable(report.live_bundle(data, wasm_directory=self.build, include_live=include_live))

    def test_missing_build_artifacts_have_an_explicit_unavailable_reason(self):
        for artifact in (self.module, self.binary):
            with self.subTest(artifact=artifact.name):
                content = artifact.read_bytes()
                artifact.unlink()
                try:
                    self.assert_unavailable(report.live_bundle(self.data, wasm_directory=self.build))
                finally:
                    artifact.write_bytes(content)

    def test_bad_magic_truncated_header_and_wrong_binary_version_are_unavailable(self):
        for content in (b"", b"\0asm", b"not wasm", b"\0asm\x02\0\0\0"):
            with self.subTest(content=content):
                self.binary.write_bytes(content)
                bundle = report.live_bundle(self.data, wasm_directory=self.build)
                self.assert_unavailable(bundle)
                self.assertRegex(bundle["reason"].lower(), r"wasm|webassembly|binary|header")

    def test_source_read_errors_and_non_utf8_sources_degrade_gracefully(self):
        original_read = Path.read_text

        def blocked_worker(path, *args, **kwargs):
            if path == self.worker:
                raise PermissionError("worker source unreadable")
            return original_read(path, *args, **kwargs)

        with patch.object(Path, "read_text", new=blocked_worker):
            self.assert_unavailable(report.live_bundle(self.data, wasm_directory=self.build))
        self.wrapper.write_bytes(b"\xff")
        self.assert_unavailable(report.live_bundle(self.data, wasm_directory=self.build))

    def test_renderer_embeds_sources_without_allowing_them_to_close_the_json_script(self):
        hostile = '\n// </script><script id="intruder">boom & wow</script> __REPORT_DATA__ __LIVE_DATA__\n'
        for source in (self.wrapper, self.worker, self.module):
            source.write_text(source.read_text() + hostile)
        before = {path for path in self.root.rglob("*") if path.is_file()}
        destination = self.root / "report.html"
        report.render(self.data, destination, wasm_directory=self.build)
        embedded = self.embedded_live_data(destination)
        self.assertIs(embedded["available"], True)
        for key in ("wrapper_source", "worker_source", "module_source"):
            self.assertIn(hostile, embedded[key])
        html = destination.read_text()
        self.assertNotIn('<script id="intruder">', html)
        self.assertIn("\\u003c/script>", html)
        self.assertIn("\\u0026", html)
        after = {path for path in self.root.rglob("*") if path.is_file()}
        self.assertEqual(after - before, {destination})

    def test_missing_artifacts_do_not_prevent_rendering_a_historical_report(self):
        self.binary.unlink()
        destination = self.root / "historical.html"
        report.render(self.data, destination, wasm_directory=self.build)
        self.assert_unavailable(self.embedded_live_data(destination))
        self.assertIn('"implementations"', destination.read_text())


if __name__ == "__main__":
    unittest.main()
