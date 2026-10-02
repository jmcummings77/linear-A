"""Geometry reports embed separate trace artifacts without rewriting historical data."""
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
from publication import PublicSanitizer


class EmbeddedData(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.current = None
        self.parts = {}

    def handle_starttag(self, tag, attrs):
        identifier = dict(attrs).get("id")
        if tag == "script" and identifier in ("data", "live-data"):
            self.current = identifier
            self.parts[identifier] = []

    def handle_endtag(self, tag):
        if tag == "script":
            self.current = None

    def handle_data(self, value):
        if self.current is not None:
            self.parts[self.current].append(value)

    def value(self, identifier):
        return json.loads("".join(self.parts[identifier]))


class GeometryReportTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.standard = self.root / ".build/wasm"
        self.trace = self.root / ".build/wasm-trace"
        self.standard.mkdir(parents=True)
        self.trace.mkdir()
        self.header = b"\0asm\x01\0\0\0"
        self.trace_binary = self.header + b"\x01\x01\0"  # Valid empty type section.
        (self.standard / "matrix.wasm").write_bytes(self.header)
        (self.trace / "matrix.wasm").write_bytes(self.trace_binary)
        (self.standard / "matrix.mjs").write_text("export default async function ordinary() {}\n")
        (self.trace / "matrix.mjs").write_text("export default async function instrumented() {}\n")
        files = {
            "ports/wasm/matrix.mjs": "export async function createMatrixAPI() {}\n",
            "benchmarks/live-worker.mjs": "export async function runChecks() {}\n",
            "benchmarks/geometry-worker.mjs": "export async function computeGeometry() {}\n",
            "benchmarks/live-report.mjs": "void 0;\n",
            "benchmarks/geometry-report.mjs": "void 1;\n",
        }
        for name, content in files.items():
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content)
        self.manifest = {"path": "ports/c/matrix.c", "start_line": 173,
                         "lines": ["matrix_status m_multiply(const matrix *a, const matrix *b, matrix *out) {",
                                   "    double sum = 0.0;",
                                   "    m_trace_record(0, 0, 0, 1, 1, 1, sum, __LINE__);", "}"],
                         "source_sha256": "a" * 64}
        self.write_manifest()
        self.data = {"schema_version": 1, "created_at": "2026-10-02T12:34:56+00:00",
                     "source_sha256": "b" * 64, "revision": "c" * 40, "seed": 17, "suite": "quick",
                     "implementations": [{"id": "wasm", "name": "WebAssembly", "status": "passed"}],
                     "results": [{"implementation": "wasm", "operation": "multiply", "size": 16,
                                  "status": "passed", "iterations": 3, "median_ns": 1234.5,
                                  "samples": [{"elapsed_ns": 3703.5, "ns_per_op": 1234.5,
                                               "iterations": 3, "checksum": -12.75}]}],
                     "profiles": [{"implementation": "wasm", "status": "available", "chronological": True,
                                   "stacks": [["m_multiply (matrix.wasm:1)"]], "weights": [1.25]}]}
        self.root_patch = patch.object(report, "ROOT", self.root)
        self.root_patch.start()
        self.addCleanup(self.root_patch.stop)

    def write_manifest(self):
        (self.trace / "trace-source.json").write_text(json.dumps(self.manifest))

    def test_geometry_and_trace_are_separate_from_the_ordinary_module(self):
        before = copy.deepcopy(self.data)
        bundle = report.live_bundle(self.data)
        self.assertIs(bundle["available"], True)
        self.assertEqual(bundle["geometry_worker_source"], (self.root / "benchmarks/geometry-worker.mjs").read_text())
        self.assertIn("ordinary", bundle["module_source"])
        self.assertIn("instrumented", bundle["trace_module_source"])
        self.assertEqual(base64.b64decode(bundle["wasm_base64"]), self.header)
        self.assertEqual(base64.b64decode(bundle["trace_wasm_base64"]), self.trace_binary)
        self.assertEqual(bundle["geometry_source"], self.manifest)
        self.assertEqual(bundle["geometry_source"]["start_line"] + 2, 175)
        self.assertIn("m_trace_record", bundle["geometry_source"]["lines"][2])
        self.assertEqual(len(bundle["fixtures"]), 39 + len(report.eigen_fixtures()) + len(report.general_eigen_fixtures()) + len(report.vector_fixtures()))
        self.assertEqual(self.data, before)

    def test_hash_tracks_worker_trace_module_binary_and_build_source_manifest(self):
        hashes = [report.live_bundle(self.data)["sha256"]]
        for path, suffix in ((self.root / "benchmarks/geometry-worker.mjs", "// geometry update\n"),
                             (self.trace / "matrix.mjs", "// trace update\n")):
            path.write_text(path.read_text() + suffix)
            hashes.append(report.live_bundle(self.data)["sha256"])
        (self.trace / "matrix.wasm").write_bytes(self.trace_binary + b"\x03\x01\0")
        hashes.append(report.live_bundle(self.data)["sha256"])
        self.manifest["start_line"] += 1
        self.manifest["source_sha256"] = "d" * 64
        self.write_manifest()
        hashes.append(report.live_bundle(self.data)["sha256"])
        self.assertEqual(len(set(hashes)), len(hashes))
        self.assertEqual(self.data["source_sha256"], "b" * 64)

    def test_source_lines_come_from_the_trace_build_not_the_current_checkout(self):
        original = report.live_bundle(self.data)
        current = self.root / "ports/c/matrix.c"
        current.parent.mkdir(parents=True)
        current.write_text("// a later edit with different source line numbers\n" * 40
                           + "matrix_status m_multiply() {\n    return M_OK;\n}\n")
        bundle = report.live_bundle(self.data)
        self.assertEqual(bundle["geometry_source"], self.manifest)
        self.assertEqual(bundle["sha256"], original["sha256"])

    def test_missing_or_invalid_optional_trace_leaves_geometry_available(self):
        for name in ("matrix.mjs", "matrix.wasm"):
            path = self.trace / name
            before = path.read_bytes()
            for content in (None, b"\xff"):
                with self.subTest(name=name, content=content):
                    if content is None:
                        path.unlink()
                    else:
                        path.write_bytes(content)
                    bundle = report.live_bundle(self.data)
                    self.assertIs(bundle["available"], True)
                    self.assertIn("geometry_worker_source", bundle)
                    self.assertIn("wasm_base64", bundle)
                    self.assertNotIn("trace_wasm_base64", bundle)
                    self.assertTrue(bundle["trace_reason"])
                    path.write_bytes(before)

    def test_unreadable_optional_geometry_preserves_the_standard_live_bundle(self):
        geometry = self.root / "benchmarks/geometry-worker.mjs"
        original_read = Path.read_text

        def unreadable(path, *args, **kwargs):
            if path == geometry:
                raise PermissionError("optional geometry worker is unreadable")
            return original_read(path, *args, **kwargs)

        with patch.object(Path, "read_text", new=unreadable):
            bundle = report.live_bundle(self.data)
        self.assertIs(bundle["available"], True)
        self.assertIn("worker_source", bundle)
        self.assertNotIn("geometry_worker_source", bundle)
        geometry.write_bytes(b"\xff")
        bundle = report.live_bundle(self.data)
        self.assertIs(bundle["available"], True)
        self.assertNotIn("geometry_worker_source", bundle)

    def test_missing_source_manifest_does_not_substitute_unverified_current_source(self):
        (self.trace / "trace-source.json").unlink()
        current = self.root / "ports/c/matrix.c"
        current.parent.mkdir(parents=True)
        current.write_text("matrix_status m_multiply() {\n    return M_OK;\n}\n")
        bundle = report.live_bundle(self.data)
        self.assertIs(bundle["available"], True)
        self.assertIn("geometry_worker_source", bundle)
        self.assertIn("trace_wasm_base64", bundle)
        self.assertNotIn("geometry_source", bundle)

    def test_explicit_trace_directory_and_public_source_metadata(self):
        external = self.root / "custom-trace"
        external.mkdir()
        for name in ("matrix.mjs", "matrix.wasm", "trace-source.json"):
            (external / name).write_bytes((self.trace / name).read_bytes())
            (self.trace / name).unlink()
        # Only strings/comments and optional binary metadata carry these markers.
        (external / "matrix.mjs").write_text('const source = "/Users/fixtureperson/project/ports/wasm/trace.c:27";\nexport default async function instrumented() {}\n')
        (self.root / "benchmarks/geometry-worker.mjs").write_text('const source = "/Users/fixtureperson/project/benchmarks/geometry-worker.mjs:10";\nexport async function computeGeometry() {}\n')
        metadata = bytes([5]) + b"debug" + b"/Users/fixtureperson/private-build"
        (external / "matrix.wasm").write_bytes(self.trace_binary + b"\0" + bytes([len(metadata)]) + metadata)
        manifest = copy.deepcopy(self.manifest)
        manifest["lines"][1] = '    /* /Users/fixtureperson/project/ports/c/matrix.c:174 */ double sum = 0.0;'
        manifest["build_host"] = "private-host-marker"
        (external / "trace-source.json").write_text(json.dumps(manifest))
        sanitizer = PublicSanitizer(root="/Users/fixtureperson/project", home="/Users/fixtureperson")
        bundle = report.live_bundle(self.data, trace_directory=external, sanitizer=sanitizer)
        self.assertIs(bundle["available"], True)
        self.assertIn("ports/wasm/trace.c:27", bundle["trace_module_source"])
        self.assertIn("benchmarks/geometry-worker.mjs:10", bundle["geometry_worker_source"])
        self.assertIn("ports/c/matrix.c:174", bundle["geometry_source"]["lines"][1])
        self.assertEqual(set(bundle["geometry_source"]), set(self.manifest))
        self.assertEqual(bundle["geometry_source"]["start_line"], self.manifest["start_line"])
        self.assertEqual(bundle["geometry_source"]["source_sha256"], self.manifest["source_sha256"])
        self.assertEqual(len(bundle["geometry_source"]["lines"]), len(manifest["lines"]))
        self.assertNotIn("fixtureperson", json.dumps(bundle))
        self.assertNotIn("private-host-marker", json.dumps(bundle))
        self.assertNotIn(b"fixtureperson", base64.b64decode(bundle["trace_wasm_base64"]))

    def test_render_retains_historical_measurements_and_escapes_geometry_sources(self):
        hostile = '\n// </script><script id="injected">bad & metadata</script>\n'
        for path in (self.root / "benchmarks/geometry-worker.mjs", self.trace / "matrix.mjs"):
            path.write_text(path.read_text() + hostile)
        before = copy.deepcopy(self.data)
        destination = self.root / "geometry.html"
        report.render(self.data, destination, self.standard, True, self.trace)
        parsed = EmbeddedData()
        parsed.feed(destination.read_text())
        recorded, bundle = parsed.value("data"), parsed.value("live-data")
        self.assertEqual(recorded["results"], self.data["results"])
        self.assertEqual(recorded["profiles"], self.data["profiles"])
        self.assertEqual(recorded["created_at"], self.data["created_at"])
        self.assertEqual(recorded["source_sha256"], self.data["source_sha256"])
        self.assertIn(hostile, bundle["geometry_worker_source"])
        self.assertIn(hostile, bundle["trace_module_source"])
        self.assertNotIn('<script id="injected">', destination.read_text())
        self.assertIn("\\u003c/script>", destination.read_text())
        self.assertEqual(self.data, before)

    def test_cli_forwards_separate_trace_directory_without_measuring(self):
        source = self.root / "results.json"
        source.write_text(json.dumps(self.data))
        destination = self.root / "public/index.html"
        args = ["report.py", str(source), str(destination), "--wasm-dir", str(self.standard),
                "--trace-wasm-dir", str(self.trace)]
        with patch.object(sys, "argv", args), patch.object(report, "publish") as publish:
            report.main()
        publish.assert_called_once_with(self.data, destination, source.parent, self.standard, True, self.trace)


if __name__ == "__main__":
    unittest.main()
