"""Publishing raw profiles keeps measurements and rejects unsafe artifact paths."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "benchmarks"))
from profile_publication import (
    sanitize_existing_profiles, sanitize_go_raw, sanitize_julia,
    sanitize_macos_sample, sanitize_node, sanitize_python,
)
from profiles import parse_go_raw, parse_macos_sample, parse_node, parse_python
from publication import PublicSanitizer


class ProfilePublicationTests(unittest.TestCase):
    def setUp(self):
        self.sanitizer = PublicSanitizer(root=Path("/Users/privateperson/project"), home=Path("/Users/privateperson"))

    def test_node_allowlist_keeps_sample_order_positions_and_clocks(self):
        source = {"nodes": [
            {"id": 1, "hitCount": 0, "callFrame": {"functionName": "(root)", "scriptId": "0"}, "children": [2]},
            {"id": 2, "hitCount": 3, "callFrame": {"functionName": "multiply", "scriptId": "12",
             "url": "file:///Users/privateperson/project/ports/typescript/matrix.ts", "lineNumber": 42, "columnNumber": 7,
             "machine": "private-machine"}, "positionTicks": [{"line": 43, "ticks": 3, "private": "secret"}],
             "unrecognized": {"private": "secret"}},
        ], "startTime": 123456789, "endTime": 123462789,
            "samples": [2, 1, 2], "timeDeltas": [1000, 2000, 3000], "machine": "private-machine"}
        original = copy.deepcopy(source)
        clean = sanitize_node(source, self.sanitizer)
        self.assertEqual(source, original)
        self.assertEqual(set(clean), {"nodes", "startTime", "endTime", "samples", "timeDeltas"})
        for key in ("startTime", "endTime", "samples", "timeDeltas"):
            self.assertEqual(clean[key], source[key])
        self.assertEqual(clean["nodes"][1]["positionTicks"], [{"line": 43, "ticks": 3}])
        self.assertEqual(clean["nodes"][1]["hitCount"], 3)
        self.assertEqual(clean["nodes"][1]["callFrame"]["url"], "ports/typescript/matrix.ts")
        self.assertEqual(parse_node(clean), parse_node(source))
        self.assertNotIn("private", json.dumps(clean))

    def test_python_and_julia_preserve_every_interval_weight_and_checksum(self):
        python = {"samples": [{"at_ns": 1234000,
                   "stack": ["main (/Users/privateperson/project/ports/python/runner.py:90)", "multiply (matrix.py:97)"],
                   "thread_id": 12345}, {"at_ns": 4567000, "stack": ["cleanup (matrix.py:5)"]}],
                  "end_ns": 6789000, "host": "private-machine"}
        clean = sanitize_python(python, self.sanitizer)
        self.assertEqual(clean["end_ns"], python["end_ns"])
        self.assertEqual([s["at_ns"] for s in clean["samples"]], [s["at_ns"] for s in python["samples"]])
        self.assertEqual(parse_python(clean)[1], parse_python(python)[1])
        self.assertEqual(clean["samples"][0]["stack"], ["main (ports/python/runner.py:90)", "multiply (matrix.py:97)"])
        self.assertNotIn("thread_id", clean["samples"][0])
        julia = {"implementation": "julia", "kind": "sampled", "chronological": False,
                 "unit": "milliseconds", "checksum": -12.875, "stacks": [["multiply (/tmp/build/LinearA.jl:42)"]],
                 "weights": [1.0], "hostname": "private-machine"}
        result = sanitize_julia(julia, self.sanitizer)
        self.assertEqual(result["weights"], julia["weights"])
        self.assertEqual(result["checksum"], julia["checksum"])
        self.assertEqual(result["stacks"], [["multiply (LinearA.jl:42)"]])
        self.assertNotIn("hostname", result)

    def test_go_keeps_timestamps_counts_inlining_and_source_lines(self):
        source = """PeriodType: cpu nanoseconds
Period: 10000000
Time: 2026-10-02 02:41:59.480309 -0700 PDT
Duration: 4.65
Machine: private-machine
Samples:
samples/count cpu/nanoseconds
  3 30000000: 1 2
  2 20000000: 2
                hostname:[private-machine]
Locations
  1: 0x100 M=1 math.IsInf /opt/toolchain/src/math/isinf.go:1:0 s=1
          matrix.Multiply /Users/privateperson/project/ports/go/matrix.go:10:0 s=2
  2: 0x200 M=1 main.main /Users/privateperson/project/ports/go/main.go:20:0 s=3
Mappings
1: 0x100 0x500 /Users/privateperson/project/.build/runner secret-build-id
"""
        clean = sanitize_go_raw(source, self.sanitizer)
        self.assertIn("Time: 2026-10-02 02:41:59.480309 -0700 PDT", clean)
        self.assertIn("Duration: 4.65", clean)
        self.assertIn("matrix.Multiply ports/go/matrix.go:10:0 s=2", clean)
        self.assertNotIn("private", clean)
        self.assertNotIn("secret-build-id", clean)
        before, weights = parse_go_raw(source)
        after, actual = parse_go_raw(clean)
        self.assertEqual(weights, actual)
        self.assertEqual(after, [[self.sanitizer.text(frame) for frame in stack] for stack in before])

    def test_native_drops_process_machine_and_thread_identity_without_changing_tree(self):
        source = """Process: private-runner [12345]
Path: /Users/privateperson/project/runner
Hardware Model: private-machine
Analysis of sampling private-runner (pid 12345) every 1.5 millisecond
Call graph:
    10 Thread_8888 DispatchQueue_999: com.apple.main-thread (serial)
      10 main (in runner) + 32 [0x10] /Users/privateperson/project/ports/c/runner.c:139
        7 multiply (in runner) + 20 [0x20] /Users/privateperson/project/ports/c/matrix.c:177
        + 5 dot (in runner) + 12 [0x30]
        2 allocate (in runner) + 4 [0x40]
    10 Thread_9999 private-thread-name
      10 unrelated (in system) + 0 [0x50]
Total number in stack:
Binary Images:
  0x10 private-machine /Users/privateperson/project/runner
"""
        clean = sanitize_macos_sample(source, self.sanitizer)
        self.assertNotIn("private", clean)
        self.assertNotIn("12345", clean)
        self.assertNotIn("8888", clean)
        self.assertNotIn("9999", clean)
        self.assertIn("ports/c/matrix.c:177", clean)
        self.assertEqual(parse_macos_sample(clean), parse_macos_sample(source))
        self.assertEqual(sum(parse_macos_sample(clean)[1]), 15)

    @staticmethod
    def _profile():
        return {"nodes": [{"id": 1, "callFrame": {"functionName": "multiply", "url": "/tmp/private/matrix.mjs", "lineNumber": 20}}],
                "startTime": 100, "endTime": 200, "samples": [1], "timeDeltas": [100], "hostname": "private-machine"}

    def test_copy_scans_ignored_unreferenced_profiles_and_keeps_source_and_report(self):
        with tempfile.TemporaryDirectory() as temporary:
            source, destination = Path(temporary) / "source", Path(temporary) / "destination"
            raw = source / "profiles/wasm/profile.cpuprofile"
            raw.parent.mkdir(parents=True)
            raw.write_text(json.dumps(self._profile()))
            ignored = source / "profiles/typescript/profile.cpuprofile"
            ignored.parent.mkdir()
            ignored.write_bytes(raw.read_bytes())
            before = raw.read_bytes()
            report = {"profiles": [{"raw_file": "profiles/wasm/profile.cpuprofile", "weights": [0.1]}], "results": [{"elapsed_ns": 9}]}
            snapshot = copy.deepcopy(report)
            paths = sanitize_existing_profiles(report, source, self.sanitizer, destination)
            self.assertEqual(paths, ["profiles/typescript/profile.cpuprofile", "profiles/wasm/profile.cpuprofile"])
            self.assertEqual(raw.read_bytes(), before)
            self.assertEqual(report, snapshot)
            published = destination / paths[1]
            self.assertNotIn("private", published.read_text())
            self.assertEqual(parse_node(json.loads(published.read_text())), parse_node(json.loads(before)))
            once = published.read_bytes()
            sanitize_existing_profiles(report, destination, self.sanitizer)
            self.assertEqual(published.read_bytes(), once)

    def test_rejects_traversal_missing_unsupported_and_symlink_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            outside = root / "outside.cpuprofile"
            outside.write_text(json.dumps(self._profile()))
            source = root / "report"
            source.mkdir()
            for reference in ("../outside.cpuprofile", str(outside), "profiles/../outside.cpuprofile", "profiles\\outside.cpuprofile", "profiles/missing.cpuprofile"):
                with self.subTest(reference=reference), self.assertRaises((ValueError, OSError)):
                    sanitize_existing_profiles({"profiles": [{"raw_file": reference}]}, source, self.sanitizer)
            files = source / "profiles/wasm"
            files.mkdir(parents=True)
            link = files / "profile.cpuprofile"
            link.symlink_to(outside)
            with self.assertRaises(ValueError):
                sanitize_existing_profiles({}, source, self.sanitizer)
            link.unlink()
            (files / "cpu.pprof").write_bytes(b"unsupported binary")
            with self.assertRaises(ValueError):
                sanitize_existing_profiles({}, source, self.sanitizer)
            self.assertIn("private-machine", outside.read_text())

    def test_rejects_destination_symlink_without_modifying_external_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, destination, outside = root / "source", root / "destination", root / "outside"
            raw = source / "profiles/wasm/profile.cpuprofile"
            raw.parent.mkdir(parents=True)
            raw.write_text(json.dumps(self._profile()))
            destination.mkdir()
            outside.mkdir()
            (destination / "profiles").symlink_to(outside, target_is_directory=True)
            with self.assertRaises((ValueError, OSError)):
                sanitize_existing_profiles({}, source, self.sanitizer, destination)
            self.assertEqual(list(outside.iterdir()), [])

    def test_copy_also_sanitizes_stale_destination_profiles_and_source_wins_collisions(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source, destination = root / "source", root / "destination"
            raw = source / "profiles/wasm/profile.cpuprofile"
            raw.parent.mkdir(parents=True)
            raw.write_text(json.dumps(self._profile()))
            old = self._profile()
            old["timeDeltas"] = [20]
            for name in ("wasm", "typescript"):
                stale = destination / ("profiles/%s/profile.cpuprofile" % name)
                stale.parent.mkdir(parents=True)
                stale.write_text(json.dumps(old))
            unrelated = destination / "unrelated.json"
            unrelated.write_text("unrelated private data")
            original = raw.read_bytes()
            paths = sanitize_existing_profiles({}, source, self.sanitizer, destination)
            self.assertEqual(len(paths), 2)
            for relative in paths:
                self.assertNotIn("private", (destination / relative).read_text())
            self.assertEqual(json.loads((destination / "profiles/wasm/profile.cpuprofile").read_text())["timeDeltas"], [100])
            self.assertEqual(json.loads((destination / "profiles/typescript/profile.cpuprofile").read_text())["timeDeltas"], [20])
            self.assertEqual(raw.read_bytes(), original)
            self.assertEqual(unrelated.read_text(), "unrelated private data")
            unsupported = destination / "profiles/wasm/cpu.pprof"
            unsupported.write_bytes(b"unsupported binary")
            with self.assertRaises(ValueError):
                sanitize_existing_profiles({}, source, self.sanitizer, destination)


if __name__ == "__main__":
    unittest.main()
