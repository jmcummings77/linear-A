"""WebAssembly participates in toolchain discovery, portable profiling and provenance."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "benchmarks"))
import profiles
import report
import run
from publication import PublicSanitizer


class WasmIntegrationTests(unittest.TestCase):
    def test_missing_emcc_is_reported_as_unavailable(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "report"
            argv = ["run.py", "--implementations", "wasm", "--verify-only",
                    "--require-all", "--output", str(destination)]
            with patch.object(sys, "argv", argv), \
                 patch.object(run, "BUILD", Path(directory) / "build"), \
                 patch.dict(run.os.environ, {"EMCC": "/missing-toolchain/emcc"}), \
                 patch.object(run.shutil, "which", return_value=None), \
                 patch.object(run, "require", side_effect=["revision\n", ""]), \
                 patch.object(run, "fingerprint", return_value="source-fingerprint"), \
                 patch.object(report, "render"), patch("builtins.print"):
                self.assertEqual(run.main(), 1)
            data = json.loads((destination / "results.json").read_text())
            self.assertEqual(len(data["implementations"]), 1)
            implementation = data["implementations"][0]
            self.assertEqual(implementation["id"], "wasm")
            self.assertEqual(implementation["status"], "unavailable")
            self.assertIn("EMCC", implementation["error"])
            self.assertEqual(data["results"], [])

    def test_linux_wasm_profile_uses_node_and_preserves_sample_chronology(self):
        # Repeated WASM samples surrounding a JS sample must remain ordered,
        # rather than being reduced to aggregate flame-graph weights.
        raw_profile = {"nodes": [
            {"id": 1, "callFrame": {"functionName": "(root)"}, "children": [2]},
            {"id": 2, "callFrame": {"functionName": "benchMultiply",
                "url": "file:///fixture/runner.mjs", "lineNumber": 12}, "children": [3]},
            {"id": 3, "callFrame": {"functionName": "m_multiply",
                "url": "wasm://wasm/module-hash", "lineNumber": 0}},
        ], "samples": [3, 2, 3], "timeDeltas": [1000, 2000, 3000]}
        implementation = {"id": "wasm", "runner": ["/fixture/node", "/fixture/runner.mjs"],
                          "env": {"LINEAR_A_WASM_MODULE": "/fixture/matrix.mjs"}}
        iterations, checksum = 12, 7.0

        def node_process(command, **kwargs):
            self.assertEqual(command[0], implementation["runner"][0])
            self.assertIn("--cpu-prof", command)
            self.assertIn(implementation["runner"][1], command)
            self.assertEqual(kwargs["env"]["LINEAR_A_WASM_MODULE"], "/fixture/matrix.mjs")
            profile_dir = next(arg.split("=", 1)[1] for arg in command if arg.startswith("--cpu-prof-dir="))
            profile_name = next(arg.split("=", 1)[1] for arg in command if arg.startswith("--cpu-prof-name="))
            (Path(profile_dir) / profile_name).write_text(json.dumps(raw_profile))
            return subprocess.CompletedProcess(command, 0, stdout=json.dumps({
                "elapsed_ns": 1000, "iterations": iterations, "checksum": checksum * iterations,
            }), stderr="")

        sanitizer = PublicSanitizer(root="/fixture", home="/home/profile-test-user")
        with tempfile.TemporaryDirectory() as directory:
            with patch("publication.PublicSanitizer", return_value=sanitizer), \
                 patch.object(profiles.platform, "system", return_value="Linux"), \
                 patch.object(profiles, "_calibrate", return_value=(iterations, checksum)), \
                 patch.object(profiles.subprocess, "run", side_effect=node_process) as process, \
                 patch.object(profiles, "_native") as native:
                collected = profiles.collect_profiles([implementation], directory, seed=17)
            self.assertEqual(len(collected), 1)
            profile = collected[0]
            self.assertEqual(profile["status"], "available", profile.get("note"))
            self.assertTrue(profile["chronological"])
            js_stack = ["(root)", "benchMultiply (runner.mjs:13)"]
            wasm_stack = js_stack + ["m_multiply (module-hash:1)"]
            self.assertEqual(profile["stacks"], [wasm_stack, js_stack, wasm_stack])
            self.assertEqual(profile["weights"], [1.0, 2.0, 3.0])
            self.assertIn("WebAssembly", profile["note"])
            self.assertEqual(profile["workload"]["iterations"], iterations)
            self.assertTrue((Path(directory) / profile["raw_file"]).is_file())
            process.assert_called_once()
            native.assert_not_called()

    def test_wasm_source_changes_affect_fingerprint_without_generated_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "ports/wasm/runner.mjs"
            source.parent.mkdir(parents=True)
            source.write_text("export const revision = 1;\n")
            with patch.object(run, "ROOT", root):
                before = run.fingerprint()
                source.write_text("export const revision = 2;\n")
                after = run.fingerprint()
                self.assertNotEqual(before, after)
                generated = root / "ports/wasm/dist/matrix.mjs"
                generated.parent.mkdir()
                generated.write_text("generated module content\n")
                self.assertEqual(after, run.fingerprint())


if __name__ == "__main__":
    unittest.main()
