"""Published measurements retain a recoverable commit and unchanged artifacts."""
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "benchmarks"))
from check_provenance import artifact_hashes, check_reports, validate_report


class ProvenanceGateTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.git("init", "--quiet")
        (self.root / "source.py").write_text("print('measured source')\n")
        self.git("add", "source.py")
        self.git("-c", "user.name=Benchmark Test", "-c", "user.email=benchmark@example.test",
                 "commit", "--quiet", "-m", "Measured source")
        self.revision = self.git("rev-parse", "HEAD")
        self.tree = self.git("rev-parse", "HEAD^{tree}")
        self.report = self.root / "benchmarks/reports/test"
        self.report.mkdir(parents=True)
        self.data = {"revision": self.revision, "dirty": False, "results": [{"elapsed_ns": 123}]}
        self.write_results()
        (self.report / "profiles").mkdir()
        (self.report / "profiles/raw.json").write_text('{"samples": [1, 2, 3]}\n')
        self.provenance = {
            "schema_version": 1, "capture_method": "clean-checkout-v1",
            "revision": self.revision, "source_tree": self.tree,
            "command": ["python3", "benchmarks/run.py", "--implementations", "python",
                        "--require-all", "--output", "<output>"],
            "python_version": "3.12.0", "artifacts": artifact_hashes(self.report),
        }
        self.write_provenance()

    def git(self, *arguments):
        return subprocess.check_output(["git", "-C", str(self.root), *arguments],
                                       text=True, stderr=subprocess.PIPE).strip()

    def write_results(self):
        (self.report / "results.json").write_text(json.dumps(self.data) + "\n")
        (self.report / "index.html").write_text(
            '<html><script type="application/json" id="data">' + json.dumps(self.data) +
            "</script></html>")

    def write_provenance(self):
        (self.report / "provenance.json").write_text(json.dumps(self.provenance) + "\n")

    def make_legacy(self):
        (self.report / "provenance.json").unlink()
        self.data["dirty"] = True
        self.write_results()
        registry = {"schema_version": 1, "reason": "Original dirty source is unrecoverable.",
                    "reports": {"benchmarks/reports/test/results.json": hashlib.sha256(
                        (self.report / "results.json").read_bytes()).hexdigest()}}
        (self.root / "benchmarks/legacy-results.json").write_text(json.dumps(registry))

    def assert_rejected(self, message):
        errors, _ = check_reports(self.root)
        self.assertTrue(errors)
        self.assertIn(message, "\n".join(errors))

    def test_clean_report_validates_from_an_available_commit(self):
        validate_report(self.report, repo_root=self.root)
        self.assertEqual(check_reports(self.root), ([], []))

    def test_dirty_or_mismatched_revision_is_rejected_even_with_matching_hashes(self):
        for field, value in (("dirty", True), ("dirty", 0), ("revision", "a" * 40)):
            with self.subTest(field=field, value=value):
                self.data = {"revision": self.revision, "dirty": False, field: value}
                self.write_results()
                self.provenance["artifacts"] = artifact_hashes(self.report)
                self.write_provenance()
                self.assert_rejected("dirty: false")

    def test_result_and_profile_changes_invalidate_digests(self):
        for relative in ("results.json", "profiles/raw.json"):
            with self.subTest(relative=relative):
                path = self.report / relative
                original = path.read_bytes()
                path.write_bytes(original + b" ")
                self.assert_rejected("artifact set or SHA-256 digest")
                path.write_bytes(original)

    def test_added_or_missing_artifacts_are_rejected(self):
        extra = self.report / "profiles/extra.json"
        extra.write_text("{}")
        self.assert_rejected("artifact set or SHA-256 digest")
        extra.unlink()
        (self.report / "profiles/raw.json").unlink()
        self.assert_rejected("artifact set or SHA-256 digest")

    def test_report_restyling_preserves_provenance_but_measurement_changes_do_not(self):
        html = self.report / "index.html"
        html.write_text('<h1>New presentation</h1><script id="data" type="application/json">' +
                        json.dumps(self.data) + "</script>")
        self.assertEqual(check_reports(self.root), ([], []))
        html.write_text('<script id="data" type="application/json">{}</script>')
        self.assert_rejected("embedded measurements differ")

    def test_missing_index_or_duplicate_data_is_rejected(self):
        html = self.report / "index.html"
        html.write_text(html.read_text() * 2)
        self.assert_rejected("exactly one")
        html.unlink()
        self.assert_rejected("index.html")

    def test_missing_provenance_unknown_or_incomplete_report_is_rejected(self):
        (self.report / "provenance.json").unlink()
        self.assert_rejected("missing valid provenance")
        (self.report / "results.json").unlink()
        self.assert_rejected("missing valid provenance")

    def test_legacy_is_accepted_only_at_the_exact_registered_digest(self):
        self.make_legacy()
        errors, warnings = check_reports(self.root)
        self.assertEqual(errors, [])
        self.assertEqual(len(warnings), 1)
        self.assertIn("unrecoverable", warnings[0])
        self.data["results"][0]["elapsed_ns"] = 456
        self.write_results()
        self.assert_rejected("missing valid provenance")

    def test_legacy_html_cannot_present_different_measurements(self):
        self.make_legacy()
        (self.report / "index.html").write_text(
            '<script id="data" type="application/json">{}</script>')
        self.assert_rejected("embedded measurements differ")

    def test_invalid_new_sidecar_cannot_fall_back_to_legacy(self):
        self.make_legacy()
        self.write_provenance()
        self.assert_rejected("dirty: false")

    def test_missing_commit_wrong_tree_and_unreachable_commit_are_rejected(self):
        self.provenance["source_tree"] = "0" * 40
        self.write_provenance()
        self.assert_rejected("source_tree does not match")
        self.provenance["source_tree"] = self.tree
        self.provenance["revision"] = self.data["revision"] = "0" * 40
        self.write_results()
        self.provenance["artifacts"] = artifact_hashes(self.report)
        self.write_provenance()
        self.assert_rejected("fetch full Git history")
        orphan = self.git("-c", "user.name=Benchmark Test", "-c", "user.email=benchmark@example.test",
                          "commit-tree", self.tree, "-m", "Unreachable measured source")
        self.provenance["revision"] = self.data["revision"] = orphan
        self.write_results()
        self.provenance["artifacts"] = artifact_hashes(self.report)
        self.write_provenance()
        self.assert_rejected("ancestor of the published checkout")

    def test_unsafe_commands_and_abbreviations_are_rejected(self):
        original = copy.deepcopy(self.provenance["command"])
        for option in ("--no-build", "--no-b", "--render-only", "--render", "--out=elsewhere",
                       "--help", "-h", "--", "--amd-library=external.so", "--amd-l"):
            with self.subTest(option=option):
                self.provenance["command"] = original[:-2] + [option] + original[-2:]
                self.write_provenance()
                self.assert_rejected("forbidden or abbreviated option")
        self.provenance["command"] = ["python3", "untracked.py", "--output", "<output>"]
        self.write_provenance()
        self.assert_rejected("approved Python runner")

    def test_invalid_schema_and_artifact_paths_are_rejected(self):
        self.provenance["schema_version"] = True
        self.write_provenance()
        self.assert_rejected("schema_version")
        self.provenance["schema_version"] = 1
        self.provenance["artifacts"]["../elsewhere"] = "0" * 64
        self.write_provenance()
        self.assert_rejected("relative measurement paths")

    def test_symlink_artifacts_are_rejected(self):
        (self.report / "profiles/link.json").symlink_to(self.root / "source.py")
        self.assert_rejected("symlink artifact")


if __name__ == "__main__":
    unittest.main()
