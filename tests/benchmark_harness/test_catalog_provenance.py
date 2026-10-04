import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "benchmarks"))
import compare
from check_catalog import validate_catalog
from compare import select_public


class CatalogProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = self.root / "benchmarks/reports/run/results.json"
        self.source.parent.mkdir(parents=True)
        self.report = {
            "schema_version": 1, "revision": "a" * 40, "dirty": False,
            "machine": {}, "methodology": {},
            "implementations": [{"id": "c", "status": "passed"}],
            "results": [{"implementation": "c", "operation": "multiply",
                         "size": 2, "status": "passed",
                         "samples": [{"elapsed_ns": 100, "iterations": 2, "ns_per_op": 50}]}],
        }
        self.source.write_text(json.dumps(self.report))
        self.catalog = self.root / "benchmarks/comparison/index.html"
        self.catalog.parent.mkdir(parents=True)
        self.snapshots = [{"label": label, "data": select_public(self.report)}
                          for label in ("Baseline", "Candidate")]
        self.write_catalog()

    def write_catalog(self):
        self.catalog.write_text('<script id="snapshots" type="application/json">' +
                                json.dumps(self.snapshots) + '</script>')

    def test_saved_projection_allows_restyling_and_label_changes(self):
        self.snapshots[0]["label"] = "An explanatory label"
        self.write_catalog()
        self.catalog.write_text('<h1>Restyled</h1>' + self.catalog.read_text())
        validate_catalog(self.root)

    def test_published_paths_are_independent_of_validator_username_and_checkout(self):
        command = ["dotnet", "build", "benchmarks/dotnet/Runner.csproj",
                   "--output", ".build/c/runner"]
        self.report["implementations"][0]["build_commands"] = [command]
        self.source.write_text(json.dumps(self.report))
        for snapshot in self.snapshots:
            snapshot["data"]["implementations"][0]["build_commands"] = [command]
        self.write_catalog()

        # These usernames collide with a public source filename or language ID.
        # Verification must not reinterpret either as the current host's identity.
        hosts = (("/Users/runner", "/Users/runner/projects/linear-A"),
                 ("/home/runner", "/home/runner/work/linear-A/linear-A"),
                 ("/home/c", "/home/c/linear-A"))
        for home, checkout in hosts:
            with self.subTest(home=home), \
                 patch("publication.Path.home", return_value=Path(home)), \
                 patch.object(compare, "ROOT", Path(checkout)):
                validate_catalog(self.root)
                compare.render([("Baseline", self.report), ("Candidate", self.report)],
                               self.catalog, published=True)
                validate_catalog(self.root)
                encoded = self.catalog.read_text().split(
                    '<script id="snapshots" type="application/json">', 1)[1].split('</script>', 1)[0]
                rendered = json.loads(encoded)[0]["data"]
                self.assertEqual(rendered["implementations"][0]["id"], "c")
                self.assertEqual(rendered["implementations"][0]["build_commands"], [command])
                self.assertEqual(rendered["results"], self.snapshots[0]["data"]["results"])

    def test_altered_measurements_or_provenance_cannot_bypass_report_gate(self):
        original = copy.deepcopy(self.snapshots)
        for field in ("samples", "revision", "dirty"):
            with self.subTest(field=field):
                self.snapshots = copy.deepcopy(original)
                data = self.snapshots[1]["data"]
                if field == "samples":
                    data["results"][0]["samples"][0]["ns_per_op"] = 1
                else:
                    data[field] = "b" * 40 if field == "revision" else True
                self.write_catalog()
                with self.assertRaisesRegex(ValueError, "does not match a recorded report"):
                    validate_catalog(self.root)

    def test_catalog_must_be_refreshed_after_replacing_its_source_report(self):
        self.report["results"][0]["samples"][0]["ns_per_op"] = 60
        self.source.write_text(json.dumps(self.report))
        with self.assertRaisesRegex(ValueError, "regenerate the public catalog"):
            validate_catalog(self.root)

    def test_missing_duplicate_or_unclosed_snapshots_are_rejected(self):
        original = self.catalog.read_text()
        for content in ("<html></html>", original + original, original.removesuffix("</script>")):
            with self.subTest(content=content):
                self.catalog.write_text(content)
                with self.assertRaisesRegex(ValueError, "exactly one completed"):
                    validate_catalog(self.root)

    def test_malformed_catalog_and_unknown_sources_are_rejected(self):
        for snapshots in ({}, [self.snapshots[0]], [None, None], [{"data": {}}, {"data": {}}]):
            with self.subTest(snapshots=snapshots):
                self.snapshots = snapshots
                self.write_catalog()
                with self.assertRaises(ValueError):
                    validate_catalog(self.root)

    def test_symlink_catalog_is_rejected(self):
        target = self.root / "catalog.html"
        self.catalog.rename(target)
        self.catalog.symlink_to(target)
        with self.assertRaisesRegex(ValueError, "symlink"):
            validate_catalog(self.root)

    def test_repository_without_a_catalog_can_validate_its_reports(self):
        self.catalog.unlink()
        validate_catalog(self.root)


if __name__ == "__main__":
    unittest.main()
