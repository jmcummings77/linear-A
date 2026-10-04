"""Check that the published comparison catalog uses recorded report data."""
from html.parser import HTMLParser
import json
from pathlib import Path

from compare import select_public


class _CatalogData(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.values = []
        self.active = False

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "script" and attrs.get("id") == "snapshots":
            if attrs.get("type") != "application/json":
                raise ValueError("catalog snapshots must have type application/json")
            self.values.append("")
            self.active = True

    def handle_endtag(self, tag):
        if tag == "script":
            self.active = False

    def handle_data(self, data):
        if self.active:
            self.values[-1] += data


def validate_catalog(repo_root):
    """Check the public catalog after its source reports pass the provenance gate.

    Local comparison exports may use arbitrary input files. The catalog deployed
    by Pages must instead be a projection of the repository's verified reports.
    """
    repo_root = Path(repo_root)
    path = repo_root / "benchmarks/comparison/index.html"
    if path.is_symlink():
        raise ValueError("comparison catalog must not be a symlink")
    if not path.exists():
        return

    parser = _CatalogData()
    parser.feed(path.read_text(encoding="utf-8"))
    if parser.active or len(parser.values) != 1:
        raise ValueError("comparison catalog must embed exactly one completed snapshots script")
    snapshots = json.loads(parser.values[0])
    if not isinstance(snapshots, list) or not 2 <= len(snapshots) <= 30:
        raise ValueError("comparison catalog must contain 2–30 snapshots")

    reports = [select_public(json.loads(source.read_text(encoding="utf-8")))
               for source in sorted((repo_root / "benchmarks/reports").glob("*/results.json"))]
    for index, snapshot in enumerate(snapshots):
        if not isinstance(snapshot, dict) or snapshot.get("data") not in reports:
            raise ValueError(
                "comparison snapshot %d does not match a recorded report; "
                "regenerate the public catalog with benchmarks/compare.py --catalog" % (index + 1))
