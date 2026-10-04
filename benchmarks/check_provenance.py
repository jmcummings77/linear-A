#!/usr/bin/env python3
"""Reject published benchmark data without recoverable source provenance.

Historical measurements are accepted only at their explicitly recorded digest.
They remain disclosed as unreproducible; a digest is not a source snapshot.
"""
import argparse
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys

from check_catalog import validate_catalog

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT_RUNNERS = {
    "machine-code-dot": "experiments/machine-code-dot/run.py",
    "matmul-locality": "experiments/matmul-locality/run.py",
}
RUNNERS = frozenset({
    "benchmarks/run.py", "benchmarks/determinants.py", "benchmarks/sparse.py",
    "benchmarks/gmres_bench.py", "benchmarks/ilu_bench.py",
    "benchmarks/ordering_bench.py", "benchmarks/cholesky_bench.py",
    "benchmarks/amd_bench.py", "benchmarks/ic0_bench.py", "benchmarks/multigrid_bench.py",
}) | frozenset(EXPERIMENT_RUNNERS.values())


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def artifact_hashes(directory):
    """Hash measurement artifacts; HTML may be restyled without changing data."""
    directory = Path(directory)
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError("report directory must be a real directory")
    artifacts = {}
    for path in sorted(directory.rglob("*")):
        if path.is_symlink():
            raise ValueError("symlink artifact is not allowed: %s" % path)
        if not path.is_file():
            continue
        relative = path.relative_to(directory).as_posix()
        if relative == "provenance.json" or path.suffix.lower() == ".html":
            continue
        artifacts[relative] = sha256(path)
    return artifacts


def _json(path):
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("%s must contain a JSON object" % path)
    return value


class _ReportData(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=False)
        self.values = []
        self.active = False

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "script" and attrs.get("id") == "data":
            if attrs.get("type") != "application/json":
                raise ValueError("report data script must have type application/json")
            self.values.append("")
            self.active = True

    def handle_endtag(self, tag):
        if tag == "script":
            self.active = False

    def handle_data(self, data):
        if self.active:
            self.values[-1] += data


def _check_html(directory, data):
    path = directory / "index.html"
    if path.is_symlink():
        raise ValueError("index.html must not be a symlink")
    parser = _ReportData()
    parser.feed(path.read_text(encoding="utf-8"))
    if parser.active or len(parser.values) != 1:
        raise ValueError("index.html must embed exactly one completed report data script")
    if json.loads(parser.values[0]) != data:
        raise ValueError("index.html embedded measurements differ from results.json")


def _hex(value, length):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{%d}" % length, value)


def _relative_path(value):
    return (isinstance(value, str) and bool(value) and "\\" not in value
            and not PurePosixPath(value).is_absolute()
            and all(part not in ("", ".", "..") for part in value.split("/")))


def forbidden_options(runner):
    flags = ("--no-build", "--render-only", "--output", "--help", "--amd-library")
    if runner == EXPERIMENT_RUNNERS["machine-code-dot"]:
        # Published commands must not retain host-specific executable paths.
        # Toolchains are selected through PATH and their versions are recorded.
        flags += ("--dotnet", "--julia")
    elif runner == EXPERIMENT_RUNNERS["matmul-locality"]:
        flags += ("--verify-only", "--sanitize")
    return flags


def _check_command(command):
    if (not isinstance(command, list) or not all(isinstance(arg, str) for arg in command)
            or len(command) < 4 or command[0] != "python3" or command[1] not in RUNNERS
            or command[-2:] != ["--output", "<output>"]):
        raise ValueError("command must record an approved Python runner and --output <output>")
    for arg in command[2:-2]:
        option = arg.split("=", 1)[0]
        if (option in ("--", "-h") or (option.startswith("--")
                and any(flag.startswith(option) for flag in forbidden_options(command[1])))):
            raise ValueError("command contains a forbidden or abbreviated option: %s" % arg)


def validate_report(directory, repo_root=ROOT):
    """Validate one new report, raising ValueError/OSError on invalid provenance."""
    directory, repo_root = Path(directory), Path(repo_root)
    actual = artifact_hashes(directory)
    data = _json(directory / "results.json")
    provenance = _json(directory / "provenance.json")
    if type(provenance.get("schema_version")) is not int or provenance["schema_version"] != 1:
        raise ValueError("unsupported provenance schema_version")
    if provenance.get("capture_method") != "clean-checkout-v1":
        raise ValueError("capture_method must be clean-checkout-v1")
    revision, tree = provenance.get("revision"), provenance.get("source_tree")
    if not _hex(revision, 40) or not _hex(tree, 40):
        raise ValueError("revision and source_tree must be full Git SHA-1 identifiers")
    if data.get("dirty") is not False or data.get("revision") != revision:
        raise ValueError("results must record dirty: false and the provenance revision")
    if not isinstance(provenance.get("python_version"), str) or not provenance["python_version"].strip():
        raise ValueError("python_version is required")
    _check_command(provenance.get("command"))
    artifacts = provenance.get("artifacts")
    if (not isinstance(artifacts, dict) or "results.json" not in artifacts
            or any(not _relative_path(path) or not _hex(digest, 64)
                   for path, digest in artifacts.items())):
        raise ValueError("artifacts must map relative measurement paths to SHA-256 digests")
    if artifacts != actual:
        raise ValueError("measurement artifact set or SHA-256 digest differs from provenance")
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(repo_root), "rev-parse", "--verify", revision + "^{commit}"],
            text=True, stderr=subprocess.PIPE).strip()
        recorded_tree = subprocess.check_output(
            ["git", "-C", str(repo_root), "rev-parse", "--verify", revision + "^{tree}"],
            text=True, stderr=subprocess.PIPE).strip()
    except subprocess.CalledProcessError as error:
        raise ValueError("recorded revision is unavailable; fetch full Git history") from error
    if commit != revision or recorded_tree != tree:
        raise ValueError("source_tree does not match the recorded Git commit")
    ancestor = subprocess.run(
        ["git", "-C", str(repo_root), "merge-base", "--is-ancestor", revision, "HEAD"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ancestor.returncode:
        raise ValueError("recorded revision must be an ancestor of the published checkout")
    _check_html(directory, data)


def _legacy_reports(repo_root):
    path = repo_root / "benchmarks/legacy-results.json"
    if not path.exists():
        return {}
    registry = _json(path)
    reports, artifacts = registry.get("reports"), registry.get("artifacts")
    if (type(registry.get("schema_version")) is not int or registry["schema_version"] != 2
            or not isinstance(registry.get("reason"), str) or not registry["reason"].strip()
            or not isinstance(reports, dict)
            or any(not re.fullmatch(r"(?:benchmarks/reports/[^/]+|experiments/(?:machine-code-dot|matmul-locality)/results)/results\.json", name)
                   or not _relative_path(name) or not _hex(digest, 64)
                   for name, digest in reports.items())
            or not isinstance(artifacts, dict) or artifacts.keys() != reports.keys()
            or any(not isinstance(files, dict) or files.get("results.json") != reports[name]
                   or any(not _relative_path(path) or not _hex(digest, 64)
                          for path, digest in files.items())
                   for name, files in artifacts.items())):
        raise ValueError("invalid legacy-results.json registry")
    return artifacts


def check_reports(repo_root=ROOT, report_paths=None):
    """Return (errors, warnings) for published reports or selected directories."""
    repo_root = Path(repo_root).resolve()
    errors, warnings = [], []
    check_catalog = report_paths is None
    try:
        legacy = _legacy_reports(repo_root)
    except (ValueError, OSError) as error:
        return [str(error)], warnings
    if report_paths is None:
        # Include incomplete report directories so deleting results cannot evade the gate.
        report_paths = sorted(path for path in (repo_root / "benchmarks/reports").glob("*")
                              if path.is_dir())
        report_paths += [repo_root / "experiments" / name / "results"
                         for name in EXPERIMENT_RUNNERS
                         if (repo_root / "experiments" / name / "results").exists()]
    for directory in report_paths:
        directory = Path(directory)
        if not directory.is_absolute():
            directory = repo_root / directory
        try:
            relative = (directory / "results.json").relative_to(repo_root).as_posix()
            if (directory / "provenance.json").exists() or (directory / "provenance.json").is_symlink():
                validate_report(directory, repo_root)
            elif relative in legacy and artifact_hashes(directory) == legacy[relative]:
                _check_html(directory, _json(directory / "results.json"))
                warnings.append("%s: unchanged legacy measurement; original measured source is "
                                "unrecoverable, so this is not a reproducible source snapshot" % relative)
            else:
                raise ValueError("missing valid provenance; new or changed measurements must "
                                 "be captured with benchmarks/reproduce.py")
        except (ValueError, OSError) as error:
            errors.append("%s: %s" % (directory, error))
    if check_catalog and not errors:
        try:
            validate_catalog(repo_root)
        except (ValueError, OSError, TypeError, AttributeError) as error:
            errors.append("benchmarks/comparison/index.html: %s" % error)
    return errors, warnings


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs="*", type=Path, help="optional report directories")
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    errors, warnings = check_reports(args.repo_root, args.reports or None)
    for warning in warnings:
        print("WARNING: " + warning, file=sys.stderr)
    for error in errors:
        print("ERROR: " + error, file=sys.stderr)
    if not errors:
        print("Benchmark provenance checks passed (%d legacy reports)." % len(warnings))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
