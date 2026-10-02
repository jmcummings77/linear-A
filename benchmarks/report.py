"""Render a portable HTML report, optionally embedding live WebAssembly tests."""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import re

from reference import fixtures as arithmetic_fixtures
from eigen_reference import eigen_fixtures
from general_eigen_reference import general_eigen_fixtures
from solve_reference import solve_fixtures
from vector_reference import vector_fixtures
from publication import PublicSanitizer
from wasm_publication import sanitize_live_sources
from report_art import load_report_art

ROOT = Path(__file__).resolve().parents[1]


def live_bundle(data, wasm_directory=None, include_live=True, sanitizer=None, trace_directory=None):
    """Embed a separately identified executable bundle without changing recorded results."""
    if not include_live:
        return {"available": False, "reason": "Live testing was disabled when this report was generated."}
    if not any(item.get("id") == "wasm" and item.get("status") == "passed"
               for item in data.get("implementations", [])):
        return {"available": False, "reason": "This report has no verified WebAssembly implementation. Run the harness with WebAssembly enabled to add live tests."}
    directory = Path(wasm_directory) if wasm_directory is not None else ROOT / ".build/wasm"
    try:
        binary = (directory / "matrix.wasm").read_bytes()
        if not binary.startswith(b"\0asm\x01\0\0\0"):
            return {"available": False, "reason": "The WebAssembly build is invalid. Rebuild the WASM port and regenerate this report."}
        sources = {
            "module_source": (directory / "matrix.mjs").read_text(encoding="utf-8"),
            "wrapper_source": (ROOT / "ports/wasm/matrix.mjs").read_text(encoding="utf-8"),
            "worker_source": (ROOT / "benchmarks/live-worker.mjs").read_text(encoding="utf-8"),
            "accuracy_worker_source": (ROOT / "benchmarks/accuracy-worker.mjs").read_text(encoding="utf-8"),
        }
    except (OSError, UnicodeError):
        return {"available": False, "reason": "The WebAssembly build is unavailable. Build the WASM port and regenerate this report to enable live tests."}
    sanitizer = sanitizer or PublicSanitizer(root=ROOT)
    geometry = ROOT / "benchmarks/geometry-worker.mjs"
    try:
        sources["geometry_worker_source"] = geometry.read_text(encoding="utf-8")
    except (OSError, UnicodeError):
        pass
    try:
        sources, binary = sanitize_live_sources(sources, binary, sanitizer)
    except ValueError:
        return {"available": False, "reason": "The WebAssembly bundle could not be prepared for publication. Rebuild the WASM port and regenerate this report."}
    digest = hashlib.sha256(binary)
    for name, source in sorted(sources.items()):
        digest.update(name.encode("utf-8"))
        digest.update(source.encode("utf-8"))
    extra = {}
    if "geometry_worker_source" in sources:
        trace = Path(trace_directory) if trace_directory is not None else directory.parent / "wasm-trace"
        try:
            trace_sources, trace_binary = sanitize_live_sources(
                {"trace_module_source": (trace / "matrix.mjs").read_text(encoding="utf-8")},
                (trace / "matrix.wasm").read_bytes(), sanitizer)
            extra.update(trace_sources, trace_wasm_base64=base64.b64encode(trace_binary).decode("ascii"))
            digest.update(b"trace\0" + trace_binary + trace_sources["trace_module_source"].encode("utf-8"))
        except (OSError, UnicodeError, ValueError):
            extra["trace_reason"] = "Build the optional trace bundle to inspect captured calculation steps."
        try:
            source = json.loads((trace / "trace-source.json").read_text(encoding="utf-8"))
            if ("trace_wasm_base64" not in extra or source.get("path") != "ports/c/matrix.c"
                    or type(source.get("start_line")) is not int or source["start_line"] < 1
                    or not isinstance(source.get("lines"), list) or not 1 <= len(source["lines"]) <= 200
                    or not all(isinstance(line, str) for line in source["lines"])
                    or not re.fullmatch(r"[0-9a-f]{64}", source.get("source_sha256", ""))):
                raise ValueError("invalid trace source manifest")
            extra["geometry_source"] = {"path": "ports/c/matrix.c", "start_line": source["start_line"],
                                        "lines": [sanitizer.text(line) for line in source["lines"]],
                                        "source_sha256": source["source_sha256"]}
            digest.update(json.dumps(extra["geometry_source"], sort_keys=True).encode("utf-8"))
        except (OSError, UnicodeError, ValueError, TypeError, AttributeError):
            pass
    cases = json.loads(json.dumps(arithmetic_fixtures() + eigen_fixtures() + general_eigen_fixtures() + vector_fixtures() + solve_fixtures(), default=float, allow_nan=False))
    digest.update(b"fixtures\0")
    digest.update(json.dumps(cases, sort_keys=True, allow_nan=False).encode("utf-8"))
    return {"available": True, **sources, **extra, "wasm_base64": base64.b64encode(binary).decode("ascii"),
            "fixtures": cases, "sha256": digest.hexdigest(), "byte_length": len(binary)}


def json_for_html(value):
    return json.dumps(value, allow_nan=False).replace("<", "\\u003c").replace("&", "\\u0026")


def render(data, destination, wasm_directory=None, include_live=True, trace_directory=None):
    sanitizer = PublicSanitizer(root=ROOT)
    data = sanitizer.report(data)
    template = Path(__file__).with_name("report.html").read_text()
    live = live_bundle(data, wasm_directory, include_live, sanitizer, trace_directory)
    script = (ROOT / "benchmarks/live-report.mjs").read_text(encoding="utf-8")
    geometry_script = ROOT / "benchmarks/geometry-report.mjs"
    # Replace template markers in one pass so embedded source cannot introduce a
    # marker that would be interpreted as another template instruction.
    replacements = {"__REPORT_DATA__": json_for_html(data), "__LIVE_DATA__": json_for_html(live),
                    "__REPORT_ART__": load_report_art(),
                    "__TIMING_SCRIPT__": Path(__file__).with_name("timing-report.js").read_text(encoding="utf-8"),
                    "__THEME_SCRIPT__": Path(__file__).with_name("report-theme.js").read_text(encoding="utf-8"),
                    "__LIVE_SCRIPT__": script,
                    "__ACCURACY_SCRIPT__": (ROOT / "benchmarks/accuracy-report.mjs").read_text(encoding="utf-8"),
                    "__GEOMETRY_SCRIPT__": geometry_script.read_text(encoding="utf-8") if geometry_script.is_file() else ""}
    html = re.sub("|".join(replacements), lambda match: replacements[match[0]], template)
    Path(destination).write_text(html, encoding="utf-8")


def publish(data, destination, source_directory=None, wasm_directory=None, include_live=True, trace_directory=None):
    """Publish HTML, JSON, and raw profiles without running any measurements.

    The destination directory receives a public copy of the results and profiles.
    When source and destination match, existing artifacts are sanitized in place.
    """
    from profile_publication import sanitize_existing_profiles
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    sanitizer = PublicSanitizer(root=ROOT)
    sanitize_existing_profiles(data, source_directory or destination.parent, sanitizer,
                               destination_dir=destination.parent)
    public = sanitizer.report(data)
    render(public, destination, wasm_directory, include_live, trace_directory)
    (destination.parent / "results.json").write_text(
        json.dumps(public, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return public


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Existing results.json; recorded timings are retained")
    parser.add_argument("output", type=Path, help="Destination HTML file; sanitized results.json and profiles are written alongside it")
    parser.add_argument("--wasm-dir", type=Path, help="Directory containing matrix.mjs and matrix.wasm")
    parser.add_argument("--trace-wasm-dir", type=Path, help="Optional instrumented WASM build for calculation steps")
    parser.add_argument("--without-live", action="store_true", help="Omit the executable WebAssembly bundle")
    args = parser.parse_args()
    publish(json.loads(args.input.read_text()), args.output, args.input.parent,
            args.wasm_dir, not args.without_live, args.trace_wasm_dir)


if __name__ == "__main__":
    main()
