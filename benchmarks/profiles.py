"""Real stack profiles from separate, bounded invocations of known runners.

Sampling never alters ordinary timing runs. Native sampling attaches only to the
process this module launches; it never discovers or attaches to user processes.
"""
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess

from profile_publication import (
    sanitize_go_raw, sanitize_julia, sanitize_macos_sample, sanitize_node,
    sanitize_python, sanitized_profile_text,
)

ROOT = Path(__file__).resolve().parents[1]
SIZE = 48
TARGET_NS = 3_000_000_000


def _portable_paths(text):
    return text.replace(str(ROOT) + "/", "").replace(str(Path.home()) + "/", "~/")


def parse_node(profile):
    nodes = {node["id"]: node for node in profile.get("nodes", [])}
    parents = {}
    for node in nodes.values():
        for child in node.get("children", []):
            parents[child] = node["id"]

    def stack(node_id):
        result, seen = [], set()
        while node_id in nodes:
            if node_id in seen:
                raise ValueError("cyclic Node profile stack")
            seen.add(node_id)
            frame = nodes[node_id].get("callFrame", {})
            name = frame.get("functionName") or "(anonymous)"
            if frame.get("url"):
                name += " (%s:%s)" % (Path(frame["url"]).name, frame.get("lineNumber", 0) + 1)
            result.append(name)
            node_id = parents.get(node_id)
        return result[::-1]

    samples, deltas = profile.get("samples", []), profile.get("timeDeltas", [])
    if len(samples) != len(deltas):
        raise ValueError("Node profile sample/delta counts disagree")
    stacks, weights = [], []
    for sample, delta in zip(samples, deltas):
        frames = stack(sample)
        if frames and delta > 0:
            stacks.append(frames)
            weights.append(delta / 1000.0)  # V8 timeDeltas are microseconds.
    return stacks, weights


def parse_python(profile):
    samples = profile.get("samples", [])
    stacks, weights = [], []
    for index, sample in enumerate(samples):
        end = samples[index + 1]["at_ns"] if index + 1 < len(samples) else profile["end_ns"]
        delta = end - sample["at_ns"]
        if sample["stack"] and delta > 0:
            stacks.append(sample["stack"])
            weights.append(delta / 1_000_000.0)
    return stacks, weights


def parse_macos_sample(text):
    """Flatten inclusive call-tree counts into non-overlapping stack weights."""
    interval = re.search(r"every ([\d.]+) millisecond", text)
    if not interval:
        raise ValueError("sample report has no sampling interval")
    interval_ms = float(interval.group(1))
    in_graph, selected, root = False, False, None
    ancestors = []
    for line in text.splitlines():
        if line.strip() == "Call graph:":
            in_graph = True
            continue
        if not in_graph:
            continue
        if line.startswith("Total number") or line.startswith("Binary Images:"):
            break
        match = re.match(r"^([ |+!:]*)(\d+) (.+)$", line)
        if not match:
            continue
        depth, count, name = len(match[1]), int(match[2]), match[3]
        if re.match(r"Thread[_ ]", name):
            if selected:
                break
            selected = "com.apple.main-thread" in name
            if selected:
                root = {"name": "main thread", "count": count, "children": []}
                ancestors = [(depth, root)]
            continue
        if not selected:
            continue
        # Retain actual symbol + image, dropping instruction addresses/offsets.
        name = re.sub(r"\s+\+\s+[\d,\.]+.*$", "", name)
        name = re.sub(r"\s+\[0x[\da-fA-F,\.]+\].*$", "", name)
        node = {"name": name, "count": count, "children": []}
        while ancestors and depth <= ancestors[-1][0]:
            ancestors.pop()
        if not ancestors:
            raise ValueError("invalid sample call-tree indentation")
        ancestors[-1][1]["children"].append(node)
        ancestors.append((depth, node))
    if root is None:
        raise ValueError("sample did not identify the main thread")
    stacks, weights = [], []

    def visit(node, prefix):
        path = prefix + [node["name"]]
        residual = node["count"] - sum(child["count"] for child in node["children"])
        if residual < 0:
            raise ValueError("sample child counts exceed their parent")
        if residual:
            stacks.append(path)
            weights.append(residual * interval_ms)
        for child in node["children"]:
            visit(child, path)

    visit(root, [])
    return stacks, weights


def parse_go_raw(text):
    section, locations, samples = None, {}, []
    current = None
    for line in text.splitlines():
        if line.startswith("Samples:"):
            section = "samples"
            continue
        if line.startswith("Locations"):
            section = "locations"
            continue
        if line.startswith("Mappings"):
            section = None
        if section == "samples":
            match = re.match(r"^\s*(\d+)\s+(\d+):\s*([\d ]+)\s*$", line)
            if match:
                samples.append((int(match[2]) / 1_000_000.0, [int(value) for value in match[3].split()]))
        elif section == "locations":
            match = re.match(r"^\s*(\d+):\s+(0x[\da-fA-F]+)(?:\s+M=\d+)?\s*(.*)$", line)
            if match:
                current = int(match[1])
                description = match[3].strip()
                locations[current] = [description or match[2]]
            elif current is not None and line.strip():
                locations[current].append(line.strip())
    stacks, weights = [], []
    for weight, ids in samples:
        if weight <= 0:
            continue
        path = []
        for location in reversed(ids):
            if location not in locations:
                raise ValueError("pprof sample refers to an unknown location")
            path.extend(_portable_paths(frame) for frame in reversed(locations[location]))
        if path:
            stacks.append(path)
            weights.append(weight)
    return stacks, weights


def _command(command, env, timeout=45):
    result = subprocess.run([str(part) for part in command], cwd=ROOT, env=env,
                            capture_output=True, text=True, timeout=timeout)
    if result.returncode:
        raise RuntimeError((result.stderr or result.stdout or "profiling command failed")[-2500:])
    return result


def _check_output(output, iterations, checksum):
    result = json.loads(output)
    value = result.get("checksum")
    if type(result.get("iterations")) is not int or result["iterations"] != iterations or isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("invalid profiling runner output")
    if not math.isclose(value, checksum * iterations, rel_tol=1e-7, abs_tol=1e-8 * iterations):
        raise ValueError("profiling checksum differs from the independent reference")


def _calibrate(implementation, seed, operation="multiply", size=SIZE):
    from run import measure
    from reference import benchmark_checksum
    expected = benchmark_checksum(operation, size, seed)
    trial = measure(implementation, operation, size, 1, seed, expected, timeout=90)
    count = max(1, min(1000, math.ceil(50_000_000 / trial["ns_per_op"])))
    trial = measure(implementation, operation, size, count, seed, expected, timeout=90)
    iterations = max(1, min(500_000_000, math.ceil(TARGET_NS / trial["ns_per_op"])))
    return iterations, expected


def _native(implementation, arguments, raw, iterations, expected):
    sampler = "/usr/bin/sample"
    if platform.system() != "Darwin" or not Path(sampler).is_file():
        raise FileNotFoundError("Native stack sampling is currently supported only by macOS /usr/bin/sample.")
    process = subprocess.Popen(implementation["runner"] + arguments, cwd=ROOT, env=implementation["env"],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        if process.poll() is not None:
            raise RuntimeError("profiled runner exited before sampling began")
        result = subprocess.run([sampler, str(process.pid), "2", "1", "-file", str(raw)],
                                text=True, capture_output=True, timeout=15)
        if result.returncode:
            raise FileNotFoundError("macOS sampling unavailable: " + (result.stderr or result.stdout)[-1500:])
        stdout, stderr = process.communicate(timeout=35)
        if process.returncode:
            raise RuntimeError(stderr[-2500:])
        _check_output(stdout, iterations, expected)
    finally:
        if process.poll() is None:
            process.terminate()
            try:
                process.communicate(timeout=3)
            except subprocess.TimeoutExpired:
                process.kill()
                process.communicate(timeout=3)
        if raw.exists():
            try:
                _sanitize_native_file(raw)
            except ValueError:
                # A partial/unsupported report can contain irrelevant machine
                # headers; do not retain it even when profiling failed.
                raw.unlink()
    if not raw.exists():
        raise ValueError("sample output has no supported call graph")
    return parse_macos_sample(raw.read_text(encoding="utf-8"))


def _sanitize_native_file(raw):
    from publication import PublicSanitizer
    text = sanitize_macos_sample(raw.read_text(encoding="utf-8", errors="replace"), PublicSanitizer())
    raw.write_text(text, encoding="utf-8")


def collect_profiles(implementations, output_dir, seed, operation="multiply", size=SIZE):
    from publication import PublicSanitizer
    sanitizer = PublicSanitizer()
    output_dir = Path(output_dir).resolve()
    results = []
    for implementation in implementations:
        name = implementation["id"]
        profile = {"implementation": name, "status": "failed", "chronological": False, "stacks": [], "weights": []}
        results.append(profile)
        directory = output_dir / "profiles" / name
        raw = None
        try:
            if name not in ("csharp", "fsharp", "rust", "go", "typescript", "python", "cpp", "c", "assembly", "julia", "wasm"):
                raise ValueError("unknown profiling implementation")
            if (output_dir / "profiles").is_symlink() or directory.is_symlink():
                raise ValueError("raw profile directories must not be symlinks")
            directory.mkdir(parents=True, exist_ok=True)
            if any(path.is_symlink() or (path.is_file() and path.stat().st_nlink != 1) for path in directory.iterdir()):
                raise ValueError("raw profile files must not be links")
            if name not in ("python", "typescript", "wasm", "go", "julia") and platform.system() != "Darwin":
                raise FileNotFoundError("No native stack profiler is configured on this platform; benchmark timings remain available.")
            iterations, expected = _calibrate(implementation, seed, operation, size)
            arguments = ["bench", operation, str(size), str(iterations), str(seed)]
            runner, env = implementation["runner"], implementation["env"].copy()
            if name == "python":
                if len(runner) != 2:
                    raise FileNotFoundError("Python profiling requires the standard interpreter + script runner.")
                raw = directory / "samples.json"
                completed = _command([runner[0], ROOT / "benchmarks/profile_python.py", raw, runner[1]] + arguments, env)
                _check_output(completed.stdout, iterations, expected)
                data = sanitize_python(json.loads(raw.read_text()), sanitizer)
                raw.write_text(json.dumps(data, allow_nan=False), encoding="utf-8")
                stacks, weights = parse_python(data)
                profile.update(chronological=True, note="Python main-thread stacks, sampled by a Python thread; weights use actual sample intervals. GIL scheduling biases sampling. Includes startup and warmup in a separate run.")
            elif name in ("typescript", "wasm"):
                raw = directory / "profile.cpuprofile"
                completed = _command([runner[0], "--cpu-prof", "--cpu-prof-interval=1000", "--cpu-prof-dir=" + str(directory), "--cpu-prof-name=" + raw.name] + runner[1:] + arguments, env)
                _check_output(completed.stdout, iterations, expected)
                data = sanitize_node(json.loads(raw.read_text()), sanitizer)
                raw.write_text(json.dumps(data, allow_nan=False), encoding="utf-8")
                stacks, weights = parse_node(data)
                profile.update(chronological=True, note="V8 CPU profiler samples with recorded time deltas; chronological sample timeline. Includes Node startup, warmup, and measured workload in a separate run.")
                if name == "wasm":
                    profile["note"] += " Includes JavaScript wrapper and WebAssembly frames; optimized/inlined WASM functions may be merged or omitted."
            elif name == "go":
                raw = directory / "cpu.pprof"
                env["LINEAR_A_CPU_PROFILE"] = str(raw)
                completed = _command(runner + arguments, env)
                _check_output(completed.stdout, iterations, expected)
                go = env.get("GO") or shutil.which("go", path=env.get("PATH"))
                if not go:
                    raise FileNotFoundError("Go toolchain is needed to decode runtime/pprof CPU samples.")
                decoded = _command([go, "tool", "pprof", "-raw", raw], env)
                # Retain decoded real samples, without machine-specific source or
                # executable paths embedded in the binary protobuf profile.
                raw.unlink()
                raw = directory / "cpu.raw.txt"
                text = sanitize_go_raw(decoded.stdout, sanitizer)
                raw.write_text(text, encoding="utf-8")
                stacks, weights = parse_go_raw(text)
                profile["note"] = "Go runtime/pprof CPU samples across goroutines; aggregate CPU-time weights in milliseconds, without chronology. Includes warmup and workload in a separate run."
            elif name == "julia":
                raw = directory / "samples.json"
                wrapper = ROOT / "ports/julia/profile.jl"
                if not wrapper.is_file():
                    raise FileNotFoundError("Julia Profile wrapper is unavailable.")
                completed = _command(runner[:-1] + [str(wrapper), operation, str(size), str(iterations), str(seed), str(raw)], env, timeout=90)
                # Julia wrapper profiles the same workload after its own JIT warmup.
                data = sanitize_julia(json.loads(raw.read_text()), sanitizer)
                raw.write_text(json.dumps(data, allow_nan=False), encoding="utf-8")
                _check_output(json.dumps({"iterations": iterations, "checksum": data.get("checksum")}), iterations, expected)
                stacks, weights = data["stacks"], data["weights"]
                profile["note"] = data.get("note", "Julia Profile stack samples; aggregate nominal sample-interval weights, without chronology. Collected separately after JIT warmup.")
            else:
                raw = directory / "sample.txt"
                stacks, weights = _native(implementation, arguments, raw, iterations, expected)
                profile["note"] = "macOS main-thread stack samples: aggregate counts × nominal 1 ms sampling interval, without chronology. Native/JIT symbols may be incomplete; includes allocation and runtime work."
                if name in ("csharp", "fsharp"):
                    profile["note"] += " This is an external native-process profile; managed .NET method names may remain unresolved, so it is not a managed method-level profile."
            if not stacks or len(stacks) != len(weights) or any(isinstance(weight, bool) or not isinstance(weight, (int, float)) or not math.isfinite(weight) or weight <= 0 for weight in weights):
                raise ValueError("profiler produced no usable stack samples")
            profile.update(status="available", stacks=stacks, weights=weights, raw_file=str(raw.relative_to(output_dir)),
                           workload={"operation": operation, "size": size, "iterations": iterations, "seed": seed})
        except FileNotFoundError as error:
            profile.update(status="unavailable", note=str(error))
        except Exception as error:
            profile.update(status="failed", note=str(error)[-2500:])
        finally:
            # Profilers can write files before exiting with an error. Retain only
            # supported, sanitized text; never leave the Go protobuf or a partial
            # private artifact behind when a collection fails.
            if raw is not None and raw.is_file() and not raw.is_symlink():
                try:
                    clean = sanitized_profile_text(raw.name, raw.read_text(encoding="utf-8"), sanitizer)
                    raw.write_text(clean, encoding="utf-8")
                except (ValueError, KeyError, TypeError, UnicodeError):
                    raw.unlink()
    return results
