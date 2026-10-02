"""Publish raw stack samples without machine metadata or private source paths."""
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat


def _number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("profile contains an invalid numeric measurement")
    return value


def _numbers(values):
    if not isinstance(values, list):
        raise ValueError("profile measurements must be an array")
    return [_number(value) for value in values]


def _label(value, sanitizer):
    if not isinstance(value, str):
        raise ValueError("profile source labels must be strings")
    return sanitizer.text(value)


def sanitize_node(profile, sanitizer):
    """Retain the V8 CPU-profile fields needed by consumers and sample parsers."""
    result = {}
    for key in ("startTime", "endTime"):
        if key in profile:
            result[key] = _number(profile[key])
    for key in ("samples", "timeDeltas"):
        if key in profile:
            result[key] = _numbers(profile[key])
    nodes = []
    for node in profile.get("nodes", []):
        clean = {key: _number(node[key]) for key in ("id", "hitCount") if key in node}
        if "children" in node:
            clean["children"] = _numbers(node["children"])
        frame = node.get("callFrame", {})
        clean["callFrame"] = {
            key: _label(frame[key], sanitizer)
            for key in ("functionName", "scriptId", "url") if key in frame
        }
        clean["callFrame"].update({
            key: _number(frame[key]) for key in ("lineNumber", "columnNumber") if key in frame
        })
        if "positionTicks" in node:
            clean["positionTicks"] = [
                {key: _number(tick[key]) for key in ("line", "ticks") if key in tick}
                for tick in node["positionTicks"]
            ]
        nodes.append(clean)
    result["nodes"] = nodes
    return result


def sanitize_python(profile, sanitizer):
    return {"samples": [
        {"at_ns": _number(sample["at_ns"]),
         "stack": [_label(frame, sanitizer) for frame in sample["stack"]]}
        for sample in profile.get("samples", [])
    ], "end_ns": _number(profile["end_ns"])}


def sanitize_julia(profile, sanitizer):
    result = {key: _label(profile[key], sanitizer)
              for key in ("implementation", "kind", "unit", "note") if key in profile}
    if "chronological" in profile:
        if type(profile["chronological"]) is not bool:
            raise ValueError("profile chronology must be a boolean")
        result["chronological"] = profile["chronological"]
    if "checksum" in profile:
        result["checksum"] = _number(profile["checksum"])
    result["stacks"] = [[_label(frame, sanitizer) for frame in stack] for stack in profile["stacks"]]
    result["weights"] = _numbers(profile["weights"])
    return result


def sanitize_go_raw(text, sanitizer):
    """Keep pprof's sample and location sections; discard labels and image inventory."""
    output, section, location = [], None, False
    for line in text.splitlines():
        if re.fullmatch(r"PeriodType: [A-Za-z]+ [A-Za-z]+", line):
            output.append(line)
        elif re.fullmatch(r"(?:Period|Duration): [\d.eE+\-]+", line):
            output.append(line)
        elif re.fullmatch(r"Time: \d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}(?:\.\d+)?(?: [+-]\d{4})?(?: [A-Za-z]+)?", line):
            output.append(line)
        elif line == "Samples:":
            section = "samples"
            output.append(line)
        elif line.startswith("Locations"):
            section = "locations"
            output.append("Locations")
        elif line.startswith("Mappings"):
            section = None
            output.append("Mappings")
        elif section == "samples":
            if re.fullmatch(r"\s*samples/count\s+cpu/nanoseconds\s*", line) or re.fullmatch(r"\s*\d+\s+\d+:\s*[\d ]+\s*", line):
                output.append(line)
        elif section == "locations":
            match = re.match(r"^(\s*\d+:\s+0x[\da-fA-F]+(?:\s+M=\d+)?\s*)(.*)$", line)
            if match:
                location = True
                output.append(match[1] + sanitizer.text(match[2]))
            elif location and line.strip():
                # pprof prints inlined functions on continuation lines.
                output.append(line[:len(line) - len(line.lstrip())] + sanitizer.text(line.strip()))
    if "Samples:" not in output or "Locations" not in output:
        raise ValueError("unsupported Go raw profile")
    return "\n".join(output) + "\n"


def sanitize_macos_sample(text, sanitizer):
    """Retain counted call trees and interval, with anonymous thread identifiers."""
    interval = re.search(r"^Analysis of sampling .* every ([\d.]+) millisecond.*$", text, re.MULTILINE)
    graph = re.search(r"(?ms)^Call graph:\n(.*?)(?=^Total number|^Binary Images:|\Z)", text)
    if not interval or not graph:
        raise ValueError("sample output has no supported call graph")
    output = ["Analysis of sampling runner every %s millisecond" % interval[1], "", "Call graph:"]
    thread = 0
    for line in graph[1].splitlines():
        match = re.match(r"^([ |+!:]*)(\d+) (.+)$", line)
        if not match:
            continue
        prefix, count, name = match.groups()
        if re.match(r"Thread[_ ]", name):
            thread += 1
            queue = " DispatchQueue_1: com.apple.main-thread (serial)" if "com.apple.main-thread" in name else " worker"
            name = "Thread_%s%s" % (thread, queue)
        else:
            name = sanitizer.text(name)
        output.append(prefix + count + " " + name)
    return "\n".join(output) + "\n"


def sanitized_profile_text(relative, text, sanitizer):
    name = PurePosixPath(relative).name
    if name.endswith(".cpuprofile"):
        result = sanitize_node(json.loads(text), sanitizer)
    elif name == "samples.json":
        data = json.loads(text)
        result = sanitize_python(data, sanitizer) if "samples" in data else sanitize_julia(data, sanitizer)
    elif name == "cpu.raw.txt":
        return sanitize_go_raw(text, sanitizer)
    elif name == "sample.txt":
        return sanitize_macos_sample(text, sanitizer)
    else:
        raise ValueError("unsupported raw profile format")
    return json.dumps(result, allow_nan=False, separators=(",", ":")) + "\n"


def _relative(value):
    if not isinstance(value, str) or "\\" in value:
        raise ValueError("raw profile path must be a relative POSIX path")
    path = PurePosixPath(value)
    if path.is_absolute() or not path.parts or path.parts[0] != "profiles" or any(part in (".", "..") for part in value.split("/")):
        raise ValueError("raw profile path must remain within profiles/")
    return path


def _root(directory, create=False):
    path = Path(directory).absolute()
    if path.is_symlink():
        raise ValueError("profile report directory must not be a symlink")
    if create:
        path.mkdir(parents=True, exist_ok=True)
    return path.resolve(strict=True)


def _file(root, relative, write=None):
    """Open each component relative to an anchored directory, without following links."""
    parts = _relative(str(relative)).parts
    directory = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        for part in parts[:-1]:
            if write is not None:
                try:
                    os.mkdir(part, dir_fd=directory)
                except FileExistsError:
                    pass
            child = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=directory)
            os.close(directory)
            directory = child
        flags = os.O_NOFOLLOW | (os.O_WRONLY | os.O_CREAT if write is not None else os.O_RDONLY)
        descriptor = os.open(parts[-1], flags, 0o644, dir_fd=directory)
        with os.fdopen(descriptor, "w" if write is not None else "r", encoding="utf-8") as stream:
            info = os.fstat(stream.fileno())
            if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                raise ValueError("raw profiles must be regular files without hard links")
            if write is not None:
                stream.truncate(0)
                stream.write(write)
            else:
                return stream.read()
    finally:
        os.close(directory)


def _profile_paths(root):
    paths = set()
    tree = root / "profiles"
    if tree.is_symlink():
        raise ValueError("raw profile directories must not be symlinks")
    if tree.exists():
        for directory, subdirs, files in os.walk(tree, followlinks=False):
            for name in subdirs + files:
                path = Path(directory) / name
                if path.is_symlink():
                    raise ValueError("raw profile paths must not be symlinks")
            for name in files:
                paths.add((Path(directory) / name).relative_to(root).as_posix())
    return paths


def sanitize_existing_profiles(report, report_dir, sanitizer, destination_dir=None):
    """Rewrite/copy raw profiles, returning paths; recorded measurements are untouched.

    All files under profiles/ in both directories are included, even ignored or
    unreferenced profiles. Source files take precedence over destination copies.
    Unsafe paths, symlinks, missing references, and unsupported formats fail closed.
    When destination_dir differs, source bytes are left unchanged.
    """
    source = _root(report_dir)
    destination = _root(destination_dir, create=True) if destination_dir is not None else source
    paths = _profile_paths(source)
    for profile in report.get("profiles", []):
        if "raw_file" in profile:
            paths.add(str(_relative(profile["raw_file"])))
    # Validate and sanitize everything before replacing any existing files.
    prepared = {}
    if destination != source:
        for relative in sorted(_profile_paths(destination)):
            prepared[relative] = sanitized_profile_text(relative, _file(destination, relative), sanitizer)
    for relative in sorted(paths):
        prepared[relative] = sanitized_profile_text(relative, _file(source, relative), sanitizer)
    for relative in sorted(prepared):
        _file(destination, relative, write=prepared[relative])
    return sorted(prepared)
