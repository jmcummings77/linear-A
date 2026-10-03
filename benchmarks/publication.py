"""Public report metadata: retain measurements and source context, remove host details.

This module sanitizes human-readable metadata. It must not be applied to whole
JavaScript programs, WebAssembly binaries, or other executable content.
"""
import copy
import os
from pathlib import Path
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]

_PATH_START = re.compile(
    r"(?<![\w./:\\<])(?:file://|[A-Za-z]:[\\/]|~[\\/]|\\\\[^\\\s]+[\\/]|/(?!/)[A-Za-z0-9_.-])",
    re.IGNORECASE,
)
_SOURCE_END = re.compile(
    r"\.(?:c|cc|h|cpp|hpp|cs|fs|fsx|rs|go|py|pyc|jl|ts|js|mjs|json|toml|wasm|"
    r"so(?:\.\d+)*|dll|dylib|exe|csproj|fsproj|slnx?|log|txt|lock|zip|gz)"
    r"(?::\d+){0,3}(?=\s|$)", re.IGNORECASE,
)
_EXECUTABLE_END = re.compile(
    r"[\\/](?:dotnet|cargo|rustc|go|node|npm|python(?:\d+(?:\.\d+)*)?|emcc|"
    r"clang(?:\+\+)?(?:-\d+)?|gcc|g\+\+|julia|runner)(?:\.exe)?(?=\s|$)",
    re.IGNORECASE,
)
_SECRET_FLAG = re.compile(r"^--?(?:[\w-]*(?:token|secret|password|credential)|api[-_]key|access[-_]key|authorization)$", re.IGNORECASE)
_UUID = re.compile(r"\b[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}\b", re.IGNORECASE)


def _decode(value):
    # URI-encoded file names occasionally arrive through both JSON and a URL.
    for _ in range(2):
        decoded = unquote(value)
        if decoded == value:
            break
        value = decoded
    return value


def _normalized(value):
    value = _decode(str(value)).replace("\\", "/")
    if value.lower().startswith("file://"):
        value = value[7:]
        if value.lower().startswith("localhost/"):
            value = value[9:]
        if not value.startswith("/") and not re.match(r"^[A-Za-z]:/", value):
            value = "//" + value
    if re.match(r"^/[A-Za-z]:/", value):
        value = value[1:]
    return value.rstrip("/")


class PublicSanitizer:
    def __init__(self, root=ROOT, home=None):
        self.root = _normalized(root)
        self.home = _normalized(Path.home() if home is None else home)
        self._usernames = set()
        for path in (self.home, self.root):
            match = re.search(r"(?:^|/)(?:Users|home)/([^/]+)", path, re.IGNORECASE)
            if match:
                self._usernames.add(match[1])

    def is_local_path(self, value):
        """Recognize host paths; arbitrary browser URLs like /assets/x are excluded."""
        if not isinstance(value, str):
            return False
        decoded = _decode(value)
        path = _normalized(decoded)
        if decoded.lower().startswith("file://") or re.match(r"^(?:[A-Za-z]:[\\/]|\\\\|~[/\\])", decoded):
            return True
        for prefix in (self.root, self.home):
            if prefix and (path.casefold() == prefix.casefold() or path.casefold().startswith(prefix.casefold() + "/")):
                return True
        return bool(re.match(r"^/(?:Users|home|private|tmp|var|usr|opt|Applications|Library|System|nix|build|workspace|workspaces|root)(?:/|$)", path, re.IGNORECASE))

    def _path(self, value):
        path = _normalized(value)
        suffix = ""
        location = re.search(r"(?::\d+){1,3}$", path)
        if location:
            suffix, path = location[0], path[:location.start()]
        if path.casefold() == self.root.casefold():
            return "." + suffix
        if path.casefold().startswith(self.root.casefold() + "/"):
            relative = path[len(self.root) + 1:]
            if not re.search(r"(?:^|/)(?:\.codex|\.chatgpt)(?:/|$)", relative, re.IGNORECASE):
                return relative + suffix
        # Preserve useful public library source context, without retaining the
        # host checkout, package-cache path, username, or temporary directory.
        for marker in ("/rustlib/src/rust/library/", "/node_modules/"):
            if marker in path:
                tail = path.split(marker, 1)[1]
                prefix = "library/" if "rustlib" in marker else "node_modules/"
                return prefix + tail + suffix
        if "/src/" in path:
            return "src/" + path.rsplit("/src/", 1)[1] + suffix
        basename = path.rstrip("/").rsplit("/", 1)[-1]
        if not basename or basename in self._usernames or basename.lower() in ("tmp", "private", "home", "users", ".codex", ".chatgpt"):
            basename = "[local path]"
        return basename + suffix

    @staticmethod
    def _redact_credentials(value):
        value = re.sub(r"-----BEGIN [^-]*PRIVATE KEY-----.*?-----END [^-]*PRIVATE KEY-----", "[redacted credential]", value, flags=re.DOTALL)
        value = re.sub(r"(?i)\b(?:authorization\s*[:=]\s*)?(?:bearer|basic)\s+[A-Za-z0-9+/=._~-]+", "[redacted credential]", value)
        value = re.sub(r"(?i)\b(https?://)[^/@\s]+@", r"\1", value)
        value = re.sub(
            r"(?i)(?:--?[\w-]*(?:token|secret|password|credential)|--?api[-_]key)(?:=|\s+)"
            r"(?:\"[^\"]*\"|'[^']*'|[^\s,;]+)", "[redacted credential]", value)
        value = re.sub(
            r"(?i)\b(?:[\w-]*(?:token|secret|password|credential|api[_-]?key|access[_-]?key)|authorization|cookie)"
            r"\s*[:=]\s*(?:\"[^\"]*\"|'[^']*'|[^\s,;&]+)", "[redacted credential]", value)
        value = re.sub(r"\b(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|sk-(?:proj-)?[A-Za-z0-9_-]{16,}|AKIA[A-Z0-9]{16}|AIza[A-Za-z0-9_-]{30,})\b", "[redacted credential]", value)
        value = re.sub(r"\beyJ[A-Za-z0-9_-]+\.[A-Za-z0-9_-]+\.[A-Za-z0-9_-]+\b", "[redacted credential]", value)
        value = re.sub(r"(?<![\w-])(?:export\s+)?[A-Z][A-Z0-9_]{2,}\s*=\s*(?:\"[^\"]*\"|'[^']*'|[^\s,;]+)", "[redacted environment]", value)
        value = re.sub(r"(?i)\b(?:host|hostname|computer[_ -]?name|machine[_ -]?(?:name|id)|serial[_ -]?(?:number|id)|username|user|agent[_ -]?(?:id|name))\s*[:=]\s*(?:\"[^\"]*\"|'[^']*'|[^\s,;]+)", "[redacted identity]", value)
        return value

    def text(self, value):
        """Sanitize a human-readable diagnostic, toolchain string, or stack label."""
        if not isinstance(value, str):
            raise TypeError("text must be a string")
        value = self._redact_credentials(_decode(value))
        pieces, cursor = [], 0
        while True:
            match = _PATH_START.search(value, cursor)
            if not match:
                pieces.append(value[cursor:])
                break
            start = match.start()
            end = len(value)
            delimiter = re.search(r"[\r\n\t\"'<>|(){}\[\],;]", value[match.end():])
            if delimiter:
                end = match.end() + delimiter.start()
            next_path = _PATH_START.search(value, match.end())
            if next_path and next_path.start() < end:
                end = next_path.start()
            candidate = value[start:end]
            source_end, executable_end = _SOURCE_END.search(candidate), _EXECUTABLE_END.search(candidate)
            endings = [item.end() for item in (source_end, executable_end) if item]
            if endings:
                end = start + min(endings)
            while end > start and value[end - 1].isspace():
                end -= 1
            pieces.append(value[cursor:start])
            pieces.append(self._path(value[start:end]))
            cursor = end
        value = "".join(pieces)
        value = _UUID.sub("[redacted id]", value)
        value = re.sub(r"(?i)\b(?:generated|created|written|authored|assisted|produced)\s+(?:with|by)\s+(?:an?\s+)?(?:ai(?: agent)?|codex|chatgpt|openai)(?:\s+agent)?\b", "[redacted attribution]", value)
        value = re.sub(r"(?i)\b(?:codex|chatgpt)\s+(?:bundled|workspace|agent)\b", "[local tooling]", value)
        def identities(segment):
            for username in self._usernames:
                segment = re.sub(r"\b" + re.escape(username) + r"@[^\s,;<>]+", "[redacted identity]", segment, flags=re.IGNORECASE)
                segment = re.sub(r"(?<![\w-])" + re.escape(username) + r"(?![\w-])", "[user]", segment, flags=re.IGNORECASE)
            return re.sub(r"(?i)\b(?:codex|chatgpt)\b", "[local tooling]", segment)

        # Public import/module identities are not local account names, even when
        # a developer deliberately used the same name for both accounts.
        public = re.compile(r"(?:https?://)?(?:github\.com|gitlab\.com|bitbucket\.org)/[^\s<>\"']+")
        pieces, cursor = [], 0
        for match in public.finditer(value):
            pieces.extend((identities(value[cursor:match.start()]), match[0]))
            cursor = match.end()
        pieces.append(identities(value[cursor:]))
        return "".join(pieces)

    @staticmethod
    def _tool(value):
        name = _normalized(value).rsplit("/", 1)[-1]
        name = re.sub(r"\.exe$", "", name, flags=re.IGNORECASE)
        if re.fullmatch(r"python(?:\d+(?:\.\d+)*)?", name):
            return "python3"
        if name in ("emcc.py", "em++", "em++.py"):
            return "emcc" if name == "emcc.py" else "em++"
        if name == "npm-cli.js":
            return "npm"
        if name == "nodejs":
            return "node"
        return re.sub(r"^(clang(?:\+\+)?|gcc|g\+\+)-\d+$", r"\1", name)

    def command(self, arguments):
        """Retain argv structure, replacing executable locations with public names."""
        if not isinstance(arguments, (list, tuple)):
            raise TypeError("command must be an argv sequence")
        result, secret_value = [], False
        for index, argument in enumerate(arguments):
            if not isinstance(argument, (str, os.PathLike, int, float)):
                continue
            value = str(argument)
            if secret_value:
                result.append("[redacted]")
                secret_value = False
            elif index == 0:
                result.append(self.text(self._tool(value)))
            elif _SECRET_FLAG.match(value):
                result.append(value)
                secret_value = True
            else:
                result.append(self.text(value))
        return result

    def _leaves(self, value, keys):
        if not isinstance(value, dict):
            return {}
        result = {}
        for key in keys:
            if key not in value:
                continue
            item = value[key]
            if isinstance(item, str):
                result[key] = self.text(item)
            elif item is None or type(item) in (int, float, bool):
                result[key] = item
        return result

    def _checks(self, checks):
        return [self._leaves(check, ("name", "passed", "error")) for check in checks if isinstance(check, dict)] if isinstance(checks, list) else []

    def report(self, data):
        """Return a schema-allowlisted copy; never mutate the original report."""
        if not isinstance(data, dict):
            raise TypeError("report must be an object")
        result = self._leaves(data, ("schema_version", "created_at", "revision", "dirty", "source_sha256", "suite", "seed", "total_seconds"))
        machine = data.get("machine", {})
        result["machine"] = {}
        if isinstance(machine, dict):
            for key in ("os", "release", "architecture"):
                value = machine.get(key)
                if isinstance(value, str) and re.fullmatch(r"[\w .+()-]{1,128}", value):
                    result["machine"][key] = self.text(value)
            for key, maximum in (("logical_cpus", 2**20), ("memory_bytes", 2**70)):
                value = machine.get(key)
                if type(value) is int and 0 < value <= maximum:
                    result["machine"][key] = value
            model = machine.get("cpu_model")
            if isinstance(model, str) and re.fullmatch(r"[\w .+()@-]{1,200}", model) and "\n" not in model:
                result["machine"]["cpu_model"] = self.text(model)
        result["methodology"] = self._leaves(data.get("methodology"), ("numeric_type", "matrix_layout", "timing", "warmup", "sampling", "profiles", "limitations"))
        result["implementations"] = []
        for item in data.get("implementations", []):
            if not isinstance(item, dict):
                continue
            cleaned = self._leaves(item, ("id", "name", "status", "toolchain", "detail", "error"))
            if isinstance(item.get("build_commands"), list):
                cleaned["build_commands"] = [self.command(command) for command in item["build_commands"] if isinstance(command, (list, tuple))]
            if "checks" in item:
                cleaned["checks"] = self._checks(item["checks"])
            result["implementations"].append(cleaned)
        result["results"] = []
        for item in data.get("results", []):
            if not isinstance(item, dict):
                continue
            cleaned = self._leaves(item, ("implementation", "operation", "size", "status", "iterations", "median_ns", "min_ns", "max_ns", "mad_ns", "error", "nnz", "unknowns", "logical_dense_bytes", "logical_csr_bytes", "cg_iterations", "contrast"))
            if isinstance(item.get("samples"), list):
                cleaned["samples"] = [self._leaves(sample, ("elapsed_ns", "iterations", "checksum", "ns_per_op")) for sample in item["samples"] if isinstance(sample, dict)]
            result["results"].append(cleaned)
        result["profiles"] = []
        for item in data.get("profiles", []):
            if not isinstance(item, dict):
                continue
            cleaned = self._leaves(item, ("implementation", "status", "chronological", "note", "error", "raw_file", "created_at", "revision", "dirty", "source_sha256", "toolchain"))
            if isinstance(item.get("workload"), dict):
                cleaned["workload"] = self._leaves(item["workload"], ("operation", "size", "iterations", "seed"))
            if isinstance(item.get("stacks"), list):
                cleaned["stacks"] = [[self.text(frame) for frame in stack if isinstance(frame, str)] for stack in item["stacks"] if isinstance(stack, list)]
            if isinstance(item.get("weights"), list):
                cleaned["weights"] = copy.deepcopy([weight for weight in item["weights"] if type(weight) in (int, float)])
            result["profiles"].append(cleaned)
        return result
