"""Remove machine paths from the embedded static WASM bundle without recompiling.

Executable WASM sections are copied byte-for-byte, except identified NUL-ended
path strings in data segments. Those replacements retain their exact byte length.
Function names and source line counts remain useful for profiling/debugging.
"""
import json
import re


def _has_local_path(value, sanitizer):
    if value in ("/", "\\", "file:", "file:/", "file://", "file:///"):
        return False
    if sanitizer.is_local_path(value):
        return True
    # Diagnostic text can contain paths after a function name or other prose.
    candidates = re.findall(r"(?<![\w./:\\])(?:file://|[A-Za-z]:[\\/]|\\\\|/)[^\s\x00\"'`<>]+", value)
    return any(candidate not in ("/", "file://", "file:///") and sanitizer.is_local_path(candidate) for candidate in candidates)


def _sensitive_metadata(value):
    return bool(re.search(
        r"(?i)\b(?:codex|chatgpt)\b|(?:generated|created|written|authored|assisted|produced)\s+(?:with|by)\s+(?:an?\s+)?(?:ai|openai)\b"
        r"|\b(?:[\w-]*(?:token|secret|password|credential|api[_-]?key|access[_-]?key)|authorization|cookie|hostname|username|agent[_ -]?(?:id|name))\s*(?::|=(?!=))"
        r"|\b(?:bearer|basic)\s+[A-Za-z0-9+/=._~-]+|-----BEGIN [^-]*PRIVATE KEY-----"
        r"|\b(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|sk-(?:proj-)?[A-Za-z0-9_-]{16,}|AKIA[A-Z0-9]{16}|AIza[A-Za-z0-9_-]{30,})\b"
        r"|\beyJ[A-Za-z0-9_-]+\.[A-Za-z0-9_-]+\.[A-Za-z0-9_-]+\b|https?://[^/@\s]+@", value)
        or re.search(r"\b[A-Z][A-Z0-9_]{2,}\s*=(?!=)\s*\S+", value))


def _text(value, sanitizer):
    has_path = _has_local_path(value, sanitizer)
    if not has_path and not _sensitive_metadata(value):
        return value
    result = sanitizer.text(value)
    if has_path and value.endswith(("/", "\\")) and not result.endswith(value[-1]):
        result += value[-1]
    if _has_local_path(result, sanitizer):
        raise ValueError("sanitizer left a local path in publication metadata")
    return result


def _uleb(data, position, end):
    value = 0
    for shift in range(0, 35, 7):
        if position >= end:
            raise ValueError("truncated WebAssembly integer")
        byte = data[position]
        position += 1
        value |= (byte & 0x7f) << shift
        if not byte & 0x80:
            if value > 0xffffffff:
                raise ValueError("WebAssembly integer exceeds u32")
            return value, position
    raise ValueError("oversized WebAssembly integer")


def _encode_uleb(value):
    result = bytearray()
    while True:
        byte = value & 0x7f
        value >>= 7
        result.append(byte | (0x80 if value else 0))
        if not value:
            return bytes(result)


def _string(data, position, end):
    length, position = _uleb(data, position, end)
    if length > end - position:
        raise ValueError("truncated WebAssembly string")
    try:
        value = data[position:position + length].decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("invalid WebAssembly metadata UTF-8") from error
    return value, position + length


def _encode_string(value):
    encoded = value.encode("utf-8")
    return _encode_uleb(len(encoded)) + encoded


def _name_map(data, position, end, sanitizer, indirect=False):
    count, position = _uleb(data, position, end)
    output = bytearray(_encode_uleb(count))
    for _ in range(count):
        index, position = _uleb(data, position, end)
        output += _encode_uleb(index)
        if indirect:
            entries, position = _name_map(data, position, end, sanitizer)
            output += entries
        else:
            value, position = _string(data, position, end)
            output += _encode_string(_text(value, sanitizer))
    return bytes(output), position


def _name_section(payload, sanitizer):
    position, output = 0, bytearray()
    while position < len(payload):
        subsection = payload[position]
        length, start = _uleb(payload, position + 1, len(payload))
        end = start + length
        if end > len(payload):
            raise ValueError("truncated WebAssembly name subsection")
        if subsection == 0:
            name, after = _string(payload, start, end)
            result = _encode_string(_text(name, sanitizer))
        elif subsection in (1, 4, 5, 6, 7, 8, 9, 11):
            result, after = _name_map(payload, start, end, sanitizer)
        elif subsection in (2, 3, 10):
            result, after = _name_map(payload, start, end, sanitizer, indirect=True)
        else:
            # Unknown name subsections are optional debug metadata. Omitting
            # them is safer than interpreting a future binary format as text.
            position = end
            continue
        if after != end:
            raise ValueError("unexpected bytes in WebAssembly name subsection")
        output.append(subsection)
        output += _encode_uleb(len(result)) + result
        position = end
    return bytes(output)


def _skip_signed_leb(data, position, end, maximum):
    for _ in range(maximum):
        if position >= end:
            raise ValueError("truncated WebAssembly initializer")
        byte = data[position]
        position += 1
        if not byte & 0x80:
            return position
    raise ValueError("oversized WebAssembly initializer")


def _initializer(data, position, end):
    # Data offsets use integer const/global initializers, optionally extended
    # with integer arithmetic. Reject unknown forms instead of guessing bounds.
    while position < end:
        opcode = data[position]
        position += 1
        if opcode == 0x0b:
            return position
        if opcode in (0x41, 0x42):
            position = _skip_signed_leb(data, position, end, 5 if opcode == 0x41 else 10)
        elif opcode == 0x23:
            _, position = _uleb(data, position, end)
        elif opcode not in (0x6a, 0x6b, 0x6c, 0x7c, 0x7d, 0x7e):
            raise ValueError("unsupported WebAssembly data initializer")
    raise ValueError("unterminated WebAssembly data initializer")


def _data_strings(data, sanitizer):
    output = bytearray(data)
    start = 0
    while start < len(data):
        end = data.find(b"\0", start)
        if end < 0:
            break  # Nonterminated bytes are not identified C strings.
        raw = data[start:end]
        try:
            value = raw.decode("utf-8")
        except UnicodeDecodeError:
            value = ""
        if value and all(character.isprintable() for character in value) and (_has_local_path(value, sanitizer) or _sensitive_metadata(value)):
            safe = _text(value, sanitizer).encode("utf-8")
            if safe == raw:
                raise ValueError("sanitizer did not remove a local WebAssembly data path")
            if len(safe) > len(raw):
                safe = b"[path]" if len(raw) >= 6 else b"~"
            # Preserve segment length and all pointer offsets. A shorter C
            # string ends at the first padding NUL; following bytes stay put.
            output[start:end] = safe + b"\0" * (len(raw) - len(safe))
        start = end + 1
    result = bytes(output)
    remaining = result.decode("utf-8", errors="ignore")
    if _has_local_path(remaining, sanitizer) or _sensitive_metadata(remaining):
        raise ValueError("private metadata remains in opaque WebAssembly bytes; refusing to alter executable data")
    return result


def _data_section(payload, sanitizer):
    count, position = _uleb(payload, 0, len(payload))
    output = bytearray(payload)
    for _ in range(count):
        flags, position = _uleb(payload, position, len(payload))
        if flags == 2:
            _, position = _uleb(payload, position, len(payload))
        if flags in (0, 2):
            position = _initializer(payload, position, len(payload))
        elif flags != 1:
            raise ValueError("unsupported WebAssembly data segment flags")
        length, position = _uleb(payload, position, len(payload))
        end = position + length
        if end > len(payload):
            raise ValueError("truncated WebAssembly data segment")
        output[position:end] = _data_strings(payload[position:end], sanitizer)
        position = end
    if position != len(payload):
        raise ValueError("unexpected bytes after WebAssembly data segments")
    return bytes(output)


def sanitize_wasm(binary, sanitizer):
    """Sanitize a static core WASM module; preserve executable section bytes."""
    if not isinstance(binary, bytes) or binary[:8] != b"\0asm\x01\0\0\0":
        raise ValueError("expected a WebAssembly version 1 core module")
    output = bytearray(binary[:8])
    position = 8
    while position < len(binary):
        section = binary[position]
        if section > 13:
            raise ValueError("unsupported WebAssembly section identifier")
        length, start = _uleb(binary, position + 1, len(binary))
        end = start + length
        if end > len(binary):
            raise ValueError("truncated WebAssembly section")
        payload = binary[start:end]
        if section == 0:
            name, content = _string(payload, 0, len(payload))
            if name in ("dylink", "dylink.0"):
                raise ValueError("dynamic-link WebAssembly modules are not supported for embedding")
            if name == "name":
                updated = _encode_string(name) + _name_section(payload[content:], sanitizer)
                if updated == payload:
                    output += binary[position:end]
                else:
                    output += b"\0" + _encode_uleb(len(updated)) + updated
            # All other custom sections are non-executing metadata in this
            # static bundle: drop DWARF, source maps, paths, and build inventory.
        elif section == 11:
            updated = _data_section(payload, sanitizer)
            if len(updated) != len(payload):
                raise ValueError("WebAssembly data sanitization changed storage size")
            output += binary[position:start] + updated
        else:
            output += binary[position:end]
        position = end
    result = bytes(output)
    remaining = result.decode("utf-8", errors="ignore")
    if _has_local_path(remaining, sanitizer) or _sensitive_metadata(remaining):
        raise ValueError("private metadata remains in opaque WebAssembly bytes; refusing to alter executable data")
    return result


def _decode_js_string(raw):
    output, position = [], 0
    escapes = {"n": "\n", "r": "\r", "t": "\t", "b": "\b", "f": "\f", "v": "\v", "0": "\0"}
    while position < len(raw):
        char = raw[position]
        position += 1
        if char != "\\":
            output.append(char)
            continue
        if position >= len(raw):
            raise ValueError("truncated JavaScript escape")
        char = raw[position]
        position += 1
        if char == "\n":
            continue
        if char == "\r":
            if position < len(raw) and raw[position] == "\n":
                position += 1
            continue
        if char in ("x", "u"):
            if char == "u" and position < len(raw) and raw[position] == "{":
                end = raw.find("}", position + 1)
                if end < 0:
                    raise ValueError("unterminated JavaScript Unicode escape")
                digits = raw[position + 1:end]
                position = end + 1
            else:
                length = 2 if char == "x" else 4
                digits = raw[position:position + length]
                position += length
                if len(digits) != length:
                    raise ValueError("truncated JavaScript Unicode escape")
            try:
                output.append(chr(int(digits, 16)))
            except (ValueError, OverflowError) as error:
                raise ValueError("invalid JavaScript character escape") from error
        else:
            output.append(escapes.get(char, char))
    return "".join(output)


def _quoted(raw, quote, sanitizer):
    value = _decode_js_string(raw)
    safe = _text(value, sanitizer)
    if safe == value:
        return quote + raw + quote
    encoded = json.dumps(safe, ensure_ascii=True)[1:-1]
    if quote == "'":
        encoded = encoded.replace("'", "\\'")
    elif quote == "`":
        encoded = encoded.replace("`", "\\`").replace("${", "\\${")
    # Retain the following code's original source line numbers, even when a
    # path literal used escaped line continuations or a multiline template.
    encoded += "\\\n" * raw.count("\n")
    return quote + encoded + quote


def sanitize_javascript(source, sanitizer):
    """Lex JavaScript strings/comments; leave operators and regex literals alone."""
    if not isinstance(source, str):
        raise TypeError("JavaScript source must be text")
    length = len(source)

    def comment(raw):
        if re.match(r"//[#@]\s*source(?:Mapping)?URL\s*=", raw):
            return "//"
        block = raw.startswith("/*")
        body = raw[2:-2] if block else raw[2:]
        safe = _text(body, sanitizer)
        difference = body.count("\n") - safe.count("\n")
        if difference < 0:
            raise ValueError("sanitization introduced JavaScript comment lines")
        safe += "\n" * difference
        return raw[:2] + safe + ("*/" if block else "")

    def template(position):
        output, start = ["`"], position
        while position < length:
            char = source[position]
            if char == "\\":
                position += 2
            elif char == "`" or source.startswith("${", position):
                chunk = _quoted(source[start:position], "`", sanitizer)[1:-1]
                output.append(chunk)
                if char == "`":
                    return "".join(output) + "`", position + 1
                expression, position = scan(position + 2, True)
                output += ["${", expression, "}"]
                start = position
            else:
                position += 1
        raise ValueError("unterminated JavaScript template")

    def regex_end(position):
        start, in_class = position, False
        position += 1
        while position < length:
            char = source[position]
            if char in "\r\n":
                raise ValueError("unterminated JavaScript regex")
            if char == "\\":
                position += 2
                continue
            if char == "[":
                in_class = True
            elif char == "]":
                in_class = False
            elif char == "/" and not in_class:
                position += 1
                while position < length and (source[position].isalpha() or source[position].isdigit()):
                    position += 1
                raw = source[start:position]
                decoded = raw.replace("\\/", "/").replace("\\\\", "\\")
                if _has_local_path(decoded, sanitizer) or _sensitive_metadata(decoded):
                    raise ValueError("private metadata occurs in executable JavaScript regex code")
                return position
            position += 1
        raise ValueError("unterminated JavaScript regex")

    def scan(position, interpolation=False):
        output, braces, regex_allowed, previous, parentheses = [], 0, True, "", []
        while position < length:
            char = source[position]
            if char.isspace():
                output.append(char)
                position += 1
                continue
            if source.startswith("//", position):
                end = source.find("\n", position)
                if end < 0:
                    end = length
                output.append(comment(source[position:end]))
                position = end
                continue
            if source.startswith("/*", position):
                end = source.find("*/", position + 2)
                if end < 0:
                    raise ValueError("unterminated JavaScript comment")
                end += 2
                output.append(comment(source[position:end]))
                position = end
                continue
            if char in ("'", '"'):
                end = position + 1
                while end < length and source[end] != char:
                    end += 2 if source[end] == "\\" else 1
                if end >= length:
                    raise ValueError("unterminated JavaScript string")
                output.append(_quoted(source[position + 1:end], char, sanitizer))
                position, regex_allowed, previous = end + 1, False, "literal"
                continue
            if char == "`":
                text, position = template(position + 1)
                output.append(text)
                regex_allowed, previous = False, "literal"
                continue
            if char == "/" and regex_allowed:
                end = regex_end(position)
                output.append(source[position:end])
                position, regex_allowed, previous = end, False, "literal"
                continue
            if char.isalpha() or char in "_$":
                end = position + 1
                while end < length and (source[end].isalnum() or source[end] in "_$"):
                    end += 1
                token = source[position:end]
                output.append(token)
                regex_allowed = token in ("return", "throw", "case", "delete", "void", "typeof", "yield", "await", "in", "of", "instanceof")
                previous, position = token, end
                continue
            if char.isdigit():
                end = position + 1
                while end < length and (source[end].isalnum() or source[end] in "._"):
                    end += 1
                output.append(source[position:end])
                position, regex_allowed, previous = end, False, "number"
                continue
            if char == "}" and interpolation and braces == 0:
                return "".join(output), position + 1
            if char == "{":
                braces += 1
            elif char == "}":
                braces -= 1
            if char == "(":
                parentheses.append(previous in ("if", "while", "for", "with", "switch", "catch"))
                regex_allowed = True
            elif char == ")":
                regex_allowed = parentheses.pop() if parentheses else False
            elif char in ("]", "."):
                regex_allowed = False
            elif source.startswith("++", position) or source.startswith("--", position):
                output.append(source[position:position + 2])
                position += 2
                previous = "postfix"
                continue
            else:
                regex_allowed = char != "}"
            output.append(char)
            previous = char
            position += 1
        if interpolation:
            raise ValueError("unterminated JavaScript template expression")
        return "".join(output), position

    result, _ = scan(0)
    if result.count("\n") != source.count("\n"):
        raise ValueError("JavaScript sanitization changed source line counts")
    return result


def sanitize_live_sources(sources, binary, sanitizer):
    """Return independent sanitized sources/binary for report.live_bundle."""
    if not isinstance(sources, dict) or not all(isinstance(name, str) for name in sources):
        raise TypeError("sources must map names to JavaScript source text")
    return ({name: sanitize_javascript(source, sanitizer) for name, source in sources.items()},
            sanitize_wasm(binary, sanitizer))
