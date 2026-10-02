"""Public WASM bundles retain arithmetic, source lines, and storage offsets."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benchmarks"))
from publication import PublicSanitizer
from wasm_publication import sanitize_javascript, sanitize_live_sources, sanitize_wasm


def leb(value):
    result = bytearray()
    while True:
        byte = value & 127
        value >>= 7
        result.append(byte | (128 if value else 0))
        if not value:
            return bytes(result)


def string(value):
    encoded = value.encode()
    return leb(len(encoded)) + encoded


def section(index, payload):
    return bytes([index]) + leb(len(payload)) + payload


def custom(name, content):
    return section(0, string(name) + content)


def read_leb(binary, position):
    result, shift = 0, 0
    while True:
        byte = binary[position]
        position += 1
        result |= (byte & 127) << shift
        if not byte & 128:
            return result, position
        shift += 7


def sections(binary):
    result, position = [], 8
    while position < len(binary):
        length, start = read_leb(binary, position + 1)
        end = start + length
        result.append((binary[position], binary[position:end], binary[start:end]))
        position = end
    return result


def module(data, extras=b""):
    # Standalone answer() -> i32(7), memory, and an active data segment.
    return (b"\0asm\x01\0\0\0"
            + section(1, b"\x01\x60\x00\x01\x7f")
            + section(3, b"\x01\x00")
            + section(5, b"\x01\x00\x01")
            + section(7, b"\x02" + string("answer") + b"\x00\x00" + string("memory") + b"\x02\x00")
            + section(10, b"\x01\x04\x00\x41\x07\x0b")
            + section(11, b"\x01\x00\x41\x00\x0b" + leb(len(data)) + data)
            + extras)


class WasmPublicationTests(unittest.TestCase):
    def setUp(self):
        self.root = "/Users/PRIVATE_USER/Private Project"
        self.sanitizer = PublicSanitizer(root=self.root, home="/Users/PRIVATE_USER")

    def test_debug_metadata_removed_and_function_names_reencoded(self):
        function_names = b"\x01\x00" + string("answer " + self.root + "/ports/c/matrix.c:12")
        names = b"\x01" + leb(len(function_names)) + function_names
        original = module(b"safe\0", custom("name", names)
                          + custom(".debug_info", b"/Users/PRIVATE_USER")
                          + custom("sourceMappingURL", b"file:///tmp/private.map"))
        result = sanitize_wasm(original, self.sanitizer)
        self.assertIn(b"answer ports/c/matrix.c:12", result)
        self.assertNotIn(b"PRIVATE_USER", result)
        self.assertNotIn(b".debug_info", result)
        self.assertNotIn(b"sourceMappingURL", result)
        self.assertEqual([raw for kind, raw, _ in sections(original) if kind],
                         [raw for kind, raw, _ in sections(result) if kind])

    def test_c_file_string_preserves_data_length_and_adjacent_bytes(self):
        path = (self.root + "/ports/c/matrix.c").encode()
        data = b"prefix\0" + path + b"\0" + b"\xff\xfe/opaque/value\0suffix\0"
        original = module(data)
        result = sanitize_wasm(original, self.sanitizer)
        before = {kind: raw for kind, raw, _ in sections(original)}
        after = {kind: raw for kind, raw, _ in sections(result)}
        for kind in before:
            self.assertEqual(len(before[kind]), len(after[kind]))
            if kind != 11:
                self.assertEqual(before[kind], after[kind])
        self.assertIn(b"ports/c/matrix.c\0", result)
        self.assertIn(b"\xff\xfe/opaque/value\0suffix\0", result)
        self.assertEqual(len(original), len(result))

    def test_short_sensitive_data_uses_bounded_fallback_without_overrun(self):
        original = module(b"TOKEN=a\0AFTER\0")
        result = sanitize_wasm(original, self.sanitizer)
        self.assertEqual(len(original), len(result))
        self.assertIn(b"[path]\0\0AFTER\0", result)
        self.assertNotIn(b"TOKEN=", result)

    def test_opaque_path_refuses_publication_instead_of_mutating(self):
        with self.assertRaisesRegex(ValueError, "opaque"):
            sanitize_wasm(module(b"\xff\xfe/private/secret\0"), self.sanitizer)

    def test_indirect_name_maps_keep_indices(self):
        local_map = b"\x01\x00\x01\x03" + string(self.root + "/src/local.c")
        names = b"\x02" + leb(len(local_map)) + local_map
        result = sanitize_wasm(module(b"", custom("name", names)), self.sanitizer)
        self.assertIn(b"\x01\x00\x01\x03" + string("src/local.c"), result)

    def test_truncated_and_dynamic_wasm_fail_explicitly(self):
        for binary in [b"bad", b"\0asm\x01\0\0\0\x00\x80", module(b"x")[:-1],
                       module(b"", custom("dylink.0", b""))]:
            with self.assertRaises(ValueError):
                sanitize_wasm(binary, self.sanitizer)

    def test_js_literals_templates_regex_and_protocol_strings_keep_syntax(self):
        tick, expr = chr(96), "$" + "{1 + 2}"
        source = (f"// built from {self.root}/ports/c/matrix.c\n"
                  f'const metadata = "{self.root}/ports/c/matrix.c";\n'
                  f"const single = '{self.root}/source.js';\n"
                  f'const template = {tick}{self.root}/src/{expr}/matrix.c{tick};\n'
                  + r'const regex = /["\x27/]+/g;' + "\n"
                  + 'const scheme = "file://"; const slash = "/";\n'
                  + 'const resource = "/assets/matrix.wasm";\n'
                  + 'export const add = (a,b) => a / b + 7;\n'
                  + f'//# sourceMappingURL=file://{self.root}/private.map\n')
        result = sanitize_javascript(source, self.sanitizer)
        self.assertEqual(source.count("\n"), result.count("\n"))
        self.assertNotIn(self.root, result)
        self.assertNotIn("sourceMappingURL", result)
        self.assertIn('const metadata = "ports/c/matrix.c";', result)
        self.assertIn(f'const template = {tick}src/{expr}/matrix.c{tick};', result)
        self.assertIn(r'const regex = /["\x27/]+/g;', result)
        self.assertIn('const resource = "/assets/matrix.wasm";', result)
        self.assertIn('const scheme = "file://";', result)

    def test_credentials_environment_and_attribution_are_redacted_in_metadata(self):
        source = ('// TOKEN=fake-private-token; Generated by Codex\n'
                  'const diagnostic = "PRIVATE_SETTING=fake-private-value";\n'
                  'const attribution = "Generated by Codex";\n'
                  'const resource = "/assets/matrix.wasm";\n')
        result = sanitize_javascript(source, self.sanitizer)
        for forbidden in ("fake-private-token", "fake-private-value", "Generated by Codex"):
            self.assertNotIn(forbidden, result)
        self.assertEqual(source.count("\n"), result.count("\n"))
        data = module(b"Generated by Codex\0TOKEN=fake-private-token\0tail\0")
        safe = sanitize_wasm(data, self.sanitizer)
        self.assertNotIn(b"Codex", safe)
        self.assertNotIn(b"fake-private-token", safe)
        self.assertEqual(len(data), len(safe))
        self.assertIn(b"tail\0", safe)

    def test_block_comment_delimiters_survive_credential_redaction(self):
        source = '/* TOKEN=private-token*/\nexport const answer = 7;\n'
        result = sanitize_javascript(source, self.sanitizer)
        self.assertIn('*/\nexport const answer = 7;', result)
        self.assertNotIn('private-token', result)

    def test_escaped_paths_and_continuations_keep_source_line_counts(self):
        source = 'const file = "\\x2fUsers/PRIVATE_USER/Private Project/ports/\\\nwasm/bridge.c";\n'
        result = sanitize_javascript(source, self.sanitizer)
        self.assertNotIn("PRIVATE_USER", result)
        self.assertEqual(source.count("\n"), result.count("\n"))
        self.assertIn("ports/wasm/bridge.c", result)

    @unittest.skipUnless(shutil.which("node"), "Node is needed to execute sanitized WASM")
    def test_sanitized_fixture_executes_arithmetic_and_safe_c_string(self):
        original = module((self.root + "/ports/c/matrix.c").encode() + b"\0")
        result = sanitize_wasm(original, self.sanitizer)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "safe.wasm"
            path.write_bytes(result)
            code = """const fs=require('node:fs');
                WebAssembly.instantiate(fs.readFileSync(process.argv[1])).then(({instance})=>{
                  const bytes=new Uint8Array(instance.exports.memory.buffer);
                  const end=bytes.indexOf(0);
                  console.log(JSON.stringify({answer:instance.exports.answer(),path:new TextDecoder().decode(bytes.slice(0,end))}));
                });"""
            done = subprocess.run([shutil.which("node"), "-e", code, str(path)], text=True, capture_output=True, timeout=10, check=True)
            self.assertEqual(json.loads(done.stdout), {"answer": 7, "path": "ports/c/matrix.c"})

    @unittest.skipUnless(shutil.which("node") and (ROOT / ".build/wasm/matrix.wasm").exists(), "built WASM and Node are needed")
    def test_actual_bundle_preserves_code_and_runs_matrix_operations(self):
        sources = {
            "wrapper_source": (ROOT / "ports/wasm/matrix.mjs").read_text(),
            "module_source": (ROOT / ".build/wasm/matrix.mjs").read_text(),
            "worker_source": (ROOT / "benchmarks/live-worker.mjs").read_text(),
        }
        binary = (ROOT / ".build/wasm/matrix.wasm").read_bytes()
        clean, safe = sanitize_live_sources(sources, binary, PublicSanitizer())
        self.assertEqual(sources, clean)
        self.assertEqual([raw for kind, raw, _ in sections(binary) if kind != 0],
                         [raw for kind, raw, _ in sections(safe) if kind != 0])
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)
            (target / "wrapper.mjs").write_text(clean["wrapper_source"])
            (target / "matrix.mjs").write_text(clean["module_source"])
            (target / "matrix.wasm").write_bytes(safe)
            (target / "check.mjs").write_text("""import {readFileSync} from 'node:fs';
                import {createMatrixAPI} from './wrapper.mjs';
                const {Matrix}=await createMatrixAPI({moduleUrl:new URL('./matrix.mjs',import.meta.url),
                  wasmBinary:readFileSync(new URL('./matrix.wasm',import.meta.url)),locateFile:()=> 'not-fetched.wasm'});
                const a=new Matrix(2,3,[1,2,3,4,5,6]), b=new Matrix(3,2,[7,8,9,10,11,12]);
                const product=a.multiply(b), square=new Matrix(2,2,[1,2,3,4]);
                console.log(JSON.stringify({values:[...product.toArray()],det:square.determinant()}));
                product.dispose();square.dispose();a.dispose();b.dispose();""")
            done = subprocess.run([shutil.which("node"), str(target / "check.mjs")], text=True, capture_output=True, timeout=10, check=True)
            self.assertEqual(json.loads(done.stdout), {"values": [58, 64, 139, 154], "det": -2})


if __name__ == "__main__":
    unittest.main()
