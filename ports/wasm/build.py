#!/usr/bin/env python3
"""Compile the C float64 implementation to an ES module and WebAssembly binary."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
EXPORTS = [
    "wm_multigrid_matrix", "wm_multigrid_apply", "wm_csr_cg_multigrid",
    "wm_ic0_create", "wm_csr_cg_preconditioned", "wm_cholesky_analyze", "wm_cholesky_symbolic_destroy", "wm_cholesky_pattern", "wm_cholesky_factorize", "wm_cholesky_destroy", "wm_cholesky_lower", "wm_cholesky_solve",
    "wm_csr_rcm", "wm_csr_amd", "wm_csr_permute", "wm_permute_vector", "wm_csr_gmres", "wm_csr_gmres_preconditioned", "wm_ilu0_create", "wm_ilu0_apply", "wm_ilu0_destroy", "wm_gmres_destroy", "wm_gmres_iterations", "wm_gmres_reason", "wm_gmres_restart_count", "wm_gmres_restarts", "wm_gmres_data", "wm_csr_create", "wm_csr_destroy", "wm_csr_matvec", "wm_csr_cg", "wm_cg_destroy", "wm_cg_iterations", "wm_cg_reason", "wm_cg_data",
    "wm_solve_ridge", "wm_svd", "wm_pseudoinverse", "wm_solve_minimum_norm", "wm_spectral_diagnostics",
    "wm_factorize", "wm_factor_destroy", "wm_factor_solve", "wm_factor_rcond",
    "wm_create", "wm_identity", "wm_destroy", "wm_rows", "wm_cols", "wm_data", "wm_copy",
    "wm_add", "wm_subtract", "wm_scale", "wm_transpose", "wm_multiply", "wm_cross", "wm_row", "wm_column",
    "wm_rotation_2d", "wm_rotation_x", "wm_rotation_y", "wm_rotation_z", "wm_rotation_axis_angle",
    "wm_set", "wm_trace", "wm_determinant", "wm_determinant_algorithm", "wm_eigen_symmetric", "wm_eigen_general", "wm_triangular", "wm_checksum", "wm_last_error", "wm_error_message",
]


def trace_source_snapshot(source):
    """Keep original line numbers using this file's top-level brace style."""
    lines = source.decode("utf-8").splitlines()
    start = next(index for index, line in enumerate(lines) if line.startswith("matrix_status m_multiply("))
    end = next(index for index in range(start + 1, len(lines)) if lines[index] == "}")
    return {
        "path": "ports/c/matrix.c",
        "start_line": start + 1,
        "lines": lines[start:end + 1],
        "source_sha256": hashlib.sha256(source).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--emcc", default=shutil.which("emcc"), help="Emscripten compiler path")
    parser.add_argument("--output", type=Path, help="output directory (default: .build/wasm or .build/wasm-trace)")
    parser.add_argument("--trace", action="store_true", help="build bounded multiplication capture for visualization")
    args = parser.parse_args()
    if not args.emcc:
        parser.error("Emscripten emcc was not found; pass --emcc or activate an emsdk environment")
    output = (args.output or ROOT / (".build/wasm-trace" if args.trace else ".build/wasm")).resolve()
    if args.trace and output == (ROOT / ".build/wasm").resolve():
        parser.error("the trace build must use a separate output directory from .build/wasm")
    output.mkdir(parents=True, exist_ok=True)
    exports = list(EXPORTS)
    trace_options = []
    source_path = ROOT / "ports/c/matrix.c"
    if args.trace:
        source = source_path.read_bytes()
        snapshot = trace_source_snapshot(source)
        manifest_path = output / "trace-source.json"
        # A failed or interrupted build must not leave an old source manifest
        # attached to a binary that the compiler might already have replaced.
        manifest_path.unlink(missing_ok=True)
        trace_options = ["-DMATRIX_TRACE", str(ROOT / "ports/wasm/trace.c")]
        exports += ["wm_multiply_trace", "wm_trace_count", "wm_trace_data", "wm_trace_stride"]
    command = [
        args.emcc, str(ROOT / "ports/c/matrix.c"), str(ROOT / "ports/wasm/bridge.c"),
        *trace_options,
        "-std=c11", "-O3", "-g2", "-ffp-contract=off", "-fno-fast-math",
        "--no-entry", "-sMODULARIZE=1", "-sEXPORT_ES6=1", "-sENVIRONMENT=web,worker,node",
        "-sALLOW_MEMORY_GROWTH=1", "-sABORTING_MALLOC=0", "-sFILESYSTEM=0",
        '-sINCOMING_MODULE_JS_API=["locateFile","wasmBinary"]',
        "-sEXPORTED_FUNCTIONS=" + json.dumps(["_" + name for name in exports]),
        '-sEXPORTED_RUNTIME_METHODS=["UTF8ToString","HEAPF64"]',
        "-o", str(output / "matrix.mjs"),
    ]
    subprocess.run(command, cwd=ROOT, check=True)
    if args.trace:
        if source_path.read_bytes() != source:
            raise SystemExit("matrix.c changed during compilation; rerun the trace build to capture its matching source")
        manifest_path.write_text(json.dumps(snapshot, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
