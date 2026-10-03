#!/usr/bin/env python3
"""Build, verify, benchmark, and report the matrix implementations."""
import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import shutil
import statistics
import subprocess
import sys
import time

from reference import OPERATIONS, fixtures as arithmetic_fixtures, assert_result, benchmark_checksum, finite_number
from eigen_reference import eigen_fixtures, generated_eigen, eigen_spectrum, assert_eigen_result
from general_eigen_reference import general_eigen_fixtures, generated_general_eigen, general_eigen_spectrum, assert_general_eigen_result
from solve_reference import solve_fixtures
from pseudoinverse_reference import pseudoinverse_fixtures, assert_inverse_result
from ridge_reference import ridge_fixtures, assert_ridge_result
from svd_reference import svd_fixtures, assert_svd_result
from vector_reference import vector_fixtures, vector_inputs, vector_expected
from publication import PublicSanitizer

OPERATIONS = (*OPERATIONS, "eigen_symmetric", "eigen_general", "cross", "rotation2d", "rotation3d")


def fixtures():
    return arithmetic_fixtures() + eigen_fixtures() + general_eigen_fixtures() + vector_fixtures() + solve_fixtures() + svd_fixtures() + pseudoinverse_fixtures() + ridge_fixtures()

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / ".build"
NAMES = ("csharp", "fsharp", "rust", "go", "typescript", "python", "cpp", "c", "assembly", "julia", "wasm")
LABELS = {"csharp":"C#", "fsharp":"F#", "rust":"Rust", "go":"Go", "typescript":"TypeScript",
          "python":"Python", "cpp":"C++", "c":"C", "assembly":"ARM64 Assembly", "julia":"Julia", "wasm":"WebAssembly"}
DETAILS = {
    "csharp":"Matrix<double>; type-aware determinant with pivoting and scaling; rectangular array storage",
    "fsharp":"F# float64 matrix; managed flat storage; partial-pivot determinant",
    "rust":"Rust f64 matrix; contiguous storage; optimized release build",
    "go":"Go float64 matrix; contiguous storage; default optimizing compiler",
    "typescript":"TypeScript compiled to JavaScript; Float64Array; Node.js JIT",
    "python":"Python float matrix; interpreted loops; no NumPy",
    "cpp":"C++ double matrix; optimized native build; no fast-math; shared C general eigen, SVD and factorization kernels",
    "c":"C double matrix; optimized native build; no fast-math",
    "assembly":"ARM64 arithmetic kernels with C allocation, validation, pivot control, shared real eigensolvers, SVD and system solvers",
    "julia":"Julia Float64 matrix; handwritten loops; JIT; no BLAS calls",
    "wasm":"C double matrix compiled with Emscripten; Node.js WebAssembly runtime; JS/WASM dispatch, allocation, checksum reads, and disposal included"
}


def execute(command, *, stdin=None, timeout=180, env=None, cwd=ROOT):
    completed = subprocess.run([str(x) for x in command], input=stdin, text=True, capture_output=True,
                               cwd=cwd, timeout=timeout, env=env)
    return completed


def require(command, **kwargs):
    result = execute(command, **kwargs)
    if result.returncode:
        raise RuntimeError("%s\n%s" % (" ".join(map(str,command)), (result.stderr or result.stdout)[-6000:]))
    return result.stdout


class ToolchainUnavailable(RuntimeError):
    """The requested compiler/runtime is missing or unsupported on this host."""


def binary(env_name, default):
    value = os.environ.get(env_name) or shutil.which(default)
    if not value or not shutil.which(value):
        raise ToolchainUnavailable("%s not found; install it or set %s" % (value or default,env_name))
    return value


def build_one(name, no_build=False):
    target = BUILD / name
    target.mkdir(parents=True,exist_ok=True)
    env = os.environ.copy()
    env.setdefault("DOTNET_CLI_TELEMETRY_OPTOUT","1")
    env.setdefault("DOTNET_NOLOGO","1")
    commands = []
    if name in ("csharp","fsharp"):
        compiler = binary("DOTNET","dotnet")
        project = ROOT / ("benchmarks/dotnet/Runner.csproj" if name=="csharp" else "ports/fsharp/Runner.fsproj")
        commands = [[compiler,"build",project,"-c","Release","-o",target,"--disable-build-servers","-m:1"]]
        runner = [compiler,str(target / "Runner.dll")]
        version = [compiler,"--version"]
    elif name == "rust":
        compiler = binary("CARGO","cargo")
        commands = [[compiler,"build","--release","--locked","--manifest-path",ROOT/"ports/rust/Cargo.toml","--target-dir",target]]
        runner = [str(target/"release/linear-a-rust")]
        version = [binary("RUSTC","rustc"),"--version"]
    elif name == "go":
        compiler = binary("GO","go")
        commands = [[compiler,"build","-o",target/"runner","./cmd/runner"]]
        runner = [str(target/"runner")]
        version = [compiler,"version"]
        env.setdefault("GOCACHE",str(BUILD/"go-cache"))
    elif name == "typescript":
        compiler = binary("NPM","npm")
        commands = [[compiler,"ci","--ignore-scripts"],[compiler,"run","build"]]
        runner = [binary("NODE","node"),str(ROOT/"ports/typescript/dist/runner.js")]
        env["PATH"] = str(Path(runner[0]).parent) + os.pathsep + env.get("PATH", "")
        version = [runner[0],"--version"]
    elif name == "python":
        compiler = binary("PYTHON","python3")
        runner = [compiler,str(ROOT/"ports/python/runner.py")]
        version = [compiler,"--version"]
    elif name == "julia":
        compiler = binary("JULIA","julia")
        runner = [compiler,"--startup-file=no","--threads=1","--project="+str(ROOT/"ports/julia"),str(ROOT/"ports/julia/runner.jl")]
        version = [compiler,"--version"]
    elif name == "wasm":
        compiler = binary("EMCC","emcc")
        node = binary("NODE","node")
        commands = [[sys.executable,ROOT/"ports/wasm/build.py","--emcc",compiler,"--output",target]]
        runner = [node,str(ROOT/"ports/wasm/runner.mjs")]
        env["LINEAR_A_WASM_MODULE"] = str(target/"matrix.mjs")
        env["PATH"] = str(Path(node).parent) + os.pathsep + env.get("PATH", "")
        version = [compiler,"--version"]
    elif name in ("c","cpp","assembly"):
        if name == "assembly" and platform.machine().lower() not in ("arm64","aarch64"):
            raise ToolchainUnavailable("ARM64 kernels require an ARM64 host")
        compiler = binary("CXX" if name=="cpp" else "CC","clang++" if name=="cpp" else "clang")
        flags = ["-std=c++17" if name=="cpp" else "-std=c11","-O3","-g","-ffp-contract=off"]
        sources = [ROOT/"ports/cpp/runner.cpp"] if name=="cpp" else [ROOT/"ports/c/matrix.c",ROOT/"ports/c/runner.c"]
        if name=="assembly":
            flags += ["-DMATRIX_USE_ASM"]
            sources += [ROOT/"ports/assembly/kernels.S"]
        commands = [[compiler]+flags+sources+["-lm","-o",target/"runner"]]
        runner = [str(target/"runner")]
        version = [compiler,"--version"]
    else:
        raise ValueError(name)
    cwd = ROOT / ("ports/go" if name=="go" else "ports/typescript" if name=="typescript" else ".")
    toolchain = require(version,env=env).strip().splitlines()[0]
    if not no_build:
        for command in commands:
            require(command,env=env,cwd=cwd,timeout=600)
    if name == "typescript" and (ROOT/"ports/typescript/node_modules/typescript/package.json").exists():
        toolchain += "; TypeScript " + json.loads((ROOT/"ports/typescript/node_modules/typescript/package.json").read_text())["version"]
    if name == "wasm":
        toolchain += "; Node.js " + require([runner[0],"--version"],env=env).strip()
    return {"id":name,"name":LABELS[name],"status":"built","toolchain":toolchain,"detail":DETAILS[name],
            "build_commands":[PublicSanitizer(root=ROOT).command(command) for command in commands],
            "runner":runner,"env":env}


def verify(implementation, cases=None):
    outcomes=[]
    for case in fixtures() if cases is None else cases:
        a,b = case["a"],case["b"]
        arguments = ["check",case["op"],str(a["rows"]),str(a["cols"])]
        values = a["values"][:]
        if b is not None:
            arguments += [str(b["rows"]),str(b["cols"])]
            values += b["values"]
        elif case["op"] in ("scale", "rotation2d", "rotation3d", "pseudoinverse", "spectral_diagnostics"): arguments.append(str(case["scalar"]))
        if case["op"] == "solve_ridge": arguments.append(str(case["scalar"]))
        outcome = {"name":case["name"],"passed":False}
        try:
            result = execute(implementation["runner"]+arguments,stdin=" ".join(format(x,".17g") for x in values),
                             env=implementation["env"],timeout=90)
            if case["invalid"]:
                if result.returncode==0: raise AssertionError("Invalid shape accepted")
                if not result.stderr.strip(): raise AssertionError("Invalid shape produced no error message")
                if result.stdout.strip(): raise AssertionError("Invalid request produced output on stdout")
            else:
                if result.returncode: raise AssertionError(result.stderr.strip()[-2000:])
                if case["op"] == "solve_ridge":
                    assert_ridge_result(json.loads(result.stdout),case)
                elif case["op"] in ("pseudoinverse", "solve_minimum_norm", "spectral_diagnostics"):
                    assert_inverse_result(json.loads(result.stdout),case)
                elif case["op"] == "svd":
                    assert_svd_result(json.loads(result.stdout), a, case.get("spectrum"))
                elif case["op"] == "eigen_general":
                    assert_general_eigen_result(json.loads(result.stdout), a, case.get("expected_complex_eigenvalues"),
                                                spectrum_scale=case.get("spectrum_scale"), componentwise=case.get("componentwise",False))
                elif case["op"] == "eigen_symmetric":
                    assert_eigen_result(json.loads(result.stdout), a, case.get("expected_eigenvalues"))
                else:
                    assert_result(json.loads(result.stdout),case["expected"])
            outcome["passed"]=True
        except (ValueError,AssertionError,subprocess.TimeoutExpired,OverflowError,OSError) as error:
            outcome["error"]=str(error)
        outcomes.append(outcome)
    implementation["checks"]=outcomes
    implementation["status"]="passed" if all(x["passed"] for x in outcomes) else "failed"
    return implementation["status"]=="passed"


def measure(implementation,op,size,iterations,seed,expected_checksum,timeout=180):
    command_op = implementation.get("operation_map", {}).get(op, op)
    output = require(implementation["runner"]+["bench",command_op,str(size),str(iterations),str(seed)],env=implementation["env"],timeout=timeout)
    data = json.loads(output)
    if not isinstance(data,dict) or set(data)!={"elapsed_ns","iterations","checksum"}:
        raise ValueError("Timing result must be an object with elapsed_ns, iterations, and checksum")
    if type(data["iterations"]) is not int or data["iterations"]!=iterations:
        raise ValueError("Runner reported wrong iteration count")
    elapsed,checksum = data["elapsed_ns"],data["checksum"]
    if not finite_number(elapsed) or elapsed <= 0:
        raise ValueError("Invalid elapsed time")
    if not finite_number(checksum) or not math.isclose(checksum,expected_checksum*iterations,rel_tol=1e-7,abs_tol=1e-8*iterations):
        raise ValueError("Benchmark checksum mismatch: %r, expected %r" % (checksum,expected_checksum*iterations))
    return {"elapsed_ns":elapsed,"iterations":iterations,"checksum":checksum,"ns_per_op":elapsed/iterations}


def benchmark(implementations,options):
    rows=[]
    rng=random.Random(options.seed)
    sizes = [16,48] if options.suite=="quick" else [64,128,256]
    det_sizes = getattr(options, "det_sizes", None) or ([8,16] if options.suite=="quick" else [8,24,48])
    eigen_sizes = [8,16] if options.suite=="quick" else [16,32,48]
    samples=options.samples or (3 if options.suite=="quick" else 5)
    target_ns=20_000_000 if options.suite=="quick" else 100_000_000
    max_iterations=4096 if options.suite=="quick" else 100_000
    for op in options.operations:
        operation_sizes = ([2] if op == "rotation2d" else [3]) if op in ("cross", "rotation2d", "rotation3d") else det_sizes if op.startswith("determinant") else eigen_sizes if op in ("eigen_symmetric", "eigen_general") else sizes
        for size in operation_sizes:
            expected_checksum=benchmark_checksum(op,size,options.seed)
            order=implementations[:]
            rng.shuffle(order)
            active=[]
            for implementation in order:
                row={"implementation":implementation["id"],"operation":op,"size":size,"samples":[],"status":"passed"}
                rows.append(row)
                try:
                    if op in ("cross", "rotation2d", "rotation3d"):
                        a, b, scalar = vector_inputs(op, size, options.seed)
                        arguments = ["check", op, str(a["rows"]), str(a["cols"])]
                        values = a["values"][:]
                        if b is not None:
                            arguments += [str(b["rows"]), str(b["cols"])]
                            values += b["values"]
                        else:
                            arguments.append(str(scalar))
                        output = require(implementation["runner"] + arguments,
                                         stdin=" ".join(format(x, ".17g") for x in values), env=implementation["env"], timeout=90)
                        assert_result(json.loads(output), vector_expected(op, a, b, scalar))
                    if op == "eigen_general":
                        a = generated_general_eigen(size, options.seed)
                        output = require(implementation["runner"] + ["check", op, str(size), str(size)],
                            stdin=" ".join(format(x,".17g") for x in a["values"]),env=implementation["env"],timeout=90)
                        assert_general_eigen_result(json.loads(output),a,general_eigen_spectrum(size,options.seed))
                    if op == "eigen_symmetric":
                        # Validate the actual timed input's complete eigenbasis outside the clock.
                        a = generated_eigen(size, options.seed)
                        values = " ".join(format(x, ".17g") for x in a["values"])
                        output = require(implementation["runner"] + ["check", op, str(size), str(size)],
                                         stdin=values, env=implementation["env"], timeout=90)
                        assert_eigen_result(json.loads(output), a, eigen_spectrum(size, options.seed))
                    trial=measure(implementation,op,size,1,options.seed,expected_checksum)
                    iterations=max(1,min(max_iterations,round(target_ns/trial["ns_per_op"])))
                    # Recalibrate after JIT warmup, without including calibration in the reported samples.
                    calibration=measure(implementation,op,size,iterations,options.seed,expected_checksum)
                    iterations=max(1,min(max_iterations,round(target_ns/calibration["ns_per_op"])))
                    row["iterations"]=iterations
                    active.append((implementation,row))
                except Exception as error:
                    row.update(status="failed",error=str(error)[-3000:])
            for _ in range(samples):
                rng.shuffle(active)
                for implementation,row in active:
                    if row["status"]!="passed": continue
                    try:
                        row["samples"].append(measure(implementation,op,size,row["iterations"],options.seed,expected_checksum))
                    except Exception as error:
                        row.update(status="failed",error=str(error)[-3000:])
            for _,row in active:
                if row["status"]=="passed":
                    values=[sample["ns_per_op"] for sample in row["samples"]]
                    median=statistics.median(values)
                    row.update(median_ns=median,min_ns=min(values),max_ns=max(values),
                               mad_ns=statistics.median(abs(value-median) for value in values))
            print("  timed %s %dx%d" % (op,size,size),flush=True)
    return rows


def fingerprint():
    digest=hashlib.sha256()
    exclusions={"node_modules","target","bin","obj","dist","__pycache__"}
    for directory in (ROOT/"ports",ROOT/"linear-A/linear-A",ROOT/"benchmarks"):
        for path in sorted(directory.rglob("*")):
            if path.is_file() and not any(part in exclusions or part=="reports" for part in path.relative_to(directory).parts) and path.suffix in (".c",".h",".cpp",".hpp",".S",".rs",".go",".ts",".mjs",".py",".jl",".fs",".cs",".json",".toml"):
                digest.update(str(path.relative_to(ROOT)).encode())
                digest.update(path.read_bytes())
    return digest.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--implementations",nargs="+",choices=NAMES,default=list(NAMES))
    parser.add_argument("--suite",choices=("quick","full"),default="quick")
    parser.add_argument("--operations",nargs="+",choices=OPERATIONS,default=list(OPERATIONS))
    parser.add_argument("--samples",type=int)
    parser.add_argument("--seed",type=int,default=17)
    parser.add_argument("--output",type=Path,default=ROOT/"benchmarks/reports/latest")
    parser.add_argument("--verify-only",action="store_true")
    parser.add_argument("--no-build",action="store_true")
    parser.add_argument("--require-all",action="store_true",help="Fail when any requested toolchain is unavailable")
    parser.add_argument("--profiles",action="store_true",help="Collect real stack profiles separately from benchmark timings")
    options=parser.parse_args()
    if not 0<=options.seed<=2147483646 or (options.samples is not None and options.samples<1): parser.error("Invalid seed or sample count")
    options.output.mkdir(parents=True,exist_ok=True)
    implementations=[]
    available=[]
    started=time.perf_counter()
    for name in options.implementations:
        print("Building and checking %s..." % LABELS[name],flush=True)
        implementation={"id":name,"name":LABELS[name]}
        implementations.append(implementation)
        try:
            implementation.update(build_one(name,options.no_build))
            if verify(implementation):
                available.append(implementation)
                print("  %d/%d checks passed" % (len(implementation["checks"]),len(implementation["checks"])),flush=True)
            else:
                for check in implementation["checks"]:
                    if not check["passed"]: print("  FAIL: %s: %s" % (check["name"],check["error"]),flush=True)
        except ToolchainUnavailable as error:
            implementation.update(status="unavailable",error=str(error))
            print("  unavailable: "+str(error),flush=True)
        except Exception as error:
            implementation.update(status="failed",error=str(error)[-5000:])
            print("  build/check failed: "+str(error)[-1500:],flush=True)
    results=[] if options.verify_only else benchmark(available,options)
    profiles=[]
    if options.profiles and not options.verify_only:
        from profiles import collect_profiles
        profiles=collect_profiles(available,options.output,options.seed)
    git=require(["git","rev-parse","HEAD"]).strip()
    dirty=bool(require(["git","status","--porcelain"]).strip())
    report={"schema_version":1,"created_at":datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "machine":{"os":platform.system(),"release":platform.release(),"architecture":platform.machine(),"logical_cpus":os.cpu_count()},
            "revision":git,"dirty":dirty,"source_sha256":fingerprint(),"suite":options.suite,"seed":options.seed,
            "total_seconds":time.perf_counter()-started,"methodology":{
                "workload_version":"matrix-protocol-v1",
                "sample_unit":"fresh_process_per_batch",
                "numeric_type":"IEEE 754 binary64 (double precision)","matrix_layout":"Row-major protocol; implementation storage documented separately",
                "timing":"In-process monotonic timer; allocation and checksum included; startup, input creation, compilation and JSON excluded",
                "warmup":"max(5, min(iterations, 100)) operations before each timing sample; two calibration batches excluded",
                "sampling":"Serial runs, deterministic randomized implementation order each sample round; median and MAD reported",
                "profiles":"Separate instrumented or sampled runs; never used for comparative timings",
                "limitations":"Symmetric eigenpairs use cyclic Jacobi at relative Frobenius tolerance 1e-12 and 50 sweeps; vectors are columns and nonsymmetric inputs are rejected. General real eigenpairs use balanced Hessenberg reduction and bounded double-shift QR, with normalized complex right columns that may be dependent for defective inputs. Eigen benchmark inputs have analytic Toeplitz or block-triangular spectra and their complete eigenbasis is checked outside timing. These compare implementations and algorithms, not languages in isolation. Runtime/JIT/GC, compiler age, OS scheduling, thermals and allocation affect results. Small quick workloads are a smoke test, not a universal ranking."},
            "implementations":[{k:v for k,v in item.items() if k not in ("runner","env")} for item in implementations],
            "results":results,"profiles":profiles}
    from report import publish
    publish(report,options.output/"index.html")
    failures=any(item["status"]=="failed" or (options.require_all and item["status"]=="unavailable") for item in implementations) or any(row["status"]!="passed" for row in results)
    print("Report: %s" % (options.output/"index.html"),flush=True)
    return 1 if failures or not available else 0


if __name__=="__main__":
    raise SystemExit(main())
