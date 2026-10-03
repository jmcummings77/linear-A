#!/usr/bin/env python3
"""Compare determinant algorithms on identical general and SPD matrices."""
import argparse
import copy
import datetime
import os
from pathlib import Path
import platform
import time
from types import SimpleNamespace

from reference import fixtures, matrix, determinant, generated_spd
from run import ROOT, NAMES, LABELS, build_one, verify, benchmark, fingerprint, require, ToolchainUnavailable
from report import publish


def algorithm_cases(algorithm):
    if algorithm != "cholesky":
        cases = [copy.deepcopy(case) for case in fixtures() if case["op"] == "determinant"]
        for case in cases:
            case["op"] = "determinant" if algorithm == "auto" else "determinant_" + algorithm
        return cases
    cases = []
    inputs = [("empty SPD", matrix(0, 0, []), False)]
    inputs += [("generated SPD %d" % n, generated_spd(n, 17), False) for n in range(1, 6)]
    inputs += [("nonsymmetric", matrix(2, 2, [1, 2, 0, 1]), True),
               ("indefinite", matrix(2, 2, [1, 2, 2, 1]), True),
               ("semidefinite", matrix(2, 2, [1, 1, 1, 1]), True),
               ("nonsquare", matrix(1, 2, [1, 2]), True)]
    for name, a, invalid in inputs:
        case = {"name":name, "op":"determinant_cholesky", "a":a, "b":None, "invalid":invalid}
        if not invalid:
            case["expected"] = {"value":determinant(a)}
        cases.append(case)
    return cases


def variant(base, algorithm):
    item = dict(base)
    item["id"] = base["id"] if algorithm == "auto" else base["id"] + "-" + algorithm
    item["name"] = base["name"] + " / " + {"auto":"Auto", "lu":"LU", "cholesky":"Cholesky", "cofactor":"Cofactor"}[algorithm]
    item["detail"] = "Determinant algorithm: " + algorithm + ". " + base["detail"]
    item["operation_map"] = {"determinant":"determinant" if algorithm == "auto" else "determinant_" + algorithm,
                             "determinant_spd":"determinant_cholesky" if algorithm == "cholesky" else "determinant_spd_lu",
                             "determinant_small":"determinant" if algorithm == "auto" else "determinant_" + algorithm}
    verify(item, algorithm_cases(algorithm))
    return item


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--implementations", nargs="+", choices=NAMES, default=list(NAMES))
    parser.add_argument("--suite", choices=("quick", "full"), default="quick")
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--require-all", action="store_true")
    parser.add_argument("--output", type=Path, default=ROOT/"benchmarks/reports/determinants")
    options = parser.parse_args()
    if options.samples < 1 or not 0 <= options.seed <= 2147483646:
        parser.error("Invalid sample count or seed")
    # The normal saved comparison is intentionally a different default destination.
    options.output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    implementations, general, spd, small = [], [], [], []
    for name in options.implementations:
        print("Building and checking %s determinant algorithms..." % LABELS[name], flush=True)
        try:
            base = build_one(name, options.no_build)
            if not verify(base):
                implementations.append(base)
                continue
            for algorithm in ("auto", "lu", "cholesky") + (("cofactor",) if name == "csharp" else ()):
                item = variant(base, algorithm)
                implementations.append(item)
                if item["status"] != "passed":
                    for check in item["checks"]:
                        if not check["passed"]: print("  FAIL %s: %s" % (check["name"], check["error"]), flush=True)
                    continue
                if algorithm in ("auto", "lu"): general.append(item)
                if algorithm in ("lu", "cholesky"): spd.append(item)
                if name == "csharp" and algorithm != "cholesky": small.append(item)
        except Exception as error:
            unavailable = isinstance(error, ToolchainUnavailable)
            implementations.append({"id":name, "name":LABELS[name], "status":"unavailable" if unavailable else "failed", "error":str(error)})
            print("  %s: %s" % ("unavailable" if unavailable else "failed", error), flush=True)
    results = []
    if not options.verify_only:
        for cohort, operation, sizes in [(general, "determinant", None), (spd, "determinant_spd", None), (small, "determinant_small", [3,5,7])]:
            if cohort:
                config = SimpleNamespace(**vars(options), operations=[operation], det_sizes=sizes)
                results.extend(benchmark(cohort, config))
    data = {"schema_version":1, "created_at":datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "machine":{"os":platform.system(), "release":platform.release(), "architecture":platform.machine(), "logical_cpus":os.cpu_count()},
        "revision":require(["git","rev-parse","HEAD"]).strip(), "dirty":bool(require(["git","status","--porcelain"]).strip()),
        "source_sha256":fingerprint(), "suite":options.suite, "seed":options.seed, "total_seconds":time.perf_counter()-started,
        "methodology":{
            "workload_version":"determinant-study-v1",
            "sample_unit":"fresh_process_per_batch",
            "numeric_type":"IEEE 754 binary64; exact rational determinant references",
            "matrix_layout":"Identical row-major inputs within each workload; determinant_spd is symmetric positive definite; determinant_small compares the C# cofactor baseline",
            "timing":"In-process monotonic timer; allocation, factorization, checksum, and disposal included; setup, startup and JSON excluded",
            "warmup":"max(5, min(iterations, 100)) operations before each sample; two calibration batches excluded",
            "sampling":"Serial deterministic randomized order, median and MAD; independent checksum validation for every sample",
            "profiles":"No profiles collected in this determinant comparison",
            "limitations":"Compare algorithms within the same workload and size. General and SPD inputs differ; do not take speed ratios across them. Cofactor is restricted to sizes 3, 5 and 7. Cholesky requires exact symmetry and positive computed pivots. These are implementation comparisons, not language rankings."},
        "implementations":[{k:v for k,v in item.items() if k not in ("runner","env","operation_map")} for item in implementations],
        "results":results, "profiles":[]}
    publish(data, options.output/"index.html")
    print("Report: %s" % (options.output/"index.html"), flush=True)
    failed = any(item["status"] == "failed" or (options.require_all and item["status"] == "unavailable") for item in implementations)
    return int(failed or any(row["status"] != "passed" for row in results) or not general)


if __name__ == "__main__":
    raise SystemExit(main())
