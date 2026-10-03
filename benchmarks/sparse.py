#!/usr/bin/env python3
"""Cross-port CSR/CG conformance, scaling measurements, and portable live report."""
import argparse
import base64
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import statistics
import time

import run
from publication import PublicSanitizer
from report import json_for_html
from report_design import apply_report_design
from sparse_reference import fixtures as csr_fixtures, protocol, check_result as check_csr, diffusion, multiply
from gmres_reference import fixtures as gmres_fixtures, check_result as check_gmres
from ilu_reference import fixtures as ilu_fixtures
from cholesky_reference import fixtures as cholesky_fixtures
from ordering_reference import fixtures as ordering_fixtures, check_solve as check_ordering_solve
from wasm_publication import sanitize_live_sources

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent


def fixtures():
    return csr_fixtures()+gmres_fixtures()+ilu_fixtures()+ordering_fixtures()+cholesky_fixtures()

def check_result(actual,case):
    if case["op"].startswith("chol_"):
        from cholesky_reference import check_result as check_cholesky
        return check_cholesky(actual,case)
    if case["op"] in ("rcm_solve","ilu_solve"):return check_ordering_solve(actual,case)
    return check_gmres(actual,case) if case["op"]=="gmres" else check_csr(actual,case)

def verify(implementation, cases=None):
    checks = []
    for case in fixtures() if cases is None else cases:
        try:
            args, data = protocol(case)
            result = run.execute(implementation['runner'] + args, stdin=data, env=implementation['env'], timeout=90)
            if case['invalid']:
                assert result.returncode and result.stderr and not result.stdout, 'invalid input accepted'
            else:
                assert result.returncode == 0, result.stderr
                check_result(json.loads(result.stdout), case)
            checks.append(dict(name=case['name'], passed=True))
        except Exception as error:
            checks.append(dict(name=case['name'], passed=False, error=str(error) or 'incorrect sparse result'))
    implementation.update(checks=checks, status='passed' if all(c['passed'] for c in checks) else 'failed')
    return implementation['status'] == 'passed'


def measure(implementation, case, iterations):
    args, data = protocol(case, iterations)
    result = run.execute(implementation['runner'] + args, stdin=data, env=implementation['env'], timeout=120)
    if result.returncode:
        raise ValueError(result.stderr)
    sample = json.loads(result.stdout)
    expected = sum(case['expected']) * iterations
    if sample['iterations'] != iterations or not math.isfinite(sample['elapsed_ns']) or sample['elapsed_ns'] <= 0:
        raise ValueError('invalid timer result')
    if not math.isfinite(sample['checksum']) or abs(sample['checksum'] - expected) > 1e-6 * max(1, abs(expected)):
        raise ValueError('incorrect timed checksum')
    sample['ns_per_op'] = sample['elapsed_ns'] / iterations
    return sample


def benchmark(implementations, sizes, samples):
    results = []
    rng = random.Random(2026)
    for size in sizes:
        a = diffusion(size, 2)
        n, nnz = a['rows'], len(a['values'])
        for operation, op, jacobi in [('dense_spmv', 'dense', 0), ('csr_spmv', 'spmv', 0), ('cg', 'cg', 0), ('cg_jacobi', 'cg', 1)]:
            case = dict(a=a, b=multiply(a, [1.] * n) if op == 'cg' else [1.] * n,
                        op=op, expected=[1.] * n if op == 'cg' else multiply(a, [1.] * n), reason='converged',
                        options=dict(rtol=1e-10, atol=0, limit=2000, jacobi=jacobi, capture=0))
            rows = []
            for impl in implementations:
                row = dict(implementation=impl['id'], operation=operation, size=size, unknowns=n, nnz=nnz,
                           logical_dense_bytes=8*n*n, logical_csr_bytes=16*nnz+8*(n+1), contrast=2, samples=[], status='passed')
                try:
                    args, data = protocol(case)
                    checked = run.execute(impl['runner'] + args, stdin=data, env=impl['env'], timeout=120)
                    if checked.returncode: raise ValueError(checked.stderr)
                    actual = json.loads(checked.stdout)
                    check_result(actual, case)
                    if op == 'cg': row['cg_iterations'] = actual['values'][1]
                    calibration = measure(impl, case, 1)
                    row['iterations'] = max(1, min(10000, math.ceil(20e6 / calibration['elapsed_ns'])))
                except Exception as error:
                    row.update(status='failed', error=str(error))
                rows.append((impl, row))
            for _ in range(samples):
                rng.shuffle(rows)
                for impl, row in rows:
                    if row['status'] != 'passed': continue
                    try: row['samples'].append(measure(impl, case, row['iterations']))
                    except Exception as error: row.update(status='failed', error=str(error))
            for _, row in rows:
                if row['status'] == 'passed':
                    values = [s['ns_per_op'] for s in row['samples']]
                    median = statistics.median(values)
                    row.update(median_ns=median, min_ns=min(values), max_ns=max(values), mad_ns=statistics.median(abs(v-median) for v in values))
                results.append(row)
                print('%s %s %s: %s' % (row['implementation'], operation, size, row['status']), flush=True)
    return results


def live_bundle(data, extra_sources=None):
    sanitizer = PublicSanitizer()
    data = sanitizer.report(data)
    live = dict(available=False, reason='No verified WebAssembly build is available.')
    if any(i['id'] == 'wasm' and i['status'] == 'passed' for i in data['implementations']):
        sources = {'module_source': (ROOT / '.build/wasm/matrix.mjs').read_text(),
                   'wrapper_source': (ROOT / 'ports/wasm/matrix.mjs').read_text(),
                   'heat_source': (HERE / 'sparse-heat.mjs').read_text(),
                   'checks_source': (HERE / 'sparse-checks.mjs').read_text(),
                   'worker_source': (HERE / 'sparse-worker.mjs').read_text()}
        sources.update(extra_sources or {})
        sources, binary = sanitize_live_sources(sources, (ROOT / '.build/wasm/matrix.wasm').read_bytes(), sanitizer)
        digest = hashlib.sha256(binary)
        for key, value in sorted(sources.items()): digest.update((key + value).encode())
        cases = fixtures()
        digest.update(json.dumps(cases, sort_keys=True, allow_nan=False).encode())
        live = dict(available=True, capabilities=["gmres","ilu0","rhs_cache","rcm","cholesky"], **sources, fixtures=cases, wasm_base64=base64.b64encode(binary).decode(), sha256=digest.hexdigest())
    return live


def render(data, destination, *, live_override=None):
    data = PublicSanitizer().report(data)
    live = live_override if live_override is not None else live_bundle(data)
    template = (HERE / 'sparse-report.html').read_text()
    # One pass: embedded source/data cannot introduce replacement markers.
    import re
    replacements = {'DATA': json_for_html(data), 'LIVE': json_for_html(live), 'SCRIPT': (HERE / 'sparse-report.mjs').read_text()}
    if data.get('suite') == 'gmres-transport-v1':
        template = template.replace('<title>Sparse solvers · linear-A</title>', '<title>Restarted GMRES · linear-A</title>').replace('Saved benchmark / Sparse systems', 'Saved benchmark / restarted GMRES').replace('<h1>Sparse systems</h1>', '<h1>Restarted GMRES</h1>').replace('Fewer entries. Less work. Compare sparse matrix operations across languages, then explore how an iterative solver reaches equilibrium.', 'Compare restart lengths, runtime and workspace across eleven implementations. Explore how GMRES solves a nonsymmetric flow and diffusion system.')
    if data.get('suite') == 'ilu-reuse-v1':
        template=template.replace('<title>Sparse solvers · linear-A</title>','<title>ILU(0) and solver reuse · linear-A</title>').replace('Saved benchmark / Sparse systems','Saved benchmark / ILU(0)').replace('<h1>Sparse systems</h1>','<h1>Prepare once. Solve again.</h1>').replace('Fewer entries. Less work. Compare sparse matrix operations across languages, then explore how an iterative solver reaches equilibrium.','Compare ILU(0) setup, repeated GMRES solves and one-shot costs across eleven ports. Explore when stronger preconditioning pays off.')
    html = re.sub(r'@@(DATA|LIVE|SCRIPT)@@', lambda m: replacements[m[1]], template)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(apply_report_design(html))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-build', action='store_true')
    parser.add_argument('--verify-only', action='store_true')
    parser.add_argument('--require-all', action='store_true')
    parser.add_argument('--languages', nargs='+', choices=list(run.NAMES), default=list(run.NAMES))
    parser.add_argument('--sizes', nargs='+', type=int, default=[8, 16, 24])
    parser.add_argument('--samples', type=int, default=3)
    parser.add_argument('--output', type=Path, default=HERE / 'reports/sparse')
    parser.add_argument('--render-only', type=Path)
    args = parser.parse_args()
    if args.render_only:
        render(json.loads(args.render_only.read_text()), args.render_only.with_name('index.html'))
        return 0
    if not 1 <= args.samples <= 100 or any(s < 2 or s > 64 for s in args.sizes): parser.error('samples must be 1–100; grid sizes 2–64')
    started = time.perf_counter()
    implementations = []
    for name in args.languages:
        impl = dict(id=name, name=run.LABELS[name])
        try:
            impl = run.build_one(name, no_build=args.no_build)
            verify(impl)
        except run.ToolchainUnavailable as error: impl.update(status='unavailable', error=str(error))
        except Exception as error: impl.update(status='failed', error=str(error))
        implementations.append(impl)
        print(name, impl['status'], flush=True)
        for check in impl.get('checks', []):
            if not check['passed']: print(check, flush=True)
    available = [i for i in implementations if i['status'] == 'passed']
    results = [] if args.verify_only else benchmark(available, args.sizes, args.samples)
    failed = not available or any(i['status'] == 'failed' or (args.require_all and i['status'] == 'unavailable') for i in implementations) or any(r['status'] != 'passed' for r in results)
    if not args.verify_only:
        data = dict(schema_version=1, created_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    revision=run.require(['git','rev-parse','HEAD']).strip(), dirty=bool(run.require(['git','status','--porcelain']).strip()),
                    source_sha256=run.fingerprint(), suite='sparse-diffusion-v1', seed=2026, total_seconds=time.perf_counter()-started,
                    machine=dict(os=platform.system(), release=platform.release(), architecture=platform.machine(), logical_cpus=os.cpu_count()),
                    methodology=dict(numeric_type='IEEE 754 binary64', matrix_layout='Canonical CSR; 5-point SPD diffusion, conductivity 1–100, grids with zero Dirichlet boundaries.',
                        timing='In-process timers include validation, allocation and checksum. Exclude process startup, parsing, CSR construction, dense conversion, compilation and serialization.',
                        warmup='Three calls before every timed batch; one calibration batch excluded. Target 20 ms per batch, capped at 10,000 calls. Julia warms the identical checksum loop before timing.',
                        sampling='Fresh process per batch; serial, seeded shuffled port order per sample round; median and MAD.',
                        profiles='No profiles collected for this suite.',
                        limitations='Small single-threaded workloads compare these implementations, not languages in isolation. JIT, GC, compiler, scheduling and copying costs differ. Logical storage uses float64 values and 64-bit indices; excludes objects, allocator overhead, vectors and solver workspace; it is not measured process memory. CG starts at zero, rtol=1e-10, atol=0, limit=2000; true residual recomputed each step (two sparse products). Both CG variants solve the same system; iteration counts may differ. No claim of monotone residuals.'),
                    implementations=implementations, results=results, profiles=[])
        data = PublicSanitizer().report(data)
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output / 'results.json').write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
        render(data, args.output / 'index.html')
    return int(failed)


if __name__ == '__main__': raise SystemExit(main())
