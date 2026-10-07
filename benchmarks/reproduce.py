#!/usr/bin/env python3
"""Capture publishable benchmarks from a clean commit in a fresh local clone."""
import argparse
import json
import math
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile

from check_provenance import EXPERIMENT_RUNNERS, artifact_hashes, forbidden_options, validate_report

ROOT = Path(__file__).resolve().parents[1]
RUNNERS = (
    'run.py', 'determinants.py', 'sparse.py', 'gmres_bench.py', 'ilu_bench.py',
    'ordering_bench.py', 'cholesky_bench.py', 'amd_bench.py', 'ic0_bench.py', 'multigrid_bench.py',
) + tuple(EXPERIMENT_RUNNERS)


def git(root, *arguments):
    return subprocess.check_output(['git', *arguments], cwd=root, text=True).strip()


def checked_arguments(runner, arguments):
    if runner not in RUNNERS:
        raise ValueError('Unsupported benchmark runner: ' + runner)
    arguments = list(arguments)
    # The harness parsers accept abbreviated long options. Reject prefixes too,
    # so --no-b cannot accidentally benchmark a stale binary.
    forbidden = forbidden_options(EXPERIMENT_RUNNERS.get(runner, 'benchmarks/' + runner))
    for argument in arguments:
        option = argument.split('=', 1)[0]
        if option in ('--', '-h') or (option.startswith('--') and
                any(flag.startswith(option) for flag in forbidden)):
            raise ValueError('Publication capture does not accept ' + option)
    # IC(0), multigrid and the experiments already fail unless every selected
    # implementation passes; their parsers do not expose --require-all.
    if runner not in ('ic0_bench.py', 'multigrid_bench.py', *EXPERIMENT_RUNNERS) and '--require-all' not in arguments:
        arguments.append('--require-all')
    return arguments


def validate_nonnormal_completion(data):
    """Require the complete bounded diagnostic family and its captured evidence."""
    def require(condition, detail):
        if not condition:
            raise ValueError('Nonnormal GMRES experiment: ' + detail)

    def finite(value):
        try:
            return type(value) in (int, float) and math.isfinite(value)
        except OverflowError:
            return False

    def vector(value, length):
        return isinstance(value, list) and len(value) == length and all(map(finite, value))

    def close(actual, expected, tolerance=1e-12):
        return finite(actual) and abs(actual - expected) <= tolerance

    couplings = {0.0, 0.25, 0.75, 2.0, 4.0}
    limit = data.get('max_iterations')
    require(type(data.get('schema_version')) is int and data['schema_version'] == 1
            and data.get('suite') == 'nonnormal-gmres-v1'
            and data.get('implementation') == 'python', 'unexpected result schema')
    require(type(limit) is int and 2 <= limit <= 100 and data.get('rtol') == 1e-10,
            'invalid solver limits')
    polynomials, runs = data.get('polynomials'), data.get('runs')
    require(isinstance(polynomials, list) and len(polynomials) == 5,
            'all five polynomial checks are required')
    seen = set()
    for row in polynomials:
        require(isinstance(row, dict) and finite(row.get('s')) and row['s'] in couplings,
                'invalid polynomial coupling')
        strength = row['s']
        require(strength not in seen, 'duplicate polynomial coupling')
        seen.add(strength)
        require(row.get('effective_matrix') == [1.0, 2*strength, 0.0, 1.0]
                and row.get('eigenvalues') == [1.0, 1.0]
                and row.get('numerical_range_center') == [1.0, 0.0]
                and row.get('numerical_range_radius') == abs(strength)
                and row.get('polynomial') == '1-z'
                and close(row.get('polynomial_norm'), 2*abs(strength))
                and row.get('polynomial_supremum') == abs(strength)
                and row.get('crouzeix_bound') == 2*abs(strength)
                and row.get('annihilating_polynomial') == '(1-z)^2'
                and row.get('squared_polynomial') == [0.0]*4,
                'incomplete or inconsistent polynomial evidence')
    require(isinstance(runs, list) and len(runs) == 20, 'all twenty trajectories are required')
    seen = set()
    for row in runs:
        require(isinstance(row, dict) and finite(row.get('s')) and row['s'] in couplings
                and type(row.get('restart')) is int and row['restart'] in (1, 2)
                and type(row.get('right_jacobi')) is bool, 'invalid trajectory combination')
        strength, restart, jacobi = row['s'], row['restart'], row['right_jacobi']
        key = strength, restart, jacobi
        require(key not in seen, 'duplicate trajectory combination')
        seen.add(key)
        diagonal = 8.0 if jacobi else 1.0
        require(row.get('matrix') == [1.0, 2*strength*diagonal, 0.0, diagonal]
                and row.get('rhs') == [0.0, 1.0]
                and row.get('effective_matrix') == [1.0, 2*strength, 0.0, 1.0]
                and row.get('exact_solution') == [-2*strength, 1/diagonal],
                'inconsistent original or preconditioned system')
        iterations, reason, history = row.get('iterations'), row.get('reason'), row.get('history')
        require(type(iterations) is int and 1 <= iterations <= limit
                and reason in ('converged', 'iteration_limit', 'stagnation'), 'invalid termination record')
        require(reason != 'iteration_limit' or iterations == limit, 'incomplete limited trajectory')
        require(restart != 2 or (reason == 'converged' and iterations <= 2),
                'restart two must solve the bounded system')
        require(row.get('restarts') == list(range(restart, iterations, restart)),
                'incomplete restart boundaries')
        require(isinstance(history, list) and len(history) == iterations + 1,
                'incomplete trajectory history')
        for index, frame in enumerate(history):
            require(isinstance(frame, dict) and type(frame.get('iteration')) is int
                    and frame['iteration'] == index and vector(frame.get('x'), 2),
                    'invalid captured iterate')
            cycle = ((index - 1)//restart)*restart if index else 0
            require(type(frame.get('cycle_start')) is int and frame['cycle_start'] == cycle
                    and type(frame.get('local_step')) is int and frame['local_step'] == index - cycle,
                    'invalid cycle position')
            for name in ('true_residual', 'estimated_residual', 'cycle_envelope', 'rounding_allowance'):
                require(finite(frame.get(name)) and frame[name] >= 0, 'invalid ' + name)
            x, y = frame['x']
            scale = 1 + abs(x) + abs(2*strength*diagonal*y) + abs(diagonal*y)
            allowance = 256*sys.float_info.epsilon*scale
            require(finite(allowance) and close(frame['rounding_allowance'], allowance,
                                               8*sys.float_info.epsilon*allowance),
                    'incorrect rounding allowance')
            residual = math.hypot(-x - 2*strength*diagonal*y, 1 - diagonal*y)
            factor = min(1.0, 2*abs(strength)**(index - cycle)) if index else 1.0
            envelope = history[cycle]['true_residual']*factor
            require(close(frame['true_residual'], residual, allowance)
                    and close(frame['cycle_envelope'], envelope, allowance)
                    and residual <= envelope + allowance, 'inconsistent residual or cycle envelope')
        require(history[0]['x'] == [0.0, 0.0]
                and history[0]['true_residual'] == history[0]['estimated_residual'] == 1.0,
                'invalid initial state')
        final = history[-1]
        require(reason != 'converged' or final['true_residual'] <= 1e-10 + final['rounding_allowance'],
                'false convergence')


    graphs = data.get('graphs')
    require(isinstance(graphs, dict) and type(graphs.get('schema_version')) is int
            and graphs['schema_version'] == 1 and graphs.get('suite') == 'shifted-graph-diagnostics-v1'
            and type(graphs.get('max_iterations')) is int and graphs['max_iterations'] == 2000,
            'missing or invalid graph diagnostics')
    expected_graphs = {'petersen': 10, 'triangular-prism': 6, 'prism-8': 16,
                       'prism-16': 32, 'prism-32': 64}
    cases = graphs.get('cases')
    require(isinstance(cases, list) and len(cases) == len(expected_graphs),
            'all five graph cases are required')
    seen = set()
    for case in cases:
        require(isinstance(case, dict) and isinstance(case.get('id'), str)
                and case['id'] in expected_graphs and case['id'] not in seen,
                'invalid or duplicate graph case')
        seen.add(case['id'])
        n = expected_graphs[case['id']]
        require(type(case.get('n')) is int and case['n'] == n, 'invalid graph size')
        solves = case.get('solves')
        require(isinstance(solves, list) and len(solves) == 2, 'both graph solves are required')
        solvers = set()
        for solve in solves:
            require(isinstance(solve, dict) and isinstance(solve.get('solver'), str)
                    and solve['solver'] in ('cg', 'cg_jacobi') and solve['solver'] not in solvers,
                    'invalid or duplicate graph solver')
            solvers.add(solve['solver'])
            require(type(solve.get('jacobi')) is bool
                    and solve['jacobi'] == (solve['solver'] == 'cg_jacobi')
                    and solve.get('converged') is True and solve.get('reason') == 'converged',
                    'graph solves must converge')
            iterations, history = solve.get('iterations'), solve.get('history')
            require(type(iterations) is int and 0 <= iterations <= graphs['max_iterations']
                    and isinstance(history, list) and len(history) == iterations + 1
                    and vector(solve.get('residuals'), iterations + 1) and vector(solve.get('x'), n),
                    'incomplete graph solve history')
            for name in ('true_residual', 'relative_true_residual', 'solution_error',
                         'relative_solution_error', 'stopping_threshold', 'rounding_allowance'):
                require(finite(solve.get(name)) and solve[name] >= 0, 'invalid graph ' + name)
            for index, frame in enumerate(history):
                require(isinstance(frame, dict) and type(frame.get('iteration')) is int
                        and frame['iteration'] == index, 'invalid graph history index')
                for name in ('true_residual', 'solver_residual', 'solution_error', 'rounding_allowance'):
                    require(finite(frame.get(name)) and frame[name] >= 0, 'invalid graph history ' + name)
                require(solve['residuals'][index] == frame['true_residual'],
                        'inconsistent graph residual history')
            final = history[-1]
            require(all(solve[name] == final[name] for name in
                        ('true_residual', 'solution_error', 'rounding_allowance')),
                    'inconsistent final graph diagnostics')
            require(solve['true_residual'] <= solve['stopping_threshold'] + solve['rounding_allowance'],
                    'graph solve reports false convergence')


def validate_completion(data, runner):
    """Each harness fails differently; require its completed verification record."""
    if runner == 'nonnormal-gmres':
        validate_nonnormal_completion(data)
        return
    if runner == 'machine-code-dot':
        # This experiment raises on any build, checksum or timing failure. It
        # writes the final report only after every implementation is measured.
        if (type(data.get('correctness_checks')) is not int or data['correctness_checks'] <= 0
                or not data.get('timings') or not data.get('samples')):
            raise ValueError('Dot experiment must complete correctness checks and timings')
        if any(len(row.get('ns_per_call', [])) != data['samples'] for row in data['timings']):
            raise ValueError('Dot experiment has incomplete timing samples')
        return
    if runner == 'matmul-locality':
        variants = {item['id'] for item in data.get('variants', [])}
        checks, results = data.get('checks', []), data.get('results', [])
        if (not variants or not checks or any(item.get('passed') is not True for item in checks)
                or {item.get('variant') for item in checks} != variants
                or not results or {item.get('variant') for item in results} != variants):
            raise ValueError('Locality experiment must verify and measure every requested variant')
        # Missing optional profiles are disclosed by this harness and do not
        # invalidate completed timings; failed correctness checks always do.
        return
    implementations = data.get('implementations', [])
    if not implementations or any(item.get('status') != 'passed' for item in implementations):
        raise ValueError('Every requested implementation must pass before publication')
    if any(row.get('status') != 'passed' for row in data.get('results', [])):
        raise ValueError('Failed measurements cannot be published')


def capture(output, runner='run.py', arguments=(), root=ROOT):
    """Do not modify the caller's checkout or reuse any of its build artifacts."""
    root, output = Path(root).resolve(), Path(output).resolve()
    arguments = checked_arguments(runner, arguments)
    if output.exists():
        raise ValueError('Choose a new output directory; existing reports are never overwritten')
    if git(root, 'status', '--porcelain', '--untracked-files=all'):
        raise ValueError('Publication requires a clean checkout. Commit source changes first; '
                         'use the ordinary benchmark runner for local development results.')
    revision = git(root, 'rev-parse', 'HEAD')
    source_tree = git(root, 'rev-parse', 'HEAD^{tree}')
    with tempfile.TemporaryDirectory(prefix='linear-a-benchmark-') as temporary:
        temporary = Path(temporary)
        checkout, staged = temporary / 'source', temporary / 'report'
        # A clone has no ignored build products, untracked sources, or local
        # worktree edits. --no-hardlinks also keeps its object storage independent.
        subprocess.run(['git', 'clone', '--quiet', '--no-hardlinks', '--no-checkout',
                        str(root), str(checkout)], check=True)
        subprocess.run(['git', 'checkout', '--quiet', '--detach', revision],
                       cwd=checkout, check=True)
        if git(checkout, 'status', '--porcelain', '--untracked-files=all'):
            raise ValueError('Fresh benchmark checkout is not clean')
        runner_path = EXPERIMENT_RUNNERS.get(runner, 'benchmarks/' + runner)
        command = [runner_path, *arguments, '--output', str(staged)]
        subprocess.run([sys.executable, *command], cwd=checkout, check=True)
        if (git(checkout, 'rev-parse', 'HEAD') != revision or
                git(checkout, 'status', '--porcelain', '--untracked-files=all')):
            raise ValueError('Benchmark changed its source checkout; report was not published')
        data = json.loads((staged / 'results.json').read_text(encoding='utf-8'))
        if data.get('revision') != revision or data.get('dirty') is not False:
            raise ValueError('Benchmark did not record the clean source revision')
        validate_completion(data, runner)
        provenance = {
            'schema_version': 1,
            'capture_method': 'clean-checkout-v1',
            'revision': revision,
            'source_tree': source_tree,
            'command': ['python3', *command[:-1], '<output>'],
            'python_version': platform.python_version(),
            'artifacts': artifact_hashes(staged),
        }
        (staged / 'provenance.json').write_text(
            json.dumps(provenance, indent=2, allow_nan=False) + '\n', encoding='utf-8')
        validate_report(staged, repo_root=root)
        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(staged, output)
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument('--runner', choices=RUNNERS, default='run.py')
    parser.add_argument('--output', type=Path, required=True, help='New report directory')
    parser.add_argument('arguments', nargs=argparse.REMAINDER,
                        help='Benchmark options after --; fresh builds and --require-all are enforced')
    args = parser.parse_args()
    arguments = args.arguments[1:] if args.arguments[:1] == ['--'] else args.arguments
    try:
        provenance = capture(args.output, args.runner, arguments)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'Capture failed: {error}\n')
    print(f"Captured {provenance['revision']} in {args.output}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
