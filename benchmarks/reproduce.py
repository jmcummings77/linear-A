#!/usr/bin/env python3
"""Capture publishable benchmarks from a clean commit in a fresh local clone."""
import argparse
import json
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
    'ordering_bench.py', 'cholesky_bench.py', 'amd_bench.py', 'ic0_bench.py',
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
    # IC(0) and the experiments already fail unless every selected
    # implementation passes; their parsers do not expose --require-all.
    if runner not in ('ic0_bench.py', *EXPERIMENT_RUNNERS) and '--require-all' not in arguments:
        arguments.append('--require-all')
    return arguments


def validate_completion(data, runner):
    """Each harness fails differently; require its completed verification record."""
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
