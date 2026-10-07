#!/usr/bin/env python3
"""Check nonnormal GMRES trajectories against small independent analytic oracles."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'ports/python'))
from gmres import gmres
from matrix import Matrix
from sparse import CSRMatrix

COUPLINGS = (0.0, 0.25, 0.75, 2.0, 4.0)
EPSILON = sys.float_info.epsilon


def require_close(actual, expected, tolerance, label):
    if not math.isfinite(actual) or abs(actual - expected) > tolerance:
        raise AssertionError('%s: got %r, expected %r ± %r' %
                             (label, actual, expected, tolerance))


def exact_restarted_steps(s, count=2):
    """Exact minimal-residual line searches, independent of Arnoldi and Givens.

    Only two steps are needed: repeated rational steps grow integers rapidly.
    The Krylov line for GMRES(1) is x + alpha*r for the effective operator B.
    """
    t = 2 * Fraction(s)
    x = [Fraction(0), Fraction(0)]
    r = [Fraction(0), Fraction(1)]
    frames = [x.copy()]
    for _ in range(count):
        br = [r[0] + t*r[1], r[1]]
        denominator = sum(v*v for v in br)
        if not denominator:
            break
        alpha = sum(a*b for a, b in zip(r, br)) / denominator
        x = [a + alpha*b for a, b in zip(x, r)]
        r = [a - alpha*b for a, b in zip(r, br)]
        frames.append(x.copy())
    return [[float(v) for v in frame] for frame in frames]


def polynomial_checks(s):
    """Check the spectrum and sharp polynomial norm using existing dense APIs."""
    b = Matrix(2, 2, [1.0, 2*s, 0.0, 1.0])
    p = Matrix.identity(2).subtract(b)
    measured = p.svd().values[0]
    require_close(measured, 2*abs(s), 64*EPSILON*max(1, 2*abs(s)),
                  'SVD norm of I-B')
    square = p.multiply(p)
    if square.values != [0.0]*4:
        raise AssertionError('(I-B)^2 must vanish for this dyadic family')
    eigen = b.eigen_general()
    for value in eigen.values_real:
        require_close(value, 1.0, 64*EPSILON, 'eigenvalue real part')
    for value in eigen.values_imag:
        require_close(value, 0.0, 64*EPSILON, 'eigenvalue imaginary part')
    return dict(s=s, effective_matrix=b.values, eigenvalues=[1.0, 1.0],
                numerical_range_center=[1.0, 0.0], numerical_range_radius=abs(s),
                polynomial='1-z', polynomial_norm=measured,
                polynomial_supremum=abs(s), crouzeix_bound=2*abs(s),
                annihilating_polynomial='(1-z)^2', squared_polynomial=square.values)


def trajectory(s, restart, jacobi, limit=24):
    """Run the shipped solver; check physical residuals and first projections."""
    if s not in COUPLINGS or restart not in (1, 2) or type(jacobi) is not bool:
        raise ValueError('use a documented coupling, restart 1 or 2, and Boolean Jacobi')
    if type(limit) is not int or not 2 <= limit <= 100:
        raise ValueError('iteration limit must be an integer from 2 through 100')
    diagonal = 8.0 if jacobi else 1.0
    # A=B*M; the right-Jacobi effective operator A*M^-1 is always B.
    a = CSRMatrix(2, 2, [0, 2, 3], [0, 1, 1], [1.0, 2*s*diagonal, diagonal])
    result = gmres(a, [0.0, 1.0], restart=restart, rtol=1e-10,
                   max_iterations=limit, jacobi=jacobi, capture=True)
    if result.reason not in ('converged', 'iteration_limit', 'stagnation'):
        raise AssertionError('unexpected GMRES termination: ' + result.reason)
    if len(result.iterates) != result.iterations + 1:
        raise AssertionError('missing captured iterates')
    expected_starts = list(range(restart, result.iterations, restart))
    if result.restarts != expected_starts:
        raise AssertionError('incorrect restart boundaries')
    exact = exact_restarted_steps(s)
    # Both methods use the same first Krylov line; only GMRES(1) restarts after it.
    for k in range(1, min(len(exact), len(result.iterates), 3 if restart == 1 else 2)):
        expected = [exact[k][0], exact[k][1]/diagonal]
        for actual, target in zip(result.iterates[k], expected):
            require_close(actual, target, 256*EPSILON*max(1, abs(target)),
                          'exact minimal-residual projection at step %d' % k)
    if restart == 2 and (not result.converged or result.iterations > 2):
        raise AssertionError('degree-two annihilator must solve this bounded fixture')

    history = []
    for k, x in enumerate(result.iterates):
        # Independent two-entry expression, not CSR matvec or solver histories.
        residual = math.hypot(-x[0] - 2*s*diagonal*x[1], 1 - diagonal*x[1])
        scale = 1 + abs(x[0]) + abs(2*s*diagonal*x[1]) + abs(diagonal*x[1])
        tolerance = 256*EPSILON*scale
        require_close(result.residuals[k], residual, tolerance, 'captured true residual')
        cycle_start = ((k-1)//restart)*restart if k else 0
        local_step = k-cycle_start
        # Candidate p(z)=(1-z)^j, plus the no-progress bound p=1. This is
        # a per-cycle upper bound, not a global unrestarted minimization.
        factor = min(1.0, 2*abs(s)**local_step) if k else 1.0
        envelope = result.residuals[cycle_start]*factor
        if residual > envelope + tolerance:
            raise AssertionError('cycle polynomial envelope exceeded beyond rounding allowance')
        history.append(dict(iteration=k, cycle_start=cycle_start, local_step=local_step,
                            x=x, true_residual=residual,
                            estimated_residual=result.estimated_residuals[k],
                            cycle_envelope=envelope, rounding_allowance=tolerance))
    if result.converged and history[-1]['true_residual'] > 1e-10 + history[-1]['rounding_allowance']:
        raise AssertionError('false convergence')
    return dict(s=s, restart=restart, right_jacobi=jacobi,
                matrix=[1.0, 2*s*diagonal, 0.0, diagonal], rhs=[0.0, 1.0],
                effective_matrix=[1.0, 2*s, 0.0, 1.0],
                exact_solution=[-2*s, 1/diagonal], iterations=result.iterations,
                reason=result.reason, restarts=result.restarts, history=history)


def source_record():
    paths = ['experiments/nonnormal-gmres/run.py', 'ports/python/gmres.py',
             'ports/python/sparse.py', 'ports/python/matrix.py',
             'ports/python/svd.py', 'ports/python/general_eigen.py']
    hashes = {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}
    def git(*args):
        try:
            result = subprocess.run(['git', *args], cwd=ROOT, text=True, capture_output=True)
        except OSError:
            return None
        return result.stdout.strip() if result.returncode == 0 else None
    status = git('status', '--porcelain')
    return dict(revision=git('rev-parse', 'HEAD'), dirty=None if status is None else bool(status),
                files_sha256=hashes,
                note='Local diagnostic run; hashes do not archive uncommitted source. No timings or publication capture.')


def experiment(limit=24):
    polynomials = [polynomial_checks(s) for s in COUPLINGS]
    runs = [trajectory(s, restart, jacobi, limit)
            for s in COUPLINGS for restart in (1, 2) for jacobi in (False, True)]
    # Powers-of-two column scaling is exact for this family. Compare residual
    # histories only after each run independently passes its analytic checks.
    for plain, scaled in zip(runs[::2], runs[1::2]):
        if plain['iterations'] != scaled['iterations'] or plain['reason'] != scaled['reason']:
            raise AssertionError('right Jacobi did not recover the same effective problem')
        for a, b in zip(plain['history'], scaled['history']):
            require_close(a['true_residual'], b['true_residual'],
                          max(a['rounding_allowance'], b['rounding_allowance']),
                          'equivalent preconditioned residual history')
    return dict(schema_version=1, suite='nonnormal-gmres-v1', implementation='python',
                numeric_type='IEEE 754 binary64', max_iterations=limit,
                rtol=1e-10, source=source_record(), polynomials=polynomials, runs=runs)


def render_markdown(data):
    lines = ['# Nonnormal GMRES experiment', '',
             'Local correctness and convergence diagnostics; no performance measurements.', '',
             'All effective matrices B = [[1, 2s], [0, 1]] have eigenvalues 1, 1. '
             'Their numerical ranges are disks centered at 1 with radius |s|. '
             'The SVD independently checks ||I-B||₂ = 2|s|, attaining the Crouzeix constant for nonzero s.', '',
             '| s | SVD norm of I−B | Supremum of |1−z| on W(B) |',
             '| --- | ---: | ---: |']
    for row in data['polynomials']:
        lines.append('| %g | %.6g | %.6g |' %
                     (row['s'], row['polynomial_norm'], row['polynomial_supremum']))
    lines += ['', '| s | Right Jacobi | Restart | Steps | Stop | Final true residual / initial |',
              '| --- | --- | ---: | ---: | --- | ---: |']
    for row in data['runs']:
        lines.append('| %g | %s | %d | %d | %s | %.6g |' %
                     (row['s'], 'yes' if row['right_jacobi'] else 'no', row['restart'],
                      row['iterations'], row['reason'], row['history'][-1]['true_residual']))
    lines += ['', '## Interpretation', '',
              '- Restart 2 retains the degree-two annihilating polynomial (1−z)²; '
              'the bounded test systems solve in at most two steps.',
              '- Restart 1 discards that Krylov space. A common spectrum does not imply '
              'common convergence. Iteration-limit and stagnation results are retained.',
              '- Right Jacobi uses A = B·diag(1,8), so the effective operator is '
              'A·diag(1,8)⁻¹ = B. Residuals are checked in the original system b−Ax.',
              '- Each captured step is checked against the per-cycle candidate bound '
              '||r_start||·min(1, 2|s|ʲ), where j is the step within that cycle. '
              'This is an exact-arithmetic comparison with an explicit scale-dependent '
              'rounding allowance, not a certified error bound or stopping criterion. '
              'For |s| ≥ 1 the disk contains zero and this envelope is uninformative.',
              '- The first projections are independently checked by exact rational '
              'minimal-residual line searches. Full trajectories, bounds, allowances, '
              'inputs and source hashes are in results.json.', '',
              'Source state: `%s`, dirty: `%s`. %s' %
              (data['source']['revision'], data['source']['dirty'], data['source']['note']), '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'.build/nonnormal-gmres')
    parser.add_argument('--limit', type=int, default=24)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    if not 2 <= args.limit <= 100:
        parser.error('--limit must be from 2 through 100')
    data = experiment(args.limit)
    if not args.verify_only:
        args.output.mkdir(parents=True, exist_ok=True)
        (args.output/'results.json').write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
        (args.output/'RESULTS.md').write_text(render_markdown(data))
    print('Passed %d polynomial checks and %d GMRES trajectories.' %
          (len(data['polynomials']), len(data['runs'])))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
