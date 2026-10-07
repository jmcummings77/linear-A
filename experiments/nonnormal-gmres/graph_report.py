"""Measured Python CG diagnostics on classical graphs with known analytic spectra.

This module produces no timings and does not construct a new expander family.
Its file-based imports keep benchmark and library modules named ``sparse`` apart.
"""
import importlib.util
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
EPSILON = sys.float_info.epsilon
RTOL = 1e-10
MAX_ITERATIONS = 2000


def _load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


graph_reference = _load_module('nonnormal_graph_reference', 'benchmarks/graph_reference.py')
CSRMatrix = _load_module('nonnormal_graph_python_sparse', 'ports/python/sparse.py').CSRMatrix


def _norm(vector):
    return math.hypot(*vector)


def _require_close(actual, expected, allowance, label):
    if not math.isfinite(actual) or not math.isfinite(expected) or abs(actual - expected) > allowance:
        raise AssertionError('%s: got %r, expected %r within %r' %
                             (label, actual, expected, allowance))


def _solve_case(matrix, n, edges, truth, rhs, shift, jacobi):
    result = matrix.conjugate_gradient(rhs, rtol=RTOL, atol=0.0,
                                       max_iterations=MAX_ITERATIONS,
                                       jacobi=jacobi, capture=True)
    if not result.converged or result.reason != 'converged':
        raise AssertionError('graph CG unexpectedly failed to converge: ' + result.reason)
    if not 0 <= result.iterations <= MAX_ITERATIONS:
        raise AssertionError('invalid graph CG iteration count')
    if len(result.iterates) != result.iterations + 1 or len(result.residuals) != len(result.iterates):
        raise AssertionError('missing graph CG history')
    if result.iterates[0] != [0.0] * n or result.iterates[-1] != result.x:
        raise AssertionError('incorrect graph CG captured endpoints')

    rhs_norm = _norm(rhs)
    truth_norm = _norm(truth)
    threshold = RTOL * rhs_norm
    history = []
    for iteration, (x, reported) in enumerate(zip(result.iterates, result.residuals)):
        if len(x) != n or not all(math.isfinite(value) for value in x):
            raise AssertionError('invalid graph CG iterate')
        # Edge differences use no CSR multiplication or solver residual history.
        product = graph_reference.graph_product(n, edges, x, shift)
        residual = _norm([b - ax for b, ax in zip(rhs, product)])
        error = _norm([value - expected for value, expected in zip(x, truth)])
        # A bounded diagnostic allowance, not a certified floating-point theorem.
        scale = max(1.0, rhs_norm + (shift + 6.0) * _norm(x))
        allowance = 256 * EPSILON * scale
        _require_close(reported, residual, allowance, 'independent graph true residual')
        # The full shifted Laplacian has lambda_min = shift; this also checks
        # the known solution without relying on a second numerical solver.
        error_allowance = (residual + allowance) / shift + 256 * EPSILON * truth_norm
        if not math.isfinite(error) or error > error_allowance:
            raise AssertionError('graph solution disagrees with the manufactured answer')
        history.append(dict(iteration=iteration, true_residual=residual,
                            solver_residual=reported, solution_error=error,
                            rounding_allowance=allowance))

    final = history[-1]
    if final['true_residual'] > threshold + final['rounding_allowance']:
        raise AssertionError('graph CG reported false convergence')
    return dict(solver='cg_jacobi' if jacobi else 'cg', jacobi=jacobi,
                reason=result.reason, converged=result.converged,
                iterations=result.iterations, x=list(result.x),
                true_residual=final['true_residual'],
                relative_true_residual=final['true_residual'] / rhs_norm,
                solution_error=final['solution_error'],
                relative_solution_error=final['solution_error'] / truth_norm,
                stopping_threshold=threshold,
                rounding_allowance=final['rounding_allowance'],
                residuals=[row['true_residual'] for row in history], history=history)


def _graph_case(identifier, name, cycle_size=None):
    if cycle_size is None:
        n, edges = graph_reference.petersen_graph()
        spectrum = [-2.0] * 4 + [1.0] * 5 + [3.0]
        spectrum_source = 'Exact Petersen spectrum: -2 (multiplicity 4), 1 (5), 3 (1).'
        bipartite = False
    else:
        n, edges = graph_reference.prism_graph(cycle_size)
        spectrum = graph_reference.prism_spectrum(cycle_size)
        spectrum_source = 'Known analytic prism spectrum: 2*cos(2*pi*k/m) +/- 1, k=0,...,m-1; trig evaluated in binary64.'
        bipartite = cycle_size % 2 == 0
    degree = 3
    degrees = [0] * n
    for u, v in edges:
        degrees[u] += 1
        degrees[v] += 1
    if degrees != [degree] * n:
        raise AssertionError('graph fixture is not cubic')
    if len(spectrum) != n or spectrum != sorted(spectrum):
        raise AssertionError('invalid analytic adjacency spectrum')

    shift = graph_reference.SHIFT
    csr = graph_reference.shifted_laplacian(n, edges, shift)
    matrix = CSRMatrix(n, n, csr['offsets'], csr['indices'], csr['values'])
    if matrix.nnz != 4 * n:
        raise AssertionError('unexpected shifted cubic graph storage')
    truth = graph_reference.manufactured_solution(n)
    rhs = graph_reference.graph_product(n, edges, truth, shift)
    if matrix.matvec(truth) != rhs:
        raise AssertionError('CSR action disagrees with independent graph edge action')
    energy = graph_reference.graph_energy(n, edges, truth, shift)
    energy_from_rhs = math.fsum(x * b for x, b in zip(truth, rhs))
    energy_allowance = 256 * EPSILON * max(1.0, abs(energy))
    _require_close(energy_from_rhs, energy, energy_allowance, 'graph edge energy identity')
    if energy <= 0:
        raise AssertionError('nonpositive graph energy for a nonzero solution')
    constant_action = graph_reference.graph_product(n, edges, [1.0] * n, shift)
    if constant_action != [shift] * n:
        raise AssertionError('incorrect shifted graph constant mode')

    # Connectivity gives lambda_max(H)=3 and lambda_min(shift*I+L)=shift.
    # The algebraic gap is lambda_2(L), not a gap that removes both +/-3.
    spectral_gap = degree - spectrum[-2]
    shifted_min = shift
    shifted_max = degree + shift - spectrum[0]
    # Bipartite Ramanujan graphs allow both trivial eigenvalues +3 and -3.
    nontrivial = spectrum[1:-1] if bipartite else spectrum[:-1]
    nontrivial_abs_max = max(abs(value) for value in nontrivial)
    ramanujan_threshold = 2 * math.sqrt(degree - 1)
    solves = [_solve_case(matrix, n, edges, truth, rhs, shift, jacobi)
              for jacobi in (False, True)]
    return dict(id=identifier, name=name, cycle_size=cycle_size, n=n,
                edge_count=len(edges), nnz=matrix.nnz, degree=degree, shift=shift,
                edges=[list(edge) for edge in edges], csr=csr,
                adjacency_spectrum=spectrum, spectrum_source=spectrum_source,
                adjacency_min=spectrum[0], adjacency_max=spectrum[-1],
                adjacency_second_largest=spectrum[-2], spectral_gap=spectral_gap,
                spectral_gap_definition='3 - second-largest adjacency eigenvalue = first positive unshifted Laplacian eigenvalue',
                bipartite=bipartite, nontrivial_adjacency_abs_max=nontrivial_abs_max,
                ramanujan_threshold=ramanujan_threshold,
                is_ramanujan=nontrivial_abs_max <= ramanujan_threshold,
                shifted_lambda_min=shifted_min, shifted_lambda_max=shifted_max,
                kappa2=shifted_max / shifted_min,
                conditioning_scope='Full shifted Laplacian, including the constant mode; not restricted to mean-zero vectors.',
                jacobi_diagonal=degree + shift,
                manufactured_solution=truth, rhs=rhs,
                energy=energy, energy_from_rhs=energy_from_rhs,
                energy_rounding_allowance=energy_allowance,
                spd_identity='x^T A x = shift*sum(x_i^2) + sum_edges((x_u-x_v)^2); shift > 0 makes A SPD.',
                solves=solves)


def graph_diagnostics():
    """Return JSON-serializable graph facts and checked, untimed CG outcomes."""
    cases = [_graph_case('petersen', 'Petersen'),
             _graph_case('triangular-prism', 'Triangular prism', 3)]
    cases.extend(_graph_case('prism-%d' % m, 'Prism C_%d x K_2' % m, m)
                 for m in (8, 16, 32))
    return dict(schema_version=1, suite='shifted-graph-diagnostics-v1',
                numeric_type='IEEE 754 binary64', shift=graph_reference.SHIFT,
                rtol=RTOL, atol=0.0, max_iterations=MAX_ITERATIONS,
                scope=dict(
                    spectrum='Known analytic spectra, not measured eigensolver outputs or a new spectral certification.',
                    graphs='Classical controls. Petersen and triangular prism are nonbipartite Ramanujan; the even prisms are bipartite. Prism 8 is also Ramanujan, but arbitrarily large prisms are not a Ramanujan family.',
                    conditioning='The full shifted-Laplacian condition number includes the constant eigenvector with eigenvalue shift.',
                    jacobi='Every diagonal equals 3 + shift. Jacobi is scalar scaling and gives no exact-arithmetic condition-number improvement; floating-point trajectories may differ.',
                    validation='Physical residuals and solution errors are checked using an independent edge-difference action and a deterministic dyadic manufactured solution. Roundoff allowances are bounded diagnostics, not certified error bounds.',
                    timing='No timings or performance conclusions.'),
                cases=cases)
