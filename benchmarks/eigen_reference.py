"""Analytic symmetric spectra and eigensystem checks independent of Jacobi iteration."""
from fractions import Fraction
import math

from reference import finite_number, matrix


def generated_eigen(size, seed):
    """Dyadic symmetric Toeplitz input, with a known sine eigenbasis."""
    if type(size) is not int or size < 0:
        raise ValueError("Size must be a nonnegative integer")
    diagonal = 2 + (seed % 17) / 16.0
    return matrix(size, size, [diagonal if i == j else -1.0 if abs(i-j) == 1 else 0.0
                              for i in range(size) for j in range(size)])


def eigen_spectrum(size, seed):
    """Ascending eigenvalues of generated_eigen; no iterative solver is used."""
    if type(size) is not int or size < 0:
        raise ValueError("Size must be a nonnegative integer")
    diagonal = 2 + (seed % 17) / 16.0
    return [diagonal - 2 * math.cos(k * math.pi / (size + 1)) for k in range(1, size + 1)]


def eigen_checksum(size, seed):
    """Weighted spectrum plus the squared Frobenius norm of orthonormal columns."""
    return math.fsum((i + 1) * value for i, value in enumerate(eigen_spectrum(size, seed))) + size


def _householder_matrix(spectrum, direction):
    """Construct Q D Q^T exactly over rationals before converting the input to double."""
    n = len(spectrum)
    norm = sum(value * value for value in direction)
    q = [[Fraction(int(i == j)) - Fraction(2 * direction[i] * direction[j], norm)
          for j in range(n)] for i in range(n)]
    values = [float(sum(q[i][k] * Fraction(spectrum[k]) * q[j][k] for k in range(n)))
              for i in range(n) for j in range(n)]
    return matrix(n, n, values)


def eigen_fixtures():
    """Shared valid and invalid eigen_symmetric inputs with independently known spectra."""
    cases = []

    def add(name, a, spectrum=None, invalid=False):
        case = {"name": name, "op": "eigen_symmetric", "a": a, "b": None, "invalid": invalid}
        if not invalid:
            case["expected_eigenvalues"] = sorted(spectrum)
        cases.append(case)

    add("empty symmetric eigensystem", matrix(0, 0, []), [])
    add("scalar symmetric eigensystem", matrix(1, 1, [-7.25]), [-7.25])
    add("zero symmetric eigensystem", matrix(3, 3, [0] * 9), [0, 0, 0])
    add("unsorted repeated diagonal spectrum", matrix(4, 4,
        [3, 0, 0, 0, 0, -2, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0]), [-2, 0, 3, 3])
    add("equal diagonal rotation", matrix(2, 2, [2, 1, 1, 2]), [1, 3])
    add("negative off diagonal rotation", matrix(2, 2, [2, -1, -1, 2]), [1, 3])
    add("indefinite symmetric eigensystem", matrix(2, 2, [0, 1, 1, 0]), [-1, 1])
    add("fractional symmetric eigensystem", matrix(2, 2, [.5, .125, .125, .5]), [.375, .625])
    add("singular symmetric eigensystem", matrix(2, 2, [1, 1, 1, 1]), [0, 2])
    near = 1 + 1e-12
    center = (1 + near) / 2
    radius = math.hypot((near - 1) / 2, 1e-8)
    add("nearly equal diagonal entries", matrix(2, 2, [1, 1e-8, 1e-8, near]),
        [center - radius, center + radius])
    add("dense repeated eigenspace", _householder_matrix([-2, -2, 7], [1, 1, 1]), [-2, -2, 7])
    add("dense mixed spectrum", _householder_matrix([-5, -.5, 0, 4], [1, 2, -2, 1]), [-5, -.5, 0, 4])
    for name, scale in (("tiny", 1e-200), ("huge", 1e200)):
        add(name + " symmetric eigensystem", matrix(2, 2, [2*scale, scale, scale, 2*scale]),
            [scale, 3*scale])
    add("mixed scale exact diagonal", matrix(2, 2, [1e300, 0, 0, 1e-300]), [1e-300, 1e300])
    add("subnormal scalar eigensystem", matrix(1, 1, [5e-324]), [5e-324])
    for size, seed in ((3, 0), (5, 9), (8, 17), (16, 9)):
        add("analytic Toeplitz spectrum %d" % size, generated_eigen(size, seed), eigen_spectrum(size, seed))
    add("nonsquare eigensystem", matrix(2, 3, [1, 2, 3, 4, 5, 6]), invalid=True)
    add("nonsymmetric eigensystem", matrix(2, 2, [1, 2, 0, 1]), invalid=True)
    add("tiny asymmetry is rejected", matrix(2, 2, [1, 5e-324, 0, 1]), invalid=True)
    return cases


def assert_eigen_result(actual, a, expected_eigenvalues=None, tolerance=1e-9):
    """Check A Q = Q Λ and Q^T Q = I without selecting signs or an eigenspace basis.

    Residuals are evaluated after scaling by the largest input/output eigenvalue,
    avoiding overflow and preventing an absolute tolerance from hiding tiny errors.
    The supplied spectrum, when present, is an additional independent check.
    """
    n = a["rows"]
    if a["cols"] != n:
        raise AssertionError("Eigenvalue reference input must be square")
    if not isinstance(actual, dict) or set(actual) != {"eigenvalues", "eigenvectors"}:
        raise AssertionError("Expected eigenvalues and eigenvectors result keys")
    eigenvalues = actual["eigenvalues"]
    vectors = actual["eigenvectors"]
    if not isinstance(eigenvalues, list) or len(eigenvalues) != n or not all(map(finite_number, eigenvalues)):
        raise AssertionError("Eigenvalues must be a finite array of length n")
    if any(left > right for left, right in zip(eigenvalues, eigenvalues[1:])):
        raise AssertionError("Eigenvalues must be sorted in ascending order")
    if not isinstance(vectors, dict) or set(vectors) != {"rows", "cols", "values"}:
        raise AssertionError("Eigenvectors must be a matrix object")
    if type(vectors["rows"]) is not int or type(vectors["cols"]) is not int or vectors["rows"] != n or vectors["cols"] != n:
        raise AssertionError("Eigenvectors must have shape n by n")
    q = vectors["values"]
    if not isinstance(q, list) or len(q) != n*n or not all(map(finite_number, q)):
        raise AssertionError("Eigenvectors must contain n squared finite entries")
    if any(abs(value) > 1 + tolerance*n for value in q):
        raise AssertionError("Eigenvectors are not normalized")
    orthogonality = math.sqrt(math.fsum(
        (math.fsum(q[k*n+i]*q[k*n+j] for k in range(n)) - int(i == j))**2
        for i in range(n) for j in range(n)))
    if orthogonality > tolerance * n:
        raise AssertionError("Eigenvectors are not orthonormal: Frobenius error %g" % orthogonality)
    scale = max([abs(value) for value in a["values"]] + [abs(value) for value in eigenvalues] + [0.0])
    if scale:
        scaled = [value/scale for value in a["values"]]
        lambdas = [value/scale for value in eigenvalues]
        norm = math.sqrt(math.fsum(value*value for value in scaled))
        residual = math.sqrt(math.fsum(
            (math.fsum(scaled[i*n+k]*q[k*n+j] for k in range(n)) - q[i*n+j]*lambdas[j])**2
            for i in range(n) for j in range(n)))
        if residual > tolerance*norm:
            raise AssertionError("Eigenpair residual exceeds tolerance: %g / %g" % (residual, norm))
        relative_residual = residual/norm if norm else 0.0
    else:
        norm, relative_residual = 0.0, 0.0
    if expected_eigenvalues is not None:
        if len(expected_eigenvalues) != n:
            raise AssertionError("Reference spectrum has incorrect length")
        for index, (value, wanted) in enumerate(zip(eigenvalues, expected_eigenvalues)):
            # Use a normwise bound: near-zero eigenvalues need not have small relative error.
            bound = tolerance * norm if scale else 0.0
            error = abs(value/scale - wanted/scale) if scale else abs(value-wanted)
            if error > bound:
                raise AssertionError("Eigenvalue %d differs from the independent spectrum" % index)
    return {"relative_residual": relative_residual, "orthogonality_error": orthogonality}
