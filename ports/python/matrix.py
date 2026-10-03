"""A small, dependency-free float64 matrix with row-major storage."""

import math
from typing import NamedTuple
from general_eigen import general_eigen


class GeneralEigenDecomposition(NamedTuple):
    """Complex eigenpairs as real/imaginary parts; right vectors are unit columns."""
    values_real: list
    values_imag: list
    vectors_real: "Matrix"
    vectors_imag: "Matrix"


class SymmetricEigenDecomposition(NamedTuple):
    """Ascending eigenvalues and a matrix whose corresponding columns are eigenvectors."""
    values: list
    vectors: "Matrix"


def finite(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("matrix arithmetic requires finite float64 values")
    return value


def _diagonal_product(values, size, exponent=0, sign=1.0, power=1):
    """Accumulate binary exponents separately; round to float64 only at the end."""
    mantissa = sign
    for index in range(size):
        value = values[index * size + index]
        if value == 0.0:
            return 0.0
        part, shift = math.frexp(value)
        for _ in range(power):
            mantissa, carry = math.frexp(mantissa * part)
            exponent += shift + carry
    try:
        return finite(math.ldexp(mantissa, exponent))
    except OverflowError as error:
        raise ValueError("determinant exceeds the finite float64 range") from error


def _quotient_product(numerator, denominator, value):
    """Avoid rounding an extremely small elimination factor before multiplying."""
    left, left_exponent = math.frexp(numerator)
    right, right_exponent = math.frexp(denominator)
    item, item_exponent = math.frexp(value)
    return math.ldexp((left / right) * item, left_exponent - right_exponent + item_exponent)


class Matrix:
    def __init__(self, rows=0, cols=0, values=None):
        if type(rows) is not int or type(cols) is not int or rows < 0 or cols < 0:
            raise ValueError("dimensions must be nonnegative integers")
        self._rows, self._cols = rows, cols
        self._values = [0.0] * (rows * cols) if values is None else [finite(x) for x in values]
        if len(self._values) != rows * cols:
            raise ValueError("value count does not match dimensions")

    @property
    def rows(self):
        return self._rows

    @property
    def cols(self):
        return self._cols

    @property
    def values(self):
        return self._values.copy()

    @classmethod
    def identity(cls, size):
        result = cls(size, size)
        for i in range(size):
            result._values[i * size + i] = 1.0
        return result

    def copy(self):
        return Matrix(self.rows, self.cols, self._values)

    def _index(self, row, col):
        if type(row) is not int or type(col) is not int or not 0 <= row < self.rows or not 0 <= col < self.cols:
            raise IndexError("matrix index out of range")
        return row * self.cols + col

    def __getitem__(self, index):
        return self._values[self._index(*index)]

    def __setitem__(self, index, value):
        position = self._index(*index)
        self._values[position] = finite(value)

    def row(self, index):
        if type(index) is not int or not 0 <= index < self.rows:
            raise IndexError("row index out of range")
        return self._values[index * self.cols:(index + 1) * self.cols]

    def column(self, index):
        if type(index) is not int or not 0 <= index < self.cols:
            raise IndexError("column index out of range")
        return [self._values[row * self.cols + index] for row in range(self.rows)]

    def _same_shape(self, other):
        if not isinstance(other, Matrix) or (self.rows, self.cols) != (other.rows, other.cols):
            raise ValueError("matrix dimensions must match")

    def add(self, other):
        self._same_shape(other)
        return Matrix(self.rows, self.cols, [a + b for a, b in zip(self._values, other._values)])

    def subtract(self, other):
        self._same_shape(other)
        return Matrix(self.rows, self.cols, [a - b for a, b in zip(self._values, other._values)])

    def scale(self, scalar):
        scalar = finite(scalar)
        return Matrix(self.rows, self.cols, [value * scalar for value in self._values])

    def transpose(self):
        values = [0.0] * len(self._values)
        for row in range(self.rows):
            for col in range(self.cols):
                values[col * self.rows + row] = self._values[row * self.cols + col]
        return Matrix(self.cols, self.rows, values)

    def multiply(self, other):
        if not isinstance(other, Matrix) or self.cols != other.rows:
            raise ValueError("left columns must equal right rows")
        values = [0.0] * (self.rows * other.cols)
        for row in range(self.rows):
            for col in range(other.cols):
                total = 0.0
                for k in range(self.cols):
                    total += self._values[row * self.cols + k] * other._values[k * other.cols + col]
                values[row * other.cols + col] = total
        return Matrix(self.rows, other.cols, values)

    def _vector3(self):
        if (self.rows, self.cols) not in ((3, 1), (1, 3)):
            raise ValueError("expected a three-component row or column vector")
        return [finite(value) for value in self._values]

    def cross(self, other):
        """Return the right-handed vector cross product, preserving the left shape."""
        if not isinstance(other, Matrix):
            raise ValueError("cross product requires two three-component vectors")
        left, right = self._vector3(), other._vector3()
        return Matrix(self.rows, self.cols,
                      [left[1]*right[2]-left[2]*right[1],
                       left[2]*right[0]-left[0]*right[2],
                       left[0]*right[1]-left[1]*right[0]])

    @staticmethod
    def rotation2d(radians):
        """Active counterclockwise 2D rotation acting on column vectors; angle in radians."""
        angle = finite(radians)
        cosine, sine = math.cos(angle), math.sin(angle)
        return Matrix(2, 2, [cosine, -sine, sine, cosine])

    @staticmethod
    def rotation_x(radians):
        """Active right-handed 3D rotation about X, acting on column vectors."""
        angle = finite(radians)
        cosine, sine = math.cos(angle), math.sin(angle)
        return Matrix(3, 3, [1,0,0,0,cosine,-sine,0,sine,cosine])

    @staticmethod
    def rotation_y(radians):
        """Active right-handed 3D rotation about Y, acting on column vectors."""
        angle = finite(radians)
        cosine, sine = math.cos(angle), math.sin(angle)
        return Matrix(3, 3, [cosine,0,sine,0,1,0,-sine,0,cosine])

    @staticmethod
    def rotation_z(radians):
        """Active right-handed 3D rotation about Z, acting on column vectors."""
        angle = finite(radians)
        cosine, sine = math.cos(angle), math.sin(angle)
        return Matrix(3, 3, [cosine,-sine,0,sine,cosine,0,0,0,1])

    @staticmethod
    def rotation_axis_angle(axis, radians):
        """Rodrigues rotation about a finite nonzero three-component row or column axis.

        The axis is normalized without modifying it. Positive angles follow the
        right-hand rule, and the resulting matrix acts on column vectors.
        """
        angle = finite(radians)
        if not isinstance(axis, Matrix):
            raise ValueError("rotation axis must be a three-component vector")
        values = axis._vector3()
        largest = max(map(abs, values))
        if largest == 0:
            raise ValueError("rotation axis must be nonzero")
        scaled = [value/largest for value in values]
        length = math.sqrt(sum(value*value for value in scaled))
        x, y, z = [value/length for value in scaled]
        cosine, sine = math.cos(angle), math.sin(angle)
        complement = 2*math.sin(angle/2)**2 if abs(angle) < 1 else 1-cosine
        return Matrix(3, 3, [cosine+x*x*complement, x*y*complement-z*sine, x*z*complement+y*sine,
                              x*y*complement+z*sine, cosine+y*y*complement, y*z*complement-x*sine,
                              x*z*complement-y*sine, y*z*complement+x*sine, cosine+z*z*complement])

    def _square(self):
        if self.rows != self.cols:
            raise ValueError("operation requires a square matrix")

    def trace(self):
        self._square()
        result = 0.0
        for i in range(self.rows):
            result += self._values[i * self.cols + i]
        return finite(result)

    def determinant(self, algorithm="auto"):
        """Return a nonmutating determinant using auto, lu, or SPD-only cholesky.

        Auto uses exact tiny dyadic formulas and triangular diagonal products,
        otherwise partial-pivot LU. Cholesky requires exact symmetry and positive
        computed pivots; rounding can reject a nearly singular SPD matrix. Final
        overflow raises ValueError; final underflow can return zero.
        """
        self._square()
        if algorithm not in ("auto", "lu", "cholesky"):
            raise ValueError("determinant algorithm must be auto, lu, or cholesky")
        size = self.rows
        if algorithm == "cholesky":
            return self._cholesky_determinant()
        if algorithm == "auto":
            if size == 0:
                return 1.0
            if size == 1:
                return self._values[0]
            # Products and sums of at most six /16 terms remain exact within
            # this bound. Other tiny inputs use LU, including subnormals.
            if size <= 3 and all((value * 16).is_integer() and abs(value * 16) <= 65536 for value in self._values):
                a = self._values
                if size == 2:
                    return a[0] * a[3] - a[1] * a[2]
                return (a[0] * a[4] * a[8] + a[1] * a[5] * a[6] + a[2] * a[3] * a[7]
                        - a[2] * a[4] * a[6] - a[1] * a[3] * a[8] - a[0] * a[5] * a[7])
            if any(self.triangular()):
                return _diagonal_product(self._values, size)
        return self._lu_determinant()

    def _lu_determinant(self):
        size = self.rows
        work = self._values.copy()
        exponent = 0
        for row in range(size):
            start = row * size
            largest = max(abs(value) for value in work[start:start + size])
            if largest == 0.0:
                return 0.0
            shift = math.frexp(largest)[1]
            scaled = [math.ldexp(value, -shift) for value in work[start:start + size]]
            if all(math.ldexp(value, shift) == work[start + index] for index, value in enumerate(scaled)):
                work[start:start + size] = scaled
                exponent += shift
        sign = 1.0
        for k in range(size):
            pivot_row = k
            for row in range(k + 1, size):
                if abs(work[row * size + k]) > abs(work[pivot_row * size + k]):
                    pivot_row = row
            pivot = work[pivot_row * size + k]
            if pivot == 0.0:
                return 0.0
            if pivot_row != k:
                for col in range(k, size):
                    a, b = k * size + col, pivot_row * size + col
                    work[a], work[b] = work[b], work[a]
                sign = -sign
            for row in range(k + 1, size):
                numerator = work[row * size + k]
                factor = numerator / pivot
                tiny_factor = numerator != 0.0 and abs(factor) < 2.0 ** -1022
                work[row * size + k] = 0.0
                for col in range(k + 1, size):
                    index = row * size + col
                    entry = work[k * size + col]
                    product = _quotient_product(numerator, pivot, entry) if tiny_factor else factor * entry
                    work[index] = finite(work[index] - product)
        return _diagonal_product(work, size, exponent, sign)

    def _cholesky_determinant(self):
        size = self.rows
        for row in range(size):
            for col in range(row):
                if self._values[row * size + col] != self._values[col * size + row]:
                    raise ValueError("Cholesky requires an exactly symmetric positive-definite matrix")
        lower = [0.0] * (size * size)
        for row in range(size):
            for col in range(row + 1):
                total = self._values[row * size + col]
                for k in range(col):
                    total -= lower[row * size + k] * lower[col * size + k]
                total = finite(total)
                if row == col:
                    if total <= 0.0:
                        raise ValueError("Cholesky requires positive computed pivots")
                    lower[row * size + col] = math.sqrt(total)
                else:
                    lower[row * size + col] = finite(total / lower[col * size + col])
        return _diagonal_product(lower, size, power=2)

    def triangular(self):
        if self.rows != self.cols:
            return False, False
        upper, lower = True, True
        for row in range(self.rows):
            for col in range(self.cols):
                if self._values[row * self.cols + col] != 0.0:
                    if row > col:
                        upper = False
                    elif row < col:
                        lower = False
        return upper, lower

    def eigen_general(self, max_iterations=1000):
        """Eigenpairs for finite real square inputs, including complex results.

        Columns need not be orthogonal or independent for defective matrices.
        QR iterations are bounded between deflations; failures raise ValueError.
        The result owns its storage and the input is never modified.
        """
        if self.rows != self.cols:
            raise ValueError("eigendecomposition requires a square matrix")
        vr, vi, qr, qi = general_eigen(self._values, self.rows, max_iterations)
        return GeneralEigenDecomposition(vr, vi, Matrix(self.rows, self.cols, qr),
                                          Matrix(self.rows, self.cols, qi))

    def eigen_symmetric(self, tolerance=1e-12, max_sweeps=50):
        """Return real eigenpairs of an exactly symmetric matrix without modifying it.

        Cyclic Jacobi rotations return ascending eigenvalues and orthonormal
        eigenvector columns. The tolerance bounds the off-diagonal Frobenius
        norm relative to the original matrix norm. Global scaling gives normwise
        accuracy; tiny eigenvalues in a mixed-scale input may have large relative
        error. Nonconvergence and nonfinite eigenvalues raise ValueError.
        """
        self._square()
        if type(tolerance) not in (int, float) or not math.isfinite(tolerance) or not 0 < tolerance < 1:
            raise ValueError("eigenvalue tolerance must be finite and between zero and one")
        if type(max_sweeps) is not int or max_sweeps < 1:
            raise ValueError("maximum eigenvalue sweeps must be a positive integer")
        n = self.rows
        diagonal = True
        for i in range(n):
            for j in range(n):
                finite(self._values[i*n+j])
                if self._values[i*n+j] != self._values[j*n+i]:
                    raise ValueError("eigendecomposition requires an exactly symmetric matrix")
                if i != j and self._values[i*n+j] != 0:
                    diagonal = False
        q = Matrix.identity(n)
        if diagonal:
            # Do not globally scale an already diagonal matrix: preserve mixed magnitudes.
            values = [self._values[i*n+i] for i in range(n)]
        else:
            scale = max(abs(value) for value in self._values)
            work = [value/scale for value in self._values]
            norm = math.sqrt(sum(value*value for value in work))
            target = tolerance*norm
            threshold = target/(2*n)

            def converged():
                off_norm = 0.0
                for i in range(n):
                    for j in range(i+1, n):
                        off_norm = math.hypot(off_norm, work[i*n+j], work[i*n+j])
                return off_norm <= target

            for _ in range(max_sweeps):
                if converged():
                    break
                for p in range(n-1):
                    for r in range(p+1, n):
                        off = work[p*n+r]
                        if abs(off) <= threshold:
                            continue
                        delta = (work[r*n+r]-work[p*n+p])/2
                        tangent = off/(delta+math.copysign(math.hypot(delta, off), delta))
                        cosine = 1/math.sqrt(1+tangent*tangent)
                        sine = tangent*cosine
                        work[p*n+p] -= tangent*off
                        work[r*n+r] += tangent*off
                        work[p*n+r] = work[r*n+p] = 0.0
                        for k in range(n):
                            if k != p and k != r:
                                left, right = work[k*n+p], work[k*n+r]
                                work[k*n+p] = work[p*n+k] = cosine*left-sine*right
                                work[k*n+r] = work[r*n+k] = sine*left+cosine*right
                            left, right = q._values[k*n+p], q._values[k*n+r]
                            q._values[k*n+p] = cosine*left-sine*right
                            q._values[k*n+r] = sine*left+cosine*right
            if not converged():
                raise ValueError("symmetric eigendecomposition did not converge within maximum sweeps")
            values = [finite(work[i*n+i]*scale) for i in range(n)]
        order = sorted(range(n), key=lambda index: values[index])
        vectors = Matrix(n, n)
        for col, original in enumerate(order):
            length = math.sqrt(sum(q._values[row*n+original]**2 for row in range(n)))
            pivot = max(range(n), key=lambda row: abs(q._values[row*n+original]))
            sign = -1 if q._values[pivot*n+original] < 0 else 1
            for row in range(n):
                vectors._values[row*n+col] = q._values[row*n+original]*sign/length
        return SymmetricEigenDecomposition([values[index] for index in order], vectors)

    def factor_lu(self):
        from solve import Factorization
        return Factorization(self, "lu")

    def factor_cholesky(self):
        from solve import Factorization
        return Factorization(self, "cholesky")

    def factor_qr(self):
        from solve import Factorization
        return Factorization(self, "qr")

    def solve_ridge(self, rhs, regularization):
        from svd import apply_inverse
        if not isinstance(rhs, Matrix): raise TypeError('right-hand side must be a Matrix')
        if not math.isfinite(regularization) or regularization < 0: raise ValueError('regularization must be finite and nonnegative')
        return apply_inverse(self, rhs, 0 if regularization else None, regularization)

    def pseudoinverse(self, relative_cutoff=None):
        from svd import apply_inverse
        return apply_inverse(self, relative_cutoff=relative_cutoff)

    def solve_minimum_norm(self, rhs, relative_cutoff=None):
        if not isinstance(rhs, Matrix):
            raise TypeError("right-hand side must be a Matrix")
        from svd import apply_inverse
        return apply_inverse(self, rhs, relative_cutoff)

    def spectral_diagnostics(self, relative_cutoff=None):
        from svd import spectral_diagnostics
        return spectral_diagnostics(self, relative_cutoff)

    def svd(self, tolerance=1e-12, max_sweeps=100):
        from svd import decompose
        return decompose(self, tolerance, max_sweeps)

    def solve(self, rhs):
        return self.factor_lu().solve(rhs)

    def least_squares(self, rhs):
        return self.factor_qr().solve(rhs)

    def checksum(self):
        size = len(self._values)
        return finite(self._values[0] + self._values[size // 2] + self._values[-1]) if size else 0.0

    __add__ = add
    __sub__ = subtract
    __matmul__ = multiply

from sparse import CSRMatrix, CGResult
from gmres import GMRESResult
