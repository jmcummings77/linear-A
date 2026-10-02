"""Reusable float64 LU, Cholesky and column-pivoted Householder QR factors."""
import math

EPSILON = 2.220446049250313e-16


def checked(value):
    if not math.isfinite(value):
        raise ValueError("solver arithmetic exceeds the finite float64 range")
    return value


def scaled(value, scale):
    result = checked(value / scale)
    if value and not result:
        raise ValueError("solver scaling would discard a nonzero value")
    return result


class Factorization:
    """An independent snapshot. solve accepts one or more right-hand-side columns.

    LU rejects zero computed pivots. QR requires rows >= columns and rejects
    rank below eps * max(rows, cols) * largest initial column norm. Cholesky
    requires exact symmetry and positive computed pivots. Scaling may fail on
    inputs spanning more than the representable float64 range.
    """
    def __init__(self, matrix, algorithm):
        if algorithm not in ("lu", "cholesky", "qr"):
            raise ValueError("unknown factorization")
        m, n = matrix.rows, matrix.cols
        if (algorithm != "qr" and m != n) or m < n:
            raise ValueError("factorization requires square input, or rows >= columns for QR")
        original = matrix.values
        maximum = max(map(abs, original), default=0.0)
        self._scale = math.ldexp(1.0, math.frexp(maximum)[1]-1) if maximum else 1.0
        a = [scaled(value, self._scale) for value in original]
        self._rows, self._cols, self._algorithm = m, n, algorithm
        self._norm = max((sum(abs(a[i*n+j]) for j in range(n)) for i in range(m)), default=0.0)
        self._permutation, self._tau = list(range(m if algorithm == "lu" else n)), [0.0]*n
        p = self._permutation
        if algorithm == "lu":
            for k in range(n):
                pivot = max(range(k,n), key=lambda i: abs(a[i*n+k]))
                if a[pivot*n+k] == 0:
                    raise ValueError("singular matrix: zero computed LU pivot")
                a[k*n:(k+1)*n], a[pivot*n:(pivot+1)*n] = a[pivot*n:(pivot+1)*n], a[k*n:(k+1)*n]
                p[k], p[pivot] = p[pivot], p[k]
                for i in range(k+1,n):
                    a[i*n+k] = checked(a[i*n+k]/a[k*n+k])
                    for j in range(k+1,n):
                        a[i*n+j] = checked(a[i*n+j]-a[i*n+k]*a[k*n+j])
        elif algorithm == "cholesky":
            if any(original[i*n+j] != original[j*n+i] for i in range(n) for j in range(i)):
                raise ValueError("Cholesky requires exact symmetry")
            for i in range(n):
                for j in range(i+1):
                    value = a[i*n+j]
                    for k in range(j):
                        value = checked(value-a[i*n+k]*a[j*n+k])
                    if i == j:
                        if value <= 0:
                            raise ValueError("Cholesky requires positive computed pivots")
                        a[i*n+j] = math.sqrt(value)
                    else:
                        a[i*n+j] = checked(value/a[j*n+j])
        else:
            largest = max((math.hypot(*(a[i*n+j] for i in range(m))) for j in range(n)), default=0.0)
            threshold = EPSILON*max(m,n)*largest
            for k in range(n):
                norms = [math.hypot(*(a[i*n+j] for i in range(k,m))) for j in range(k,n)]
                pivot = k+max(range(len(norms)), key=norms.__getitem__)
                norm = norms[pivot-k]
                if norm <= threshold:
                    raise ValueError("QR input is numerically rank deficient")
                for i in range(m):
                    a[i*n+k], a[i*n+pivot] = a[i*n+pivot], a[i*n+k]
                p[k],p[pivot] = p[pivot],p[k]
                old = a[k*n+k]
                alpha = -math.copysign(norm, old)
                divisor = old-alpha
                self._tau[k] = (alpha-old)/alpha
                for i in range(k+1,m):
                    a[i*n+k] /= divisor
                a[k*n+k] = alpha
                for j in range(k+1,n):
                    dot = a[k*n+j]
                    for i in range(k+1,m):
                        dot += a[i*n+k]*a[i*n+j]
                    dot *= self._tau[k]
                    a[k*n+j] = checked(a[k*n+j]-dot)
                    for i in range(k+1,m):
                        a[i*n+j] = checked(a[i*n+j]-a[i*n+k]*dot)
        self._data = a

    def _solve(self, rhs, rescale):
        from matrix import Matrix
        m,n,p = self._rows,self._cols,rhs.cols
        if rhs.rows != m:
            raise ValueError("right-hand side row count must match the factorization")
        a = self._data
        b = rhs.values
        work = [scaled(b[(self._permutation[i] if self._algorithm == "lu" else i)*p+j], self._scale) if rescale
                else b[(self._permutation[i] if self._algorithm == "lu" else i)*p+j] for i in range(m) for j in range(p)]
        if self._algorithm == "qr":
            for k in range(n):
                for j in range(p):
                    dot = work[k*p+j]
                    for i in range(k+1,m):
                        dot = checked(dot+a[i*n+k]*work[i*p+j])
                    dot = checked(dot*self._tau[k])
                    work[k*p+j] = checked(work[k*p+j]-dot)
                    for i in range(k+1,m):
                        work[i*p+j] = checked(work[i*p+j]-a[i*n+k]*dot)
        else:
            for i in range(n):
                for j in range(p):
                    value = work[i*p+j]
                    for k in range(i):
                        value = checked(value-a[i*n+k]*work[k*p+j])
                    work[i*p+j] = checked(value/a[i*n+i]) if self._algorithm == "cholesky" else value
        for i in range(n-1,-1,-1):
            for j in range(p):
                value = work[i*p+j]
                for k in range(i+1,n):
                    coefficient = a[k*n+i] if self._algorithm == "cholesky" else a[i*n+k]
                    value = checked(value-coefficient*work[k*p+j])
                work[i*p+j] = checked(value/a[i*n+i])
        result = [0.0]*(n*p)
        for i in range(n):
            row = self._permutation[i] if self._algorithm == "qr" else i
            result[row*p:(row+1)*p] = work[i*p:(i+1)*p]
        return Matrix(n,p,result)

    def solve(self, rhs):
        return self._solve(rhs, True)

    def reciprocal_condition(self):
        """Computed infinity-norm rcond, using n solves; O(n^3), not a bound.

        Returns zero when inverse arithmetic overflows. Empty input returns one.
        Only square factors have this diagnostic.
        """
        from matrix import Matrix
        n = self._cols
        if self._rows != n:
            raise ValueError("condition diagnostic requires a square matrix")
        if not n:
            return 1.0
        try:
            inverse = self._solve(Matrix.identity(n), False).values
            norm = max(sum(abs(inverse[i*n+j]) for j in range(n)) for i in range(n))
            return min(1.0, (1.0/self._norm)/norm)
        except ValueError:
            return 0.0
