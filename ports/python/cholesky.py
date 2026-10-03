"""Sparse Cholesky with reusable symbolic analysis and explicit fill history."""
import math
from sparse import CSRMatrix


class SparseCholeskySymbolic:
    """Analyze stored structure, including zeros; order the input before analysis."""
    def __init__(self, a):
        if a.rows != a.cols:
            raise ValueError('Cholesky requires a square matrix')
        self._n = a.rows
        self._source = (tuple(a.row_offsets), tuple(a.column_indices))
        graph = [set() for _ in range(self._n)]
        birth = {}
        for i in range(self._n):
            for p in range(a._rp[i], a._rp[i + 1]):
                j = a._ci[p]
                if i != j:
                    graph[i].add(j)
                    graph[j].add(i)
                    birth[max(i, j), min(i, j)] = -1
        for k in range(self._n):
            neighbors = sorted(j for j in graph[k] if j > k)
            for pos, i in enumerate(neighbors):
                for j in neighbors[:pos]:
                    if j not in graph[i]:
                        graph[i].add(j)
                        graph[j].add(i)
                        birth[i, j] = k
        self._rp, self._ci, self._steps = [0], [], []
        for i in range(self._n):
            for j in sorted(j for j in graph[i] if j < i):
                self._ci.append(j)
                self._steps.append(birth[i, j])
            self._ci.append(i)
            self._steps.append(-1)
            self._rp.append(len(self._ci))

    @property
    def size(self): return self._n
    @property
    def nnz(self): return len(self._ci)
    @property
    def fill_count(self): return sum(k >= 0 for k in self._steps)
    @property
    def row_offsets(self): return list(self._rp)
    @property
    def column_indices(self): return list(self._ci)
    @property
    def fill_steps(self): return list(self._steps)

    def factorize(self, a):
        if a.rows != self._n or a.cols != self._n or (tuple(a.row_offsets), tuple(a.column_indices)) != self._source:
            raise ValueError('Cholesky symbolic pattern mismatch')
        rows = [dict(zip(a._ci[a._rp[i]:a._rp[i+1]], a._v[a._rp[i]:a._rp[i+1]])) for i in range(self._n)]
        for i, row in enumerate(rows):
            for j, v in row.items():
                if v != rows[j].get(i, 0.):
                    raise ValueError('Cholesky requires symmetric values')
        rp, ci = self._rp, self._ci
        v = [0.] * len(ci)
        for i in range(self._n):
            for p in range(rp[i], rp[i+1]):
                j = ci[p]
                s = rows[i].get(j, 0.)
                u, w = rp[i], rp[j]
                while u < p and w < rp[j+1] - 1:
                    if ci[u] == ci[w]:
                        s -= v[u] * v[w]
                        u += 1
                        w += 1
                    elif ci[u] < ci[w]: u += 1
                    else: w += 1
                if not math.isfinite(s):
                    raise ArithmeticError('nonfinite Cholesky factor')
                if i == j:
                    if s <= 0: raise ArithmeticError('nonpositive Cholesky pivot')
                    v[p] = math.sqrt(s)
                else:
                    v[p] = s / v[rp[j+1]-1]
                    if not math.isfinite(v[p]): raise ArithmeticError('nonfinite Cholesky factor')
        return SparseCholesky(CSRMatrix(self._n, self._n, rp, ci, v))


class SparseCholesky:
    """Owned lower factor. Construct with SparseCholeskySymbolic.factorize."""
    def __init__(self, lower):
        # Keep this constructor checked: callers cannot supply an invalid factor.
        if lower.rows != lower.cols: raise ValueError('invalid lower factor')
        for i in range(lower.rows):
            start, end = lower._rp[i:i+2]
            if start == end or lower._ci[end-1] != i or lower._v[end-1] <= 0:
                raise ValueError('invalid lower factor')
        self._lower = CSRMatrix(lower.rows, lower.cols, lower.row_offsets, lower.column_indices, lower.values)

    @property
    def size(self): return self._lower.rows
    @property
    def nnz(self): return self._lower.nnz
    @property
    def lower(self):
        a = self._lower
        return CSRMatrix(a.rows, a.cols, a.row_offsets, a.column_indices, a.values)

    def solve(self, b):
        x = list(b)
        if len(x) != self.size or any(not math.isfinite(z) for z in x):
            raise ValueError('invalid Cholesky right-hand side')
        rp, ci, v = self._lower._rp, self._lower._ci, self._lower._v
        for i in range(self.size):
            for p in range(rp[i], rp[i+1]-1): x[i] -= v[p]*x[ci[p]]
            x[i] /= v[rp[i+1]-1]
            if not math.isfinite(x[i]): raise ArithmeticError('nonfinite Cholesky solve')
        for i in range(self.size-1, -1, -1):
            x[i] /= v[rp[i+1]-1]
            if not math.isfinite(x[i]): raise ArithmeticError('nonfinite Cholesky solve')
            for p in range(rp[i], rp[i+1]-1):
                x[ci[p]] -= v[p]*x[i]
                if not math.isfinite(x[ci[p]]): raise ArithmeticError('nonfinite Cholesky solve')
        return x
