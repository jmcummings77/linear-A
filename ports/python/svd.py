"""Economy SVD using cyclic one-sided Jacobi rotations (no Gram matrix)."""
import math
from typing import NamedTuple


class SingularValueDecomposition(NamedTuple):
    u: object
    values: list
    vt: object


def decompose(source, tolerance=1e-12, max_sweeps=100):
    from matrix import Matrix
    if not math.isfinite(tolerance) or not 0 < tolerance < 1 or type(max_sweeps) is not int or not 1 <= max_sweeps <= 10000:
        raise ValueError('invalid SVD tolerance or sweep limit')
    rows, cols = source.rows, source.cols
    if rows < cols:
        result = decompose(source.transpose(), tolerance, max_sweeps)
        return SingularValueDecomposition(result.vt.transpose(), result.values, result.u.transpose())
    m, n = rows, cols
    raw = source.values
    if any(not math.isfinite(x) for x in raw):
        raise ValueError('SVD requires finite entries')
    scale = max(map(abs, raw), default=0) or 1.0
    b = [x / scale for x in raw]
    if any(x and not y for x, y in zip(raw, b)):
        raise ArithmeticError('SVD scaling would discard a nonzero entry')
    v = [float(i == j) for i in range(n) for j in range(n)]
    def norm(j):
        total = 0.0
        for i in range(m): total = math.hypot(total, b[i*n+j])
        return total
    converged = False
    for sweep in range(max_sweeps + 1):
        changed = False
        for p in range(n):
            for q in range(p+1, n):
                np, nq = norm(p), norm(q)
                if not np or not nq: continue
                corr = sum((b[i*n+p]/np)*(b[i*n+q]/nq) for i in range(m))
                if abs(corr) <= tolerance: continue
                changed = True
                if sweep == max_sweeps: continue
                pair = max(np, nq)
                ap, aq = np/pair, nq/pair
                delta, g = aq*aq-ap*ap, 2*ap*aq*corr
                t = math.copysign(1.0, g) if delta == 0 else g/(delta+math.copysign(math.hypot(delta,g),delta))
                if abs(t) < 2.2250738585072014e-308:
                    # The required angle is below the normal float64 range.
                    # Deflate the smaller column; its norm is negligible at this scale.
                    small = p if np < nq else q
                    for i in range(m): b[i*n+small] = 0.0
                    continue
                c = 1/math.hypot(1,t); s = c*t
                for data, count in ((b,m),(v,n)):
                    for i in range(count):
                        x,y=data[i*n+p],data[i*n+q]
                        data[i*n+p],data[i*n+q]=c*x-s*y,s*x+c*y
        if not changed:
            converged = True; break
    if not converged: raise ArithmeticError('SVD did not converge')
    norms = [norm(j) for j in range(n)]
    order = sorted(range(n), key=lambda j: -norms[j])
    u, vt, values = [0.0]*(m*n), [0.0]*(n*n), []
    for j,k in enumerate(order):
        value = norms[k]*scale
        if not math.isfinite(value) or (norms[k] and not value):
            raise ArithmeticError('singular value is outside the representable range')
        values.append(value)
        for i in range(n): vt[j*n+i]=v[i*n+k]
        if norms[k]:
            for i in range(m): u[i*n+j]=b[i*n+k]/norms[k]
        else:
            # Complete exact null directions using twice-reorthogonalized unit axes.
            for axis in range(m):
                candidate=[float(i==axis) for i in range(m)]
                for _ in range(2):
                    for k2 in range(j):
                        dot=sum(candidate[i]*u[i*n+k2] for i in range(m))
                        for i in range(m): candidate[i]-=dot*u[i*n+k2]
                length=math.sqrt(sum(x*x for x in candidate))
                if length > 0.5/math.sqrt(m):
                    for i in range(m): u[i*n+j]=candidate[i]/length
                    break
            else: raise ArithmeticError('cannot complete SVD null basis')
    return SingularValueDecomposition(Matrix(m,n,u),values,Matrix(n,n,vt))
