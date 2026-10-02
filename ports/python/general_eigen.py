"""Real Hessenberg/double-shift QR adapted from public-domain NIST/MathWorks
JAMA 1.0.3 orthes/hqr2 (EISPACK/Martin-Wilkinson).
https://math.nist.gov/javanumerics/jama/
"""
import math


def _scale(value, exponent):
    return math.ldexp(value, exponent)


def _balance(work):
    """Power-of-two similarity balancing, preserving every nonzero entry."""
    n, powers = len(work), [0]*len(work)
    def log_norm(values):
        largest = max(map(abs,values),default=0)
        return None if largest==0 else math.log2(largest)+math.log2(sum(abs(x)/largest for x in values))
    def log_sum(a,b):
        return max(a,b)+math.log2(1+2**(-abs(a-b)))
    for _ in range(64):
        changed=False
        for i in range(n):
            row=[work[i][j] for j in range(n) if j!=i]
            col=[work[j][i] for j in range(n) if j!=i]
            r,c=log_norm(row),log_norm(col)
            if r is None or c is None: continue
            shift=max(-512,min(512,round((r-c)/2)))
            if not shift or log_sum(r-shift,c+shift)>=log_sum(r,c)+math.log2(.95): continue
            try:
                scaled_row=[_scale(x,-shift) for x in row]
                scaled_col=[_scale(x,shift) for x in col]
                if any(not math.isfinite(y) or x!=0 and y==0 or _scale(y,shift)!=x for x,y in zip(row,scaled_row)): continue
                if any(not math.isfinite(y) or x!=0 and y==0 or _scale(y,-shift)!=x for x,y in zip(col,scaled_col)): continue
            except OverflowError:
                continue
            for j in range(n):
                if j!=i:
                    work[i][j]=_scale(work[i][j],-shift)
                    work[j][i]=_scale(work[j][i],shift)
            powers[i]+=shift
            changed=True
        if not changed: break
    return powers


def general_eigen(input_values, size, max_iterations):
    if type(max_iterations) is not int or not 1 <= max_iterations <= 100000:
        raise ValueError("max_iterations must be an integer between 1 and 100000")
    if not all(map(math.isfinite, input_values)):
        raise ValueError("eigendecomposition requires finite inputs")
    n = size
    diagonal = all(input_values[i*n+j] == 0 for i in range(n) for j in range(n) if i != j)
    upper = all(input_values[i*n+j]==0 for i in range(n) for j in range(i))
    lower = all(input_values[i*n+j]==0 for i in range(n) for j in range(i+1,n))
    reversed_basis = lower and not upper
    H = [[input_values[(n-1-i)*n+n-1-j] if reversed_basis else input_values[i*n+j]
          for j in range(n)] for i in range(n)]
    original_diagonal = [H[i][i] for i in range(n)]
    powers = [0]*n if upper or lower else _balance(H)
    scale = 1 if diagonal else max((abs(x) for row in H for x in row),default=0) or 1
    H = [[x/scale for x in row] for row in H]
    V = [[float(i==j) for j in range(n)] for i in range(n)]
    d, e, ort = [0.0]*n, [0.0]*n, [0.0]*n
    cdivr = cdivi = 0.0
    def orthes ():
        low = 0
        high = n-1
        for m in range(low+1, (high-1) + 1):
            scale = 0.0
            for i in range(m, (high) + 1):
                scale = scale + abs(H[i][m-1])
            if scale != 0.0:
                h = 0.0
                for i in range(high, (m) - 1, -1):
                    ort[i] = H[i][m-1]/scale
                    h += ort[i] * ort[i]
                g = math.sqrt(h)
                if ort[m] > 0:
                    g = -g
                h = h - ort[m] * g
                ort[m] = ort[m] - g
                for j in range(m, n):
                    f = 0.0
                    for i in range(high, (m) - 1, -1):
                        f += ort[i]*H[i][j]
                    f = f/h
                    for i in range(m, (high) + 1):
                        H[i][j] -= f*ort[i]
                for i in range(0, (high) + 1):
                    f = 0.0
                    for j in range(high, (m) - 1, -1):
                        f += ort[j]*H[i][j]
                    f = f/h
                    for j in range(m, (high) + 1):
                        H[i][j] -= f*ort[j]
                ort[m] = scale*ort[m]
                H[m][m-1] = scale*g
        for i in range(0, n):
            for j in range(0, n):
                V[i][j] = (1.0 if i == j else 0.0)
        for m in range(high-1, (low+1) - 1, -1):
            if H[m][m-1] != 0.0:
                for i in range(m+1, (high) + 1):
                    ort[i] = H[i][m-1]
                for j in range(m, (high) + 1):
                    g = 0.0
                    for i in range(m, (high) + 1):
                        g += ort[i] * V[i][j]
                    g = (g / ort[m]) / H[m][m-1]
                    for i in range(m, (high) + 1):
                        V[i][j] += g * ort[i]
    def cdiv(xr, xi, yr, yi):
        nonlocal cdivr, cdivi
        r = 0
        d = 0
        if abs(yr) > abs(yi):
            r = yi/yr
            d = yr + r*yi
            cdivr = (xr + r*xi)/d
            cdivi = (xi - r*xr)/d
        else:
            r = yr/yi
            d = yi + r*yr
            cdivr = (r*xr + xi)/d
            cdivi = (r*xi - xr)/d
    def hqr2 ():
        nn = size
        n = nn-1
        low = 0
        high = nn-1
        eps = pow(2.0,-52.0)
        exshift = 0.0
        p=0
        q=0
        r=0
        s=0
        z=0
        t = 0
        w = 0
        x = 0
        y = 0
        norm = 0.0
        for i in range(0, nn):
            if i < low or i > high:
                d[i] = H[i][i]
                e[i] = 0.0
            for j in range(max(i-1,0), nn):
                norm = norm + abs(H[i][j])
        iter = 0
        while n >= low:
            l = n
            while l > low:
                s = abs(H[l-1][l-1]) + abs(H[l][l])
                if s == 0.0:
                    s = norm
                if abs(H[l][l-1]) <= eps * s:
                    break
                l -= 1
            if l == n:
                H[n][n] = H[n][n] + exshift
                d[n] = H[n][n]
                e[n] = 0.0
                n -= 1
                iter = 0
            elif l == n-1:
                w = H[n][n-1] * H[n-1][n]
                p = (H[n-1][n-1] - H[n][n]) / 2.0
                q = p * p + w
                z = math.sqrt(abs(q))
                H[n][n] = H[n][n] + exshift
                H[n-1][n-1] = H[n-1][n-1] + exshift
                x = H[n][n]
                if q >= 0:
                    if p >= 0:
                        z = p + z
                    else:
                        z = p - z
                    d[n-1] = x + z
                    d[n] = d[n-1]
                    if z != 0.0:
                        d[n] = x - w / z
                    e[n-1] = 0.0
                    e[n] = 0.0
                    x = H[n][n-1]
                    s = abs(x) + abs(z)
                    p = x / s
                    q = z / s
                    r = math.sqrt(p * p+q * q)
                    p = p / r
                    q = q / r
                    for j in range(n-1, nn):
                        z = H[n-1][j]
                        H[n-1][j] = q * z + p * H[n][j]
                        H[n][j] = q * H[n][j] - p * z
                    for i in range(0, (n) + 1):
                        z = H[i][n-1]
                        H[i][n-1] = q * z + p * H[i][n]
                        H[i][n] = q * H[i][n] - p * z
                    for i in range(low, (high) + 1):
                        z = V[i][n-1]
                        V[i][n-1] = q * z + p * V[i][n]
                        V[i][n] = q * V[i][n] - p * z
                else:
                    d[n-1] = x + p
                    d[n] = x + p
                    e[n-1] = z
                    e[n] = -z
                n = n - 2
                iter = 0
            else:
                x = H[n][n]
                y = 0.0
                w = 0.0
                if l < n:
                    y = H[n-1][n-1]
                    w = H[n][n-1] * H[n-1][n]
                if iter == 10:
                    exshift += x
                    for i in range(low, (n) + 1):
                        H[i][i] -= x
                    s = abs(H[n][n-1]) + abs(H[n-1][n-2])
                    x = y = 0.75 * s
                    w = -0.4375 * s * s
                if iter == 30:
                    s = (y - x) / 2.0
                    s = s * s + w
                    if s > 0:
                        s = math.sqrt(s)
                        if y < x:
                            s = -s
                        s = x - w / ((y - x) / 2.0 + s)
                        for i in range(low, (n) + 1):
                            H[i][i] -= s
                        exshift += s
                        x = y = w = 0.964
                iter = iter + 1
                if iter > max_iterations:
                    raise ValueError("general eigendecomposition did not converge within maximum iterations")
                m = n-2
                while m >= l:
                    z = H[m][m]
                    r = x - z
                    s = y - z
                    p = (r * s - w) / H[m+1][m] + H[m][m+1]
                    q = H[m+1][m+1] - z - r - s
                    r = H[m+2][m+1]
                    s = abs(p) + abs(q) + abs(r)
                    p = p / s
                    q = q / s
                    r = r / s
                    if m == l:
                        break
                    if abs(H[m][m-1]) * (abs(q) + abs(r)) < eps * (abs(p) * (abs(H[m-1][m-1]) + abs(z) + abs(H[m+1][m+1]))):
                        break
                    m -= 1
                for i in range(m+2, (n) + 1):
                    H[i][i-2] = 0.0
                    if i > m+2:
                        H[i][i-3] = 0.0
                for k in range(m, (n-1) + 1):
                    notlast = (k != n-1)
                    if k != m:
                        p = H[k][k-1]
                        q = H[k+1][k-1]
                        r = (H[k+2][k-1] if notlast else 0.0)
                        x = abs(p) + abs(q) + abs(r)
                        if x == 0.0:
                            continue
                        p = p / x
                        q = q / x
                        r = r / x
                    s = math.sqrt(p * p + q * q + r * r)
                    if p < 0:
                        s = -s
                    if s != 0:
                        if k != m:
                            H[k][k-1] = -s * x
                        elif l != m:
                            H[k][k-1] = -H[k][k-1]
                        p = p + s
                        x = p / s
                        y = q / s
                        z = r / s
                        q = q / p
                        r = r / p
                        for j in range(k, nn):
                            p = H[k][j] + q * H[k+1][j]
                            if notlast:
                                p = p + r * H[k+2][j]
                                H[k+2][j] = H[k+2][j] - p * z
                            H[k][j] = H[k][j] - p * x
                            H[k+1][j] = H[k+1][j] - p * y
                        for i in range(0, (min(n,k+3)) + 1):
                            p = x * H[i][k] + y * H[i][k+1]
                            if notlast:
                                p = p + z * H[i][k+2]
                                H[i][k+2] = H[i][k+2] - p * r
                            H[i][k] = H[i][k] - p
                            H[i][k+1] = H[i][k+1] - p * q
                        for i in range(low, (high) + 1):
                            p = x * V[i][k] + y * V[i][k+1]
                            if notlast:
                                p = p + z * V[i][k+2]
                                V[i][k+2] = V[i][k+2] - p * r
                            V[i][k] = V[i][k] - p
                            V[i][k+1] = V[i][k+1] - p * q
        if norm == 0.0:
            return
        for n in range(nn-1, (0) - 1, -1):
            p = d[n]
            q = e[n]
            if q == 0:
                l = n
                H[n][n] = 1.0
                for i in range(n-1, (0) - 1, -1):
                    w = H[i][i] - p
                    r = 0.0
                    for j in range(l, (n) + 1):
                        r = r + H[i][j] * H[j][n]
                    if e[i] < 0.0:
                        z = w
                        s = r
                    else:
                        l = i
                        if e[i] == 0.0:
                            if w != 0.0:
                                H[i][n] = -r / w
                            else:
                                H[i][n] = -r / (eps * norm)
                        else:
                            x = H[i][i+1]
                            y = H[i+1][i]
                            q = (d[i] - p) * (d[i] - p) + e[i] * e[i]
                            t = (x * s - z * r) / q
                            H[i][n] = t
                            if abs(x) > abs(z):
                                H[i+1][n] = (-r - w * t) / x
                            else:
                                H[i+1][n] = (-s - y * t) / z
                        t = abs(H[i][n])
                        if (eps * t) * t > 1:
                            for j in range(i, (n) + 1):
                                H[j][n] = H[j][n] / t
            elif q < 0:
                l = n-1
                if abs(H[n][n-1]) > abs(H[n-1][n]):
                    H[n-1][n-1] = q / H[n][n-1]
                    H[n-1][n] = -(H[n][n] - p) / H[n][n-1]
                else:
                    cdiv(0.0,-H[n-1][n],H[n-1][n-1]-p,q)
                    H[n-1][n-1] = cdivr
                    H[n-1][n] = cdivi
                H[n][n-1] = 0.0
                H[n][n] = 1.0
                for i in range(n-2, (0) - 1, -1):
                    ra = 0
                    sa = 0
                    vr = 0
                    vi = 0
                    ra = 0.0
                    sa = 0.0
                    for j in range(l, (n) + 1):
                        ra = ra + H[i][j] * H[j][n-1]
                        sa = sa + H[i][j] * H[j][n]
                    w = H[i][i] - p
                    if e[i] < 0.0:
                        z = w
                        r = ra
                        s = sa
                    else:
                        l = i
                        if e[i] == 0:
                            cdiv(-ra,-sa,w,q)
                            H[i][n-1] = cdivr
                            H[i][n] = cdivi
                        else:
                            x = H[i][i+1]
                            y = H[i+1][i]
                            vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q
                            vi = (d[i] - p) * 2.0 * q
                            if vr == 0.0 and vi == 0.0:
                                vr = eps * norm * (abs(w) + abs(q) + abs(x) + abs(y) + abs(z))
                            cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi)
                            H[i][n-1] = cdivr
                            H[i][n] = cdivi
                            if abs(x) > (abs(z) + abs(q)):
                                H[i+1][n-1] = (-ra - w * H[i][n-1] + q * H[i][n]) / x
                                H[i+1][n] = (-sa - w * H[i][n] - q * H[i][n-1]) / x
                            else:
                                cdiv(-r-y*H[i][n-1],-s-y*H[i][n],z,q)
                                H[i+1][n-1] = cdivr
                                H[i+1][n] = cdivi
                        t = max(abs(H[i][n-1]),abs(H[i][n]))
                        if (eps * t) * t > 1:
                            for j in range(i, (n) + 1):
                                H[j][n-1] = H[j][n-1] / t
                                H[j][n] = H[j][n] / t
        for i in range(0, nn):
            if i < low or i > high:
                for j in range(i, nn):
                    V[i][j] = H[i][j]
        for j in range(nn-1, (low) - 1, -1):
            for i in range(low, (high) + 1):
                z = 0.0
                for k in range(low, (min(j,high)) + 1):
                    z = z + V[i][k] * H[k][j]
                V[i][j] = z

    if diagonal:
        d = [input_values[i*n+i] for i in range(n)]
    else:
        try:
            orthes()
            hqr2()
        except (ZeroDivisionError, OverflowError) as error:
            raise ValueError("eigendecomposition exceeded the finite working range") from error
    if upper or lower:
        d, e, scale = original_diagonal, [0.0]*n, 1
    if reversed_basis:
        V.reverse()
    order = sorted(range(n), key=lambda i: (d[i],e[i]))
    real_values, imag_values = [d[i]*scale for i in order], [e[i]*scale for i in order]
    real_vectors, imag_vectors = [0.0]*(n*n), [0.0]*(n*n)
    for col, source in enumerate(order):
        real_col = source-1 if e[source] < 0 else source
        re = [row[real_col] for row in V]
        im = [0.0 if e[source]==0 else row[real_col+1]*(-1 if e[source]<0 else 1) for row in V]
        # Undo similarity scaling relative to a shared exponent to avoid overflow.
        if any(powers):
            exponents=[math.frexp(max(abs(re[i]),abs(im[i])))[1]+powers[i]
                       for i in range(n) if re[i]!=0 or im[i]!=0]
            exponent=max(exponents)
            re=[_scale(x,powers[i]-exponent) for i,x in enumerate(re)]
            im=[_scale(x,powers[i]-exponent) for i,x in enumerate(im)]
        magnitude = max(abs(x) for x in re+im)
        if not math.isfinite(magnitude) or magnitude == 0:
            raise ValueError("nonfinite or zero eigenvector")
        re, im = [x/magnitude for x in re], [x/magnitude for x in im]
        norm = math.sqrt(math.fsum(x*x for x in re+im))
        pivot = max(range(n),key=lambda i: math.hypot(re[i],im[i]))
        phase = math.hypot(re[pivot],im[pivot])
        pr, pi = re[pivot]/phase, im[pivot]/phase
        for row in range(n):
            real_vectors[row*n+col] = (re[row]*pr+im[row]*pi)/norm
            imag_vectors[row*n+col] = (im[row]*pr-re[row]*pi)/norm
    if not all(map(math.isfinite,real_values+imag_values+real_vectors+imag_vectors)):
        raise ValueError("nonfinite eigendecomposition result")
    return real_values, imag_values, real_vectors, imag_vectors
