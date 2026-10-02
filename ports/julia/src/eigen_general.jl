# Hessenberg reduction and real double-shift QR follow public-domain JAMA's
# implementation of the EISPACK orthes/hqr2 algorithms.
# Reference: https://math.nist.gov/javanumerics/jama/

"""General complex right eigenpairs of a finite real square matrix.

Uses machine-epsilon deflation and bounded Hessenberg double-shift QR. Values
are sorted by (real, imaginary); matching complex columns have unit norm. No
orthogonality or independent-basis guarantee is made for defective matrices.
The input is unchanged. max_iterations must be an integer from 1 to 100000.
"""
function eigen_general(a::Matrix64; max_iterations=1000)
    a.rows == a.cols || throw(DimensionMismatch("eigendecomposition requires a square matrix"))
    (!(max_iterations isa Integer) || max_iterations isa Bool || !(1 <= max_iterations <= 100000)) && throw(ArgumentError("iteration limit must be between 1 and 100000"))
    n=a.rows
    all(isfinite,a.data) || throw(ArgumentError("eigendecomposition requires finite entries"))
    H=[a[i,j] for i in 1:n,j in 1:n]
    V=zeros(n,n); d=zeros(n); e=zeros(n)
    for i in 1:n; V[i,i]=1; d[i]=H[i,i]; end
    diagonal=all(i==j || H[i,j]==0 for i in 1:n,j in 1:n)
    upper=all(i<=j || H[i,j]==0 for i in 1:n,j in 1:n)
    lower=all(i>=j || H[i,j]==0 for i in 1:n,j in 1:n)
    reverse=lower&&!upper
    scale=1.0
    balance=zeros(Int,n)
    if !diagonal
        reverse && (H=H[end:-1:1,end:-1:1])
        if !upper&&!lower;balance=general_balance!(H);end
        scale=maximum(abs,H)
        bits=(reinterpret(UInt64,scale)>>52)&0x7ff
        bits!=0 && (scale=reinterpret(Float64,bits<<52))
        H ./= scale
        general_hessenberg!(H,V,zeros(n),n)
        general_qr!(H,V,d,e,max_iterations)
    end
    imag_signs=copy(e)
    if upper||lower
        for i in 1:n;original=reverse ? n+1-i : i;d[i]=a[original,original];end
    else
        d .*= scale
    end
    e .*= scale
    all(isfinite,d) && all(isfinite,e) || throw(ErrorException("nonfinite eigenvalue"))
    order=sortperm(1:n,by=i->(d[i],e[i]))
    values=ComplexF64[complex(d[i],e[i]) for i in order]
    vectors=zeros(ComplexF64,n,n)
    for (col,original) in enumerate(order)
        rc=original; ic=0; sign=0.0
        if imag_signs[original]>0; ic=original+1;sign=1.0
        elseif imag_signs[original]<0;rc=original-1;ic=original;sign=-1.0;end
        length=0.0
        for row in 1:n
            input_row=reverse ? n+1-row : row
            r=V[input_row,rc];im=ic==0 ? 0.0 : sign*V[input_row,ic]
            isfinite(r)&&isfinite(im) || throw(ErrorException("nonfinite eigenvector"))
            vectors[row,col]=complex(r,im)

        end
        general_unbalance!(vectors,col,balance)
        for row in 1:n;length=hypot(length,abs(vectors[row,col]));end
        isfinite(length)&&length>0 || throw(ErrorException("invalid eigenvector norm"))
        vectors[:,col] ./= length
    end
    (values=values,vectors=vectors)
end

function general_cdiv(xr,xi,yr,yi)
    if abs(yr)>abs(yi)
        r=yi/yr; den=yr+r*yi
        return ((xr+r*xi)/den,(xi-r*xr)/den)
    end
    r=yr/yi;den=yi+r*yr
    ((r*xr+xi)/den,(r*xi-xr)/den)
end

function general_hessenberg!(H,V,ort,n)
    low = 0
    high = n-1
    for m in (low+1):(high-1)
        scale = 0.0
        for i in (m):(high)
            scale = scale + abs(H[(i)+1,(m-1)+1])
        end
        if (scale != 0.0)
            h = 0.0
            for i in (high):-1:(m)
                ort[(i)+1] = H[(i)+1,(m-1)+1]/scale
                h += ort[(i)+1] * ort[(i)+1]
            end
            g = sqrt(h)
            if (ort[(m)+1] > 0)
                g = -g
            end
            h = h - ort[(m)+1] * g
            ort[(m)+1] = ort[(m)+1] - g
            for j in (m):(n-1)
                f = 0.0
                for i in (high):-1:(m)
                    f += ort[(i)+1]*H[(i)+1,(j)+1]
                end
                f = f/h
                for i in (m):(high)
                    H[(i)+1,(j)+1] -= f*ort[(i)+1]
                end
            end
            for i in (0):(high)
                f = 0.0
                for j in (high):-1:(m)
                    f += ort[(j)+1]*H[(i)+1,(j)+1]
                end
                f = f/h
                for j in (m):(high)
                    H[(i)+1,(j)+1] -= f*ort[(j)+1]
                end
            end
            ort[(m)+1] = scale*ort[(m)+1]
            H[(m)+1,(m-1)+1] = scale*g
        end
    end
    for i in (0):(n-1)
        for j in (0):(n-1)
            V[(i)+1,(j)+1] = 0.0
            if (i == j)
                V[(i)+1,(j)+1] = 1.0
            end
        end
    end
    for m in (high-1):-1:(low+1)
        if (H[(m)+1,(m-1)+1] != 0.0)
            for i in (m+1):(high)
                ort[(i)+1] = H[(i)+1,(m-1)+1]
            end
            for j in (m):(high)
                g = 0.0
                for i in (m):(high)
                    g += ort[(i)+1] * V[(i)+1,(j)+1]
                end
                g = (g / ort[(m)+1]) / H[(m)+1,(m-1)+1]
                for i in (m):(high)
                    V[(i)+1,(j)+1] += g * ort[(i)+1]
                end
            end
        end
    end
end

function general_qr!(H,V,d,e,maxIterations)
    nn = size(H,1)
    n = nn-1
    low = 0
    high = nn-1
    eps = Base.eps(Float64)
    exshift = 0.0
    p=0
    q=0
    r=0
    s=0
    z=0
    t=0.0
    w=0.0
    x=0.0
    y=0.0
    norm = 0.0
    for i in (0):(nn-1)
        if (i < low || i > high)
            d[(i)+1] = H[(i)+1,(i)+1]
            e[(i)+1] = 0.0
        end
        for j in (max(i-1,0)):(nn-1)
            norm = norm + abs(H[(i)+1,(j)+1])
        end
    end
    iter = 0
    while (n >= low)
        l = n
        while (l > low)
            s = abs(H[(l-1)+1,(l-1)+1]) + abs(H[(l)+1,(l)+1])
            if (s == 0.0)
                s = norm
            end
            if (abs(H[(l)+1,(l-1)+1]) <= eps * s)
                break
            end
            l -= 1
        end
        if (l == n)
            H[(n)+1,(n)+1] = H[(n)+1,(n)+1] + exshift
            d[(n)+1] = H[(n)+1,(n)+1]
            e[(n)+1] = 0.0
            n -= 1
            iter = 0
        elseif (l == n-1)
            w = H[(n)+1,(n-1)+1] * H[(n-1)+1,(n)+1]
            p = (H[(n-1)+1,(n-1)+1] - H[(n)+1,(n)+1]) / 2.0
            q = p * p + w
            z = sqrt(abs(q))
            H[(n)+1,(n)+1] = H[(n)+1,(n)+1] + exshift
            H[(n-1)+1,(n-1)+1] = H[(n-1)+1,(n-1)+1] + exshift
            x = H[(n)+1,(n)+1]
            if (q >= 0)
                if (p >= 0)
                    z = p + z
                else
                    z = p - z
                end
                d[(n-1)+1] = x + z
                d[(n)+1] = d[(n-1)+1]
                if (z != 0.0)
                    d[(n)+1] = x - w / z
                end
                e[(n-1)+1] = 0.0
                e[(n)+1] = 0.0
                x = H[(n)+1,(n-1)+1]
                s = abs(x) + abs(z)
                p = x / s
                q = z / s
                r = sqrt(p * p+q * q)
                p = p / r
                q = q / r
                for j in (n-1):(nn-1)
                    z = H[(n-1)+1,(j)+1]
                    H[(n-1)+1,(j)+1] = q * z + p * H[(n)+1,(j)+1]
                    H[(n)+1,(j)+1] = q * H[(n)+1,(j)+1] - p * z
                end
                for i in (0):(n)
                    z = H[(i)+1,(n-1)+1]
                    H[(i)+1,(n-1)+1] = q * z + p * H[(i)+1,(n)+1]
                    H[(i)+1,(n)+1] = q * H[(i)+1,(n)+1] - p * z
                end
                for i in (low):(high)
                    z = V[(i)+1,(n-1)+1]
                    V[(i)+1,(n-1)+1] = q * z + p * V[(i)+1,(n)+1]
                    V[(i)+1,(n)+1] = q * V[(i)+1,(n)+1] - p * z
                end
            else
                d[(n-1)+1] = x + p
                d[(n)+1] = x + p
                e[(n-1)+1] = z
                e[(n)+1] = -z
            end
            n = n - 2
            iter = 0
        else
            x = H[(n)+1,(n)+1]
            y = 0.0
            w = 0.0
            if (l < n)
                y = H[(n-1)+1,(n-1)+1]
                w = H[(n)+1,(n-1)+1] * H[(n-1)+1,(n)+1]
            end
            if (iter == 10)
                exshift += x
                for i in (low):(n)
                    H[(i)+1,(i)+1] -= x
                end
                s = abs(H[(n)+1,(n-1)+1]) + abs(H[(n-1)+1,(n-2)+1])
                x = 0.75 * s
                y = x
                w = -0.4375 * s * s
            end
            if (iter == 30)
                s = (y - x) / 2.0
                s = s * s + w
                if (s > 0)
                    s = sqrt(s)
                    if (y < x)
                        s = -s
                    end
                    s = x - w / ((y - x) / 2.0 + s)
                    for i in (low):(n)
                        H[(i)+1,(i)+1] -= s
                    end
                    exshift += s
                    x = 0.964
                    y = x
                    w = x
                end
            end
            iter = iter + 1
            if (iter > maxIterations)
                throw(ErrorException("general eigendecomposition did not converge within the iteration limit"))
            end
            m = n-2
            while (m >= l)
                z = H[(m)+1,(m)+1]
                r = x - z
                s = y - z
                p = (r * s - w) / H[(m+1)+1,(m)+1] + H[(m)+1,(m+1)+1]
                q = H[(m+1)+1,(m+1)+1] - z - r - s
                r = H[(m+2)+1,(m+1)+1]
                s = abs(p) + abs(q) + abs(r)
                p = p / s
                q = q / s
                r = r / s
                if (m == l)
                    break
                end
                if (abs(H[(m)+1,(m-1)+1]) * (abs(q) + abs(r)) <
                    eps * (abs(p) * (abs(H[(m-1)+1,(m-1)+1]) + abs(z) +
                    abs(H[(m+1)+1,(m+1)+1]))))
                    break
                end
                m -= 1
            end
            for i in (m+2):(n)
                H[(i)+1,(i-2)+1] = 0.0
                if (i > m+2)
                    H[(i)+1,(i-3)+1] = 0.0
                end
            end
            for k in (m):(n-1)
                notlast = (k != n-1)
                if (k != m)
                    p = H[(k)+1,(k-1)+1]
                    q = H[(k+1)+1,(k-1)+1]
                    r = 0.0
                    if (notlast)
                        r = H[(k+2)+1,(k-1)+1]
                    end
                    x = abs(p) + abs(q) + abs(r)
                    if (x == 0.0)
                        continue
                    end
                    p = p / x
                    q = q / x
                    r = r / x
                end
                s = sqrt(p * p + q * q + r * r)
                if (p < 0)
                    s = -s
                end
                if (s != 0)
                    if (k != m)
                        H[(k)+1,(k-1)+1] = -s * x
                    elseif (l != m)
                        H[(k)+1,(k-1)+1] = -H[(k)+1,(k-1)+1]
                    end
                    p = p + s
                    x = p / s
                    y = q / s
                    z = r / s
                    q = q / p
                    r = r / p
                    for j in (k):(nn-1)
                        p = H[(k)+1,(j)+1] + q * H[(k+1)+1,(j)+1]
                        if (notlast)
                            p = p + r * H[(k+2)+1,(j)+1]
                            H[(k+2)+1,(j)+1] = H[(k+2)+1,(j)+1] - p * z
                        end
                        H[(k)+1,(j)+1] = H[(k)+1,(j)+1] - p * x
                        H[(k+1)+1,(j)+1] = H[(k+1)+1,(j)+1] - p * y
                    end
                    for i in (0):(min(n,k+3))
                        p = x * H[(i)+1,(k)+1] + y * H[(i)+1,(k+1)+1]
                        if (notlast)
                            p = p + z * H[(i)+1,(k+2)+1]
                            H[(i)+1,(k+2)+1] = H[(i)+1,(k+2)+1] - p * r
                        end
                        H[(i)+1,(k)+1] = H[(i)+1,(k)+1] - p
                        H[(i)+1,(k+1)+1] = H[(i)+1,(k+1)+1] - p * q
                    end
                    for i in (low):(high)
                        p = x * V[(i)+1,(k)+1] + y * V[(i)+1,(k+1)+1]
                        if (notlast)
                            p = p + z * V[(i)+1,(k+2)+1]
                            V[(i)+1,(k+2)+1] = V[(i)+1,(k+2)+1] - p * r
                        end
                        V[(i)+1,(k)+1] = V[(i)+1,(k)+1] - p
                        V[(i)+1,(k+1)+1] = V[(i)+1,(k+1)+1] - p * q
                    end
                end
            end
        end
    end
    if (norm == 0.0)
        return nothing
    end
    for n in (nn-1):-1:(0)
        p = d[(n)+1]
        q = e[(n)+1]
        if (q == 0)
            l = n
            H[(n)+1,(n)+1] = 1.0
            for i in (n-1):-1:(0)
                w = H[(i)+1,(i)+1] - p
                r = 0.0
                for j in (l):(n)
                    r = r + H[(i)+1,(j)+1] * H[(j)+1,(n)+1]
                end
                if (e[(i)+1] < 0.0)
                    z = w
                    s = r
                else
                    l = i
                    if (e[(i)+1] == 0.0)
                        if (w != 0.0)
                            H[(i)+1,(n)+1] = -r / w
                        else
                            H[(i)+1,(n)+1] = -r / (eps * norm)
                        end
                    else
                        x = H[(i)+1,(i+1)+1]
                        y = H[(i+1)+1,(i)+1]
                        q = (d[(i)+1] - p) * (d[(i)+1] - p) + e[(i)+1] * e[(i)+1]
                        t = (x * s - z * r) / q
                        H[(i)+1,(n)+1] = t
                        if (abs(x) > abs(z))
                            H[(i+1)+1,(n)+1] = (-r - w * t) / x
                        else
                            H[(i+1)+1,(n)+1] = (-s - y * t) / z
                        end
                    end
                    t = abs(H[(i)+1,(n)+1])
                    if ((eps * t) * t > 1)
                        for j in (i):(n)
                            H[(j)+1,(n)+1] = H[(j)+1,(n)+1] / t
                        end
                    end
                end
            end
        elseif (q < 0)
            l = n-1
            if (abs(H[(n)+1,(n-1)+1]) > abs(H[(n-1)+1,(n)+1]))
                H[(n-1)+1,(n-1)+1] = q / H[(n)+1,(n-1)+1]
                H[(n-1)+1,(n)+1] = -(H[(n)+1,(n)+1] - p) / H[(n)+1,(n-1)+1]
            else
                cdivr,cdivi=general_cdiv(0.0,-H[(n-1)+1,(n)+1],H[(n-1)+1,(n-1)+1]-p,q)
                H[(n-1)+1,(n-1)+1] = cdivr
                H[(n-1)+1,(n)+1] = cdivi
            end
            H[(n)+1,(n-1)+1] = 0.0
            H[(n)+1,(n)+1] = 1.0
            for i in (n-2):-1:(0)
                ra=0.0
                sa=0.0
                vr=0.0
                vi=0.0
                ra = 0.0
                sa = 0.0
                for j in (l):(n)
                    ra = ra + H[(i)+1,(j)+1] * H[(j)+1,(n-1)+1]
                    sa = sa + H[(i)+1,(j)+1] * H[(j)+1,(n)+1]
                end
                w = H[(i)+1,(i)+1] - p
                if (e[(i)+1] < 0.0)
                    z = w
                    r = ra
                    s = sa
                else
                    l = i
                    if (e[(i)+1] == 0)
                        cdivr,cdivi=general_cdiv(-ra,-sa,w,q)
                        H[(i)+1,(n-1)+1] = cdivr
                        H[(i)+1,(n)+1] = cdivi
                    else
                        x = H[(i)+1,(i+1)+1]
                        y = H[(i+1)+1,(i)+1]
                        vr = (d[(i)+1] - p) * (d[(i)+1] - p) + e[(i)+1] * e[(i)+1] - q * q
                        vi = (d[(i)+1] - p) * 2.0 * q
                        if (vr == 0.0 && vi == 0.0)
                            vr = eps * norm * (abs(w) + abs(q) +
                            abs(x) + abs(y) + abs(z))
                        end
                        cdivr,cdivi=general_cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi)
                        H[(i)+1,(n-1)+1] = cdivr
                        H[(i)+1,(n)+1] = cdivi
                        if (abs(x) > (abs(z) + abs(q)))
                            H[(i+1)+1,(n-1)+1] = (-ra - w * H[(i)+1,(n-1)+1] + q * H[(i)+1,(n)+1]) / x
                            H[(i+1)+1,(n)+1] = (-sa - w * H[(i)+1,(n)+1] - q * H[(i)+1,(n-1)+1]) / x
                        else
                            cdivr,cdivi=general_cdiv(-r-y*H[(i)+1,(n-1)+1],-s-y*H[(i)+1,(n)+1],z,q)
                            H[(i+1)+1,(n-1)+1] = cdivr
                            H[(i+1)+1,(n)+1] = cdivi
                        end
                    end
                    t = max(abs(H[(i)+1,(n-1)+1]),abs(H[(i)+1,(n)+1]))
                    if ((eps * t) * t > 1)
                        for j in (i):(n)
                            H[(j)+1,(n-1)+1] = H[(j)+1,(n-1)+1] / t
                            H[(j)+1,(n)+1] = H[(j)+1,(n)+1] / t
                        end
                    end
                end
            end
        end
    end
    for i in (0):(nn-1)
        if (i < low || i > high)
            for j in (i):(nn-1)
                V[(i)+1,(j)+1] = H[(i)+1,(j)+1]
            end
        end
    end
    for j in (nn-1):-1:(low)
        for i in (low):(high)
            z = 0.0
            for k in (low):(min(j,high))
                z = z + V[(i)+1,(k)+1] * H[(k)+1,(j)+1]
            end
            V[(i)+1,(j)+1] = z
        end
    end
end

function general_balance!(H)
    n=size(H,1);exponents=zeros(Int,n)
    function lognorm(i,row)
        largest=0.0
        for j in 1:n
            i==j && continue
            largest=max(largest,abs(row ? H[i,j] : H[j,i]))
        end
        largest==0 && return -Inf
        total=0.0
        for j in 1:n
            i==j && continue
            total+=abs(row ? H[i,j] : H[j,i])/largest
        end
        log2(largest)+log2(total)
    end
    logadd(a,b)=max(a,b)+log2(exp2(a-max(a,b))+exp2(b-max(a,b)))
    for _ in 1:64
        changed=false
        for i in 1:n
            r,c=lognorm(i,true),lognorm(i,false)
            isfinite(r)&&isfinite(c) || continue
            shift=clamp(floor(Int,(r-c)/2+.5),-512,512)
            (shift==0 || logadd(r-shift,c+shift)>=logadd(r,c)+log2(.95)) && continue
            safe=true
            for j in 1:n
                i==j && continue
                for (v,power) in ((H[i,j],-shift),(H[j,i],shift))
                    scaled=ldexp(v,power)
                    if !isfinite(scaled)||(v!=0&&(scaled==0||ldexp(scaled,-power)!=v));safe=false;end
                end
            end
            safe || continue
            for j in 1:n
                i==j && continue
                H[i,j]=ldexp(H[i,j],-shift);H[j,i]=ldexp(H[j,i],shift)
            end
            exponents[i]+=shift;changed=true
        end
        changed || break
    end
    exponents
end

function general_unbalance!(vectors,col,exponents)
    common=typemin(Int)
    for row in axes(vectors,1)
        v=max(abs(real(vectors[row,col])),abs(imag(vectors[row,col])))
        v==0 && continue
        common=max(common,frexp(v)[2]-1+exponents[row])
    end
    common==typemin(Int) && return
    for row in axes(vectors,1)
        power=exponents[row]-common
        vectors[row,col]=complex(ldexp(real(vectors[row,col]),power),ldexp(imag(vectors[row,col]),power))
    end
end
