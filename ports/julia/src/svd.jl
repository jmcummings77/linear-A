"""Economy one-sided Jacobi SVD, returning `(u, values, vt)`; see ports/SVD.md."""
function svd(a::Matrix64; tolerance=1e-12, max_sweeps=100)
    isfinite(tolerance) && 0<tolerance<1 && max_sweeps isa Integer && !(max_sweeps isa Bool) && 1<=max_sweeps<=10000 || throw(ArgumentError("invalid SVD options"))
    m,n=size(a)
    if m<n
        r=svd(transpose(a); tolerance=tolerance,max_sweeps=max_sweeps)
        return (u=transpose(r.vt),values=r.values,vt=transpose(r.u))
    end
    all(isfinite,a.data)||throw(ArgumentError("SVD requires finite input"))
    scaling=isempty(a.data) ? 1.0 : maximum(abs,a.data)
    scaling==0 && (scaling=1.0)
    b=a.data ./ scaling; v=rowmajor(identitymatrix(n))
    any(i->a.data[i]!=0&&b[i]==0,eachindex(b)) && throw(ErrorException("SVD scaling discards an entry"))
    norm(j)=foldl(hypot,(b[(i-1)*n+j] for i in 1:m);init=0.0)
    converged=false
    for sweep in 0:max_sweeps
        changed=false
        for p in 1:n, q in p+1:n
            np,nq=norm(p),norm(q)
            (np==0||nq==0)&&continue
            corr=0.0
            for i in 1:m;corr+=(b[(i-1)*n+p]/np)*(b[(i-1)*n+q]/nq);end
            abs(corr)<=tolerance&&continue
            changed=true;sweep==max_sweeps&&continue
            pair=max(np,nq);ap,aq=np/pair,nq/pair;delta,g=aq*aq-ap*ap,2*ap*aq*corr
            t=delta==0 ? copysign(1.0,g) : g/(delta+copysign(hypot(delta,g),delta))
            if abs(t)<floatmin(Float64)
                small=np<nq ? p : q
                for i in 1:m;b[(i-1)*n+small]=0;end
                continue
            end
            c=1/hypot(1,t);s=c*t
            for (data,count) in ((b,m),(v,n)), i in 1:count
                x,y=data[(i-1)*n+p],data[(i-1)*n+q]
                data[(i-1)*n+p],data[(i-1)*n+q]=c*x-s*y,s*x+c*y
            end
        end
        if !changed;converged=true;break;end
    end
    converged||throw(ErrorException("SVD did not converge"))
    norms=[norm(j) for j in 1:n];order=sortperm(norms;rev=true)
    u=zeros(m*n);vt=zeros(n*n);values=zeros(n)
    for (j,k) in enumerate(order)
        values[j]=norms[k]*scaling
        (isfinite(values[j])&&!(norms[k]!=0&&values[j]==0))||throw(OverflowError("singular value outside float64 range"))
        for i in 1:n;vt[(j-1)*n+i]=v[(i-1)*n+k];end
        if norms[k]!=0
            for i in 1:m;u[(i-1)*n+j]=b[(i-1)*n+k]/norms[k];end
        else
            found=false
            for axis in 1:m
                candidate=zeros(m);candidate[axis]=1
                for pass in 1:2,col in 1:j-1
                    dot=sum(candidate[i]*u[(i-1)*n+col] for i in 1:m)
                    for i in 1:m;candidate[i]-=dot*u[(i-1)*n+col];end
                end
                len=foldl(hypot,candidate;init=0.0)
                if len>0.5/sqrt(m)
                    for i in 1:m;u[(i-1)*n+j]=candidate[i]/len;end
                    found=true;break
                end
            end
            found||throw(ErrorException("cannot complete SVD null basis"))
        end
    end
    (u=Matrix64(m,n,u),values=values,vt=Matrix64(n,n,vt))
end

function spectral_cutoff(a::Matrix64, cutoff)
    value = cutoff === nothing ? max(a.rows,a.cols)*eps(Float64) : Float64(cutoff)
    isfinite(value) && 0 <= value <= 1 || throw(ArgumentError("invalid relative cutoff"))
    value
end
spectral_rank(s,cutoff) = count(x -> x>0 && (cutoff==0 || x/s[1]>cutoff), s)
function spectral_diagnostics(a::Matrix64; relative_cutoff=nothing)
    cutoff=spectral_cutoff(a,relative_cutoff);s=svd(a).values;rank=spectral_rank(s,cutoff)
    (rank=rank,reciprocal_condition=isempty(s)||s[1]==0 ? 0.0 : s[end]/s[1],
     retained_reciprocal_condition=rank==0 ? 0.0 : s[rank]/s[1])
end
function inverse_product(a,b,c)
    fa,ea=frexp(a);fb,eb=frexp(b);fc,ec=frexp(c)
    finite(ldexp((fa/fb)*fc,ea-eb+ec))
end
function apply_inverse(a::Matrix64,rhs,cutoff,lambda=0.0)
    cutoff=spectral_cutoff(a,cutoff)
    rhs === nothing || rhs.rows==a.rows || throw(DimensionMismatch("incompatible right-hand side"))
    rhs === nothing || all(isfinite,rhs.data) || throw(ArgumentError("expected finite right-hand side"))
    r=svd(a);rank=spectral_rank(r.values,cutoff);m,n=a.rows,a.cols
    cols=rhs === nothing ? m : rhs.cols;out=Matrix64(n,cols)
    for j in 1:cols
        scaling=rhs === nothing ? 1.0 : maximum((abs(rhs[i,j]) for i in 1:m);init=0.0)
        if rhs !== nothing && rank>0 && scaling!=0 && any(rhs[i,j]!=0 && rhs[i,j]/scaling==0 for i in 1:m)
            error("right-hand side scaling discards an entry")
        end
        for p in 1:rank
            projection=rhs === nothing ? r.u[j,p] : scaling==0 ? 0.0 : sum((r.u[i,p]*(rhs[i,j]/scaling) for i in 1:m);init=0.0)
            coefficient=rhs === nothing ? 0.0 : lambda>0 ? ridge_product(projection,r.values[p],scaling,lambda) : inverse_product(projection,r.values[p],scaling)
            for i in 1:n
                term=rhs === nothing ? inverse_product(projection,r.values[p],r.vt[p,i]) : r.vt[p,i]*coefficient
                out[i,j]=finite(out[i,j]+term)
            end
        end
    end
    out
end
pseudoinverse(a::Matrix64;relative_cutoff=nothing)=apply_inverse(a,nothing,relative_cutoff)
solve_minimum_norm(a::Matrix64,b::Matrix64;relative_cutoff=nothing)=apply_inverse(a,b,relative_cutoff)

function ridge_product(p,s,b,lambda)
    (p==0 || s==0 || b==0) && return 0.0
    root=sqrt(lambda);d=max(s,root);h=(s/d)^2+(root/d)^2
    pf,pe=frexp(p);sf,se=frexp(s);bf,be=frexp(b);df,de=frexp(d)
    finite(ldexp(pf*sf*bf/(df*df*h),pe+se+be-2*de))
end
function solve_ridge(a::Matrix64,b::Matrix64,lambda::Real)
    lambda=Float64(lambda)
    isfinite(lambda) && lambda>=0 || throw(ArgumentError("lambda must be finite and nonnegative"))
    apply_inverse(a,b,lambda==0 ? nothing : 0.0,lambda)
end
