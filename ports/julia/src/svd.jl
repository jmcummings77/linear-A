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
