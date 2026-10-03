"""Restarted GMRES, right Jacobi preconditioning, and true residual stopping."""
function gmres(a::CSRMatrix,b;restart=30,rtol=1e-10,atol=0.0,max_iterations=1000,jacobi=false,capture=false,preconditioner=nothing)
    validate_csr(a.rows,a.cols,a.offsets,a.indices,a.values)
    n=a.rows;b=Float64.(b)
    a.cols==n && length(b)==n && all(isfinite,b) || throw(ArgumentError("GMRES requires square matrix and finite matching vector"))
    restart isa Integer && 1<=restart<=1024 && max_iterations isa Integer && 0<=max_iterations<=100000 && isfinite(rtol) && 0<=rtol<1 && isfinite(atol) && atol>=0 && jacobi isa Bool && capture isa Bool || throw(ArgumentError("invalid GMRES options"))
    preconditioner===nothing || (preconditioner isa ILU0 && !jacobi && preconditioner.factors.rows==n) || throw(ArgumentError("incompatible preconditioner"))
    diagonal=ones(n)
    apply(v)=preconditioner===nothing ? v./diagonal : ilu_apply(preconditioner,v)
    if jacobi
        for i in 1:n
            found=false
            for p in a.offsets[i]+1:a.offsets[i+1]
                if a.indices[p]==i-1;diagonal[i]=a.values[p];found=true;break;end
            end
            found && diagonal[i]!=0 || throw(ArgumentError("Jacobi requires nonzero diagonal"))
        end
    end
    norm(v)=foldl(hypot,v;init=0.0)
    x=zeros(n);r=copy(b);history=[norm(r)];estimates=copy(history);frames=capture ? [copy(x)] : Vector{Float64}[];restarts=Int[]
    result(reason)=(x=copy(x),converged=reason=="converged",iterations=length(history)-1,reason=reason,residuals=copy(history),estimated_residuals=copy(estimates),iterates=copy(frames),restarts=copy(restarts))
    threshold=max(atol,rtol*history[1]);m=min(restart,n,max_iterations)
    isfinite(history[1]) || return result("nonfinite")
    history[1]<=threshold && return result("converged")
    while length(history)-1<max_iterations
        length(history)>1 && push!(restarts,length(history)-1)
        base=copy(x);beta=norm(r);basis=[r./beta];h=zeros(m+1,m);cs=zeros(m);sn=zeros(m);g=zeros(m+1);g[1]=beta
        steps=min(m,max_iterations-(length(history)-1))
        for j in 1:steps
            w=try matvec(a,apply(basis[j])) catch e;e isa ArgumentError || e isa OverflowError || (e isa ErrorException && e.msg=="sparse multiplication outside float64 range") || rethrow();return result("nonfinite");end
            original=norm(w)
            for _ in 1:2,k in 1:j
                dot=sum(basis[k][i]*w[i] for i in 1:n);h[k,j]+=dot
                for i in 1:n;w[i]-=dot*basis[k][i];end
            end
            tail=norm(w)
            isfinite(original) && isfinite(tail) && all(isfinite,h[1:j,j]) || return result("nonfinite")
            happy=tail<=8*eps(Float64)*original;h[j+1,j]=happy ? 0.0 : tail
            !happy && push!(basis,w./tail)
            for k in 1:j-1
                top=cs[k]*h[k,j]+sn[k]*h[k+1,j];h[k+1,j]=-sn[k]*h[k,j]+cs[k]*h[k+1,j];h[k,j]=top
            end
            pivot=hypot(h[j,j],h[j+1,j]);isfinite(pivot) || return result("nonfinite");pivot!=0 || return result("breakdown")
            cs[j]=h[j,j]/pivot;sn[j]=h[j+1,j]/pivot;h[j,j]=pivot;h[j+1,j]=0
            g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];y=copy(g[1:j])
            for k in j:-1:1
                h[k,k]!=0 || return result("breakdown")
                total=0.0
                for q in k+1:j;total+=h[k,q]*y[q];end
                y[k]=(y[k]-total)/h[k,k]
            end
            correction=try apply([sum(basis[k][i]*y[k] for k in 1:j) for i in 1:n]) catch e;e isa ArgumentError || e isa OverflowError || rethrow();return result("nonfinite");end
            candidate=base.+correction
            all(isfinite,y) && all(isfinite,candidate) && isfinite(g[j+1]) || return result("nonfinite")
            ax=try matvec(a,candidate) catch e;e isa ArgumentError || e isa OverflowError || (e isa ErrorException && e.msg=="sparse multiplication outside float64 range") || rethrow();return result("nonfinite");end
            residual=b.-ax;resnorm=norm(residual);isfinite(resnorm) || return result("nonfinite")
            x,r=candidate,residual;push!(history,resnorm);push!(estimates,abs(g[j+1]));capture && push!(frames,copy(x))
            resnorm<=threshold && return result("converged")
            happy && return result("breakdown")
        end
        x==base && return result("stagnation")
    end
    result("iteration_limit")
end
