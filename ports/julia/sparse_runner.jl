function sparse_run(args)
    length(args)==(length(args)>0 && args[1]=="gmres" ? 11 : 10) || throw(ArgumentError("invalid sparse protocol"))
    op=args[1];rows,cols,nnz,iterations=parse.(Int,args[2:5]);rtol,atol=parse.(Float64,args[6:7]);limit,jacobi,capture=parse.(Int,args[8:10])
    op in ("mg_setup","mg_matrix","mg_apply","spmv","dense","cg","gmres","ic0_factor","ic0_apply","ilu_setup","ilu_apply","rcm","amd","permute","permutation_check","rcm_solve","ilu_solve","chol_symbolic","chol_factor","chol_solve","chol_total","chol_rcm_total","chol_amd_total") && min(rows,cols,nnz,iterations)>=0 && jacobi in (op=="cg" ? (0,1,2,3,4,5) : op=="gmres" ? (0,1,2,3) : (0,1)) && capture in (0,1) || throw(ArgumentError("invalid sparse options"))
    tokens=split(read(stdin,String));count=op in ("cg","gmres") ? rows : cols
    length(tokens)==rows+1+2*nnz+count || throw(ArgumentError("incorrect sparse input count"))
    rp=parse.(Int,tokens[1:rows+1]);ci=parse.(Int,tokens[rows+2:rows+1+nnz]);v=parse.(Float64,tokens[rows+2+nnz:rows+1+2*nnz]);b=parse.(Float64,tokens[rows+2+2*nnz:end])
    a=CSRMatrix(rows,cols,rp,ci,v);dense=nothing;right=nothing
    function new_mg()
        w=isqrt(rows);rows==cols && w*w==rows || throw(ArgumentError("multigrid requires a square grid"));GeometricMultigrid(w)
    end
    mg=op in ("mg_apply","mg_matrix") || (op=="cg" && jacobi==4) ? new_mg() : nothing
    factor=op=="ilu_apply" || (op=="gmres" && jacobi==2) ? ILU0(a) : nothing
    if op=="dense"
        dense=Matrix64(rows,cols);right=Matrix64(cols,1,b)
        for i in 1:rows,p in rp[i]+1:rp[i+1];dense[i,ci[p]+1]=v[p];end
    end
    ic=op=="ic0_apply" || (op=="cg" && jacobi==2) ? IC0(a) : nothing
    plan=op in ("chol_factor","chol_solve") ? SparseCholeskySymbolic(a) : nothing
    chol=op=="chol_solve" ? cholesky_factorize(plan,a) : nothing
    function compute()
        if op=="mg_setup";m=new_mg();return Float64[multigrid_size(m),multigrid_levels(m)];end
        if op=="mg_matrix";m=multigrid_matrix(mg);return vcat(Float64.(m.offsets),Float64.(m.indices),m.values);end
        op=="mg_apply" && return multigrid_apply(mg,b)
        if op=="ic0_factor";l=ic0_lower(IC0(a));return vcat(Float64.(l.offsets),Float64.(l.indices),l.values);end
        op=="ic0_apply" && return ic0_apply(ic,b)
        if op=="chol_symbolic";s=SparseCholeskySymbolic(a);return Float64.(vcat(s.offsets,s.indices,s.fill_steps));end
        if op=="chol_factor";l=cholesky_factorize(plan,a).lower;return vcat(Float64.(l.offsets),Float64.(l.indices),l.values);end
        op=="chol_solve" && return cholesky_solve(chol,b)
        if op in ("chol_total","chol_rcm_total","chol_amd_total")
            p=op=="chol_amd_total" ? approximate_minimum_degree(a) : op=="chol_rcm_total" ? reverse_cuthill_mckee(a) : nothing
            q=isnothing(p) ? a : permute_symmetric(a,p);rhs=isnothing(p) ? b : permute_vector(p,b)
            x=cholesky_solve(cholesky_factorize(SparseCholeskySymbolic(q),q),rhs)
            return isnothing(p) ? x : permute_vector(p,x,inverse=true)
        end
        if op in ("rcm_solve","ilu_solve")
            p=op=="rcm_solve" ? reverse_cuthill_mckee(a) : collect(0:rows-1)
            q=op=="rcm_solve" ? permute_symmetric(a,p) : a;rhs=op=="rcm_solve" ? permute_vector(p,b) : b
            result=gmres(q,rhs;restart=20,rtol=rtol,atol=atol,max_iterations=limit,preconditioner=ILU0(q))
            result.converged || error("ordering solve failed: "*result.reason)
            x=op=="rcm_solve" ? permute_vector(p,result.x,inverse=true) : result.x
            return vcat([Float64(result.iterations)],x)
        end
        op=="amd" && return Float64.(approximate_minimum_degree(a))
        op=="rcm" && return Float64.(reverse_cuthill_mckee(a))
        if op in ("permute","permutation_check")
            all(x->isfinite(x)&&0<=x<rows&&isinteger(x),b) || throw(ArgumentError("invalid permutation"))
            p=Int.(b);q=permute_symmetric(a,p)
            op=="permute" && return vcat(Float64.(q.offsets),Float64.(q.indices),q.values)
            x=Float64.(1:rows);y=permute_vector(p,x)
            return vcat(y,permute_vector(p,y,inverse=true),permute_vector(p,matvec(q,y),inverse=true))
        end
        if op=="ilu_setup";f=ILU0(a);return Float64[f.factors.rows,length(f.factors.values)];end
        op=="ilu_apply" && return ilu_apply(factor,b)
        op=="spmv" && return matvec(a,b)
        op=="dense" && return rowmajor(dense*right)
        r=op=="gmres" ? gmres(a,b;restart=parse(Int,args[11]),rtol=rtol,atol=atol,max_iterations=limit,jacobi=jacobi==1,capture=capture!=0,preconditioner=jacobi==3 ? ILU0(a) : factor) : conjugate_gradient(a,b;rtol=rtol,atol=atol,max_iterations=limit,jacobi=jacobi==1,capture=capture!=0,preconditioner=jacobi==5 ? new_mg() : jacobi==4 ? mg : jacobi==3 ? IC0(a) : ic)
        if iterations>0
            r.converged || error("benchmark CG did not converge: "*r.reason)
            return r.x
        end
        reason=findfirst(==(r.reason),["converged","iteration_limit","breakdown","nonfinite","stagnation"])-1
        vcat([Float64(reason),Float64(r.iterations),Float64(length(r.residuals))],r.x,r.residuals,(op=="gmres" ? vcat(r.estimated_residuals,[Float64(length(r.restarts))],Float64.(r.restarts)) : Float64[]),r.iterates...)
    end
    if iterations==0
        values=compute();return Matrix64(1,length(values),values)
    end
    sparse_batch(compute, 3) # Compile and warm the same checksum loop before timing.
    sparse_batch(compute, iterations)
end

function sparse_batch(compute::F, iterations) where F
    checksum=0.0;start=time_ns()
    for _ in 1:iterations
        for x in compute();checksum+=x;end
    end
    elapsed=time_ns()-start
    isfinite(checksum)||error("nonfinite checksum")
    (elapsed_ns=elapsed,iterations=iterations,checksum=checksum)
end
