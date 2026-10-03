function sparse_run(args)
    length(args)==10 || throw(ArgumentError("invalid sparse protocol"))
    op=args[1];rows,cols,nnz,iterations=parse.(Int,args[2:5]);rtol,atol=parse.(Float64,args[6:7]);limit,jacobi,capture=parse.(Int,args[8:10])
    op in ("spmv","dense","cg") && min(rows,cols,nnz,iterations)>=0 && jacobi in (0,1) && capture in (0,1) || throw(ArgumentError("invalid sparse options"))
    tokens=split(read(stdin,String));count=op=="cg" ? rows : cols
    length(tokens)==rows+1+2*nnz+count || throw(ArgumentError("incorrect sparse input count"))
    rp=parse.(Int,tokens[1:rows+1]);ci=parse.(Int,tokens[rows+2:rows+1+nnz]);v=parse.(Float64,tokens[rows+2+nnz:rows+1+2*nnz]);b=parse.(Float64,tokens[rows+2+2*nnz:end])
    a=CSRMatrix(rows,cols,rp,ci,v);dense=nothing;right=nothing
    if op=="dense"
        dense=Matrix64(rows,cols);right=Matrix64(cols,1,b)
        for i in 1:rows,p in rp[i]+1:rp[i+1];dense[i,ci[p]+1]=v[p];end
    end
    function compute()
        op=="spmv" && return matvec(a,b)
        op=="dense" && return rowmajor(dense*right)
        r=conjugate_gradient(a,b;rtol=rtol,atol=atol,max_iterations=limit,jacobi=jacobi!=0,capture=capture!=0)
        if iterations>0
            r.converged || error("benchmark CG did not converge: "*r.reason)
            return r.x
        end
        reason=findfirst(==(r.reason),["converged","iteration_limit","breakdown","nonfinite"])-1
        vcat([Float64(reason),Float64(r.iterations),Float64(length(r.residuals))],r.x,r.residuals,r.iterates...)
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
