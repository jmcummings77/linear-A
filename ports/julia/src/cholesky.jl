"""Owned symbolic pattern; fields are read-only by contract. Steps use -1 for originals."""
struct SparseCholeskySymbolic
    size::Int
    source_offsets::Vector{Int}
    source_indices::Vector{Int}
    offsets::Vector{Int}
    indices::Vector{Int}
    fill_steps::Vector{Int}
    function SparseCholeskySymbolic(a::CSRMatrix)
        validate_csr(a.rows,a.cols,a.offsets,a.indices,a.values)
        a.rows==a.cols || throw(ArgumentError("Cholesky requires square matrix"))
        n=a.rows;g=[Dict{Int,Int}() for _ in 1:n]
        for i in 1:n, p in a.offsets[i]+1:a.offsets[i+1]
            j=a.indices[p]+1
            if i!=j;g[i][j]=-1;g[j][i]=-1;end
        end
        for k in 1:n
            ns=sort([j for j in keys(g[k]) if j>k])
            for u in eachindex(ns), w in 1:u-1
                i,j=ns[u],ns[w]
                if !haskey(g[i],j);g[i][j]=k-1;g[j][i]=k-1;end
            end
        end
        rp=[0];ci=Int[];steps=Int[]
        for i in 1:n
            for j in sort([j for j in keys(g[i]) if j<i]);push!(ci,j-1);push!(steps,g[i][j]);end
            push!(ci,i-1);push!(steps,-1);push!(rp,length(ci))
        end
        new(n,copy(a.offsets),copy(a.indices),rp,ci,steps)
    end
end
"""Owned lower factor. Fields are read-only by contract."""
struct SparseCholesky
    lower::CSRMatrix
    function SparseCholesky(l::CSRMatrix)
        validate_csr(l.rows,l.cols,l.offsets,l.indices,l.values)
        l.rows==l.cols || throw(ArgumentError("invalid lower factor"))
        for i in 1:l.rows
            d=l.offsets[i+1]
            d>l.offsets[i] && l.indices[d]==i-1 && l.values[d]>0 || throw(ArgumentError("invalid lower factor"))
        end
        new(CSRMatrix(l.rows,l.cols,l.offsets,l.indices,l.values))
    end
end
function cholesky_factorize(s::SparseCholeskySymbolic,a::CSRMatrix)
    validate_csr(a.rows,a.cols,a.offsets,a.indices,a.values)
    n=s.size
    a.rows==a.cols==n && a.offsets==s.source_offsets && a.indices==s.source_indices || throw(ArgumentError("Cholesky symbolic pattern mismatch"))
    rows=[Dict(a.indices[p]=>a.values[p] for p in a.offsets[i]+1:a.offsets[i+1]) for i in 1:n]
    for i in 1:n, (j,v) in rows[i]
        v==get(rows[j+1],i-1,0.) || throw(ArgumentError("Cholesky requires symmetric values"))
    end
    rp,ci=s.offsets,s.indices
    validate_csr(n,n,rp,ci,zeros(length(ci)))
    for i in 1:n
        rp[i+1]>rp[i] && ci[rp[i+1]]==i-1 || throw(ArgumentError("invalid symbolic diagonal"))
    end
    v=zeros(length(ci))
    for i in 1:n, p in rp[i]+1:rp[i+1]
        j=ci[p]+1;t=get(rows[i],j-1,0.);u,w=rp[i]+1,rp[j]+1
        while u<p && w<rp[j+1]
            if ci[u]==ci[w];t-=v[u]*v[w];u+=1;w+=1
            elseif ci[u]<ci[w];u+=1
            else;w+=1;end
        end
        isfinite(t) || throw(OverflowError("nonfinite Cholesky factor"))
        if i==j
            t>0 || throw(ArgumentError("nonpositive Cholesky pivot"))
            v[p]=sqrt(t)
        else;v[p]=t/v[rp[j+1]];end
        isfinite(v[p]) || throw(OverflowError("nonfinite Cholesky factor"))
    end
    SparseCholesky(CSRMatrix(n,n,rp,ci,v))
end
function cholesky_solve(f::SparseCholesky,b)
    l=f.lower;n=l.rows
    validate_csr(l.rows,l.cols,l.offsets,l.indices,l.values)
    length(b)==n && all(isfinite,b) && l.cols==n || throw(ArgumentError("invalid Cholesky right-hand side"))
    rp,ci,v=l.offsets,l.indices,l.values
    for i in 1:n
        rp[i+1]>rp[i] && ci[rp[i+1]]==i-1 && v[rp[i+1]]>0 || throw(ArgumentError("invalid lower factor"))
    end
    x=Float64.(b)
    for i in 1:n
        for p in rp[i]+1:rp[i+1]-1;x[i]-=v[p]*x[ci[p]+1];end
        x[i]/=v[rp[i+1]]
        isfinite(x[i]) || throw(OverflowError("nonfinite Cholesky solve"))
    end
    for i in n:-1:1
        x[i]/=v[rp[i+1]]
        isfinite(x[i]) || throw(OverflowError("nonfinite Cholesky solve"))
        for p in rp[i]+1:rp[i+1]-1
            x[ci[p]+1]-=v[p]*x[i]
            isfinite(x[ci[p]+1]) || throw(OverflowError("nonfinite Cholesky solve"))
        end
    end
    x
end
