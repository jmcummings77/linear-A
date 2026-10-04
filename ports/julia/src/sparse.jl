"""Canonical CSR. Offsets and stored column indices are zero-based in every port."""
struct CSRMatrix
    rows::Int
    cols::Int
    offsets::Vector{Int}
    indices::Vector{Int}
    values::Vector{Float64}
    function CSRMatrix(rows::Integer,cols::Integer,offsets,indices,values)
        rp,ci,v=Int.(offsets),Int.(indices),Float64.(values)
        validate_csr(rows,cols,rp,ci,v)
        new(rows,cols,copy(rp),copy(ci),copy(v))
    end
end
function validate_csr(rows,cols,rp,ci,v)
    rows>=0 && rows<typemax(Int) && cols>=0 && length(rp)==rows+1 && length(ci)==length(v) && rp[1]==0 && rp[end]==length(v) || throw(ArgumentError("invalid CSR dimensions or arrays"))
    all(x->0<=x<=length(v),rp) && all(isfinite,v) || throw(ArgumentError("invalid CSR offsets or values"))
    for i in 1:rows
        rp[i]<=rp[i+1] || throw(ArgumentError("CSR offsets must be monotone"))
        previous=-1
        for p in rp[i]+1:rp[i+1]
            previous<ci[p]<cols || throw(ArgumentError("CSR columns must be sorted unique and in range"))
            previous=ci[p]
        end
    end
end
function csr_from_dense(a::Matrix64)
    rp,ci,v=Int[0],Int[],Float64[]
    for i in 1:a.rows
        for j in 1:a.cols
            if a[i,j]!=0;push!(ci,j-1);push!(v,a[i,j]);end
        end
        push!(rp,length(v))
    end
    CSRMatrix(a.rows,a.cols,rp,ci,v)
end
function matvec(a::CSRMatrix,x)
    validate_csr(a.rows,a.cols,a.offsets,a.indices,a.values)
    length(x)==a.cols && all(isfinite,x) || throw(ArgumentError("invalid vector"))
    out=zeros(a.rows)
    for i in 1:a.rows
        for p in a.offsets[i]+1:a.offsets[i+1];out[i]+=a.values[p]*x[a.indices[p]+1];end
        isfinite(out[i]) || error("sparse multiplication outside float64 range")
    end
    out
end
function csr_find(a,row,col)
    lo,hi=a.offsets[row]+1,a.offsets[row+1]+1
    while lo<hi
        mid=lo+(hi-lo)÷2
        if a.indices[mid]<col;lo=mid+1;else;hi=mid;end
    end
    lo
end
function conjugate_gradient(a::CSRMatrix,b;rtol=1e-10,atol=0.0,max_iterations=1000,jacobi=false,capture=false)
    validate_csr(a.rows,a.cols,a.offsets,a.indices,a.values)
    n=a.rows
    a.cols==n && length(b)==n && all(isfinite,b) || throw(ArgumentError("CG requires square matrix and finite matching vector"))
    isfinite(rtol) && 0<=rtol<1 && isfinite(atol) && atol>=0 && max_iterations isa Integer && 0<=max_iterations<=100000 && jacobi isa Bool && capture isa Bool || throw(ArgumentError("invalid CG options"))
    diagonal=ones(n)
    for i in 1:n
        for p in a.offsets[i]+1:a.offsets[i+1]
            j=a.indices[p]+1;q=csr_find(a,j,i-1)
            other=q<=a.offsets[j+1] && a.indices[q]==i-1 ? a.values[q] : 0.0
            a.values[p]==other || throw(ArgumentError("CG requires exact symmetry"))
        end
        if jacobi
            q=csr_find(a,i,i-1)
            q<=a.offsets[i+1] && a.indices[q]==i-1 && a.values[q]>0 || throw(ArgumentError("Jacobi requires positive diagonal"))
            diagonal[i]=a.values[q]
        end
    end
    norm(v)=foldl(hypot,v;init=0.0)
    dot(a,b)=sum((x*y for (x,y) in zip(a,b));init=0.0)
    x=zeros(n);r=Float64.(b);history=[norm(r)];frames=capture ? [copy(x)] : Vector{Float64}[]
    threshold=max(atol,rtol*history[1])
    result(reason)=(x=copy(x),converged=reason=="converged",iterations=length(history)-1,reason=reason,residuals=copy(history),iterates=copy(frames))
    !isfinite(history[1]) && return result("nonfinite")
    history[1]<=threshold && return result("converged")
    z=r./diagonal;p=copy(z);rho=dot(r,z)
    for _ in 1:max_iterations
        q=try matvec(a,p) catch;return result("nonfinite");end
        curvature=dot(p,q)
        (!isfinite(rho)||!isfinite(curvature)) && return result("nonfinite")
        (rho<=0||curvature<=0) && return result("breakdown")
        alpha=rho/curvature;candidate=x.+alpha.*p
        ax=try matvec(a,candidate) catch;return result("nonfinite");end
        residual=b.-ax;length_residual=norm(residual)
        !isfinite(length_residual) && return result("nonfinite")
        x,r=candidate,residual;push!(history,length_residual)
        capture && push!(frames,copy(x))
        length_residual<=threshold && return result("converged")
        z=r./diagonal;next=dot(r,z)
        !isfinite(next) && return result("nonfinite")
        next<=0 && return result("breakdown")
        beta=next/rho;p=z.+beta.*p;rho=next
    end
    result("iteration_limit")
end
export CSRMatrix,csr_from_dense,matvec,conjugate_gradient

"""Quotient-graph AMD, returning deterministic zero-based new-to-old indices."""
function approximate_minimum_degree(a::CSRMatrix)
 validate_csr(a.rows,a.cols,a.offsets,a.indices,a.values)
 a.rows==a.cols || throw(ArgumentError("AMD requires square matrix"))
 n=a.rows;direct=[Set{Int}() for _ in 1:n];elements=[Set{Int}() for _ in 1:n]
 for i in 1:n, k in a.offsets[i]+1:a.offsets[i+1]
  j=a.indices[k]+1
  if i!=j;push!(direct[i],j);push!(direct[j],i);end
 end
 active=Set(1:n);degree=length.(direct);order=Int[]
 while !isempty(active)
  pivot=first(sort(collect(active),by=i->(min(degree[i],length(active)-1),i)))
  neighbors=copy(direct[pivot])
  for e in elements
   if pivot in e;union!(neighbors,e);empty!(e);end
  end
  delete!(neighbors,pivot);delete!(active,pivot);push!(order,pivot-1);empty!(direct[pivot])
  for i in neighbors;delete!(direct[i],pivot);setdiff!(direct[i],neighbors);end
  elements[pivot]=neighbors
  for i in neighbors
   bound=length(neighbors)-1+length(direct[i])
   for e in 1:n
    if e!=pivot && i in elements[e];bound+=length(setdiff(elements[e],neighbors));end
   end
   degree[i]=min(length(active)-1,bound)
  end
 end
 order
end
export approximate_minimum_degree

"""Deterministic new-to-old RCM indices (zero-based, as in CSR)."""
function reverse_cuthill_mckee(a::CSRMatrix)
 validate_csr(a.rows,a.cols,a.offsets,a.indices,a.values)
 a.rows==a.cols || throw(ArgumentError("RCM requires square matrix"))
 n=a.rows;graph=[Set{Int}() for _ in 1:n]
 for i in 1:n, k in a.offsets[i]+1:a.offsets[i+1]
  j=a.indices[k]+1
  if i!=j;push!(graph[i],j);push!(graph[j],i);end
 end
 key(i)=(length(graph[i]),i);seen=falses(n);order=Int[]
 for start in sort(collect(1:n),by=key)
  seen[start] && continue
  q=[start];seen[start]=true;h=1
  while h<=length(q)
   for j in sort([j for j in graph[q[h]] if !seen[j]],by=key);seen[j]=true;push!(q,j);end
   h+=1
  end
  append!(order,q)
 end
 reverse(order).-1
end
function permute_vector(order,x;inverse=false)
 n=length(x)
 length(order)==n && all(i->i isa Integer && !(i isa Bool) && 0<=i<n,order) && length(Set(order))==n && all(isfinite,x) || throw(ArgumentError("invalid permutation or vector"))
 out=zeros(n)
 for i in 1:n
  if inverse;out[order[i]+1]=x[i];else;out[i]=x[order[i]+1];end
 end
 out
end
function permute_symmetric(a::CSRMatrix,order)
 validate_csr(a.rows,a.cols,a.offsets,a.indices,a.values)
 a.rows==a.cols || throw(ArgumentError("permutation requires square matrix"))
 inv=Int.(permute_vector(order,collect(0:a.rows-1),inverse=true));rp=[0];ci=Int[];v=Float64[]
 for old in order
  i=old+1
  for k in sort(collect(a.offsets[i]+1:a.offsets[i+1]),by=k->inv[a.indices[k]+1]);push!(ci,inv[a.indices[k]+1]);push!(v,a.values[k]);end
  push!(rp,length(v))
 end
 CSRMatrix(a.rows,a.cols,rp,ci,v)
end
export reverse_cuthill_mckee,permute_symmetric,permute_vector
