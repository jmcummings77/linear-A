"""Symmetric V-cycle for the unit five-point Dirichlet Laplacian."""
struct GeometricMultigrid
    width::Int
    function GeometricMultigrid(width::Integer)
        1<=width<=255 && (width&(width+1))==0 || throw(ArgumentError("grid width must be 2^k-1 in 1..255"))
        new(Int(width))
    end
end
multigrid_size(m::GeometricMultigrid)=m.width^2
multigrid_levels(m::GeometricMultigrid)=trailing_zeros(m.width+1)
function multigrid_matrix(m::GeometricMultigrid)
    w=m.width;n=w*w;rp=[0];ci=Int[];v=Float64[]
    for i in 0:n-1
        cols=[i]
        i>=w && push!(cols,i-w);i%w>0 && push!(cols,i-1);i%w+1<w && push!(cols,i+1);i+w<n && push!(cols,i+w)
        for j in sort(cols);push!(ci,j);push!(v,j==i ? 4. : -1.);end
        push!(rp,length(v))
    end
    CSRMatrix(n,n,rp,ci,v)
end
function mg_multiply(w,x)
    n=length(x)
    [4*x[i+1]-(i%w>0 ? x[i] : 0.)-(i%w+1<w ? x[i+2] : 0.)-(i>=w ? x[i-w+1] : 0.)-(i+w<n ? x[i+w+1] : 0.) for i in 0:n-1]
end
function mg_cycle(w,b)
    w==1 && return [b[1]/4]
    x=zeros(length(b))
    function smooth!(x)
        for _ in 1:2;x.+= (b.-mg_multiply(w,x))./6;end
    end
    smooth!(x);r=b.-mg_multiply(w,x);c=w÷2;bc=zeros(c*c)
    for y in 0:c-1,j in 0:c-1,dy in -1:1,dx in -1:1
        bc[y*c+j+1]+=(dy==0 ? 1. : .5)*(dx==0 ? 1. : .5)*r[(2*y+1+dy)*w+2*j+2+dx]
    end
    ec=mg_cycle(c,bc)
    for y in 0:c-1,j in 0:c-1,dy in -1:1,dx in -1:1
        x[(2*y+1+dy)*w+2*j+2+dx]+=(dy==0 ? 1. : .5)*(dx==0 ? 1. : .5)*ec[y*c+j+1]
    end
    smooth!(x);all(isfinite,x) || throw(OverflowError("nonfinite multigrid cycle"));x
end
function multigrid_apply(m::GeometricMultigrid,b)
    length(b)==multigrid_size(m) && all(isfinite,b) || throw(ArgumentError("invalid multigrid RHS"))
    mg_cycle(m.width,Float64.(b))
end
export GeometricMultigrid,multigrid_size,multigrid_levels,multigrid_matrix,multigrid_apply
