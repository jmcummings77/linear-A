using LinearAMatrices
a=Matrix64(2,2,[4.,1.,2.,3.]); b=Matrix64(2,1,[6.,8.])
x=solve(a,b)
@assert maximum(abs.(rowmajor(x) .- [1.,2.])) < 1e-12
@assert maximum(abs.(rowmajor(a*x) .- rowmajor(b))) < 1e-12
println("solution: 1, 2")

r=svd(a)
@assert length(r.values)==2 && r.values[2]>0

p=pseudoinverse(a); minimum=solve_minimum_norm(a,b); d=spectral_diagnostics(a)
@assert abs(p[1,1]-.3)<1e-12 && abs(minimum[2,1]-2)<1e-12 && d.rank==2

ridge=solve_ridge(a,b,1)
@assert abs(ridge[1,1]-140/131)<1e-12 && abs(ridge[2,1]-230/131)<1e-12

sparse=CSRMatrix(2,2,[0,2,4],[0,1,0,1],[4.,1.,1.,3.])
@assert matvec(sparse,[1.,2.])==[6.,7.]
cg=conjugate_gradient(sparse,[6.,7.];jacobi=true,capture=true)
@assert cg.converged && maximum(abs.(cg.x .- [1.,2.]))<1e-12

gm=gmres(sparse,[6.,7.];restart=2,jacobi=true,capture=true)
@assert gm.converged && maximum(abs.(gm.x .- [1.,2.]))<1e-12

f=ILU0(sparse)
for (rhs,expected) in [([6.,7.],[1.,2.]),([11.,13.],[20/11,41/11])]
    @assert maximum(abs.(ilu_apply(f,rhs).-expected))<1e-12
    @assert gmres(sparse,rhs;preconditioner=f).converged
end

order=reverse_cuthill_mckee(sparse)
reordered=permute_symmetric(sparse,order)
@assert permute_vector(order,matvec(reordered,permute_vector(order,[1.,2.])),inverse=true)==[6.,7.]
