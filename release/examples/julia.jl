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
