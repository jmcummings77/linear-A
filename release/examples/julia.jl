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
