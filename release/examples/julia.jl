using LinearAMatrices
a=Matrix64(2,2,[4.,1.,2.,3.]); b=Matrix64(2,1,[6.,8.])
x=solve(a,b)
@assert maximum(abs.(rowmajor(x) .- [1.,2.])) < 1e-12
@assert maximum(abs.(rowmajor(a*x) .- rowmajor(b))) < 1e-12
println("solution: 1, 2")
