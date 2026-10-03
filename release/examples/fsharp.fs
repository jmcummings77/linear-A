open System
open LinearA
let a = Matrix.FromArray(2,2,[|4.;1.;2.;3.|])
let b = Matrix.FromArray(2,1,[|6.;8.|])
let x = a.Solve(b)
let r = a.Multiply(x)
for i in 0..1 do
    if not (Double.IsFinite x[i,0]) || abs(x[i,0]-float(i+1))>1e-12 || not (Double.IsFinite r[i,0]) || abs(r[i,0]-b[i,0])>1e-12 then
        failwith "incorrect solution"
printfn "solution: 1, 2"

let decomposition = a.Svd()
if decomposition.Values.Length <> 2 || decomposition.Values[1] <= 0. then failwith "incorrect SVD"

let inverse = a.Pseudoinverse()
let minimum = a.SolveMinimumNorm(b)
if abs(inverse[0,0]-0.3)>1e-12 || abs(minimum[1,0]-2.)>1e-12 || (a.SpectralDiagnostics()).Rank<>2 then failwith "incorrect SVD inverse"
