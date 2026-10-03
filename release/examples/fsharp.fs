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

let ridge = a.SolveRidge(b,1.)
if abs(ridge[0,0]-140./131.)>1e-12 || abs(ridge[1,0]-230./131.)>1e-12 then failwith "incorrect ridge solution"

let sparse = CSRMatrix(2,2,[|0;2;4|],[|0;1;0;1|],[|4.;1.;1.;3.|])
if sparse.Matvec([|1.;2.|]) <> [|6.;7.|] then failwith "incorrect sparse product"
let cg = sparse.ConjugateGradient([|6.;7.|],jacobi=true,capture=true)
if not cg.Converged || abs(cg.X[0]-1.)>1e-12 || abs(cg.X[1]-2.)>1e-12 then failwith "incorrect CG"

let gm = sparse.Gmres([|6.;7.|],restart=2,jacobi=true,capture=true)
if not gm.Converged || abs(gm.X[0]-1.)>1e-12 || abs(gm.X[1]-2.)>1e-12 then failwith "incorrect GMRES"

let ilu=ILU0(sparse)
if abs(ilu.Apply([|6.;7.|])[0]-1.)>1e-12 || abs(ilu.Apply([|11.;13.|])[0]-20./11.)>1e-12 || not(sparse.Gmres([|11.;13.|],preconditioner=ilu).Converged) then failwith "incorrect ILU reuse"

let order=sparse.ReverseCuthillMcKee()
let reordered=sparse.PermuteSymmetric(order)
let permuted=CSRMatrix.PermuteVector(order,[|1.;2.|])
if CSRMatrix.PermuteVector(order,reordered.Matvec(permuted),true) <> [|6.;7.|] then failwith "incorrect permutation"

let plan=SparseCholeskySymbolic(sparse)
let chol=plan.Factorize(sparse)
if abs(chol.Solve([|6.;7.|])[0]-1.)>1e-12 || abs(chol.Solve([|11.;13.|])[0]-20./11.)>1e-12 || chol.Lower.NNZ<>3 then failwith "incorrect Cholesky reuse"
