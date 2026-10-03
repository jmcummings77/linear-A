using System;
using System.Linq;
using linear_A;
var a = new Matrix<double>(new double[,] {{4,1},{2,3}});
var b = new Matrix<double>(new double[,] {{6},{8}});
var x = a.Solve(b);
var r = a * x;
for (var i=0;i<2;i++)
    if (!double.IsFinite(x[i,0]) || Math.Abs(x[i,0]-(i+1))>1e-12
        || !double.IsFinite(r[i,0]) || Math.Abs(r[i,0]-b[i,0])>1e-12) throw new Exception("incorrect solution");
Console.WriteLine("solution: 1, 2");

var decomposition = a.Svd();
if (decomposition.Values.Length != 2 || decomposition.Values[1] <= 0) throw new Exception("incorrect SVD");

var inverse = a.Pseudoinverse(); var minimum = a.SolveMinimumNorm(b); var diagnostics = a.SpectralDiagnostics();
if (Math.Abs(inverse[0,0]-.3)>1e-12 || Math.Abs(minimum[1,0]-2)>1e-12 || diagnostics.Rank!=2) throw new Exception("incorrect SVD inverse");

var ridge = a.SolveRidge(b, 1);
if (Math.Abs(ridge[0,0]-140.0/131)>1e-12 || Math.Abs(ridge[1,0]-230.0/131)>1e-12) throw new Exception("incorrect ridge solution");

var sparse = new CSRMatrix(2,2,new[]{0,2,4},new[]{0,1,0,1},new[]{4.0,1,1,3});
var product = sparse.Matvec(new[]{1.0,2});
var cg = sparse.ConjugateGradient(new[]{6.0,7},jacobi:true,capture:true);
if(product[0]!=6 || product[1]!=7 || !cg.Converged || Math.Abs(cg.X[0]-1)>1e-12 || Math.Abs(cg.X[1]-2)>1e-12)throw new Exception("incorrect sparse solver");

var gm = sparse.Gmres(new[]{6.0,7}, restart:2, jacobi:true, capture:true);
if(!gm.Converged || Math.Abs(gm.X[0]-1)>1e-12 || Math.Abs(gm.X[1]-2)>1e-12)throw new Exception("incorrect GMRES");

var ilu=new ILU0(sparse);
var second=ilu.Apply(new[]{11.0,13});
if(Math.Abs(ilu.Apply(new[]{6.0,7})[0]-1)>1e-12 || Math.Abs(second[0]-20.0/11)>1e-12 || !sparse.Gmres(new[]{11.0,13},preconditioner:ilu).Converged)throw new Exception("incorrect ILU reuse");

var order=sparse.ReverseCuthillMcKee();var reordered=sparse.PermuteSymmetric(order);var permuted=CSRMatrix.PermuteVector(order,new[]{1.0,2});
if(!CSRMatrix.PermuteVector(order,reordered.Matvec(permuted),true).SequenceEqual(new[]{6.0,7}))throw new Exception("incorrect permutation");

var plan=new SparseCholeskySymbolic(sparse);var chol=plan.Factorize(sparse);
if(Math.Abs(chol.Solve(new[]{6.0,7})[0]-1)>1e-12||Math.Abs(chol.Solve(new[]{11.0,13})[0]-20.0/11)>1e-12||chol.Lower.NNZ!=3)throw new Exception("incorrect Cholesky reuse");
