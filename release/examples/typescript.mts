import { ILU0, Matrix } from "linear-a-typescript";
const a = new Matrix(2, 2, [4, 1, 2, 3]);
const b = new Matrix(2, 1, [6, 8]);
const x: Matrix = a.solve(b);
for (const [actual, expected] of [[x.values, [1, 2]], [a.multiply(x).values, b.values]]) {
  if (actual.some((v, i) => Math.abs(v - expected[i]) > 1e-12 || !Number.isFinite(v))) throw Error("incorrect solution");
}
console.log("solution: 1, 2");

const decomposition = a.svd();
if (decomposition.values.length !== 2 || decomposition.values[1] <= 0) throw Error("incorrect SVD");

const inverse=a.pseudoinverse(), minimum=a.solveMinimumNorm(b), diagnostics=a.spectralDiagnostics();
if(Math.abs(inverse.get(0,0)-.3)>1e-12 || Math.abs(minimum.get(1,0)-2)>1e-12 || diagnostics.rank!==2) throw Error("incorrect SVD inverse");

const ridge=a.solveRidge(b,1);
if(Math.abs(ridge.get(0,0)-140/131)>1e-12 || Math.abs(ridge.get(1,0)-230/131)>1e-12)throw Error("incorrect ridge solution");

import {CSRMatrix} from "linear-a-typescript";
const sparse=new CSRMatrix(2,2,[0,2,4],[0,1,0,1],[4,1,1,3]);
const product=sparse.matvec([1,2]),cg=sparse.conjugateGradient([6,7],{jacobi:true,capture:true});
if(product[0]!==6||product[1]!==7||!cg.converged||Math.abs(cg.x[0]-1)>1e-12||Math.abs(cg.x[1]-2)>1e-12)throw Error("incorrect sparse solver");

const gm=sparse.gmres([6,7],{restart:2,jacobi:true,capture:true});
if(!gm.converged||Math.abs(gm.x[0]-1)>1e-12||Math.abs(gm.x[1]-2)>1e-12)throw Error("incorrect GMRES");

const ilu=new ILU0(sparse);
if(Math.abs(ilu.apply([6,7])[0]-1)>1e-12||Math.abs(ilu.apply([11,13])[0]-20/11)>1e-12||!sparse.gmres([11,13],{preconditioner:ilu}).converged)throw Error("incorrect ILU reuse");

const order=sparse.reverseCuthillMcKee(),reordered=sparse.permuteSymmetric(order),permuted=CSRMatrix.permuteVector(order,[1,2]);
if(CSRMatrix.permuteVector(order,reordered.matvec(permuted),true).join()!=="6,7")throw Error("incorrect permutation");

import {SparseCholeskySymbolic} from "linear-a-typescript";
if(JSON.stringify(sparse.approximateMinimumDegree())!=="[0,1]")throw Error("incorrect AMD ordering");
const plan=new SparseCholeskySymbolic(sparse),chol=plan.factorize(sparse);
if(Math.abs(chol.solve([6,7])[0]-1)>1e-12||Math.abs(chol.solve([11,13])[0]-20/11)>1e-12||plan.fillCount!==0||chol.lower.nnz!==3)throw Error("incorrect Cholesky reuse");

import {IC0} from "linear-a-typescript";
const ic=new IC0(sparse);
for(const rhs of [[6,7],[11,13]]){const pcg=sparse.conjugateGradient(rhs,{preconditioner:ic});if(!pcg.converged||pcg.iterations!==1)throw Error("incorrect IC0 reuse");}
