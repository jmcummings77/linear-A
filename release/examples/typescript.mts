import { Matrix } from "linear-a-typescript";
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
