import { createMatrixAPI } from "linear-a-wasm";
const { Matrix, CSRMatrix, ILU0 } = await createMatrixAPI();
const a = new Matrix(2, 2, [4, 1, 2, 3]);
const b = new Matrix(2, 1, [6, 8]);
const x = a.solve(b);
const reconstructed = a.multiply(x);
try {
  for (const [actual, expected] of [[x.toArray(), [1, 2]], [reconstructed.toArray(), b.toArray()]]) {
    if (actual.some((v, i) => Math.abs(v - expected[i]) > 1e-12 || !Number.isFinite(v))) throw Error("incorrect solution");
  }
  const decomposition=a.svd();
  try { if (decomposition.values.length!==2 || decomposition.values[1]<=0) throw Error("incorrect SVD"); }
  finally { decomposition.u.dispose(); decomposition.vt.dispose(); }
  const inverse=a.pseudoinverse(), minimum=a.solveMinimumNorm(b);
  try {if(Math.abs(inverse.get(0,0)-.3)>1e-12 || Math.abs(minimum.get(1,0)-2)>1e-12 || a.spectralDiagnostics().rank!==2) throw Error("incorrect SVD inverse");}
  finally {inverse.dispose();minimum.dispose();}
  const ridge=a.solveRidge(b,1);
  try {if(Math.abs(ridge.get(0,0)-140/131)>1e-12 || Math.abs(ridge.get(1,0)-230/131)>1e-12)throw Error("incorrect ridge solution");}finally{ridge.dispose();}
  console.log("solution: 1, 2");
} finally { reconstructed.dispose(); x.dispose(); b.dispose(); a.dispose(); }

const sparse=new CSRMatrix(2,2,[0,2,4],[0,1,0,1],[4,1,1,3]);
try{const ilu=new ILU0(sparse);
try{if(Math.abs(ilu.apply([6,7])[0]-1)>1e-12||Math.abs(ilu.apply([11,13])[0]-20/11)>1e-12||!sparse.gmres([11,13],{preconditioner:ilu}).converged)throw Error("incorrect ILU reuse");}finally{ilu.dispose();}
const gm=sparse.gmres([6,7],{restart:2,jacobi:true,capture:true});
if(!gm.converged||Math.abs(gm.x[0]-1)>1e-12||Math.abs(gm.x[1]-2)>1e-12)throw Error("incorrect GMRES");
const product=sparse.matvec([1,2]),cg=sparse.conjugateGradient([6,7],{jacobi:true,capture:true});
if(product[0]!==6||product[1]!==7||!cg.converged||Math.abs(cg.x[0]-1)>1e-12||Math.abs(cg.x[1]-2)>1e-12)throw Error("incorrect sparse solver");}finally{sparse.dispose();}
