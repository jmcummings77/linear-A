import { createMatrixAPI } from "linear-a-wasm";
const { Matrix } = await createMatrixAPI();
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
  console.log("solution: 1, 2");
} finally { reconstructed.dispose(); x.dispose(); b.dispose(); a.dispose(); }
