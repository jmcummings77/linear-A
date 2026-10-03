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
