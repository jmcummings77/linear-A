/** Node check/benchmark runner for the browser-compatible WebAssembly Matrix API. */
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { createMatrixAPI } from "./matrix.mjs";

const determinants = new Set(["determinant", "determinant_lu", "determinant_cholesky", "determinant_spd_lu"]);
const operations = new Set(["pseudoinverse", "spectral_diagnostics", "solve_minimum_norm", "svd", "svd_one_sweep","add", "subtract", "multiply", "solve", "solve_cholesky", "least_squares", "rcond", "cross", "rotation2d", "rotation3d", "scale", "transpose", "trace", ...determinants, "triangular", "eigen_symmetric", "eigen_general"]);
const binary = new Set(["add", "subtract", "multiply", "cross", "solve", "solve_cholesky", "least_squares", "solve_minimum_norm"]);
const rotations = new Set(["rotation2d", "rotation3d"]);

function integer(text, minimum = 0, maximum = Number.MAX_SAFE_INTEGER) {
  if (typeof text !== "string" || !/^[+-]?\d+$/.test(text)) throw new RangeError("expected an integer argument");
  const value = Number(text);
  if (!Number.isSafeInteger(value) || value < minimum || value > maximum) throw new RangeError("integer argument is outside the allowed range");
  return value;
}

function number(text) {
  if (typeof text !== "string" || !/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/.test(text)) throw new RangeError("expected a finite number");
  const value = Number(text);
  if (!Number.isFinite(value)) throw new RangeError("expected a finite number");
  return value;
}

function count(rows, cols) {
  const result = rows * cols;
  if (rows > 0xffffffff || cols > 0xffffffff || !Number.isSafeInteger(result) || result > Math.floor(0xffffffff / 8)) {
    throw new RangeError("matrix dimensions overflow wasm32 storage");
  }
  return result;
}

function execute(op, a, b, scalar = 1.25) {
  switch (op) {
    case "add": return a.add(b);
    case "subtract": return a.subtract(b);
    case "multiply": return a.multiply(b);
    case "pseudoinverse": return a.pseudoinverse(scalar);
    case "solve_minimum_norm": return a.solveMinimumNorm(b);
    case "spectral_diagnostics": {
      const d=a.spectralDiagnostics(scalar);
      return new a.constructor(1,3,[d.rank,d.reciprocalCondition,d.retainedReciprocalCondition]);
    }
    case "svd": case "svd_one_sweep": {
      const r=a.svd(1e-12,op==="svd"?100:1);
      const k=Math.min(a.rows,a.cols),values=new Float64Array((a.rows+1+a.cols)*k);
      try {
        values.set(r.u.toArray()); values.set(r.values,a.rows*k);
        const vt=r.vt.toArray();
        for(let i=0;i<a.cols;i++)for(let j=0;j<k;j++)values[(a.rows+1+i)*k+j]=vt[j*a.cols+i];
        return new a.constructor(a.rows+1+a.cols,k,values);
      } finally { r.u.dispose(); r.vt.dispose(); }
    }
    case "solve": return a.solve(b);
    case "least_squares": return a.leastSquares(b);
    case "solve_cholesky": { const f = a.factorCholesky(); try { return f.solve(b); } finally { f.dispose(); } }
    case "rcond": { const f = a.factorLU(); try { return f.reciprocalCondition(); } finally { f.dispose(); } }
    case "cross": return a.cross(b);
    case "rotation2d":
      if (a.rows || a.cols) throw new RangeError("rotation2d requires empty 0 by 0 input");
      return a.constructor.rotation2D(scalar);
    case "rotation3d": return a.constructor.rotationAxisAngle(a, scalar);
    case "scale": return a.scale(scalar);
    case "transpose": return a.transpose();
    case "trace": return a.trace();
    case "determinant": return a.determinant();
    case "determinant_lu":
    case "determinant_spd_lu": return a.determinant("lu");
    case "determinant_cholesky": return a.determinant("cholesky");
    case "triangular": return a.triangular();
    case "eigen_general": return a.eigenGeneral();
    case "eigen_symmetric": return a.eigenSymmetric();
    default: throw new RangeError("unknown operation");
  }
}

function check(args, Matrix) {
  const op = args[0];
  if (!operations.has(op)) throw new RangeError("unknown operation");
  if (args.length !== (binary.has(op) ? 5 : op === "scale" || rotations.has(op) || op === "pseudoinverse" || op === "spectral_diagnostics" ? 4 : 3)) throw new RangeError("incorrect argument count for operation");
  const rows = integer(args[1]), cols = integer(args[2]);
  const brows = binary.has(op) ? integer(args[3]) : 0;
  const bcols = binary.has(op) ? integer(args[4]) : 0;
  const scalar = op === "scale" || rotations.has(op) || op === "pseudoinverse" || op === "spectral_diagnostics" ? number(args[3]) : 1.25;
  const firstCount = count(rows, cols), secondCount = count(brows, bcols);
  const input = readFileSync(0, "utf8").trim();
  const values = input ? input.split(/\s+/).map(number) : [];
  if (values.length !== firstCount + secondCount) throw new RangeError("input value count does not match matrix dimensions");
  let a, b, result;
  try {
    a = new Matrix(rows, cols, values.slice(0, firstCount));
    if (binary.has(op)) b = new Matrix(brows, bcols, values.slice(firstCount));
    result = execute(op, a, b, scalar);
    if (result instanceof Matrix) return { rows: result.rows, cols: result.cols, values: Array.from(result.toArray()) };
    if (typeof result === "number") return { value: result };
    if (op === "eigen_general") return {
      eigenvalues_real: Array.from(result.valuesReal), eigenvalues_imag: Array.from(result.valuesImag),
      eigenvectors_real: { rows: result.vectorsReal.rows, cols: result.vectorsReal.cols, values: Array.from(result.vectorsReal.toArray()) },
      eigenvectors_imag: { rows: result.vectorsImag.rows, cols: result.vectorsImag.cols, values: Array.from(result.vectorsImag.toArray()) },
    };
    if (op === "eigen_symmetric") return {
      eigenvalues: Array.from(result.values),
      eigenvectors: { rows: result.vectors.rows, cols: result.vectors.cols, values: Array.from(result.vectors.toArray()) },
    };
    return result;
  } finally {
    if (result instanceof Matrix) result.dispose();
    if (op === "eigen_symmetric") result?.vectors?.dispose();
    if (op === "eigen_general") { result?.vectorsReal?.dispose(); result?.vectorsImag?.dispose(); }
    a?.dispose(); b?.dispose();
  }
}

function bench(args, Matrix) {
  if (args.length !== 4 || !operations.has(args[0]) || args[0] === "triangular") throw new RangeError("usage: bench OP SIZE ITERATIONS SEED");
  const op = args[0], size = integer(args[1], 1), iterations = integer(args[2], 1), seed = integer(args[3], 0, 2147483646);
  if ((op === "cross" || op === "rotation3d") && size !== 3) throw new RangeError("cross and rotation3d require benchmark size 3");
  if (op === "rotation2d" && size !== 2) throw new RangeError("rotation2d requires benchmark size 2");
  const rows = op === "rotation2d" ? 0 : size;
  const cols = op === "rotation2d" ? 0 : op === "cross" || op === "rotation3d" ? 1 : size;
  const length = count(rows, cols);
  const leftValues = new Float64Array(length), rightValues = new Float64Array(length);
  for (let i = 0; i < length; i++) {
    leftValues[i] = (((i % 101) * 17 + (seed % 101) * 13) % 101 - 50) / 16;
    rightValues[i] = (((i % 101) * 17 + ((seed + 1) % 101) * 13) % 101 - 50) / 16;
  }
  if (op === "rotation3d") leftValues.set([1, 2, 3]);
  if (op === "cross") rightValues[2] = -rightValues[2];
  if (op === "determinant_cholesky" || op === "determinant_spd_lu") {
    for (let row = 0; row < size; row++) for (let col = 0; col < size; col++) {
      const forward = ((((row * size + col) % 101) * 17 + (seed % 101) * 13) % 101) - 50;
      const reverse = ((((col * size + row) % 101) * 17 + (seed % 101) * 13) % 101) - 50;
      leftValues[row * size + col] = (forward + reverse) / 32;
    }
  }
  if (determinants.has(op)) for (let i = 0; i < size; i++) leftValues[i * size + i] += size * 4;
  if (op === "eigen_symmetric") for (let row = 0; row < size; row++) for (let col = 0; col < size; col++) {
    leftValues[row * size + col] = row === col ? 2 + (seed % 17) / 16 : Math.abs(row - col) === 1 ? -1 : 0;
  }
  if (op === "eigen_general") for (let row = 0; row < size; row++) for (let col = 0; col < size; col++) {
    const block = Math.floor(row / 2);
    leftValues[row * size + col] = row === col ? 1 + (seed % 17) / 16 + block / 8 :
      row % 2 === 0 && col === row + 1 ? -(0.5 + block / 16) :
      row % 2 === 1 && col === row - 1 ? 0.5 + block / 16 :
      row < col ? ((row * 3 + col * 5 + seed) % 11 - 5) / 32 : 0;
  }
  let a, b;
  try {
    a = new Matrix(rows, cols, leftValues);
    b = new Matrix(rows, cols, rightValues);
    const consume = () => {
      const result = execute(op, a, b, rotations.has(op) ? 0.5 : 1.25);
      if (result instanceof Matrix) {
        try { return result.checksum(); } finally { result.dispose(); }
      }
      if (op === "eigen_general") {
        try {
          let sum = 0;
          for (let i = 0; i < result.valuesReal.length; i++) sum += (i + 1) * (result.valuesReal[i] + Math.abs(result.valuesImag[i]));
          const real = result.vectorsReal.toArray(), imag = result.vectorsImag.toArray();
          for (let i = 0; i < real.length; i++) sum += real[i] * real[i] + imag[i] * imag[i];
          return sum;
        } finally { result.vectorsReal.dispose(); result.vectorsImag.dispose(); }
      }
      if (op === "eigen_symmetric") {
        try {
          let sum = 0;
          for (let i = 0; i < result.values.length; i++) sum += (i + 1) * result.values[i];
          for (const value of result.vectors.toArray()) sum += value * value;
          return sum;
        } finally { result.vectors.dispose(); }
      }
      return result;
    };
    for (let i = 0; i < Math.max(5, Math.min(iterations, 100)); i++) consume();
    let checksum = 0;
    const start = process.hrtime.bigint();
    for (let i = 0; i < iterations; i++) checksum += consume();
    const elapsed_ns = Number(process.hrtime.bigint() - start);
    if (!Number.isFinite(checksum)) throw new RangeError("nonfinite accumulated checksum");
    return { elapsed_ns, iterations, checksum };
  } finally {
    a?.dispose(); b?.dispose();
  }
}

async function main() {
  const [command, ...args] = process.argv.slice(2);
  if (command !== "check" && command !== "bench") throw new RangeError("expected check or bench command");
  const options = process.env.LINEAR_A_WASM_MODULE
    ? { moduleUrl: pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE)) }
    : {};
  const { Matrix } = await createMatrixAPI(options);
  const result = command === "check" ? check(args, Matrix) : bench(args, Matrix);
  process.stdout.write(JSON.stringify(result) + "\n");
}

try {
  await main();
} catch (error) {
  process.stderr.write((error instanceof Error ? error.message : String(error)) + "\n");
  process.exitCode = 1;
}
