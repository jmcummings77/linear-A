import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";
import { Matrix, finite, type GeneralEigenDecomposition, type SymmetricEigenDecomposition } from "./matrix.js";

const determinants = new Set(["determinant", "determinant_lu", "determinant_cholesky", "determinant_spd_lu"]);
const operations = new Set(["pseudoinverse", "spectral_diagnostics", "solve_minimum_norm", "svd", "svd_one_sweep","add", "subtract", "scale", "transpose", "multiply", "trace", "triangular", "eigen_symmetric", "eigen_general",
  "solve", "solve_cholesky", "least_squares", "rcond", "cross", "rotation2d", "rotation3d", ...determinants]);

function integer(text: string, minimum = 0, maximum = Number.MAX_SAFE_INTEGER): number {
  if (!/^[+-]?\d+$/.test(text)) throw new RangeError("expected an integer argument");
  const value = Number(text);
  if (!Number.isSafeInteger(value) || value < minimum || value > maximum) throw new RangeError("integer argument is outside the allowed range");
  return value;
}

function number(text: string): number {
  if (!/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/.test(text)) throw new RangeError("expected a finite decimal number");
  return finite(Number(text));
}

export function operation(name: string, a: Matrix, b?: Matrix, scalar = 1.25): Matrix | number | GeneralEigenDecomposition | SymmetricEigenDecomposition | { upper: boolean; lower: boolean } {
  switch (name) {
    case "add": return a.add(b!);
    case "subtract": return a.subtract(b!);
    case "multiply": return a.multiply(b!);
    case "pseudoinverse": return a.pseudoinverse(scalar);
    case "solve_minimum_norm": return a.solveMinimumNorm(b!);
    case "spectral_diagnostics": {
      const d=a.spectralDiagnostics(scalar);
      return new Matrix(1,3,[d.rank,d.reciprocalCondition,d.retainedReciprocalCondition]);
    }
    case "svd": case "svd_one_sweep": {
      const r=a.svd(1e-12,name==="svd"?100:1);
      const k=Math.min(a.rows,a.cols),values=new Float64Array((a.rows+1+a.cols)*k);
      {
        values.set(r.u.values); values.set(r.values,a.rows*k);
        const vt=r.vt.values;
        for(let i=0;i<a.cols;i++)for(let j=0;j<k;j++)values[(a.rows+1+i)*k+j]=vt[j*a.cols+i];
        return new Matrix(a.rows+1+a.cols,k,values);
      }
    }
    case "solve": return a.solve(b!);
    case "least_squares": return a.leastSquares(b!);
    case "solve_cholesky": { return a.factorCholesky().solve(b!); }
    case "rcond": { return a.factorLU().reciprocalCondition(); }
    case "cross": return a.cross(b!);
    case "rotation2d":
      if (a.rows !== 0 || a.cols !== 0) throw new RangeError("2D rotation construction expects a 0 by 0 input");
      return Matrix.rotation2D(scalar);
    case "rotation3d": return Matrix.rotationAxisAngle(a, scalar);
    case "scale": return a.scale(scalar);
    case "transpose": return a.transpose();
    case "trace": return a.trace();
    case "determinant": return a.determinant();
    case "determinant_lu":
    case "determinant_spd_lu": return a.determinant("lu");
    case "determinant_cholesky": return a.determinant("cholesky");
    case "triangular": return a.triangular();
    case "eigen_symmetric": return a.eigenSymmetric();
    case "eigen_general": return a.eigenGeneral();
    default: throw new RangeError("unsupported operation");
  }
}

function check(args: string[]): object {
  if (args.length < 3 || !operations.has(args[0])) throw new RangeError("usage: check OP ROWS COLS [BROWS BCOLS | SCALAR]");
  const [name] = args;
  const rows = integer(args[1]), cols = integer(args[2]);
  const binary = ["add", "subtract", "multiply", "cross", "solve", "solve_cholesky", "least_squares", "solve_minimum_norm"].includes(name);
  const scalarOperation = ["scale", "rotation2d", "rotation3d", "pseudoinverse", "spectral_diagnostics"].includes(name);
  if (args.length !== (binary ? 5 : scalarOperation ? 4 : 3)) throw new RangeError("incorrect number of operation arguments");
  const brows = binary ? integer(args[3]) : 0, bcols = binary ? integer(args[4]) : 0;
  const scalar = scalarOperation ? number(args[3]) : 1.25;
  const input = readFileSync(0, "utf8").trim();
  const values = input ? input.split(/\s+/).map(number) : [];
  const count = rows * cols;
  if (values.length !== count + brows * bcols) throw new RangeError("input value count does not match matrix dimensions");
  const a = new Matrix(rows, cols, values.slice(0, count));
  const b = binary ? new Matrix(brows, bcols, values.slice(count)) : undefined;
  const result = operation(name, a, b, scalar);
  if (result instanceof Matrix) return { rows: result.rows, cols: result.cols, values: Array.from(result.values) };
  if (typeof result === "object" && "vectorsReal" in result) return {
    eigenvalues_real:Array.from(result.valuesReal),eigenvalues_imag:Array.from(result.valuesImag),
    eigenvectors_real:{rows:result.vectorsReal.rows,cols:result.vectorsReal.cols,values:Array.from(result.vectorsReal.values)},
    eigenvectors_imag:{rows:result.vectorsImag.rows,cols:result.vectorsImag.cols,values:Array.from(result.vectorsImag.values)}
  };
  if (typeof result === "object" && "vectors" in result) return {
    eigenvalues: Array.from(result.values),
    eigenvectors: { rows: result.vectors.rows, cols: result.vectors.cols, values: Array.from(result.vectors.values) }
  };
  return typeof result === "number" ? { value: finite(result) } : result;
}

export function benchmark(args: string[]): object {
  if (args.length !== 4 || !operations.has(args[0]) || args[0] === "triangular") throw new RangeError("usage: bench OP SIZE ITERATIONS SEED");
  const [name] = args;
  const size = integer(args[1], 1), iterations = integer(args[2], 1), seed = integer(args[3], 0, 2147483646);
  const requiredSize = ({ cross: 3, rotation2d: 2, rotation3d: 3 } as Record<string, number>)[name];
  if (requiredSize !== undefined && size !== requiredSize) throw new RangeError("invalid size for vector or rotation benchmark");
  let a = new Matrix(size, size), b = new Matrix(size, size);
  const spd = name === "determinant_cholesky" || name === "determinant_spd_lu";
  const numerator = (index: number, inputSeed: number): number => (index * 17 + inputSeed * 13) % 101 - 50;
  for (let row = 0; row < size; row++) for (let col = 0; col < size; col++) {
    const index = row * size + col;
    a.set(row, col, name === "eigen_symmetric" ? (row === col ? 2 + (seed % 17) / 16 : Math.abs(row - col) === 1 ? -1 : 0)
      : spd ? (numerator(index, seed) + numerator(col * size + row, seed)) / 32 : numerator(index, seed) / 16);
    b.set(row, col, numerator(index, seed + 1) / 16);
  }
  if (determinants.has(name)) for (let i = 0; i < size; i++) a.set(i, i, a.get(i, i) + size * 4);
  if (name === "eigen_general") for(let i=0;i<size;i++) for(let j=0;j<size;j++) {
    const k=Math.floor(i/2), re=1+(seed%17)/16+k/8, im=.5+k/16;
    a.set(i,j,i===j?re:i%2===0&&j===i+1?-im:i%2===1&&j===i-1?im:i<j?((i*3+j*5+seed)%11-5)/32:0);
  }
  const scalar = name === "rotation2d" || name === "rotation3d" ? .5 : 1.25;
  if (name === "cross") {
    a = new Matrix(3, 1, a.values.slice(0, 3)); b = new Matrix(3, 1, b.values.slice(0, 3));
    b.set(2, 0, -b.get(2, 0));
  }
  else if (name === "rotation2d") a = new Matrix();
  else if (name === "rotation3d") a = new Matrix(3, 1, [1, 2, 3]);
  for (let i = 0; i < Math.max(5, Math.min(iterations, 100)); i++) operation(name, a, b, scalar);
  let checksum = 0;
  const start = process.hrtime.bigint();
  for (let i = 0; i < iterations; i++) {
    const result = operation(name, a, b, scalar);
    if (result instanceof Matrix) checksum += result.checksum();
    else if (typeof result === "object" && "vectorsReal" in result) {
      for(let j=0;j<result.valuesReal.length;j++) checksum+=(j+1)*(result.valuesReal[j]+Math.abs(result.valuesImag[j]));
      for(const x of result.vectorsReal.values) checksum+=x*x;
      for(const x of result.vectorsImag.values) checksum+=x*x;
    } else if (typeof result === "object" && "vectors" in result) {
      for (let j = 0; j < result.values.length; j++) checksum += (j + 1) * result.values[j];
      for (const value of result.vectors.values) checksum += value * value;
    } else checksum += result as number;
  }
  const elapsed = Number(process.hrtime.bigint() - start);
  return { elapsed_ns: elapsed, iterations, checksum: finite(checksum) };
}

export function main(args = process.argv.slice(2)): void {
  if (args[0] !== "check" && args[0] !== "bench") throw new RangeError("usage: runner check|bench ...");
  console.log(JSON.stringify(args[0] === "check" ? check(args.slice(1)) : benchmark(args.slice(1))));
}

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  try { main(); }
  catch (error) {
    console.error(error instanceof Error ? error.message : String(error));
    process.exitCode = 1;
  }
}
