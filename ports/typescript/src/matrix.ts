import { decompose, type SingularValueDecomposition } from "./svd.js";
export type { SingularValueDecomposition } from "./svd.js";
import { Factorization } from "./solve.js";
export { Factorization } from "./solve.js";
import { generalEigen } from "./general-eigen.js";

/** Dependency-free float64 matrix with owned, contiguous row-major storage. */
export function finite(value: number): number {
  if (!Number.isFinite(value)) throw new RangeError("matrix arithmetic requires finite float64 values");
  return value;
}

export type DeterminantAlgorithm = "auto" | "lu" | "cholesky";

export interface GeneralEigenDecomposition {
  /** Lexicographic (real, imaginary) eigenvalues and matching unit right columns. */
  valuesReal: Float64Array;
  valuesImag: Float64Array;
  vectorsReal: Matrix;
  vectorsImag: Matrix;
}

export interface SymmetricEigenDecomposition {
  /** Ascending eigenvalues; column j of vectors corresponds to values[j]. */
  values: Float64Array;
  vectors: Matrix;
}

const exponentBits = new DataView(new ArrayBuffer(8));
function frexp(value: number): [number, number] {
  if (value === 0) return [value, 0];
  let correction = 0;
  if (Math.abs(value) < 2 ** -1022) { value *= 2 ** 54; correction = -54; }
  exponentBits.setFloat64(0, value);
  const high = exponentBits.getUint32(0), exponent = ((high >>> 20) & 0x7ff) - 1022 + correction;
  exponentBits.setUint32(0, (high & 0x800fffff) | (1022 << 20));
  return [exponentBits.getFloat64(0), exponent];
}

function compose(mantissa: number, exponent: number): number {
  if (mantissa === 0) return mantissa;
  if (exponent > 1024) return mantissa > 0 ? Infinity : -Infinity;
  if (exponent < -1074) return mantissa * 0;
  if (exponent === 1024) return (mantissa * 2) * 2 ** 1023;
  // Perform the final subnormal rounding once, after a normal-range multiply.
  if (exponent < -1022) return (mantissa * 2 ** (exponent + 1022)) * 2 ** -1022;
  return mantissa * 2 ** exponent;
}

function shiftExponent(value: number, shift: number): number {
  const [mantissa, exponent] = frexp(value);
  return compose(mantissa, exponent + shift);
}

function quotientProduct(numerator: number, denominator: number, value: number): number {
  const [left, leftExponent] = frexp(numerator), [right, rightExponent] = frexp(denominator);
  const [item, itemExponent] = frexp(value), [mantissa, carry] = frexp((left / right) * item);
  return compose(mantissa, leftExponent - rightExponent + itemExponent + carry);
}

function diagonalProduct(values: Float64Array, size: number, exponent = 0, sign = 1, power = 1): number {
  let mantissa = sign;
  for (let index = 0; index < size; index++) {
    const value = values[index * size + index];
    if (value === 0) return 0;
    const [part, shift] = frexp(value);
    for (let repeat = 0; repeat < power; repeat++) {
      const [normalized, carry] = frexp(mantissa * part);
      mantissa = normalized;
      exponent += shift + carry;
    }
  }
  return finite(compose(mantissa, exponent));
}

export class Matrix {
  readonly rows: number;
  readonly cols: number;
  private readonly data: Float64Array;

  constructor(rows = 0, cols = 0, values?: Iterable<number>) {
    if (!Number.isSafeInteger(rows) || !Number.isSafeInteger(cols) || rows < 0 || cols < 0 || !Number.isSafeInteger(rows * cols)) {
      throw new RangeError("dimensions must be nonnegative safe integers");
    }
    this.rows = rows;
    this.cols = cols;
    this.data = values === undefined ? new Float64Array(rows * cols) : Float64Array.from(values);
    if (this.data.length !== rows * cols) throw new RangeError("value count does not match dimensions");
    for (const value of this.data) finite(value);
  }

  static identity(size: number): Matrix {
    const result = new Matrix(size, size);
    for (let i = 0; i < size; i++) result.data[i * size + i] = 1;
    return result;
  }

  get values(): Float64Array { return this.data.slice(); }
  copy(): Matrix { return new Matrix(this.rows, this.cols, this.data); }

  private index(row: number, col: number): number {
    if (!Number.isInteger(row) || !Number.isInteger(col) || row < 0 || row >= this.rows || col < 0 || col >= this.cols) {
      throw new RangeError("matrix index out of range");
    }
    return row * this.cols + col;
  }

  get(row: number, col: number): number { return this.data[this.index(row, col)]; }
  set(row: number, col: number, value: number): void { this.data[this.index(row, col)] = finite(value); }

  row(index: number): Float64Array {
    if (!Number.isInteger(index) || index < 0 || index >= this.rows) throw new RangeError("row index out of range");
    return this.data.slice(index * this.cols, (index + 1) * this.cols);
  }

  column(index: number): Float64Array {
    if (!Number.isInteger(index) || index < 0 || index >= this.cols) throw new RangeError("column index out of range");
    const result = new Float64Array(this.rows);
    for (let row = 0; row < this.rows; row++) result[row] = this.data[row * this.cols + index];
    return result;
  }

  private combine(other: Matrix, subtract: boolean): Matrix {
    if (!(other instanceof Matrix) || this.rows !== other.rows || this.cols !== other.cols) throw new RangeError("matrix dimensions must match");
    const result = new Matrix(this.rows, this.cols);
    for (let i = 0; i < this.data.length; i++) result.data[i] = finite(subtract ? this.data[i] - other.data[i] : this.data[i] + other.data[i]);
    return result;
  }

  add(other: Matrix): Matrix { return this.combine(other, false); }
  subtract(other: Matrix): Matrix { return this.combine(other, true); }

  scale(scalar: number): Matrix {
    finite(scalar);
    const result = new Matrix(this.rows, this.cols);
    for (let i = 0; i < this.data.length; i++) result.data[i] = finite(this.data[i] * scalar);
    return result;
  }

  transpose(): Matrix {
    const result = new Matrix(this.cols, this.rows);
    for (let row = 0; row < this.rows; row++)
      for (let col = 0; col < this.cols; col++) result.data[col * this.rows + row] = this.data[row * this.cols + col];
    return result;
  }

  multiply(other: Matrix): Matrix {
    if (!(other instanceof Matrix) || this.cols !== other.rows) throw new RangeError("left columns must equal right rows");
    const result = new Matrix(this.rows, other.cols);
    for (let row = 0; row < this.rows; row++)
      for (let col = 0; col < other.cols; col++) {
        let total = 0;
        for (let k = 0; k < this.cols; k++) total += this.data[row * this.cols + k] * other.data[k * other.cols + col];
        result.data[row * other.cols + col] = finite(total);
      }
    return result;
  }

  private vector3(): Float64Array {
    if (!((this.rows === 3 && this.cols === 1) || (this.rows === 1 && this.cols === 3))) {
      throw new RangeError("expected a three-component row or column vector");
    }
    for (const value of this.data) finite(value);
    return this.data;
  }

  /** Right-handed three-component vector cross product; retains the left shape. */
  cross(other: Matrix): Matrix {
    if (!(other instanceof Matrix)) throw new RangeError("cross product requires two three-component vectors");
    const left = this.vector3(), right = other.vector3();
    return new Matrix(this.rows, this.cols, [left[1] * right[2] - left[2] * right[1],
      left[2] * right[0] - left[0] * right[2], left[0] * right[1] - left[1] * right[0]]);
  }

  /** Active counterclockwise rotation of column vectors; angle in radians. */
  static rotation2D(radians: number): Matrix {
    finite(radians);
    const cosine = Math.cos(radians), sine = Math.sin(radians);
    return new Matrix(2, 2, [cosine, -sine, sine, cosine]);
  }

  /** Active right-handed rotation about X, acting on column vectors. */
  static rotationX(radians: number): Matrix {
    finite(radians);
    const cosine = Math.cos(radians), sine = Math.sin(radians);
    return new Matrix(3, 3, [1, 0, 0, 0, cosine, -sine, 0, sine, cosine]);
  }

  /** Active right-handed rotation about Y, acting on column vectors. */
  static rotationY(radians: number): Matrix {
    finite(radians);
    const cosine = Math.cos(radians), sine = Math.sin(radians);
    return new Matrix(3, 3, [cosine, 0, sine, 0, 1, 0, -sine, 0, cosine]);
  }

  /** Active right-handed rotation about Z, acting on column vectors. */
  static rotationZ(radians: number): Matrix {
    finite(radians);
    const cosine = Math.cos(radians), sine = Math.sin(radians);
    return new Matrix(3, 3, [cosine, -sine, 0, sine, cosine, 0, 0, 0, 1]);
  }

  /** Rodrigues rotation about a finite nonzero row or column vector.
   * The axis is normalized without mutation. Positive angles in radians follow
   * the right-hand rule; the resulting matrix acts on column vectors.
   */
  static rotationAxisAngle(axis: Matrix, radians: number): Matrix {
    finite(radians);
    if (!(axis instanceof Matrix)) throw new RangeError("rotation axis must be a three-component vector");
    const values = axis.vector3();
    const largest = Math.max(...values.map(Math.abs));
    if (largest === 0) throw new RangeError("rotation axis must be nonzero");
    const scaled = values.map(value => value / largest);
    const length = Math.sqrt(scaled[0] ** 2 + scaled[1] ** 2 + scaled[2] ** 2);
    const [x, y, z] = scaled.map(value => value / length);
    const cosine = Math.cos(radians), sine = Math.sin(radians);
    const complement = Math.abs(radians) < 1 ? 2 * Math.sin(radians / 2) ** 2 : 1 - cosine;
    return new Matrix(3, 3, [cosine + x * x * complement, x * y * complement - z * sine, x * z * complement + y * sine,
      x * y * complement + z * sine, cosine + y * y * complement, y * z * complement - x * sine,
      x * z * complement - y * sine, y * z * complement + x * sine, cosine + z * z * complement]);
  }

  private square(): void {
    if (this.rows !== this.cols) throw new RangeError("operation requires a square matrix");
  }

  trace(): number {
    this.square();
    let result = 0;
    for (let i = 0; i < this.rows; i++) result += this.data[i * this.cols + i];
    return finite(result);
  }

  /** Auto uses safe tiny/triangular shortcuts; explicit Cholesky requires SPD input.
   * LU uses double precision and exact power-of-two scaling where entries survive.
   * Final overflow throws; underflow may return zero. This is not a rank test.
   */
  determinant(algorithm: DeterminantAlgorithm = "auto"): number {
    this.square();
    if (!["auto", "lu", "cholesky"].includes(algorithm)) throw new RangeError("determinant algorithm must be auto, lu, or cholesky");
    const n = this.rows;
    if (algorithm === "cholesky") return this.choleskyDeterminant();
    if (algorithm === "auto") {
      if (n === 0) return 1;
      if (n === 1) return this.data[0];
      // Six products of /16 dyadics in this range fit exactly in 53 bits.
      if (n <= 3 && this.data.every(value => Number.isInteger(value * 16) && Math.abs(value * 16) <= 65536)) {
        const a = this.data;
        if (n === 2) return a[0] * a[3] - a[1] * a[2];
        return a[0] * a[4] * a[8] + a[1] * a[5] * a[6] + a[2] * a[3] * a[7]
          - a[2] * a[4] * a[6] - a[1] * a[3] * a[8] - a[0] * a[5] * a[7];
      }
      const { upper, lower } = this.triangular();
      if (upper || lower) return diagonalProduct(this.data, n);
    }
    return this.luDeterminant();
  }

  private luDeterminant(): number {
    const n = this.rows;
    const work = this.data.slice();
    let exponent = 0;
    for (let row = 0; row < n; row++) {
      let largest = 0;
      for (let col = 0; col < n; col++) largest = Math.max(largest, Math.abs(work[row * n + col]));
      if (largest === 0) return 0;
      const shift = frexp(largest)[1];
      const scaled = new Float64Array(n);
      let preserves = true;
      for (let col = 0; col < n; col++) {
        scaled[col] = shiftExponent(work[row * n + col], -shift);
        if (shiftExponent(scaled[col], shift) !== work[row * n + col]) preserves = false;
      }
      if (preserves) { work.set(scaled, row * n); exponent += shift; }
    }
    let sign = 1;
    for (let k = 0; k < n; k++) {
      let pivotRow = k;
      for (let row = k + 1; row < n; row++) if (Math.abs(work[row * n + k]) > Math.abs(work[pivotRow * n + k])) pivotRow = row;
      const pivot = work[pivotRow * n + k];
      if (pivot === 0) return 0;
      if (pivotRow !== k) {
        for (let col = k; col < n; col++) {
          const temporary = work[k * n + col];
          work[k * n + col] = work[pivotRow * n + col];
          work[pivotRow * n + col] = temporary;
        }
        sign = -sign;
      }
      for (let row = k + 1; row < n; row++) {
        const numerator = work[row * n + k], factor = numerator / pivot;
        const tinyFactor = numerator !== 0 && Math.abs(factor) < 2 ** -1022;
        work[row * n + k] = 0;
        for (let col = k + 1; col < n; col++) {
          const entry = work[k * n + col];
          const product = tinyFactor ? quotientProduct(numerator, pivot, entry) : factor * entry;
          work[row * n + col] = finite(work[row * n + col] - product);
        }
      }
    }
    return diagonalProduct(work, n, exponent, sign);
  }

  private choleskyDeterminant(): number {
    const n = this.rows;
    for (let row = 0; row < n; row++) for (let col = 0; col < row; col++) {
      if (this.data[row * n + col] !== this.data[col * n + row]) {
        throw new RangeError("Cholesky requires an exactly symmetric positive-definite matrix");
      }
    }
    const lower = new Float64Array(n * n);
    for (let row = 0; row < n; row++) for (let col = 0; col <= row; col++) {
      let total = this.data[row * n + col];
      for (let k = 0; k < col; k++) total -= lower[row * n + k] * lower[col * n + k];
      finite(total);
      if (row === col) {
        if (total <= 0) throw new RangeError("Cholesky requires positive computed pivots");
        lower[row * n + col] = Math.sqrt(total);
      } else lower[row * n + col] = finite(total / lower[col * n + col]);
    }
    return diagonalProduct(lower, n, 0, 1, 2);
  }

  triangular(): { upper: boolean; lower: boolean } {
    if (this.rows !== this.cols) return { upper: false, lower: false };
    let upper = true, lower = true;
    for (let row = 0; row < this.rows; row++)
      for (let col = 0; col < this.cols; col++) if (this.data[row * this.cols + col] !== 0) {
        if (row > col) upper = false;
        else if (row < col) lower = false;
      }
    return { upper, lower };
  }

  /** Eigenpairs of any finite real square matrix, including complex pairs.
   * The result owns its storage. Defective matrices can have dependent columns;
   * general eigenvectors need not be orthogonal. maxIterations bounds QR steps
   * between deflations. Nonconvergence and range failures throw RangeError.
   */
  eigenGeneral(maxIterations = 1000): GeneralEigenDecomposition {
    if (this.rows !== this.cols) throw new RangeError("eigendecomposition requires a square matrix");
    const result = generalEigen(this.data, this.rows, maxIterations);
    return {valuesReal:result.valuesReal, valuesImag:result.valuesImag,
      vectorsReal:new Matrix(this.rows,this.cols,result.vectorsReal),
      vectorsImag:new Matrix(this.rows,this.cols,result.vectorsImag)};
  }

  /** Cyclic Jacobi decomposition of an exactly symmetric matrix, preserving the input.
   * The tolerance bounds the off-diagonal Frobenius norm relative to the input norm.
   * Eigenvectors are orthonormal columns matching ascending eigenvalues. Accuracy
   * is normwise; tiny eigenvalues in mixed-scale inputs can have large relative error.
   * Nonconvergence or nonfinite eigenvalues throw RangeError.
   */
  eigenSymmetric(tolerance = 1e-12, maxSweeps = 50): SymmetricEigenDecomposition {
    this.square();
    if (!Number.isFinite(tolerance) || tolerance <= 0 || tolerance >= 1) throw new RangeError("eigenvalue tolerance must be finite and between zero and one");
    if (!Number.isSafeInteger(maxSweeps) || maxSweeps < 1) throw new RangeError("maximum eigenvalue sweeps must be a positive integer");
    const n = this.rows, q = Matrix.identity(n);
    let diagonal = true;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
      finite(this.data[i * n + j]);
      if (this.data[i * n + j] !== this.data[j * n + i]) throw new RangeError("eigendecomposition requires an exactly symmetric matrix");
      if (i !== j && this.data[i * n + j] !== 0) diagonal = false;
    }
    const values = new Float64Array(n);
    if (diagonal) {
      // Preserve mixed-magnitude diagonal spectra without global scaling.
      for (let i = 0; i < n; i++) values[i] = this.data[i * n + i];
    } else {
      let scale = 0;
      for (const value of this.data) scale = Math.max(scale, Math.abs(value));
      const work = this.data.map(value => value / scale);
      let normSquared = 0;
      for (const value of work) normSquared += value * value;
      const target = tolerance * Math.sqrt(normSquared), threshold = target / (2 * n);
      const converged = (): boolean => {
        let offNorm = 0;
        for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) offNorm = Math.hypot(offNorm, work[i * n + j], work[i * n + j]);
        return offNorm <= target;
      };
      for (let sweep = 0; sweep < maxSweeps && !converged(); sweep++) {
        for (let p = 0; p < n - 1; p++) for (let r = p + 1; r < n; r++) {
          const off = work[p * n + r];
          if (Math.abs(off) <= threshold) continue;
          const delta = (work[r * n + r] - work[p * n + p]) / 2;
          const tangent = off / (delta + (delta < 0 ? -1 : 1) * Math.hypot(delta, off));
          const cosine = 1 / Math.sqrt(1 + tangent * tangent), sine = tangent * cosine;
          work[p * n + p] -= tangent * off;
          work[r * n + r] += tangent * off;
          work[p * n + r] = work[r * n + p] = 0;
          for (let k = 0; k < n; k++) {
            if (k !== p && k !== r) {
              const left = work[k * n + p], right = work[k * n + r];
              work[k * n + p] = work[p * n + k] = cosine * left - sine * right;
              work[k * n + r] = work[r * n + k] = sine * left + cosine * right;
            }
            const left = q.data[k * n + p], right = q.data[k * n + r];
            q.data[k * n + p] = cosine * left - sine * right;
            q.data[k * n + r] = sine * left + cosine * right;
          }
        }
      }
      if (!converged()) throw new RangeError("symmetric eigendecomposition did not converge within maximum sweeps");
      for (let i = 0; i < n; i++) values[i] = finite(work[i * n + i] * scale);
    }
    const order = Array.from({ length: n }, (_, i) => i).sort((left, right) => values[left] - values[right]);
    const sorted = new Float64Array(n), vectors = new Matrix(n, n);
    for (let col = 0; col < n; col++) {
      const original = order[col];
      let lengthSquared = 0, pivot = 0;
      for (let row = 0; row < n; row++) {
        const value = q.data[row * n + original];
        lengthSquared += value * value;
        if (Math.abs(value) > Math.abs(q.data[pivot * n + original])) pivot = row;
      }
      const multiplier = (q.data[pivot * n + original] < 0 ? -1 : 1) / Math.sqrt(lengthSquared);
      for (let row = 0; row < n; row++) vectors.data[row * n + col] = q.data[row * n + original] * multiplier;
      sorted[col] = values[original];
    }
    return { values: sorted, vectors };
  }

  factorLU():Factorization { return new Factorization(this,"lu"); }
  factorCholesky():Factorization { return new Factorization(this,"cholesky"); }
  factorQR():Factorization { return new Factorization(this,"qr"); }
  svd(tolerance=1e-12,maxSweeps=100): SingularValueDecomposition { return decompose(this,tolerance,maxSweeps); }
  solve(rhs:Matrix):Matrix { return this.factorLU().solve(rhs); }
  leastSquares(rhs:Matrix):Matrix { return this.factorQR().solve(rhs); }

  checksum(): number {
    const n = this.data.length;
    return n ? finite(this.data[0] + this.data[Math.floor(n / 2)] + this.data[n - 1]) : 0;
  }
}
