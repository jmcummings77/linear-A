import assert from "node:assert/strict";
import { spawnSync } from "node:child_process";
import test from "node:test";
import { Matrix } from "./matrix.js";
import "./vector.test.js";

function permutationDeterminant(values: number[], size: number): number {
  // Inputs are dyadic quarters. Integer permutation terms give an exact oracle.
  let total = 0;
  const order: number[] = [], used = Array(size).fill(false);
  function visit(): void {
    if (order.length === size) {
      let term = 1;
      for (let row = 0; row < size; row++) {
        term *= values[row * size + order[row]] * 4;
        for (let next = row + 1; next < size; next++) if (order[row] > order[next]) term = -term;
      }
      total += term;
      return;
    }
    for (let col = 0; col < size; col++) if (!used[col]) {
      used[col] = true; order.push(col); visit(); order.pop(); used[col] = false;
    }
  }
  visit();
  return total / 4 ** size;
}

function near(actual: number, expected: number, relative = 1e-12, absolute = 0): void {
  assert.ok(Number.isFinite(actual) && Math.abs(actual - expected) <= Math.max(absolute, Math.abs(expected) * relative), `${actual} != ${expected}`);
}

function assertEigenpairs(a: Matrix, expected: number[]): void {
  const before = a.values, result = a.eigenSymmetric(), n = a.rows, q = result.vectors;
  assert.deepEqual([q.rows, q.cols], [n, n]);
  assert.deepEqual([...result.values], [...result.values].sort((x, y) => x - y));
  let scale = 1e-300;
  for (const value of before) scale = Math.max(scale, Math.abs(value));
  for (let col = 0; col < n; col++) {
    near(result.values[col] / scale, expected[col] / scale, 1e-9, 1e-9);
    for (let row = 0; row < n; row++) {
      let transformed = 0, inner = 0;
      for (let k = 0; k < n; k++) {
        transformed += a.get(row, k) / scale * q.get(k, col);
        inner += q.get(k, row) * q.get(k, col);
      }
      near(transformed, result.values[col] / scale * q.get(row, col), 1e-9, 1e-9);
      near(inner, row === col ? 1 : 0, 1e-9, 1e-9);
    }
  }
  assert.deepEqual(a.values, before);
  if (n) { result.values[0] = 99; q.set(0, 0, 99); assert.deepEqual(a.values, before); }
}

test("symmetric eigenpairs have analytic spectra, orthogonal columns, and nonmutating storage", () => {
  for (const scale of [1, 1e-200, 1e200]) assertEigenpairs(new Matrix(2, 2, [2 * scale, scale, scale, 2 * scale]), [scale, 3 * scale]);
  assertEigenpairs(new Matrix(2, 2, [0, 1, 1, 0]), [-1, 1]);
  assertEigenpairs(new Matrix(2, 2, [1, 1, 1, 1]), [0, 2]);
  assertEigenpairs(new Matrix(3, 3, [6, -2, -2, -2, 3, 4, -2, 4, 3]), [-1, (13 - Math.sqrt(33)) / 2, (13 + Math.sqrt(33)) / 2]);
  assertEigenpairs(new Matrix(3, 3), [0, 0, 0]);
  assertEigenpairs(Matrix.identity(4), [1, 1, 1, 1]);
  assertEigenpairs(new Matrix(), []);
  for (const n of [3, 8, 16]) {
    const values = Array.from({ length: n * n }, (_, index) => {
      const row = Math.floor(index / n), col = index % n;
      return row === col ? 2 : Math.abs(row - col) === 1 ? -1 : 0;
    });
    assertEigenpairs(new Matrix(n, n, values), Array.from({ length: n }, (_, i) => 2 - 2 * Math.cos((i + 1) * Math.PI / (n + 1))));
  }
});

test("symmetric eigenvalues preserve extreme diagonal magnitudes and reject invalid requests", () => {
  assert.deepEqual([...new Matrix(2, 2, [1e300, 0, 0, 1e-300]).eigenSymmetric().values], [1e-300, 1e300]);
  assert.deepEqual([...new Matrix(1, 1, [Number.MIN_VALUE]).eigenSymmetric().values], [Number.MIN_VALUE]);
  assert.notEqual(new Matrix(2, 2, [1, 1e-200, 1e-200, 2]).eigenSymmetric(1e-250).vectors.get(1, 0), 0);
  for (const a of [new Matrix(2, 3), new Matrix(2, 2, [1, 2, 0, 1]), new Matrix(2, 2, [1, Number.MIN_VALUE, 0, 1])]) {
    assert.throws(() => a.eigenSymmetric(), RangeError);
  }
  const a = new Matrix(3, 3, [4, 1, 2, 1, 3, .5, 2, .5, 5]), before = a.values;
  for (const tolerance of [0, -1, 1, NaN, Infinity]) assert.throws(() => a.eigenSymmetric(tolerance), RangeError);
  for (const sweeps of [0, -1, 1.5]) assert.throws(() => a.eigenSymmetric(1e-12, sweeps), RangeError);
  assert.throws(() => a.eigenSymmetric(1e-12, 1), /converge/);
  assert.deepEqual(a.values, before);
  assert.throws(() => new Matrix(2, 2, [1e308, 1e308, 1e308, 1e308]).eigenSymmetric(), RangeError);
});

test("eigen runner exports both eigenpair arrays and consumes the analytic workload", () => {
  const runner = new URL("runner.js", import.meta.url);
  const invoke = (args: string[], input = "") => spawnSync(process.execPath, [runner.pathname, ...args], { input, encoding: "utf8" });
  let response = invoke(["check", "eigen_symmetric", "2", "2"], "2 1 1 2");
  assert.equal(response.status, 0, response.stderr);
  let result = JSON.parse(response.stdout);
  assert.deepEqual(result.eigenvalues, [1, 3]);
  assert.deepEqual([result.eigenvectors.rows, result.eigenvectors.cols], [2, 2]);
  response = invoke(["bench", "eigen_symmetric", "3", "2", "17"]);
  assert.equal(response.status, 0, response.stderr);
  result = JSON.parse(response.stdout);
  const expected = 2 * ([1, 2, 3].reduce((sum, k) => sum + k * (2 - 2 * Math.cos(k * Math.PI / 4)), 0) + 3);
  near(result.checksum, expected, 1e-9);
  assert.ok(result.elapsed_ns > 0);
});

test("rectangular arithmetic, independent storage, identity", () => {
  const values = [1, 2, 3, 4, 5, 6], a = new Matrix(2, 3, values);
  values[0] = 99;
  assert.deepEqual([...a.transpose().values], [1, 4, 2, 5, 3, 6]);
  assert.deepEqual([...a.multiply(new Matrix(3, 2, [7, 8, 9, 10, 11, 12])).values], [58, 64, 139, 154]);
  assert.deepEqual([...a.scale(0.5).values], [0.5, 1, 1.5, 2, 2.5, 3]);
  assert.deepEqual(a.add(a).subtract(a).values, a.values);
  assert.deepEqual(Matrix.identity(2).multiply(a).values, a.values);
  const copy = a.copy(); copy.set(0, 0, 9);
  const row = a.row(0), col = a.column(1), exported = a.values;
  row[0] = col[0] = exported[0] = 9;
  assert.deepEqual([...a.values], [1, 2, 3, 4, 5, 6]);
});

test("known determinants, row swaps, singularity, fractional trace and triangles", () => {
  const cases: [number, number[], number][] = [
    [1, [-7], -7], [2, [1, 2, 3, 4], -2], [3, [6, 1, 1, 4, -2, 5, 2, 8, 7], -306],
    [4, [3, 2, 0, 1, 4, 0, 1, 2, 3, 0, 2, 1, 9, 2, 3, 1], 24],
    [3, [0, 1, 0, 0, 0, 1, 1, 0, 0], 1], [2, [1, 2, 2, 4], 0]
  ];
  for (const [size, values, expected] of cases) {
    const a = new Matrix(size, size, values);
    assert.ok(Math.abs(a.determinant() - expected) < 1e-10);
    assert.deepEqual([...a.values], values);
  }
  const a = new Matrix(3, 3, [0.5, 7, -1, 0, -1.5, 4, 0, 0, 2.5]);
  assert.deepEqual(a.triangular(), { upper: true, lower: false });
  assert.deepEqual(a.transpose().triangular(), { upper: false, lower: true });
  assert.equal(a.trace(), 1.5);
  assert.equal(a.determinant(), -1.875);
  assert.deepEqual(new Matrix(2, 3).triangular(), { upper: false, lower: false });
});

test("empty dimensions and invalid requests", () => {
  assert.equal(new Matrix().determinant(), 1);
  assert.equal(new Matrix().trace(), 0);
  assert.deepEqual(new Matrix(2, 0).multiply(new Matrix(0, 3)).values, new Float64Array(6));
  assert.deepEqual([...new Matrix(2, 0).row(1)], []);
  const a = new Matrix(1, 1, [1e308]);
  assert.throws(() => a.scale(2), RangeError);
  assert.throws(() => a.set(0, 0, NaN), RangeError);
  assert.throws(() => a.get(-1, 0), RangeError);
  assert.throws(() => new Matrix(-1, 2), RangeError);
  assert.throws(() => new Matrix(2, 2, [1]), RangeError);
  assert.throws(() => a.add(new Matrix(2, 2)), RangeError);
  assert.throws(() => a.multiply(new Matrix(2, 1)), RangeError);
  assert.throws(() => new Matrix(2, 3).trace(), RangeError);
  assert.throws(() => new Matrix(2, 3).determinant(), RangeError);
  assert.equal(a.get(0, 0), 1e308);
});

test("protocol JSON, exact input count, and maximum seed checksum", () => {
  const runner = new URL("runner.js", import.meta.url);
  const invoke = (args: string[], input = "") => spawnSync(process.execPath, [runner.pathname, ...args], { input, encoding: "utf8" });
  let result = invoke(["check", "multiply", "1", "2", "2", "1"], "2 3 4 5");
  assert.equal(result.status, 0, result.stderr);
  assert.deepEqual(JSON.parse(result.stdout), { rows: 1, cols: 1, values: [23] });
  result = invoke(["bench", "scale", "2", "7", "2147483646"]);
  assert.equal(result.status, 0, result.stderr);
  const payload = JSON.parse(result.stdout);
  const values = Array.from({ length: 4 }, (_, i) => ((i * 17 + 2147483646 * 13) % 101 - 50) / 16);
  assert.equal(payload.checksum, (values[0] + values[2] + values[3]) * 1.25 * 7);
  assert.ok(payload.elapsed_ns > 0);
  result = invoke(["check", "trace", "1", "1"], "1 2");
  assert.notEqual(result.status, 0);
  assert.equal(result.stdout, "");
});

test("Auto and LU match exact permutation determinants without mutating inputs", () => {
  let state = 741;
  const next = (): number => { state = (state * 1664525 + 1013904223) >>> 0; return state; };
  for (let size = 0; size <= 6; size++) for (let trial = 0; trial < 5; trial++) {
    const values = Array.from({ length: size * size }, () => (next() % 17 - 8) / 4);
    const matrix = new Matrix(size, size, values), expected = permutationDeterminant(values, size);
    for (const algorithm of ["auto", "lu"] as const) {
      near(matrix.determinant(algorithm), expected, 1e-11, 1e-9);
      assert.deepEqual([...matrix.values], values);
    }
  }
});

test("determinants preserve mixed scales, subnormal answers, and small elimination factors", () => {
  const cases: [number, number[], number][] = [
    [1, [Number.MIN_VALUE], Number.MIN_VALUE], [2, [1e200, 1e200, 1e200, 1e200], 0],
    [2, [1e-200, 1e-200, 3e-124, 6e-124], Number.MIN_VALUE],
    [2, [1e308, 1e308, 1e-308, 2e-308], 1], [2, [1e308, 1e-308, 1e308, 2e-308], 1],
    [3, [1e308, 1e308, 1e-308, 1e-308, 2e-308, 1e308, 0, 0, 1e-308], 1e-308],
  ];
  for (const [size, values, expected] of cases) for (const algorithm of ["auto", "lu"] as const) {
    const matrix = new Matrix(size, size, values), actual = matrix.determinant(algorithm);
    if (expected === 0 || expected === Number.MIN_VALUE) assert.equal(actual, expected);
    else near(actual, expected);
    assert.deepEqual([...matrix.values], values);
  }
  for (const diagonal of [[1e-200, 1e-200, 1e200, 1e200], [1e200, 1e200, 1e-200, 1e-200]]) {
    const values = Array.from({ length: 16 }, (_, index) => Math.floor(index / 4) === index % 4 ? diagonal[index % 4] : 0);
    const matrix = new Matrix(4, 4, values);
    for (const algorithm of ["auto", "lu", "cholesky"] as const) near(matrix.determinant(algorithm), 1);
  }
});

test("Cholesky requires exact symmetry and positive pivots, and preserves input", () => {
  // L=[[2,0,0],[1,3,0],[-1,2,4]] gives det(L*L^T)=24^2.
  const values = [4, 2, -2, 2, 10, 5, -2, 5, 21], matrix = new Matrix(3, 3, values);
  near(matrix.determinant("cholesky"), 576);
  assert.deepEqual([...matrix.values], values);
  assert.equal(new Matrix().determinant("cholesky"), 1);
  assert.equal(new Matrix(1, 1, [Number.MIN_VALUE]).determinant("cholesky"), Number.MIN_VALUE);
  for (const invalid of [[2, 100, 1, 2], [1, 2, 2, 1], [1, 1, 1, 1], [-1, 0, 0, -1]]) {
    const matrix = new Matrix(2, 2, invalid);
    assert.throws(() => matrix.determinant("cholesky"), RangeError);
    assert.deepEqual([...matrix.values], invalid);
  }
  for (const algorithm of ["auto", "lu", "cholesky"] as const) {
    assert.throws(() => new Matrix(2, 3).determinant(algorithm), RangeError);
    assert.throws(() => new Matrix(2, 2, [1e308, 0, 0, 1e308]).determinant(algorithm), RangeError);
    assert.equal(new Matrix(2, 2, [1e-200, 0, 0, 1e-200]).determinant(algorithm), 0);
  }
  assert.throws(() => new Matrix().determinant("cofactor" as never), RangeError);
});

test("algorithm CLI names use identical SPD workloads for Cholesky and LU", () => {
  const runner = new URL("runner.js", import.meta.url);
  const invoke = (args: string[], input = "") => spawnSync(process.execPath, [runner.pathname, ...args], { input, encoding: "utf8" });
  for (const name of ["determinant_lu", "determinant_cholesky", "determinant_spd_lu"]) {
    const result = invoke(["check", name, "2", "2"], "4 2 2 3");
    assert.equal(result.status, 0, result.stderr); near(JSON.parse(result.stdout).value, 8);
  }
  for (const name of ["determinant_cholesky", "determinant_spd_lu"]) {
    const result = invoke(["bench", name, "3", "2", "17"]);
    assert.equal(result.status, 0, result.stderr);
    const data = JSON.parse(result.stdout);
    assert.equal(data.iterations, 2); assert.ok(data.elapsed_ns > 0);
    near(data.checksum, 2 * 28161419 / 16384);
  }
  const general = [161/16, -14/16, 3/16, 20/16, 229/16, -47/16, -30/16, -13/16, 196/16];
  const expected = general[0]*(general[4]*general[8]-general[5]*general[7])
    - general[1]*(general[3]*general[8]-general[5]*general[6]) + general[2]*(general[3]*general[7]-general[4]*general[6]);
  for (const name of ["determinant", "determinant_lu"]) {
    const result = invoke(["bench", name, "3", "2", "17"]);
    assert.equal(result.status, 0, result.stderr); near(JSON.parse(result.stdout).checksum, 2 * expected);
  }
  for (const [args, input] of [[["check", "determinant_cholesky", "2", "2"], "1 2 2 1"],
                             [["bench", "determinant_cofactor", "2", "1", "17"], ""]] as [string[], string][]) {
    const result = invoke(args, input); assert.notEqual(result.status, 0); assert.equal(result.stdout, "");
  }
});
