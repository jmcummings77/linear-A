import test from "node:test";
import assert from "node:assert/strict";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { createMatrixAPI } from "./matrix.mjs";

const options = process.env.LINEAR_A_WASM_MODULE
  ? { moduleUrl: pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE)) } : {};
const { Matrix } = await createMatrixAPI(options);

function invariants(input, size, result) {
  assert.equal(result.values.length, size);
  assert.deepEqual([result.vectors.rows, result.vectors.cols], [size, size]);
  const vectors = result.vectors.toArray();
  const scale = Math.max(...input.map(Math.abs), 0) || 1;
  for (let col = 0; col < size; col++) {
    assert.ok(Number.isFinite(result.values[col]));
    if (col) assert.ok(result.values[col - 1] <= result.values[col]);
    let largest = 0;
    for (let row = 0; row < size; row++) {
      let av = 0, dot = 0;
      for (let k = 0; k < size; k++) {
        av += (input[row * size + k] / scale) * vectors[k * size + col];
        dot += vectors[k * size + row] * vectors[k * size + col];
      }
      assert.ok(Math.abs(av - vectors[row * size + col] * (result.values[col] / scale)) <= 2e-11 * size);
      assert.ok(Math.abs(dot - (row === col ? 1 : 0)) <= 2e-13 * size);
      if (Math.abs(vectors[row * size + col]) > Math.abs(vectors[largest * size + col])) largest = row;
    }
    assert.ok(vectors[largest * size + col] >= 0);
  }
}

function decompose(input, size, validate) {
  const a = new Matrix(size, size, input);
  let result;
  try {
    result = a.eigenSymmetric();
    invariants(input, size, result);
    validate?.(result, a);
    assert.deepEqual([...a.toArray()], input);
  } finally { result?.vectors.dispose(); a.dispose(); }
}

test("symmetric eigenpairs match analytic Toeplitz spectra and preserve inputs", () => {
  for (let n = 0; n <= 12; n++) {
    const input = Array(n * n).fill(0);
    for (let row = 0; row < n; row++) {
      input[row * n + row] = 2.25;
      if (row + 1 < n) input[row * n + row + 1] = input[(row + 1) * n + row] = -1;
    }
    decompose(input, n, result => {
      for (let i = 0; i < n; i++) assert.ok(Math.abs(result.values[i] - (2.25 - 2 * Math.cos((i + 1) * Math.PI / (n + 1)))) < 2e-12);
    });
  }
});

test("eigenpairs retain uniformly scaled inputs and exact mixed-range diagonal spectra", () => {
  for (const scale of [1, 1e300, 1e-300, Number.MIN_VALUE]) {
    decompose([2 * scale, scale, scale, 2 * scale], 2, result => {
      assert.ok(Math.abs(result.values[0] / scale - 1) < 1e-14);
      assert.ok(Math.abs(result.values[1] / scale - 3) < 1e-14);
    });
  }
  decompose([Number.MAX_VALUE, 0, 0, 0, Number.MIN_VALUE, 0, 0, 0, -Number.MAX_VALUE], 3,
    result => assert.deepEqual([...result.values], [-Number.MAX_VALUE, Number.MIN_VALUE, Number.MAX_VALUE]));
  decompose([1,0,0, 0,1,0, 0,0,1], 3, result => assert.deepEqual([...result.values], [1,1,1]));
  decompose(Array(9).fill(0), 3, result => assert.deepEqual([...result.values], [0,0,0]));
});

test("random symmetric matrices return orthonormal eigenvectors with small residuals", () => {
  let state = 31;
  for (let sample = 0; sample < 30; sample++) {
    const input = Array(49).fill(0);
    for (let row = 0; row < 7; row++) for (let col = row; col < 7; col++) {
      state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
      input[row * 7 + col] = input[col * 7 + row] = state % 101 - 50;
    }
    decompose(input, 7);
  }
});

test("eigen outputs have independent ownership across calls, disposal and memory growth", () => {
  const a = new Matrix(2, 2, [2,1,1,2]);
  const first = a.eigenSymmetric(), second = a.eigenSymmetric();
  const large = new Matrix(1, 3_000_000);
  try {
    a.dispose();
    first.values[0] = 77;
    first.vectors.set(0, 0, 99);
    assert.deepEqual([...second.values], [1,3]);
    assert.ok(Math.abs(second.vectors.get(0, 0)) < 1);
    first.vectors.dispose();
    assert.equal(second.vectors.rows, 2);
    assert.throws(() => a.eigenSymmetric(), /disposed/);
  } finally { a.dispose(); first.vectors.dispose(); second.vectors.dispose(); large.dispose(); }
});

test("symmetric eigen options, shape, convergence and range failures are explicit", () => {
  const dense = new Matrix(4, 4, [4,1,2,3, 1,5,1,2, 2,1,6,1, 3,2,1,7]);
  const asymmetric = new Matrix(2, 2, [1,2,2 + 1e-15,1]);
  const rectangle = new Matrix(2, 3);
  const overflow = new Matrix(2, 2, Array(4).fill(Number.MAX_VALUE));
  try {
    assert.throws(() => dense.eigenSymmetric({tolerance: 1e-15, maxSweeps: 1}), /converge/);
    for (const tolerance of [0, -1, 1, Infinity, NaN, "1e-12", true]) {
      assert.throws(() => dense.eigenSymmetric({tolerance}), /tolerance/);
    }
    for (const maxSweeps of [0, -1, 1.5, 10001, Infinity, NaN, "50", true]) {
      assert.throws(() => dense.eigenSymmetric({maxSweeps}), /maxSweeps/);
    }
    assert.throws(() => asymmetric.eigenSymmetric(), /symmetric/);
    assert.throws(() => rectangle.eigenSymmetric(), /dimensions/);
    assert.throws(() => overflow.eigenSymmetric(), /nonfinite/);
    const recovered = dense.eigenSymmetric(); recovered.vectors.dispose();
  } finally { dense.dispose(); asymmetric.dispose(); rectangle.dispose(); overflow.dispose(); }
});
