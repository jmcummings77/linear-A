import assert from "node:assert/strict";
import test from "node:test";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { createMatrixAPI } from "./matrix.mjs";
import { assertEigenResult, runBenchmark, runChecks } from "../../benchmarks/live-worker.mjs";
import { computeGeometry } from "../../benchmarks/geometry-worker.mjs";

const options = process.env.LINEAR_A_WASM_MODULE
  ? { moduleUrl: pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE)) } : {};
const { Matrix } = await createMatrixAPI(options);
const symmetric = [2, 1, 0, 1, 2, 0, 0, 0, .5];
const identity = [1, 0, 0, 0, 1, 0, 0, 0, 1];
const asymmetric = [1, 1, 0, 0, 1, 0, 0, 0, 1];
const input = { rows: 3, cols: 3, values: symmetric };
const eigenFixture = { name: "symmetric analytic eigenpairs", op: "eigen_symmetric", a: input,
  b: null, invalid: false, expected_eigenvalues: [.5, 1, 3] };
const geometryConfig = { operation: "multiply", m: symmetric, n: identity, scalar: 1.25 };

/** Record real WASM objects, including outputs constructed by the base API. */
function trackedAPI(changeEigen) {
  const owned = new Set();
  let eigenCalls = 0;
  function track(matrix) {
    if (owned.has(matrix)) return matrix;
    owned.add(matrix);
    const multiply = matrix.multiply;
    matrix.multiply = (...args) => track(multiply.apply(matrix, args));
    return matrix;
  }
  class TrackedMatrix extends Matrix {
    // WASM operations return the base Matrix class, whose handles are equally valid.
    static [Symbol.hasInstance](value) { return value instanceof Matrix; }
    constructor(...args) { super(...args); track(this); }
    eigenGeneral(...args) {
      const result=super.eigenGeneral(...args);
      track(result.vectorsReal);track(result.vectorsImag);
      return result;
    }
    eigenSymmetric(...args) {
      eigenCalls++;
      const result = super.eigenSymmetric(...args);
      track(result.vectors);
      changeEigen?.(result, track);
      return result;
    }
  }
  return { TrackedMatrix, owned, get eigenCalls() { return eigenCalls; } };
}

function disposed(owned) {
  assert.ok(owned.size > 0, "the operation should acquire WASM objects");
  assert.ok([...owned].every(matrix => matrix.disposed), "every acquired matrix must be disposed");
}

function replaceWithIdentity(result) {
  const before = Array.from(result.values);
  const squaredNorm = result.vectors.toArray().reduce((sum, value) => sum + value * value, 0);
  const n = result.values.length;
  for (let row = 0; row < n; row++) for (let col = 0; col < n; col++) result.vectors.set(row, col, row === col ? 1 : 0);
  assert.deepEqual(Array.from(result.values), before);
  const changedNorm = result.vectors.toArray().reduce((sum, value) => sum + value * value, 0);
  assert.ok(Math.abs(changedNorm - squaredNorm) < 1e-11,
    "corruption should preserve both the eigenvalue checksum and vector squared norm");
}

test("live eigen checks verify full eigenpairs and reject a nonsymmetric input with cleanup", async () => {
  const api = trackedAPI(), progress = [];
  const fixtures = [eigenFixture, { ...eigenFixture, name: "nonsymmetric eigenpairs", invalid: true,
    a: { rows: 3, cols: 3, values: asymmetric } }];
  const result = await runChecks(api.TrackedMatrix, fixtures, event => progress.push(event));
  assert.equal(result.passed, 2, JSON.stringify(result.checks));
  assert.equal(result.total, 2);
  assert.equal(progress.length, 2);
  assert.equal(api.eigenCalls, 2);
  assert.equal(api.owned.size, 3);
  disposed(api.owned);
});

test("live checks reject a wrong orthonormal basis with unchanged eigenvalues and norm checksum", async () => {
  const api = trackedAPI(replaceWithIdentity);
  const result = await runChecks(api.TrackedMatrix, [eigenFixture]);
  assert.equal(result.passed, 0);
  assert.match(result.checks[0].error, /eigenpair residual/);
  disposed(api.owned);
});

test("live eigen benchmarks validate the full timed input and consume both result arrays", async () => {
  const api = trackedAPI();
  const result = await runBenchmark(api.TrackedMatrix, { operation: "eigen_symmetric", size: 8, samples: 3, seed: 17 });
  assert.equal(result.operation, "eigen_symmetric");
  assert.equal(result.checksum_verified, true);
  assert.equal(result.samples.length, 3);
  assert.ok(result.iterations > 0 && result.median_ns > 0);
  const expected = Array.from({ length: 8 }, (_, index) =>
    (index + 1) * (2 - 2 * Math.cos((index + 1) * Math.PI / 9))).reduce((sum, value) => sum + value, 8);
  for (const sample of result.samples) {
    assert.ok(Math.abs(sample.checksum / result.iterations - expected) < 1e-9 * expected);
    assert.ok(sample.elapsed_ms > 0);
  }
  assert.ok(api.eigenCalls > result.iterations * 3, "validation, calibration, and warmup are separate from the samples");
  disposed(api.owned);
});

test("live benchmark rejects eigenvector corruption during validation before timing", async () => {
  const api = trackedAPI(replaceWithIdentity);
  await assert.rejects(runBenchmark(api.TrackedMatrix,
    { operation: "eigen_symmetric", size: 8, samples: 3, seed: 17 }), /eigenpair residual/);
  assert.equal(api.eigenCalls, 1);
  assert.equal(api.owned.size, 3);
  disposed(api.owned);
});

test("symmetric geometry exports the analytic eigenbasis and frees every WASM object", async () => {
  const api = trackedAPI();
  const geometry = await computeGeometry(api.TrackedMatrix, geometryConfig);
  assert.deepEqual(geometry.eigen_m.values, [.5, 1, 3]);
  assertEigenResult({ eigenvalues: geometry.eigen_m.values,
    eigenvectors: { rows: 3, cols: 3, values: geometry.eigen_m.vectors } }, input, [.5, 1, 3]);
  assert.deepEqual(geometry.result, symmetric);
  assert.equal(geometry.checks.passed, geometry.checks.total);
  assert.equal(api.eigenCalls, 1);
  assert.equal(api.owned.size, 7);
  disposed(api.owned);
  // Published geometry contains independent arrays that remain usable after disposal.
  assert.equal(geometry.eigen_m.vectors.length, 9);
  assert.ok(geometry.eigen_m.vectors.every(Number.isFinite));
});

test("nonsymmetric geometry returns general eigenpairs without invoking the symmetric solver", async () => {
  const api = trackedAPI(() => assert.fail("the symmetric solver must not handle this input"));
  const geometry = await computeGeometry(api.TrackedMatrix, { ...geometryConfig, m: asymmetric });
  assert.deepEqual(geometry.eigen_m.values,[1,1,1]);
  assert.deepEqual(geometry.eigen_m.imag_values,[0,0,0]);
  assert.equal(api.eigenCalls, 0);
  assert.equal(api.owned.size, 8);
  disposed(api.owned);
});

test("geometry rejects wrong, dependent, or malformed eigenbases and cleans up on failure", async () => {
  const corruptions = [
    { mutate: replaceWithIdentity, error: /input eigenpair residual/ },
    { mutate(result) {
      for (let row = 0; row < 3; row++) for (let col = 0; col < 3; col++) result.vectors.set(row, col, 0);
    }, error: /input eigenvector orthogonality/ },
    { mutate(result, track) {
      result.vectors.dispose();
      result.vectors = track(new Matrix(1, 3, [1, 0, 0]));
    }, error: /input eigenvectors has an invalid shape/ },
  ];
  for (const { mutate, error } of corruptions) {
    const api = trackedAPI(mutate);
    await assert.rejects(computeGeometry(api.TrackedMatrix, geometryConfig), error);
    disposed(api.owned);
  }
});
