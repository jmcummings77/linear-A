import assert from "node:assert/strict";
import { test } from "node:test";
import { readFile } from "node:fs/promises";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { createMatrixAPI } from "./matrix.mjs";
import { runChecks, runBenchmark } from "../../benchmarks/live-worker.mjs";

const moduleUrl = process.env.LINEAR_A_WASM_MODULE
  ? pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE))
  : new URL("../../.build/wasm/matrix.mjs", import.meta.url);
const wasmBinary = new Uint8Array(await readFile(new URL("matrix.wasm", moduleUrl)));
// Passing bytes must avoid a second fetch/read through locateFile.
const { Matrix } = await createMatrixAPI({ moduleUrl, wasmBinary, locateFile: () => "missing-embedded-file.wasm" });
const matrix = (rows, cols, values) => ({ rows, cols, values });
const a = matrix(2, 3, [1, 2, 3, 4, 5, 6]);
const fixtures = [
  { name: "rectangular transpose", op: "transpose", a, b: null, scalar: 1.25, invalid: false,
    expected: matrix(3, 2, [1, 4, 2, 5, 3, 6]) },
  { name: "rectangular multiplication", op: "multiply", a, b: matrix(3, 2, [7, 8, 9, 10, 11, 12]), scalar: 1.25, invalid: false,
    expected: matrix(2, 2, [58, 64, 139, 154]) },
  { name: "invalid multiplication", op: "multiply", a, b: a, scalar: 1.25, invalid: true },
  { name: "nonsquare trace", op: "trace", a, b: null, scalar: 1.25, invalid: true },
  { name: "known determinant", op: "determinant", a: matrix(3, 3, [6, 1, 1, 4, -2, 5, 2, 8, 7]), b: null, scalar: 1.25, invalid: false,
    expected: { value: -306 } },
  { name: "tiny off-diagonal", op: "triangular", a: matrix(2, 2, [1, Number.MIN_VALUE, Number.MIN_VALUE, 1]), b: null, scalar: 1.25, invalid: false,
    expected: { upper: false, lower: false } },
  { name: "empty determinant", op: "determinant", a: matrix(0, 0, []), b: null, scalar: 1.25, invalid: false,
    expected: { value: 1 } },
];

test("embedded bytes initialize WASM and independent fixtures pass with progress", async () => {
  const progress = [];
  const result = await runChecks(Matrix, fixtures, event => progress.push(event));
  assert.equal(result.passed, fixtures.length);
  assert.equal(result.total, fixtures.length);
  assert.equal(progress.length, fixtures.length);
  assert.equal(progress.at(-1).completed, fixtures.length);
  assert.ok(progress.every(event => event.message.startsWith("Checking ")));
});

test("checks expose wrong values, wrong keys, and invalid requests that succeed", async () => {
  const wrong = [
    { ...fixtures[4], expected: { value: -305 } },
    { ...fixtures[4], expected: { value: -306, extra: 1 } },
    { ...fixtures[0], invalid: true },
  ];
  const result = await runChecks(Matrix, wrong);
  assert.equal(result.passed, 0);
  assert.ok(result.checks.every(check => typeof check.error === "string" && check.error));
});

test("checks reject booleans and nonfinite numbers masquerading as scalars", async () => {
  const fixture = { name: "trace", op: "trace", a: matrix(1, 1, [1]), b: null, scalar: 1.25, invalid: false, expected: { value: 1 } };
  for (const value of [true, "1", NaN, Infinity]) {
    class WrongMatrix extends Matrix { trace() { return value; } }
    const result = await runChecks(WrongMatrix, [fixture]);
    assert.equal(result.passed, 0);
  }
});

test("floating comparison tolerances match the shared reference", async () => {
  const fixture = { name: "trace", op: "trace", a: matrix(1, 1, [1]), b: null, scalar: 1.25, invalid: false, expected: { value: 1 } };
  class CloseMatrix extends Matrix { trace() { return 1 + 5e-10; } }
  class FarMatrix extends Matrix { trace() { return 1 + 2e-9; } }
  assert.equal((await runChecks(CloseMatrix, [fixture])).passed, 1);
  assert.equal((await runChecks(FarMatrix, [fixture])).passed, 0);
});

test("live benchmarks check all seven operation checksums and report valid samples", async () => {
  for (const operation of ["add", "subtract", "scale", "transpose", "multiply", "trace", "determinant"]) {
    const result = await runBenchmark(Matrix, { operation, size: 16, samples: 3, seed: 17 });
    assert.equal(result.operation, operation);
    assert.equal(result.checksum_verified, true);
    assert.equal(result.samples.length, 3);
    assert.ok(Number.isSafeInteger(result.iterations) && result.iterations >= 1 && result.iterations <= 1_048_576);
    assert.ok(result.median_ns > 0 && result.mad_ns >= 0);
    assert.ok(result.min_ns <= result.median_ns && result.max_ns >= result.median_ns);
    assert.ok(result.samples.every(sample => sample.elapsed_ms > 0 && sample.ns_per_op > 0 && Number.isFinite(sample.checksum)));
  }
});

test("largest determinant oracle uses exact integer arithmetic before timing", async () => {
  const result = await runBenchmark(Matrix, { operation: "determinant", size: 48, samples: 3, seed: 17 });
  assert.equal(result.checksum_verified, true);
  assert.ok(result.samples.every(sample => Number.isFinite(sample.checksum)));
});

test("benchmark config is bounded and strictly typed", async () => {
  const valid = { operation: "multiply", size: 16, samples: 3, seed: 0 };
  for (const changes of [
    { operation: "triangular" }, { operation: "unknown" }, { size: 1 }, { size: 24 }, { size: 512 },
    { size: "16" }, { samples: 4 }, { samples: true }, { seed: -1 }, { seed: 2147483647 }, { seed: 0.5 }, { seed: true },
  ]) await assert.rejects(runBenchmark(Matrix, { ...valid, ...changes }));
  await assert.rejects(runBenchmark(Matrix, { ...valid, operation: "determinant", size: 256 }));
  await assert.rejects(runChecks(Matrix, []));
});

test("a false checksum is rejected before benchmark results are published", async () => {
  const owned = [];
  class WrongMatrix extends Matrix {
    constructor(...args) { super(...args); owned.push(this); }
    trace() { return 999; }
  }
  await assert.rejects(runBenchmark(WrongMatrix, { operation: "trace", size: 16, samples: 3, seed: 0 }), /checksum differs/);
  assert.equal(owned.length, 2);
  assert.ok(owned.every(matrix => matrix.disposed));
});

test("live determinant algorithms verify identical SPD inputs against the exact oracle", async () => {
  const runs = [];
  for (const operation of ["determinant_lu", "determinant_spd_lu", "determinant_cholesky"]) {
    const result = await runBenchmark(Matrix, { operation, size: 16, samples: 3, seed: 17 });
    assert.equal(result.checksum_verified, true);
    runs.push(result.samples[0].checksum / result.iterations);
  }
  assert.ok(Math.abs(runs[1] - runs[2]) <= Math.abs(runs[1]) * 1e-10);
  assert.notEqual(runs[0], runs[1]);
});
