import assert from "node:assert/strict";
import { test } from "node:test";
import { existsSync } from "node:fs";
import { readFile } from "node:fs/promises";
import { resolve } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";
import { createMatrixAPI } from "./matrix.mjs";
import { computeGeometry } from "../../benchmarks/geometry-worker.mjs";

const moduleUrl = process.env.LINEAR_A_WASM_MODULE
  ? pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE))
  : new URL("../../.build/wasm/matrix.mjs", import.meta.url);
const { Matrix } = await createMatrixAPI({ moduleUrl });
const traceUrl = process.env.LINEAR_A_WASM_TRACE_MODULE
  ? pathToFileURL(resolve(process.env.LINEAR_A_WASM_TRACE_MODULE))
  : new URL("../../.build/wasm-trace/matrix.mjs", import.meta.url);
let TraceMatrix;
if (existsSync(fileURLToPath(traceUrl)) && existsSync(fileURLToPath(new URL("matrix.wasm", traceUrl)))) {
  ({ Matrix: TraceMatrix } = await createMatrixAPI({ moduleUrl: traceUrl }));
}
const identity = [1, 0, 0, 0, 1, 0, 0, 0, 1];
const rotation = [0, -1, 0, 1, 0, 0, 0, 0, 1];
const shear = [1, 1, 0, 0, 1, 0, 0, 0, 1];
const config = { operation: "multiply", m: rotation, n: shear, scalar: 1.25 };

test("five actual WASM operations produce verified 3x3 matrices and all 125 field vectors", async () => {
  const expected = {
    multiply: { result: [0, -1, 0, 1, 1, 0, 0, 0, 1], determinant: 1 },
    add: { result: [1, 0, 0, 1, 1, 0, 0, 0, 2], determinant: 2 },
    subtract: { result: [-1, -2, 0, 1, -1, 0, 0, 0, 0], determinant: 0 },
    scale: { result: [0, -1.25, 0, 1.25, 0, 0, 0, 0, 1.25], determinant: 1.953125 },
    transpose: { result: [0, 1, 0, -1, 0, 0, 0, 0, 1], determinant: 1 },
  };
  for (const [operation, wanted] of Object.entries(expected)) {
    const progress = [];
    const geometry = await computeGeometry(Matrix, { ...config, operation }, { onProgress: event => progress.push(event) });
    assert.deepEqual(geometry.result, wanted.result);
    assert.equal(geometry.determinant_m, 1);
    assert.equal(geometry.determinant_result, wanted.determinant);
    assert.equal(geometry.points.length, 125);
    assert.equal(geometry.input_vectors.length, 125);
    assert.equal(geometry.output_vectors.length, 125);
    assert.equal(new Set(geometry.points.map(point => point.join(","))).size, 125);
    assert.ok(geometry.points.flat().every(value => [-1, -0.5, 0, 0.5, 1].includes(value)));
    assert.equal(geometry.checks.passed, geometry.checks.total);
    assert.ok(geometry.checks.total > 750);
    assert.ok(progress.length >= 2);
    assert.ok(geometry.input_vectors.flat().every(Number.isFinite));
    assert.ok(geometry.output_vectors.flat().every(Number.isFinite));
    assert.match(geometry.convention, /head x \+ A\*x/);
    assert.match(geometry.animation, /not a physical trajectory/);
    assert.equal(geometry.trace_available, false);
    assert.deepEqual(geometry.steps, []);
  }
});

test("column-vector composition preserves multiplication order and vector-field meaning", async () => {
  const qs = await computeGeometry(Matrix, config);
  const sq = await computeGeometry(Matrix, { ...config, m: shear, n: rotation });
  const index = qs.points.findIndex(point => point.every(value => value === 1));
  assert.deepEqual(qs.input_vectors[index], [-1, 1, 1]);
  assert.deepEqual(qs.output_vectors[index], [-1, 2, 1]);
  assert.deepEqual(sq.result, [1, -1, 0, 1, 0, 0, 0, 0, 1]);
  assert.deepEqual(sq.output_vectors[index], [0, 1, 1]);
  const origin = qs.points.findIndex(point => point.every(value => value === 0));
  assert.deepEqual(qs.output_vectors[origin], [0, 0, 0]);
  const unchanged = await computeGeometry(Matrix, { ...config, m: identity, n: identity });
  assert.deepEqual(unchanged.output_vectors[index], [1, 1, 1]); // The field is x, not displacement zero.
});

test("reflection and rank collapse retain determinant signs and correct field values", async () => {
  for (const [m, determinant, corner] of [
    [[-1, 0, 0, 0, 1, 0, 0, 0, 1], -1, [-1, 1, 1]],
    [[1, 1, 0, 0, 0, 0, 0, 0, 1], 0, [2, 0, 1]],
  ]) {
    const result = await computeGeometry(Matrix, { ...config, m, n: identity });
    assert.equal(result.determinant_m, determinant);
    assert.equal(result.determinant_result, determinant);
    assert.deepEqual(result.output_vectors.at(-1), corner);
  }
});

test("invalid geometry is bounded and strictly typed before WASM allocation", async () => {
  for (const changed of [
    { operation: "determinant" }, { operation: "unknown" }, { m: [] }, { n: Array(10).fill(0) },
    { m: [true, ...identity.slice(1)] }, { n: ["1", ...identity.slice(1)] },
    { m: [10.001, ...identity.slice(1)] }, { m: [NaN, ...identity.slice(1)] },
    { n: [Infinity, ...identity.slice(1)] }, { scalar: 10.001 }, { scalar: true }, { scalar: Infinity },
  ]) await assert.rejects(computeGeometry(Matrix, { ...config, ...changed }));
  await assert.rejects(computeGeometry(Matrix, null));
  const original = structuredClone(config);
  const result = await computeGeometry(Matrix, config);
  result.m[0] = 99;
  result.n[0] = 99;
  assert.deepEqual(config, original);
});

test("a corrupted WASM matrix result is rejected and all acquired handles are freed", async () => {
  const owned = [];
  class WrongMatrix extends Matrix {
    constructor(...args) { super(...args); owned.push(this); }
    multiply(other) {
      const result = super.multiply(other);
      owned.push(result);
      result.set(0, 0, result.get(0, 0) + 1);
      return result;
    }
  }
  await assert.rejects(computeGeometry(WrongMatrix, config), /matrix result.*independent reference/);
  assert.equal(owned.length, 3);
  assert.ok(owned.every(matrix => matrix.disposed));
});

test("the vector verification gate detects a bad value outside matrix result checks", async () => {
  class WrongVectorMatrix extends Matrix {
    multiply(other) {
      const result = super.multiply(other);
      if (other.cols === 125) result.set(0, 0, result.get(0, 0) + 1);
      return result;
    }
  }
  await assert.rejects(computeGeometry(WrongVectorMatrix, config), /M\*x at point 0.*independent reference/);
});

test("determinant verification tolerates roundoff but rejects real discrepancies and nonnumbers", async () => {
  class CloseMatrix extends Matrix { determinant() { return super.determinant() + 5e-10; } }
  class FarMatrix extends Matrix { determinant() { return super.determinant() + 2e-9; } }
  const result = await computeGeometry(CloseMatrix, config);
  assert.equal(result.checks.passed, result.checks.total);
  await assert.rejects(computeGeometry(FarMatrix, config), /determinant of M.*independent reference/);
  for (const bad of [true, "1", NaN, Infinity]) {
    class NonNumberMatrix extends Matrix { determinant() { return bad; } }
    await assert.rejects(computeGeometry(NonNumberMatrix, config), /determinant of M.*independent reference/);
  }
});

test("real trace WASM supplies all 27 events, source lines, and captured partial matrices", { skip: !TraceMatrix }, async () => {
  const geometry = await computeGeometry(Matrix, config, { TraceMatrix });
  assert.equal(geometry.trace_available, true);
  assert.equal(geometry.trace_source, "ports/c/matrix.c");
  assert.equal(geometry.steps.length, 27);
  const source = (await readFile(new URL("../c/matrix.c", import.meta.url), "utf8")).split("\n");
  const partial = Array(9).fill(0);
  for (const [index, step] of geometry.steps.entries()) {
    assert.equal(step.row, Math.floor(index / 9));
    assert.equal(step.col, Math.floor(index / 3) % 3);
    assert.equal(step.k, index % 3);
    partial[step.row * 3 + step.col] = step.sum;
    assert.deepEqual(step.matrix, partial);
    assert.deepEqual(step.call_path, ["Matrix.multiplyWithTrace", "wm_multiply_trace", "m_multiply"]);
    assert.match(source[step.source_line - 1], /m_trace_record/);
    assert.equal("elapsed_ns" in step, false);
  }
  assert.ok(geometry.steps.some(step => step.product === 0));
  assert.deepEqual(geometry.steps[0].matrix, Array(9).fill(0));
  assert.deepEqual(geometry.steps.at(-1).matrix, geometry.result);
  geometry.steps[0].matrix[0] = 99;
  assert.notEqual(geometry.steps[1].matrix[0], 99);
});

test("corrupted instrumented sums cannot be presented as an execution trace", { skip: !TraceMatrix }, async () => {
  const handles = [];
  class WrongTraceMatrix extends TraceMatrix {
    constructor(...args) { super(...args); handles.push(this); }
    multiplyWithTrace(other) {
      const traced = super.multiplyWithTrace(other);
      handles.push(traced.matrix);
      traced.steps[0].sum += 1;
      return traced;
    }
  }
  await assert.rejects(computeGeometry(Matrix, config, { TraceMatrix: WrongTraceMatrix }), /trace 0 accumulator.*independent reference/);
  assert.ok(handles.every(matrix => matrix.disposed));
});
