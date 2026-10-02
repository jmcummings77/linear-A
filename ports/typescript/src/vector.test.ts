import assert from "node:assert/strict";
import { spawnSync } from "node:child_process";
import test from "node:test";
import { Matrix } from "./matrix.js";

function matrixNear(actual: Matrix, expected: Matrix, tolerance = 1e-13): void {
  assert.deepEqual([actual.rows, actual.cols], [expected.rows, expected.cols]);
  const values = actual.values, wanted = expected.values;
  for (let i = 0; i < values.length; i++) assert.ok(Math.abs(values[i] - wanted[i]) <= tolerance);
}

test("cross products follow the right-hand rule, retain left orientation, and own their output", () => {
  for (const [rows, cols] of [[3, 1], [1, 3]]) for (const [br, bc] of [[3, 1], [1, 3]]) {
    const a = new Matrix(rows, cols, [1, 2, 3]), b = new Matrix(br, bc, [-2, 4, 1]), result = a.cross(b);
    assert.deepEqual([result.rows, result.cols], [rows, cols]);
    assert.deepEqual([...result.values], [-10, -7, 8]);
    assert.deepEqual([...b.cross(a).values], [10, 7, -8]);
    assert.deepEqual([...a.cross(a).values], [0, 0, 0]);
    result.set(0, 0, 99);
    assert.deepEqual([...a.values], [1, 2, 3]);
    assert.deepEqual([...b.values], [-2, 4, 1]);
  }
  assert.deepEqual([...new Matrix(3, 1, [1, 0, 0]).cross(new Matrix(3, 1, [0, 1, 0])).values], [0, 0, 1]);
  assert.deepEqual([...new Matrix(3, 1, [1e300, 0, 0]).cross(new Matrix(1, 3, [0, 1e-300, 0])).values], [0, 0, 1]);
});

test("cross products are perpendicular and satisfy bilinearity and Lagrange's identity", () => {
  const a = new Matrix(3, 1, [.5, -1.25, 2.5]), b = new Matrix(3, 1, [3, .25, -2]), c = new Matrix(3, 1, [2, -.5, .25]);
  const cross = a.cross(b).values;
  const dot = (left: Float64Array, right: Float64Array): number => left.reduce((sum, value, i) => sum + value * right[i], 0);
  assert.ok(Math.abs(dot(cross, a.values)) < 1e-13);
  assert.ok(Math.abs(dot(cross, b.values)) < 1e-13);
  assert.ok(Math.abs(dot(cross, cross) - (dot(a.values, a.values) * dot(b.values, b.values) - dot(a.values, b.values) ** 2)) < 1e-13);
  matrixNear(a.add(c).cross(b), a.cross(b).add(c.cross(b)));
  matrixNear(a.scale(-.5).cross(b), a.cross(b).scale(-.5));
});

test("principal rotation matrices actively rotate column vectors using the right-hand rule", () => {
  matrixNear(Matrix.rotation2D(Math.PI / 2).multiply(new Matrix(2, 1, [1, 0])), new Matrix(2, 1, [0, 1]));
  matrixNear(Matrix.rotation2D(-Math.PI / 2).multiply(new Matrix(2, 1, [1, 0])), new Matrix(2, 1, [0, -1]));
  const cases: [(angle: number) => Matrix, number[], number[], number[]][] = [
    [Matrix.rotationX, [0, 1, 0], [0, 0, 1], [1, 0, 0]],
    [Matrix.rotationY, [0, 0, 1], [1, 0, 0], [0, 1, 0]],
    [Matrix.rotationZ, [1, 0, 0], [0, 1, 0], [0, 0, 1]],
  ];
  for (const [rotation, basis, wanted, axis] of cases) {
    const result = rotation(Math.PI / 2);
    matrixNear(result.multiply(new Matrix(3, 1, basis)), new Matrix(3, 1, wanted));
    matrixNear(result, Matrix.rotationAxisAngle(new Matrix(1, 3, axis), Math.PI / 2));
  }
});

test("axis-angle rotations preserve the axis, are proper orthogonal, and compose without scale overflow", () => {
  for (const entries of [[1, 2, 3], [1e300, 2e300, 3e300], [1e-300, 2e-300, 3e-300], [Number.MIN_VALUE, 0, 0], [1.7e308, 1.7e308, 1.7e308]]) {
    const axis = new Matrix(3, 1, entries), before = axis.values, r = Matrix.rotationAxisAngle(axis, .5);
    assert.deepEqual(axis.values, before);
    matrixNear(r.transpose().multiply(r), Matrix.identity(3));
    assert.ok(Math.abs(r.determinant() - 1) < 1e-13);
    matrixNear(r.transpose(), Matrix.rotationAxisAngle(axis, -.5));
    matrixNear(r.transpose(), Matrix.rotationAxisAngle(axis.scale(-1), .5));
    matrixNear(r.multiply(Matrix.rotationAxisAngle(axis, .7)), Matrix.rotationAxisAngle(axis, 1.2));
    const largest = Math.max(...entries.map(Math.abs));
    const direction = new Matrix(3, 1, entries.map(value => value / largest));
    matrixNear(r.multiply(direction), direction);
  }
  matrixNear(Matrix.rotationAxisAngle(new Matrix(3, 1, [1, 2, 3]), 0), Matrix.identity(3));
  const tiny = Matrix.rotationAxisAngle(new Matrix(3, 1, [1, 1, 0]), 1e-10);
  assert.ok(tiny.get(0, 1) > 0);
  assert.ok(Math.abs(tiny.get(0, 1) / 2.5e-21 - 1) < 1e-13);
});

test("invalid vector shapes, rotation axes, angles, and overflowing cross products fail without mutation", () => {
  const a = new Matrix(3, 1, [1, 2, 3]);
  for (const wrong of [new Matrix(), new Matrix(2, 1, [1, 2]), new Matrix(2, 2, [1, 2, 3, 4])]) {
    assert.throws(() => a.cross(wrong), RangeError);
    assert.throws(() => wrong.cross(a), RangeError);
    assert.throws(() => Matrix.rotationAxisAngle(wrong, .5), RangeError);
  }
  assert.throws(() => Matrix.rotationAxisAngle(new Matrix(3, 1), .5), RangeError);
  for (const angle of [NaN, Infinity, -Infinity]) {
    for (const rotation of [Matrix.rotation2D, Matrix.rotationX, Matrix.rotationY, Matrix.rotationZ]) assert.throws(() => rotation(angle), RangeError);
    assert.throws(() => Matrix.rotationAxisAngle(a, angle), RangeError);
  }
  const large = new Matrix(3, 1, [1e308, 0, 0]);
  assert.throws(() => large.cross(new Matrix(3, 1, [0, 1e308, 0])), RangeError);
  assert.deepEqual([...large.values], [1e308, 0, 0]);
  assert.deepEqual([...a.values], [1, 2, 3]);
});

test("vector and rotation runners accept the protocol and enforce fixed benchmark sizes", () => {
  const runner = new URL("runner.js", import.meta.url);
  const invoke = (args: string[], input = "") => spawnSync(process.execPath, [runner.pathname, ...args], { input, encoding: "utf8" });
  let response = invoke(["check", "cross", "1", "3", "3", "1"], "1 2 3 -2 4 1");
  assert.equal(response.status, 0, response.stderr);
  assert.deepEqual(JSON.parse(response.stdout), { rows: 1, cols: 3, values: [-10, -7, 8] });
  response = invoke(["check", "rotation2d", "0", "0", String(Math.PI / 2)]);
  assert.equal(response.status, 0, response.stderr);
  assert.equal(JSON.parse(response.stdout).rows, 2);
  response = invoke(["check", "rotation3d", "3", "1", String(Math.PI / 2)], "0 0 1");
  assert.equal(response.status, 0, response.stderr);
  assert.equal(JSON.parse(response.stdout).rows, 3);
  for (const [op, size] of [["cross", 3], ["rotation2d", 2], ["rotation3d", 3]] as const) {
    response = invoke(["bench", op, String(size), "2", "17"]);
    assert.equal(response.status, 0, response.stderr);
    const result = JSON.parse(response.stdout);
    assert.ok(Number.isFinite(result.checksum) && result.checksum !== 0 && result.elapsed_ns > 0);
    assert.notEqual(invoke(["bench", op, String(size + 1), "2", "17"]).status, 0);
  }
  assert.notEqual(invoke(["check", "rotation2d", "1", "1", ".5"], "1").status, 0);
  assert.notEqual(invoke(["check", "rotation2d", "0", "0", ".5"], "1").status, 0);
});
