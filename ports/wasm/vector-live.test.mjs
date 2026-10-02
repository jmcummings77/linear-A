import assert from "node:assert/strict";
import test from "node:test";
import { execFileSync } from "node:child_process";
import { existsSync } from "node:fs";
import { resolve } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";
import { createMatrixAPI } from "./matrix.mjs";
import { runBenchmark, runChecks } from "../../benchmarks/live-worker.mjs";
import { computeGeometry } from "../../benchmarks/geometry-worker.mjs";

const options = process.env.LINEAR_A_WASM_MODULE
  ? { moduleUrl: pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE)) } : {};
const { Matrix } = await createMatrixAPI(options);
const traceUrl = process.env.LINEAR_A_WASM_TRACE_MODULE
  ? pathToFileURL(resolve(process.env.LINEAR_A_WASM_TRACE_MODULE))
  : new URL("../../.build/wasm-trace/matrix.mjs", import.meta.url);
let TraceMatrix;
if (existsSync(fileURLToPath(traceUrl))) ({ Matrix: TraceMatrix } = await createMatrixAPI({ moduleUrl: traceUrl }));
const root = fileURLToPath(new URL("../../", import.meta.url));
const fixtures = JSON.parse(execFileSync(process.env.PYTHON || "python3", ["-c",
  "import json,sys;sys.path.insert(0,'benchmarks');from vector_reference import vector_fixtures;print(json.dumps(vector_fixtures(),allow_nan=False))"],
{ cwd: root, encoding: "utf8" }));
const identity = [1,0,0, 0,1,0, 0,0,1];
const field = [1,2,0, 0,1,3, 0,0,2];
const base = { m: field, n: identity, scalar: 1.25, axis: [0,0,2] };

function near(actual, expected, tolerance = 2e-12) {
  assert.ok(Number.isFinite(actual) && Math.abs(actual-expected) <= tolerance, `${actual} differs from ${expected}`);
}
function nearArray(actual, expected, tolerance) {
  assert.equal(actual.length, expected.length);
  actual.forEach((value,index) => near(value,expected[index],tolerance));
}
function determinant(a) {
  return a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6]);
}
function cross(a,b) { return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]; }
function apply(a,x) { return [0,1,2].map(row => a[row*3]*x[0]+a[row*3+1]*x[1]+a[row*3+2]*x[2]); }

/** Track base-class result handles as well as subclass constructors and factories. */
function trackedAPI({ afterCross, afterRotation, afterMultiply } = {}) {
  const owned = new Set(), rotationAngles = [];
  let crossCalls = 0, multiplyCalls = 0;
  function track(matrix) {
    if (owned.has(matrix)) return matrix;
    owned.add(matrix);
    const multiply = matrix.multiply, vectorCross = matrix.cross, eigen = matrix.eigenSymmetric;
    matrix.multiply = (...args) => {
      const result = track(multiply.apply(matrix,args));
      multiplyCalls++;
      afterMultiply?.(result,{source:matrix,call:multiplyCalls,args});
      return result;
    };
    matrix.cross = (...args) => {
      const result = track(vectorCross.apply(matrix,args));
      crossCalls++;
      afterCross?.(result,{source:matrix,call:crossCalls,args});
      return result;
    };
    matrix.eigenSymmetric = (...args) => {
      const result = eigen.apply(matrix,args); track(result.vectors); return result;
    };
    return matrix;
  }
  function rotation(result,angle,kind) {
    track(result);
    rotationAngles.push(angle);
    afterRotation?.(result,{angle,kind,call:rotationAngles.length});
    return result;
  }
  class TrackedMatrix extends Matrix {
    static [Symbol.hasInstance](value) { return value instanceof Matrix; }
    constructor(...args) { super(...args); track(this); }
    static rotation2D(angle) { return rotation(Matrix.rotation2D(angle),angle,"2d"); }
    static rotationAxisAngle(axis,angle) { return rotation(Matrix.rotationAxisAngle(axis,angle),angle,"3d"); }
  }
  return { TrackedMatrix, owned, rotationAngles, get crossCalls() { return crossCalls; }, get multiplyCalls() { return multiplyCalls; } };
}
function disposed(api) {
  assert.ok(api.owned.size > 0);
  assert.ok([...api.owned].every(matrix => matrix.disposed), "all acquired WASM handles must be disposed");
}

const workloads = [["cross",3], ["rotation2d",2], ["rotation3d",3]];
test("live vector and rotation checks run all shared fixtures with handle cleanup", async () => {
  const api = trackedAPI(), progress = [];
  const result = await runChecks(api.TrackedMatrix,fixtures,event => progress.push(event));
  assert.ok(fixtures.length >= 36);
  assert.equal(result.passed,fixtures.length,JSON.stringify(result.checks.filter(check => !check.passed)));
  assert.equal(result.total,fixtures.length);
  assert.equal(progress.length,fixtures.length);
  assert.ok(api.crossCalls > 0 && api.rotationAngles.length > 0);
  disposed(api);
});

test("live cross and rotation benchmarks return verified samples at their fixed dimensions", async () => {
  for (const [operation,size] of workloads) {
    const api = trackedAPI();
    const result = await runBenchmark(api.TrackedMatrix,{operation,size,samples:3,seed:17});
    assert.equal(result.operation,operation);
    assert.equal(result.size,size);
    assert.equal(result.checksum_verified,true);
    assert.equal(result.samples.length,3);
    assert.ok(result.iterations >= 1 && result.iterations <= 4096 && result.median_ns > 0);
    assert.ok(result.samples.every(sample => Number.isFinite(sample.checksum) && sample.elapsed_ms > 0));
    if (operation === "cross") assert.ok(result.samples.every(sample => sample.checksum !== 0));
    assert.ok(api.crossCalls+api.rotationAngles.length > 3*result.iterations, "probe, warmup and calibration also execute the operation");
    disposed(api);
  }
});

test("live vector benchmark dimensions are strict and rejected before allocation", async () => {
  for (const [operation,size] of workloads) for (const bad of [0,1,4,16,String(size),true,NaN]) {
    const api = trackedAPI();
    await assert.rejects(runBenchmark(api.TrackedMatrix,{operation,size:bad,samples:3,seed:17}),/size/);
    assert.equal(api.owned.size,0);
  }
  for (const [operation,size] of workloads) {
    const api = trackedAPI();
    await assert.rejects(runBenchmark(api.TrackedMatrix,{operation,size:size === 2 ? 3 : 2,samples:3,seed:17}),/size/);
    assert.equal(api.owned.size,0);
  }
});

test("live full-result probes reject wrong vector or rotation results with unchanged checksums", async () => {
  for (const [operation,size] of workloads) {
    const corrupt = result => {
      const checksum = result.checksum();
      if (operation === "cross") {
        result.set(0,0,result.get(0,0)+.25);
        result.set(1,0,result.get(1,0)-.25);
      } else result.set(0,1,result.get(0,1)+.25); // The checksum samples diagonal/end entries.
      assert.equal(result.checksum(),checksum);
    };
    const api = trackedAPI(operation === "cross" ? {afterCross:corrupt} : {afterRotation:corrupt});
    await assert.rejects(runBenchmark(api.TrackedMatrix,{operation,size,samples:3,seed:17}),/differs from the independent fixture/);
    assert.equal(api.crossCalls+api.rotationAngles.length,1,"reject the probe before entering timed loops");
    disposed(api);
  }
});

test("cross geometry maps every field vector to (Mx) cross the unnormalized direction", async () => {
  const api = trackedAPI(), config = {...base,operation:"cross"}, before = structuredClone(config);
  const result = await computeGeometry(api.TrackedMatrix,config,{TraceMatrix});
  assert.deepEqual(result.result,[0,2,6, -2,-4,0, 0,0,0]);
  assert.equal(result.determinant_result,0);
  assert.equal(result.points.length,125);
  for (let i = 0; i < result.points.length; i++) {
    nearArray(result.input_vectors[i],apply(field,result.points[i]));
    nearArray(result.output_vectors[i],cross(result.input_vectors[i],base.axis));
  }
  assert.equal(api.crossCalls,3,"the three operator columns are crossed in WASM");
  assert.equal(result.trace_available,false);
  assert.deepEqual(result.steps,[]);
  assert.deepEqual(result.rotation_frames,[]);
  assert.equal(result.checks.passed,result.checks.total);
  assert.deepEqual(config,before);
  disposed(api);
});

test("zero cross direction is valid and collapses the field while zero rotation axis is rejected", async () => {
  const api = trackedAPI();
  const result = await computeGeometry(api.TrackedMatrix,{...base,operation:"cross",axis:[0,0,0]});
  assert.ok(result.result.every(value => value === 0));
  assert.ok(result.output_vectors.flat().every(value => value === 0));
  assert.equal(result.determinant_result,0);
  disposed(api);
  const invalid = trackedAPI();
  await assert.rejects(computeGeometry(invalid.TrackedMatrix,{...base,operation:"rotate",axis:[0,0,0],radians:0}),/nonzero/);
  assert.equal(invalid.owned.size,0);
});

test("rotation geometry samples 901 actual WASM angle frames and preserves volume at a half-turn midpoint", async () => {
  const api = trackedAPI(), config = {...base,operation:"rotate",axis:[0,0,1],radians:Math.PI};
  const before = structuredClone(config);
  const result = await computeGeometry(api.TrackedMatrix,config,{TraceMatrix});
  assert.equal(result.rotation_frames.length,901);
  assert.match(result.animation,/angle samples/);
  assert.doesNotMatch(result.animation,/Linear blend/);
  assert.equal(api.rotationAngles.length,902,"one endpoint factory call plus one per playback frame");
  nearArray(result.rotation_frames[0],field);
  nearArray(result.rotation_frames[450],[0,-1,-3, 1,2,0, 0,0,2]);
  nearArray(result.rotation_frames[900],result.result);
  nearArray(result.result,[-1,-2,0, 0,-1,-3, 0,0,2]);
  for (const [index,frame] of result.rotation_frames.entries()) {
    near(api.rotationAngles[index+1],Math.PI*index/900);
    near(determinant(frame),2);
    for (let row = 0; row < 3; row++) assert.ok(result.rotation_extent >= .55*(Math.abs(frame[row*3])+Math.abs(frame[row*3+1])+Math.abs(frame[row*3+2])));
    if (index) assert.notEqual(frame,result.rotation_frames[index-1]);
  }
  // Blending the endpoints linearly would collapse both x and y at t=1/2.
  near(determinant(field.map((value,index) => (value+result.result[index])/2)),0);
  near(result.determinant_m,2); near(result.determinant_result,2);
  assert.equal(result.trace_available,false);
  assert.deepEqual(result.steps,[]);
  assert.equal(result.checks.passed,result.checks.total);
  assert.deepEqual(config,before);
  disposed(api);
});

test("the closed angle bounds allow complete turns without changing the endpoint field", async () => {
  for (const radians of [-2*Math.PI,2*Math.PI]) {
    const api = trackedAPI();
    const result = await computeGeometry(api.TrackedMatrix,{...base,operation:"rotate",radians});
    nearArray(result.result,field);
    nearArray(result.rotation_frames[900],field);
    assert.equal(result.rotation_frames.length,901);
    disposed(api);
  }
});

test("rotation and cross geometry validate bounded directions and angles before allocating", async () => {
  for (const operation of ["cross","rotate"]) for (const axis of [undefined,[],[1,2],[1,2,3,4],[true,0,0],["1",0,0],[NaN,0,0],[Infinity,0,0],[10.001,0,0]]) {
    const api = trackedAPI();
    await assert.rejects(computeGeometry(api.TrackedMatrix,{...base,operation,axis,radians:.5}),/direction/);
    assert.equal(api.owned.size,0);
  }
  for (const radians of [undefined,NaN,Infinity,-Infinity,true,"0.5",2*Math.PI+1e-10,-2*Math.PI-1e-10]) {
    const api = trackedAPI();
    await assert.rejects(computeGeometry(api.TrackedMatrix,{...base,operation:"rotate",radians}),/angle/);
    assert.equal(api.owned.size,0);
  }
});

test("bad cross columns and a corrupted intermediate rotation frame are rejected without leaked handles", async () => {
  const brokenCross = trackedAPI({afterCross(result) { result.set(0,0,result.get(0,0)+1); }});
  await assert.rejects(computeGeometry(brokenCross.TrackedMatrix,{...base,operation:"cross"}),/matrix result.*independent reference/);
  disposed(brokenCross);
  const brokenRotation = trackedAPI({afterRotation(result,{call}) {
    if (call === 102) result.set(0,0,result.get(0,0)+.25);
  }});
  await assert.rejects(computeGeometry(brokenRotation.TrackedMatrix,{...base,operation:"rotate",radians:Math.PI}),/rotation path.*independent reference/);
  assert.equal(brokenRotation.rotationAngles.length,102);
  disposed(brokenRotation);
});

test("the existing multiplication trace remains separate from rotation path samples", {skip:!TraceMatrix}, async () => {
  const result = await computeGeometry(Matrix,{...base,operation:"multiply"},{TraceMatrix});
  assert.equal(result.trace_available,true);
  assert.equal(result.steps.length,27);
  assert.deepEqual(result.steps.at(-1).matrix,field);
  assert.deepEqual(result.rotation_frames,[]);
  assert.equal(result.trace_source,"ports/c/matrix.c");
});
