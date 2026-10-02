import test from 'node:test';
import assert from 'node:assert/strict';
import {execFileSync} from 'node:child_process';
import {resolve} from 'node:path';
import {pathToFileURL} from 'node:url';
import {createMatrixAPI} from './matrix.mjs';
import {runChecks,runBenchmark,assertGeneralEigenResult} from '../../benchmarks/live-worker.mjs';
import {computeGeometry} from '../../benchmarks/geometry-worker.mjs';
const options=process.env.LINEAR_A_WASM_MODULE?{moduleUrl:pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE))}:{};
const {Matrix}=await createMatrixAPI(options);
const fixtures=JSON.parse(execFileSync(process.env.PYTHON||'python3',['-c',
  'import sys,json;sys.path.insert(0,"benchmarks");from general_eigen_reference import general_eigen_fixtures;print(json.dumps(general_eigen_fixtures()))'],
  {cwd:resolve(import.meta.dirname,'../..'),encoding:'utf8'}));
function tracked(corrupt=false) {
  const owned=[];let calls=0;
  class Tracked extends Matrix {
    static [Symbol.hasInstance](value){return value instanceof Matrix;}
    constructor(...args){super(...args);owned.push(this);}
    eigenGeneral(...args){calls++;const r=super.eigenGeneral(...args);owned.push(r.vectorsReal,r.vectorsImag);
      if(corrupt)for(let i=0;i<r.vectorsImag.rows;i++)for(let j=0;j<r.vectorsImag.cols;j++)r.vectorsImag.set(i,j,0);
      return r;}
  }
  return {Tracked,owned,get calls(){return calls;}};
}
test('live checks validate complex spectra, defective columns and ownership',async()=>{
  const api=tracked(),r=await runChecks(api.Tracked,fixtures);
  assert.equal(r.passed,20,JSON.stringify(r.checks.filter(c=>!c.passed)));
  assert.ok(api.owned.every(m=>m.disposed));
});
test('general benchmark validates complex timed input before measurements',async()=>{
  const api=tracked(),r=await runBenchmark(api.Tracked,{operation:'eigen_general',size:4,samples:3,seed:17});
  assert.equal(r.samples.length,3);assert.ok(api.calls>3*r.iterations);assert.ok(api.owned.every(m=>m.disposed));
  const bad=tracked(true);
  await assert.rejects(runBenchmark(bad.Tracked,{operation:'eigen_general',size:4,samples:3,seed:17}),/normalized|residual/);
  assert.equal(bad.calls,1);assert.ok(bad.owned.every(m=>m.disposed));
});
test('rotation geometry exposes a verified complex invariant plane',async()=>{
  const m=[0,-1,0,1,0,0,0,0,1],n=[1,0,0,0,1,0,0,0,1];
  const scene=await computeGeometry(Matrix,{operation:'multiply',m,n,scalar:1.25});
  const e=scene.eigen_m;
  assert.deepEqual(e.imag_values,[-1,1,0]);
  assertGeneralEigenResult({eigenvalues_real:e.values,eigenvalues_imag:e.imag_values,
    eigenvectors_real:{rows:3,cols:3,values:e.vectors},eigenvectors_imag:{rows:3,cols:3,values:e.imag_vectors}},
    {rows:3,cols:3,values:m},[[0,-1],[0,1],[1,0]]);
});

test('geometry rejects duplicated valid eigenpairs that omit the complex spectrum',async()=>{
  class MissingPair extends Matrix {
    eigenGeneral(){return {valuesReal:[1,1,1],valuesImag:[0,0,0],
      vectorsReal:new Matrix(3,3,[0,0,0,0,0,0,1,1,1]),vectorsImag:new Matrix(3,3)};}
  }
  await assert.rejects(computeGeometry(MissingPair,{operation:'multiply',m:[0,-1,0,1,0,0,0,0,1],n:[1,0,0,0,1,0,0,0,1],scalar:1.25}),/complex spectrum/);
});
