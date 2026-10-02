import assert from 'node:assert/strict';
import {test} from 'node:test';
import {resolve} from 'node:path';
import {pathToFileURL} from 'node:url';
import {execFileSync} from 'node:child_process';
import {createMatrixAPI} from './matrix.mjs';
import {computeAccuracy,diagnostics} from '../../benchmarks/accuracy-worker.mjs';
import {runChecks} from '../../benchmarks/live-worker.mjs';
const {Matrix}=await createMatrixAPI(process.env.LINEAR_A_WASM_MODULE?{moduleUrl:pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE))}:{});
const stable={a:[4,1,1,3],b:[6,7],target:5,delta:1e-7,algorithm:'lu'};

test('all shared solver fixtures execute in the live worker',async()=>{
 const fixtures=JSON.parse(execFileSync(process.env.PYTHON||'python3',['-c','import sys,json;sys.path.insert(0,"benchmarks");from solve_reference import solve_fixtures;print(json.dumps(solve_fixtures()))'],{encoding:'utf8'}));
 const checks=await runChecks(Matrix,fixtures);
 assert.equal(checks.passed,fixtures.length,JSON.stringify(checks));
});
test('small residual and large sensitivity are distinct diagnostics',()=>{
 const well=computeAccuracy(Matrix,stable);
 const sensitive=computeAccuracy(Matrix,{...stable,a:[1,1,1,1.000001],b:[2,2.000001]});
 assert.deepEqual(well.base.x,[1,2]);
 assert.ok(Math.abs(well.base.reciprocal_condition-.44)<1e-14);
 assert.ok(sensitive.relative_solution_change>.09);
 assert.ok(sensitive.relative_solution_change>well.relative_solution_change*1e5);
 assert.ok(sensitive.perturbed.backward_error<1e-14);
 assert.ok(sensitive.perturbed.reciprocal_condition<1e-6);
});
test('solver choices, singular perturbation, rounded-away delta, and input rejection',()=>{
 for(const algorithm of ['lu','cholesky','qr']){
  const r=computeAccuracy(Matrix,{...stable,algorithm});
  r.base.x.forEach((x,i)=>assert.ok(Math.abs(x-[1,2][i])<1e-12));
 }
 const singular=computeAccuracy(Matrix,{a:[1,1,1,1.5],b:[2,2.5],target:3,delta:-.5,algorithm:'lu'});
 assert.match(singular.perturbed_error,/singular/i);assert.equal(singular.perturbed,undefined);
 const rounded=computeAccuracy(Matrix,{...stable,delta:1e-20});assert.equal(rounded.actual_delta,0);
 for(const update of [{delta:2},{delta:NaN},{target:6},{a:[1,2]},{b:[1,Infinity]},{algorithm:'inverse'}])assert.throws(()=>computeAccuracy(Matrix,{...stable,...update}),RangeError);
 assert.throws(()=>computeAccuracy(Matrix,{...stable,a:[1,1,2,2]}),/singular/i);
 assert.deepEqual(diagnostics([0,0,0,0],[0,0],[0,0]),{backward_error:0,residual_infinity:0});
});
