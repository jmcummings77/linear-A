import assert from 'node:assert/strict';
import {test} from 'node:test';
import {execFileSync} from 'node:child_process';
import {createMatrixAPI} from './matrix.mjs';
import {checkSparse} from '../../benchmarks/sparse-checks.mjs';
import {compareAMD,graphProblem,eliminationTrace} from '../../benchmarks/amd-live.mjs';
import {diffusionGrid} from '../../benchmarks/sparse-heat.mjs';
import {CSRMatrix as TSCSR} from '../typescript/dist/sparse.js';
import {SparseCholeskySymbolic as TSPlan} from '../typescript/dist/cholesky.js';
const wasm=await createMatrixAPI();
const fixtures=JSON.parse(execFileSync('python3',['-c','import sys,json;sys.path.insert(0,"benchmarks");from amd_reference import fixtures;print(json.dumps(fixtures()))'],{encoding:'utf8'}));
for(const [name,api] of [['WASM',wasm],['TypeScript',{CSRMatrix:TSCSR,SparseCholeskySymbolic:TSPlan}]]){
 test(`${name}: AMD fixtures, original-system solves and failures`,()=>assert.equal(checkSparse(api,fixtures),fixtures.length));
 test(`${name}: live graph frames independently match every symbolic fill edge`,()=>{
  for(const grid of [diffusionGrid(4,2,{diffusivity:1,speed:0}),graphProblem(4,'tree'),graphProblem(4,'irregular')])for(const scrambled of [false,true]){
   const r=compareAMD(api,grid,{scrambled});
   for(const key of ['natural','rcm','amd']){const f=r[key];assert.ok(f.ok,f.error);assert.ok(f.residual<1e-11&&f.reconstructionError<1e-11);assert.equal(f.trace.frames.reduce((s,x)=>s+x.fill.length,0),f.fillCount);assert.deepEqual(f.trace.frames.map(f=>f.pivot),f.order);assert.deepEqual([...f.order].sort((a,b)=>a-b),Array.from({length:r.n},(_,i)=>i));}
  }
 });
 test(`${name}: ordering ownership and graph replay rejects corrupt fill`,()=>{
  const a=new api.CSRMatrix(3,3,[0,3,5,7],[0,1,2,0,1,0,2],[4,-1,-1,-1,2,-1,2]),p=a.approximateMinimumDegree();p[0]=999;
  assert.deepEqual(a.approximateMinimumDegree(),[1,0,2]);
  assert.throws(()=>eliminationTrace(a,[0,1,2],[]),/disagrees/);a.dispose?.();
  if(a.dispose)assert.throws(()=>a.approximateMinimumDegree());
 });
}
test('live graph input validation and reproducibility',()=>{assert.throws(()=>graphProblem(1,'tree'));assert.throws(()=>graphProblem(4,'other'));assert.deepEqual(graphProblem(4,'irregular'),graphProblem(4,'irregular'));});
