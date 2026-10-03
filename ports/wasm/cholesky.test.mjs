import assert from 'node:assert/strict';
import {test} from 'node:test';
import {execFileSync} from 'node:child_process';
import {createMatrixAPI} from './matrix.mjs';
import {checkSparse} from '../../benchmarks/sparse-checks.mjs';
import {compareCholesky} from '../../benchmarks/cholesky-live.mjs';
import {diffusionGrid} from '../../benchmarks/sparse-heat.mjs';
import {CSRMatrix as TSCSR} from '../typescript/dist/sparse.js';
import {SparseCholeskySymbolic as TSPlan} from '../typescript/dist/cholesky.js';
const wasm=await createMatrixAPI();
const fixtures=JSON.parse(execFileSync('python3',['-c','import sys,json;sys.path.insert(0,"benchmarks");from cholesky_reference import fixtures;print(json.dumps(fixtures()))'],{encoding:'utf8'}));
for(const [name,api] of [['WASM',wasm],['TypeScript',{CSRMatrix:TSCSR,SparseCholeskySymbolic:TSPlan}]]){
 test(`${name}: shared Cholesky fixtures`,()=>assert.equal(checkSparse(api,fixtures),fixtures.length));
 test(`${name}: plan reuse, ownership, pattern rejection and solve failures`,()=>{
  const {CSRMatrix,SparseCholeskySymbolic}=api;
  const a=new CSRMatrix(2,2,[0,2,4],[0,1,0,1],[4,2,2,3]),s=new SparseCholeskySymbolic(a),f=s.factorize(a);
  a.dispose?.();s.rowOffsets[1]=99;s.columnIndices[0]=99;s.fillSteps[0]=99;
  const changed=new CSRMatrix(2,2,[0,2,4],[0,1,0,1],[8,4,4,6]),g=s.factorize(changed),mismatch=new CSRMatrix(2,2,[0,1,2],[0,1],[4,3]);
  assert.throws(()=>s.factorize(mismatch));
  for(const [factor,rhs] of [[f,[8,8]],[g,[16,16]]]){
   const x=factor.solve(rhs);assert.ok(Math.abs(x[0]-1)<1e-12&&Math.abs(x[1]-2)<1e-12);
   assert.throws(()=>factor.solve([1]));assert.throws(()=>factor.solve([NaN,1]));
   assert.throws(()=>factor.solve([Number.MAX_VALUE, -Number.MAX_VALUE]));
   assert.deepEqual(factor.solve([0,0]),[0,0]);const l=factor.lower;l.dispose?.();
  }
  s.dispose?.();changed.dispose?.();assert.ok(Math.abs(f.solve([8,8])[0]-1)<1e-12);
  if(s.dispose){assert.throws(()=>s.factorize(mismatch));s.dispose();}
  f.dispose?.();g.dispose?.();mismatch.dispose?.();if(f.dispose){f.dispose();assert.throws(()=>f.solve([8,8]));assert.throws(()=>f.lower);}
 });
 test(`${name}: live fill histories and reconstruction checks`,()=>{
  for(const scrambled of [false,true]){
   const r=compareCholesky(api,diffusionGrid(8,2,{diffusivity:1,speed:0}),{scrambled});
   for(const f of [r.natural,r.rcm]){assert.ok(f.ok,f.error);assert.ok(f.residual<1e-11);assert.ok(f.reconstructionError<1e-11);assert.equal(f.fillCount,f.points.filter(p=>p[2]>=0).length);for(const [i,j,k] of f.points)assert.ok(j<=i&&(k===-1||k<j));}
   if(scrambled)assert.ok(r.rcm.fillCount<r.natural.fillCount);
  }
 });
}
