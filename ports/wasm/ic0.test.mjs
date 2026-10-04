import assert from 'node:assert/strict';
import {test} from 'node:test';
import {execFileSync} from 'node:child_process';
import {createMatrixAPI} from './matrix.mjs';
import {checkSparse} from '../../benchmarks/sparse-checks.mjs';
import {CSRMatrix} from '../typescript/dist/sparse.js';
import {IC0} from '../typescript/dist/cholesky.js';
const wasm=await createMatrixAPI();
const fixtures=JSON.parse(execFileSync('python3',['-c','import sys,json;sys.path.insert(0,"benchmarks");from ic0_reference import fixtures;print(json.dumps(fixtures()))'],{encoding:'utf8'}));
for(const [name,api] of [['WASM',wasm],['TypeScript',{CSRMatrix,IC0}]]){
 test(`${name}: IC0 independent fixtures`,()=>assert.equal(checkSparse(api,fixtures),fixtures.length));
 test(`${name}: IC0 reuse, ownership and invalid preconditioners`,()=>{
  const a=new api.CSRMatrix(2,2,[0,2,4],[0,1,0,1],[4,1,1,3]),f=new api.IC0(a),small=new api.CSRMatrix(1,1,[0,1],[0],[1]),g=new api.IC0(small);
  try{
   const l=f.lower;l.values[0]=999;l.dispose?.();
   for(const b of [[6,7],[3,-2]]){const r=a.conjugateGradient(b,{preconditioner:f,capture:true});assert.equal(r.iterations,1);assert.equal(r.reason,'converged');assert.ok(Math.hypot(...a.matvec(r.x).map((v,i)=>v-b[i]))<1e-12);}
   assert.throws(()=>a.conjugateGradient([0,0],{jacobi:true,preconditioner:f}));assert.throws(()=>a.conjugateGradient([0,0],{preconditioner:g}));
   assert.throws(()=>a.conjugateGradient([0,0],{preconditioner:{}}));
   a.dispose?.();assert.ok(f.apply([6,7]).every((v,i)=>Math.abs(v-i-1)<1e-12));
   if(f.dispose){f.dispose();f.dispose();assert.throws(()=>f.apply([6,7]));assert.throws(()=>f.lower);}
  }finally{a.dispose?.();f.dispose?.();small.dispose?.();g.dispose?.();}
 });
}

test('live IC0 frames, orderings and SPD breakdown remain independently checked',async()=>{
 const {compareIC0}=await import('../../benchmarks/ic0-live.mjs'),{diffusionGrid}=await import('../../benchmarks/sparse-heat.mjs');
 for(const ordering of ['natural','rcm','amd']){
  const r=compareIC0(wasm,diffusionGrid(8,2),{ordering});
  for(const method of ['cg','jacobi','ic0','cholesky']){const m=r.runs[method];assert.equal(m.reason,'converged');assert.equal(m.frames.length,m.residualFields.length);assert.ok(m.residuals.at(-1)<m.residuals[0]*1e-9);}
  assert.ok(r.runs.ic0.nnz<r.runs.cholesky.nnz);assert.ok(r.runs.ic0.iterations<r.runs.cg.iterations);
 }
 const r=compareIC0(wasm,diffusionGrid(2),{breakdown:true});assert.equal(r.runs.ic0.ok,false);assert.equal(r.runs.cholesky.reason,'converged');
});
