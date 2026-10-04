import assert from 'node:assert/strict';
import {test} from 'node:test';
import {execFileSync} from 'node:child_process';
import {createMatrixAPI} from './matrix.mjs';
import {checkSparse} from '../../benchmarks/sparse-checks.mjs';
import {CSRMatrix,GeometricMultigrid} from '../typescript/dist/matrix.js';
const wasm=await createMatrixAPI();
const fixtures=JSON.parse(execFileSync('python3',['-c','import sys,json;sys.path.insert(0,"benchmarks");from multigrid_reference import fixtures;print(json.dumps(fixtures()))'],{encoding:'utf8'}));
for(const [name,api] of [['WASM',wasm],['TypeScript',{CSRMatrix,GeometricMultigrid}]]){
 test(`${name}: multigrid exact dense fixtures`,()=>assert.equal(checkSparse(api,fixtures),fixtures.length));
 test(`${name}: multigrid owns results, validates input, captures actual residuals`,()=>{
  const m=new api.GeometricMultigrid(7),b=Array.from({length:49},(_,i)=>i%7-2),saved=b.slice();let a;
  try{
   const t=m.trace(b);assert.deepEqual(t.x,m.apply(b));assert.deepEqual(b,saved);assert.equal(t.frames.length,10);
   for(const f of t.frames){const level=new api.GeometricMultigrid(f.width),q=level.matrix;try{const ax=q.matvec(f.x);f.residual.forEach((v,i)=>assert.ok(Math.abs(v-(f.b[i]-ax[i]))<1e-12));}finally{q.dispose?.();level.dispose?.();}}
   t.x[0]=99;assert.notEqual(m.apply(b)[0],99);
   a=m.matrix;const r=a.conjugateGradient(b,{preconditioner:m,capture:true});assert.equal(r.reason,'converged');assert.ok(r.iterations<20);
   assert.throws(()=>a.conjugateGradient(b,{jacobi:true,preconditioner:m}));
   for(const w of [0,2,4,256,-1,3.5,NaN])assert.throws(()=>new api.GeometricMultigrid(w));
   for(const rhs of [[1],Array(49).fill(NaN),Array(49).fill(Infinity)])assert.throws(()=>m.apply(rhs));
   if(m.dispose){m.dispose();assert.throws(()=>m.apply(b));assert.throws(()=>m.matrix);}
  }finally{m.dispose?.();a?.dispose?.();}
 });
}
