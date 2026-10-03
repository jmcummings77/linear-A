import assert from 'node:assert/strict';
import {test} from 'node:test';
import {createMatrixAPI} from './matrix.mjs';
import {resolve} from 'node:path';
import {pathToFileURL} from 'node:url';
const {Matrix}=await createMatrixAPI(process.env.LINEAR_A_WASM_MODULE?{moduleUrl:pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE))}:{});
test('cutoff boundary is strict; factors and results own storage',()=>{
 const a=new Matrix(2,2,[4,0,0,1]);let p,b,x;
 try {
  p=a.pseudoinverse(.25);assert.deepEqual([...p.toArray()],[.25,0,0,0]);
  assert.deepEqual(a.spectralDiagnostics(.25),{rank:1,reciprocalCondition:.25,retainedReciprocalCondition:1});
  b=new Matrix(2,2,[4,8,1,2]);x=a.solveMinimumNorm(b,.25);
  assert.deepEqual([...x.toArray()],[1,2,0,0]);a.set(0,0,8);
  assert.deepEqual([...p.toArray()],[.25,0,0,0]);assert.deepEqual([...b.toArray()],[4,8,1,2]);
  for(const cutoff of [-1,NaN,Infinity,1.1,'0']) {
   assert.throws(()=>a.pseudoinverse(cutoff));assert.throws(()=>a.solveMinimumNorm(b,cutoff));assert.throws(()=>a.spectralDiagnostics(cutoff));
  }
 } finally {a.dispose();p?.dispose();b?.dispose();x?.dispose();}
 assert.throws(()=>a.spectralDiagnostics());assert.throws(()=>a.pseudoinverse());
});
test('direct SVD solve can succeed when forming the inverse would overflow',()=>{
 const a=new Matrix(1,1,[1e-310]),b=new Matrix(1,1,[1e-310]);let x;
 try {assert.throws(()=>a.pseudoinverse(0));x=a.solveMinimumNorm(b,0);assert.ok(Math.abs(x.get(0,0)-1)<1e-14);}
 finally {a.dispose();b.dispose();x?.dispose();}
});
