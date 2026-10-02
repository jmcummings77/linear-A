import assert from 'node:assert/strict';
import {test} from 'node:test';
import {Matrix} from './matrix.js';

test('reusable factors preserve snapshots and RHS values',()=>{
 for(const algorithm of ['factorLU','factorCholesky','factorQR'] as const){
  const a=new Matrix(2,2,[4,1,1,3]), b=new Matrix(2,2,[6,5,7,4]), factor=a[algorithm]();
  a.set(0,0,99);
  for(let run=0;run<3;run++)factor.solve(b).values.forEach((v,i)=>assert.ok(Math.abs(v-[1,1,2,1][i])<1e-12));
  assert.ok(Math.abs(factor.reciprocalCondition()-.44)<1e-12);
  assert.deepEqual([...b.values],[6,5,7,4]);
  assert.throws(()=>factor.solve(new Matrix(1,1,[1])),RangeError);
 }
});
