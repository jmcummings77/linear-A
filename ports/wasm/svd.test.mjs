import assert from 'node:assert/strict';
import {test} from 'node:test';
import {createMatrixAPI} from './matrix.mjs';
const {Matrix}=await createMatrixAPI();
test('SVD options and failed decompositions do not invalidate inputs or leak output handles',()=>{
  const a=new Matrix(3,3,[1,2,3,4,5,7,6,8,9]);
  try{
    for(const args of [[0,100],[NaN,100],[1e-12,0],[1e-12,1.5],[1e-12,10001]])assert.throws(()=>a.svd(...args));
    for(let run=0;run<100;run++){
      assert.throws(()=>a.svd(1e-12,1));
      const r=a.svd();assert.equal(r.values.length,3);r.u.dispose();r.vt.dispose();
    }
    assert.equal(a.get(0,0),1);
  }finally{a.dispose();}
  assert.throws(()=>a.svd());
});
test('wide economy factors are independent and own their allocations',()=>{
  const a=new Matrix(2,3,[3,0,0,0,4,0]),r=a.svd();
  a.dispose();
  assert.deepEqual([...r.values],[4,3]);
  assert.deepEqual([r.u.rows,r.u.cols,r.vt.rows,r.vt.cols],[2,2,2,3]);
  r.u.dispose();assert.equal(r.vt.toArray().length,6);r.vt.dispose();
});
