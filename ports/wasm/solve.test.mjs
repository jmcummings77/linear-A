import assert from 'node:assert/strict';
import {test} from 'node:test';
import {resolve} from 'node:path';
import {pathToFileURL} from 'node:url';
import {createMatrixAPI} from './matrix.mjs';
const options = process.env.LINEAR_A_WASM_MODULE ? {moduleUrl:pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE))} : {};
const {Matrix, Factorization} = await createMatrixAPI(options);

test('factors own snapshots, survive source disposal, and reuse multiple RHS', () => {
  for (const method of ['factorLU','factorCholesky','factorQR']) {
    const a=new Matrix(2,2,[4,1,1,3]), b=new Matrix(2,2,[6,5,7,4]), factor=a[method]();
    a.set(0,0,99); a.dispose();
    try {
      for(let run=0;run<3;run++) {
        const x=factor.solve(b);
        try { [...x.toArray()].forEach((v,i)=>assert.ok(Math.abs(v-[1,1,2,1][i])<1e-12)); }
        finally {x.dispose();}
      }
      assert.ok(Math.abs(factor.reciprocalCondition()-.44)<1e-12);
      assert.deepEqual([...b.toArray()],[6,5,7,4]);
      assert.throws(()=>new Factorization(1),TypeError);
    } finally {factor.dispose(); factor.dispose(); b.dispose();}
    assert.equal(factor.disposed,true);
    assert.throws(()=>factor.solve(b),/disposed/);
    assert.throws(()=>factor.reciprocalCondition(),/disposed/);
  }
});

test('factor rejects RHS from another runtime and disposed RHS without damaging itself',async()=>{
  const other=await createMatrixAPI(options);
  const a=Matrix.identity(2), b=other.Matrix.identity(2), factor=a.factorLU();
  try {
    assert.throws(()=>factor.solve(b),TypeError);
    b.dispose(); a.dispose();
    assert.throws(()=>factor.solve(a),/disposed/);
    const rhs=Matrix.identity(2),x=factor.solve(rhs);
    try {assert.deepEqual([...x.toArray()],[1,0,0,1]);} finally{x.dispose();rhs.dispose();}
  } finally {factor.dispose();a.dispose();b.dispose();}
});
