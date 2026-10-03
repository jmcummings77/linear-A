import assert from "node:assert/strict";
import {test} from "node:test";
import {Matrix} from "./matrix.js";
test("minimum-norm solving rejects a missing RHS rather than returning an inverse",()=>{
  const a=new Matrix(1,2,[1,1]);
  for(const rhs of [undefined,null,[],0])assert.throws(()=>a.solveMinimumNorm(rhs as unknown as Matrix),TypeError);
  assert.deepEqual([...a.solveMinimumNorm(new Matrix(1,1,[0])).values],[0,0]);
});
test("strict cutoff, minimum-norm truncation and scaling-loss guard",()=>{
  const a=new Matrix(2,2,[4,0,0,1]), b=new Matrix(2,1,[4,1]);
  assert.deepEqual([...a.pseudoinverse(.25).values],[.25,0,0,0]);
  assert.deepEqual([...a.solveMinimumNorm(b,.25).values],[1,0]);
  assert.equal(a.spectralDiagnostics(.25).rank,1);
  assert.throws(()=>a.solveMinimumNorm(new Matrix(2,1,[1e300,1e-300])));
});
