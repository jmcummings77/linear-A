import assert from 'node:assert/strict';
import {test} from 'node:test';
import {Matrix} from './matrix.js';
test('ridge validates penalty and RHS, preserves inputs and handles extreme scales',()=>{
 const a=new Matrix(1,1,[2]),b=new Matrix(1,1,[3]);
 for(const lambda of [-1,NaN,Infinity])assert.throws(()=>a.solveRidge(b,lambda));
 assert.throws(()=>a.solveRidge(null as unknown as Matrix,1));
 assert.ok(Math.abs(a.solveRidge(b,1).get(0,0)-1.2)<1e-14);
 assert.equal(a.solveRidge(b,0).get(0,0),a.solveMinimumNorm(b).get(0,0));
 assert.equal(a.get(0,0),2);assert.equal(b.get(0,0),3);
 const x=new Matrix(1,1,[1e-300]).solveRidge(new Matrix(1,1,[1e300]),1e300);
 assert.ok(Math.abs(x.get(0,0)/1e-300-1)<1e-14);
});
