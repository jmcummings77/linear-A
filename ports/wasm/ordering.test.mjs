import assert from 'node:assert/strict';
import {test} from 'node:test';
import {execFileSync} from 'node:child_process';
import {createMatrixAPI} from './matrix.mjs';
import {checkSparse} from '../../benchmarks/sparse-checks.mjs';
const api=await createMatrixAPI();
const fixtures=JSON.parse(execFileSync('python3',['-c','import sys,json;sys.path.insert(0,"benchmarks");from ordering_reference import fixtures;print(json.dumps(fixtures()))'],{encoding:'utf8'}));
test('RCM and permutations pass the shared fixtures',()=>assert.equal(checkSparse(api,fixtures),fixtures.length));
test('permuted CSR owns storage; malformed and disposed inputs are rejected',()=>{
 const a=new api.CSRMatrix(3,3,[0,2,3,5],[0,2,1,0,2],[4,0,5,2,6]),p=[2,0,1],q=a.permuteSymmetric(p);p[0]=0;a.dispose();
 assert.equal(q.nnz,5);assert.deepEqual(q.matvec([3,1,2]),new Float64Array([20,4,10]));
 for(const p of [[0,0,2],[0,1],[0,1,3],[0,1,1.5]])assert.throws(()=>q.permuteSymmetric(p));
 assert.throws(()=>api.CSRMatrix.permuteVector([0],[Infinity]));assert.throws(()=>a.reverseCuthillMcKee());q.dispose();q.dispose();assert.throws(()=>q.permuteSymmetric([0,1,2]));
});
import {compareOrdering} from '../../benchmarks/ordering-live.mjs';
import {diffusionGrid} from '../../benchmarks/sparse-heat.mjs';
test('live comparison validates restored solutions and preserves the plotted pattern',()=>{
 for(const scrambled of [false,true]){
  const grid=diffusionGrid(8,2,{diffusivity:.2,speed:4,angle:30}),r=compareOrdering(api,grid,{scrambled});
  assert.equal(r.points.length,grid.values.length);assert.equal(new Set(r.p).size,grid.n);
  r.p.forEach((old,i)=>assert.equal(r.inverse[old],i));
  for(const solver of [r.natural,r.rcm]){assert.ok(solver.converged);assert.ok(solver.residual<=solver.threshold*1.00001);assert.equal(solver.residuals.length,solver.iterations+1);}
  if(scrambled)assert.ok(r.newWidth<r.width);
 }
});
test('Jacobi reference narrows bandwidth without improving iteration count',()=>{
 const r=compareOrdering(api,diffusionGrid(8,1,{diffusivity:.2,speed:4,angle:30}),{scrambled:true,preconditioner:'jacobi'});
 assert.ok(r.newWidth<r.width);assert.ok(r.natural.converged&&r.rcm.converged);assert.equal(r.natural.iterations,r.rcm.iterations);
 r.natural.residuals.forEach((v,i)=>assert.ok(Math.abs(v-r.rcm.residuals[i])<1e-9*Math.max(1,r.natural.residuals[0])));
});
