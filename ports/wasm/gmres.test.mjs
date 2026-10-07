import assert from 'node:assert/strict';
import {test} from 'node:test';
import {execFileSync} from 'node:child_process';
import {createMatrixAPI} from './matrix.mjs';
import {checkSparse} from '../../benchmarks/sparse-checks.mjs';
import {solveHeat,diffusionGrid} from '../../benchmarks/sparse-heat.mjs';
const api=await createMatrixAPI();
const fixtures=JSON.parse(execFileSync('python3',['-c','import sys,json;sys.path.insert(0,"benchmarks");from gmres_reference import fixtures;print(json.dumps(fixtures()))'],{encoding:'utf8'}));
test('GMRES passes independent analytic, LU, restart and failure fixtures',()=>assert.equal(checkSparse(api,fixtures),fixtures.length));
test('GMRES results own their data after matrix disposal',()=>{
 const a=new api.CSRMatrix(2,2,[0,2,3],[0,1,1],[3,1,2]);
 const r=a.gmres([5,4],{restart:1,capture:true});a.dispose();
 assert.ok(r.converged);assert.ok(r.restarts.length>0);assert.ok(Math.abs(r.x[0]-1)<1e-8);
 const x=r.x[0];r.iterates.at(-1)[0]=99;assert.equal(r.x[0],x);
 assert.throws(()=>a.gmres([5,4]));
});
test('transport demo uses true residuals with independent flow settings',()=>{
 for(const angle of [-90,0,135])for(const jacobi of [false,true]){
  const r=solveHeat(api,{solver:'gmres',size:8,contrast:1,diffusivity:.2,speed:4,angle,restart:20,jacobi});
  assert.ok(r.converged);assert.ok(r.residuals.at(-1)<=r.threshold);
  assert.equal(r.estimatedResiduals.length,r.residuals.length);
  assert.equal(r.iterates.length,r.iterations+1);
  assert.ok(r.x.every(x=>x>=-1e-7&&x<=100+1e-7));
 }
 assert.throws(()=>solveHeat(api,{solver:'cg',speed:1}));
 const r=solveHeat(api,{solver:'gmres',limit:0});assert.equal(r.reason,'iteration_limit');
});
test('browser gate rejects damaged GMRES estimates',()=>{
 class Broken extends api.CSRMatrix {gmres(...args){const r=super.gmres(...args);r.estimatedResiduals[0]=0;return r;}}
 assert.throws(()=>checkSparse({...api,CSRMatrix:Broken},fixtures));
});
test('analytic browser gate rejects wrong iterates with consistent residual histories',()=>{
 class Broken extends api.CSRMatrix {gmres(b,options){
  const r=super.gmres(b,options);r.iterates[1][1]*=2;
  const product=this.matvec(r.iterates[1]);
  const residual=b.reduce((value,entry,i)=>Math.hypot(value,entry-product[i]),0);
  r.residuals[1]=residual;r.estimatedResiduals[1]=residual;return r;
 }}
 for(const jacobi of [0,1]){
  const fixture=fixtures.find(item=>item.name===`GMRES shifted Jordan s=1 restart 2 Jacobi ${jacobi}`);
  assert.ok(fixture);assert.throws(()=>checkSparse({...api,CSRMatrix:Broken},[fixture]),/analytic iterate/);
 }
});
test('analytic browser gate rejects a damaged noninitial projected residual',()=>{
 class Broken extends api.CSRMatrix {gmres(...args){const r=super.gmres(...args);r.estimatedResiduals[1]*=.5;return r;}}
 const fixture=fixtures.find(item=>item.name==='GMRES shifted Jordan s=4 restart 1 Jacobi 0');
 assert.ok(fixture);assert.throws(()=>checkSparse({...api,CSRMatrix:Broken},[fixture]),/analytic residual history/);
});

test('live transport coefficients agree with the independent benchmark assembly',()=>{
 const references=JSON.parse(execFileSync('python3',['-c','import sys,json;sys.path.insert(0,"benchmarks");from gmres_reference import transport;print(json.dumps([transport(3,1,.2,4,a) for a in [-90,0,135]]))'],{encoding:'utf8'}));
 for(const [index,angle] of [-90,0,135].entries()){
  const live=diffusionGrid(3,1,{diffusivity:.2,speed:4,angle}),reference=references[index];
  assert.deepEqual(live.offsets,reference.offsets);assert.deepEqual(live.indices,reference.indices);
  live.values.forEach((value,i)=>assert.ok(Math.abs(value-reference.values[i])<1e-12));
 }
 const downward=diffusionGrid(3,0,{diffusivity:1,speed:4,angle:90});
 assert.ok(Math.abs(downward.b[0]-200)<1e-12);
 const upward=diffusionGrid(3,0,{diffusivity:1,speed:4,angle:-90});
 assert.equal(upward.b[0],100);
});
