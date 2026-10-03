import assert from 'node:assert/strict';
import {test} from 'node:test';
import {execFileSync} from 'node:child_process';
import {createMatrixAPI} from './matrix.mjs';
import {checkSparse} from '../../benchmarks/sparse-checks.mjs';
import {diffusionGrid,solveHeat} from '../../benchmarks/sparse-heat.mjs';
const api=await createMatrixAPI();
const fixtures=JSON.parse(execFileSync('python3',['-c','import sys,json;sys.path.insert(0,"benchmarks");from sparse_reference import fixtures;print(json.dumps(fixtures()))'],{encoding:'utf8'}));
test('WASM passes all sparse fixtures including captured true residuals',()=>assert.equal(checkSparse(api,fixtures),39));
test('CSR owns snapshots and rejects disposed use',()=>{
 const rp=[0,1],ci=[0],v=[2],a=new api.CSRMatrix(1,1,rp,ci,v);rp[1]=0;ci[0]=1;v[0]=99;
 assert.deepEqual(Array.from(a.matvec([3])),[6]);const result=a.conjugateGradient([6],{capture:true});
 assert.deepEqual(result.x,[3]);assert.deepEqual(result.iterates,[[0],[3]]);
 a.dispose();a.dispose();assert.throws(()=>a.matvec([3]));assert.throws(()=>a.conjugateGradient([6]));
 assert.deepEqual(result.x,[3]);
});
test('heat uses actual converged iterates and explicit iteration limits',()=>{
 for(const jacobi of [false,true]){const r=solveHeat(api,{size:8,contrast:2,jacobi,limit:400});assert.ok(r.converged);assert.equal(r.iterates.length,r.iterations+1);assert.ok(r.residuals.at(-1)<=r.threshold);assert.ok(r.x.every(x=>x>0&&x<100));}
 const r=solveHeat(api,{size:8,contrast:2,jacobi:true,limit:0});assert.equal(r.reason,'iteration_limit');assert.equal(r.converged,false);assert.deepEqual(r.iterates,[Array(64).fill(0)]);
});
test('heat grid is sparse, symmetric and enforces bounded settings',()=>{
 const a=diffusionGrid(32,3);assert.equal(a.values.length,5*1024-4*32);assert.equal(a.offsets.length,1025);
 for(const config of [{size:33},{contrast:4},{limit:801},{jacobi:'yes'}])assert.throws(()=>solveHeat(api,config));
});
test('browser gate detects corruption instead of presenting invalid output',()=>{
 class Broken extends api.CSRMatrix {conjugateGradient(...args){const r=super.conjugateGradient(...args);r.residuals[0]=0;return r;}}
 assert.throws(()=>checkSparse({...api,CSRMatrix:Broken},fixtures),/true residual/);
});
