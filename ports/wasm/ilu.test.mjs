import assert from 'node:assert/strict';
import {test} from 'node:test';
import {execFileSync} from 'node:child_process';
import {createMatrixAPI} from './matrix.mjs';
import {checkSparse} from '../../benchmarks/sparse-checks.mjs';
import {createHeatSolver} from '../../benchmarks/sparse-heat.mjs';
const api=await createMatrixAPI();
const fixtures=JSON.parse(execFileSync('python3',['-c','import sys,json;sys.path.insert(0,"benchmarks");from ilu_reference import fixtures;print(json.dumps(fixtures()))'],{encoding:'utf8'}));
test('ILU passes independent exact masked-Doolittle fixtures',()=>assert.equal(checkSparse(api,fixtures),21));
test('ILU owns factors, reuses them across RHS and rejects disposed use',()=>{
 const a=new api.CSRMatrix(2,2,[0,2,4],[0,1,0,1],[4,1,1,3]),f=new api.ILU0(a);a.dispose();
 for(const [b,x] of [[[6,7],[1,2]],[[11,13],[20/11,41/11]]]){
  const actual=Array.from(f.apply(b));actual.forEach((v,i)=>assert.ok(Math.abs(v-x[i])<1e-12));
 }
 assert.throws(()=>f.apply([1]));assert.throws(()=>f.apply([NaN,1]));f.dispose();f.dispose();assert.throws(()=>f.apply([6,7]));
});
test('live boundary changes reuse ILU; coefficient changes rebuild it',()=>{
 const solve=createHeatSolver(api),config={solver:'gmres',size:8,contrast:1,diffusivity:.2,speed:4,angle:30,jacobi:false,ilu:true};
 const a=solve(config),b=solve({...config,hot:50}),c=solve({...config,hot:50,restart:5}),d=solve({...config,speed:1});
 assert.ok(a.converged&&b.converged&&c.converged&&d.converged);
 assert.equal(a.reusedFactor,false);assert.equal(b.reusedFactor,true);assert.equal(c.reusedFactor,true);assert.equal(d.reusedFactor,false);
 b.x.forEach((v,i)=>assert.ok(Math.abs(v-a.x[i]/2)<1e-8));assert.equal(b.setupMilliseconds,0);
 const zero=solve({...config,hot:0});assert.equal(zero.iterations,0);assert.ok(zero.x.every(v=>v===0));solve.dispose();
});
test('preconditioner compatibility and runtime overflow remain explicit',()=>{
 const a=new api.CSRMatrix(1,1,[0,1],[0],[1]),b=new api.CSRMatrix(2,2,[0,1,2],[0,1],[1,1]),f=new api.ILU0(a);
 try{assert.throws(()=>b.gmres([1,1],{preconditioner:f}));assert.throws(()=>a.gmres([1],{jacobi:true,preconditioner:f}));}finally{f.dispose();a.dispose();b.dispose();}
});

test('browser gate cannot accept false convergence on a tiny RHS',()=>{
 class Broken extends api.CSRMatrix {gmres(b,options){const norm=Math.hypot(...b);return {x:Array(b.length).fill(0),converged:true,iterations:0,reason:'converged',residuals:[norm],estimatedResiduals:[norm],iterates:[Array(b.length).fill(0)],restarts:[]};}}
 const tiny={name:'tiny false convergence',a:{rows:1,cols:1,offsets:[0,1],indices:[0],values:[1]},b:[1e-200],op:'gmres',expected:null,reason:'converged',options:{restart:1,rtol:1e-10,atol:0,limit:10,jacobi:0,capture:1}};
 assert.throws(()=>checkSparse({...api,CSRMatrix:Broken},[tiny]),/convergence/);
});
