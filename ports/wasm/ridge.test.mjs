import assert from 'node:assert/strict';
import {test} from 'node:test';
import {execFileSync} from 'node:child_process';
import {createMatrixAPI} from './matrix.mjs';
import {runChecks} from '../../benchmarks/live-worker.mjs';
import {computeAccuracy} from '../../benchmarks/accuracy-worker.mjs';
const {Matrix}=await createMatrixAPI();
test('ridge shared analytic fixtures execute through WASM',async()=>{
 const cases=JSON.parse(execFileSync('python3',['-c','import sys,json;sys.path.insert(0,"benchmarks");from ridge_reference import ridge_fixtures;print(json.dumps(ridge_fixtures()))'],{encoding:'utf8'}));
 const result=await runChecks(Matrix,cases);assert.equal(result.passed,result.total,JSON.stringify(result));
});
test('ridge rejects invalid options and disposed inputs without damaging inputs',()=>{
 const a=new Matrix(1,1,[2]),b=new Matrix(1,1,[3]);
 try{
  for(const lambda of [-1,Infinity,NaN,undefined])assert.throws(()=>a.solveRidge(b,lambda));
  assert.throws(()=>a.solveRidge(null,1));
  const x=a.solveRidge(b,1);assert.ok(Math.abs(x.get(0,0)-1.2)<1e-14);x.dispose();
  assert.equal(a.get(0,0),2);assert.equal(b.get(0,0),3);
  b.dispose();assert.throws(()=>a.solveRidge(b,1));
 }finally{a.dispose();b.dispose();}
});
test('live ridge comparison shows smooth tradeoff and common noisy inputs',()=>{
 const config={a:[1,1,1,1.000001],b:[2,2.000001],delta:1e-7,target:5,noise:1e-4,algorithm:'ridge',cutoff:1e-5,lambda:1e-6};
 const r=computeAccuracy(Matrix,config),[ordinary,truncated,ridge]=r.compare;
 assert.equal(r.b[1],config.b[1]+config.noise);
 assert.equal(r.compare.length,3);assert.equal(r.curve.length,41);
 assert.ok(ordinary.solution_norm>ridge.solution_norm*50);
 assert.ok(ridge.sensitivity<ordinary.sensitivity/100);
 assert.ok(truncated.rank===1 && ridge.rank===2);
 assert.ok(ridge.residual_norm>ordinary.residual_norm);
 for(let i=1;i<r.curve.length;i++){
  assert.ok(r.curve[i].solution_norm<=r.curve[i-1].solution_norm*(1+1e-10));
  assert.ok(r.curve[i].residual_norm+1e-10>=r.curve[i-1].residual_norm);
 }
 const zero=computeAccuracy(Matrix,{...config,lambda:0});
 zero.base.x.forEach((x,i)=>assert.ok(Math.abs(x-zero.compare[0].x[i])<1e-12));
 for(const update of [{lambda:-1},{lambda:NaN},{noise:Infinity},{noise:2}])assert.throws(()=>computeAccuracy(Matrix,{...config,...update}));
});

test('browser ridge gate rejects a vanished tiny nonzero solution',async()=>{
 class BrokenMatrix {
  constructor(rows,cols,values){this.rows=rows;this.cols=cols;this.values=values;}
  solveRidge(){return new BrokenMatrix(1,1,[0]);}
  toArray(){return this.values;}
  dispose(){}
 }
 const fixture={name:'tiny ridge',op:'solve_ridge',a:{rows:1,cols:1,values:[1]},b:{rows:1,cols:1,values:[1e-300]},scalar:1,expected:{rows:1,cols:1,values:[5e-301]},invalid:false};
 const result=await runChecks(BrokenMatrix,[fixture]);assert.equal(result.passed,0);
 assert.match(result.checks[0].error,/ridge analytic/);
});
