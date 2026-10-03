import assert from 'node:assert/strict';
import {test} from 'node:test';
import {resolve} from 'node:path';
import {pathToFileURL} from 'node:url';
import {createMatrixAPI} from './matrix.mjs';
import {pca,fit,compress} from '../../applications/math.mjs';
const {Matrix}=await createMatrixAPI(process.env.LINEAR_A_WASM_MODULE?{moduleUrl:pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE))}:{});
const near=(a,b,t=1e-10)=>assert.ok(Math.abs(a-b)<t,`${a} != ${b}`);
test('PCA identifies an analytic line, is translation invariant, and handles degeneracy',()=>{
 const points=[[-2,-4],[-1,-2],[0,0],[1,2],[2,4]],r=pca(Matrix,points);
 near(r.values[0],12.5);near(r.values[1],0);near(r.explained,1);near(Math.abs(r.axes[0][1]/r.axes[0][0]),2);
 r.projected.flat().forEach((v,i)=>near(v,points.flat()[i]));
 const shifted=pca(Matrix,points.map(([x,y])=>[x+.5,y+.5]));shifted.values.forEach((v,i)=>near(v,r.values[i]));
 const zero=pca(Matrix,[[1,1],[1,1]]);assert.equal(zero.explained,null);assert.equal(zero.ambiguous,true);
 assert.equal(pca(Matrix,[[1,0],[-1,0],[0,1],[0,-1]]).ambiguous,true);
});
test('QR fit recovers known coefficients, gives orthogonal residuals, and rejects deficient designs',()=>{
 const points=[-4,-2,0,2,4].map(x=>[x,1+.25*x+.1*x*x]),r=fit(Matrix,points,2);
 r.coefficients.forEach((v,i)=>near(v,[1,1.25,2.5][i]));near(r.rmse,0);near(r.r2,1);
 const noisy=points.map(([x,y],i)=>[x,y+(i%2?.2:-.1)]),n=fit(Matrix,noisy,2);
 for(let j=0;j<3;j++)near(noisy.reduce((s,p,i)=>s+n.residuals[i]*(p[0]/5)**j,0),0);
 assert.throws(()=>fit(Matrix,[[1,1],[1,2],[1,3]],2));
 assert.equal(fit(Matrix,[[-1,2],[0,2],[1,2]],1).r2,null);
});
test('rank reconstruction has monotone error, exact full rank, and correct scalar counts',()=>{
 const pixels=[1,0,0,0,0,.5,0,0,0,0,.25,0,0,0,0,.125];let previous=Infinity;
 for(let rank=0;rank<=4;rank++){const r=compress(Matrix,pixels,4,rank);assert.ok(r.relativeError<=previous+1e-12);previous=r.relativeError;assert.equal(r.factorScalars,rank*9);if(rank===4)near(r.relativeError,0);}
 const rank1=compress(Matrix,Array(16).fill(.5),4,1);near(rank1.relativeError,0);
 const blank=compress(Matrix,Array(16).fill(0),4,2);assert.equal(blank.retained,null);assert.equal(blank.relativeError,0);
 const r=compress(Matrix,pixels,4,2);near(r.singularValues[0],1);near(r.singularValues[1],.5);
});
test('bounded inputs fail explicitly and repeated computations remain usable',()=>{
 assert.throws(()=>pca(Matrix,[[Infinity,0],[0,0]]));assert.throws(()=>fit(Matrix,[[0,0],[1,1]],3));
 assert.throws(()=>compress(Matrix,[1],1,1));assert.throws(()=>compress(Matrix,Array(16).fill(2),4,2));
 for(let i=0;i<50;i++){pca(Matrix,[[0,0],[1,2],[2,3]]);fit(Matrix,[[0,0],[1,1],[2,2]],1);compress(Matrix,Array(16).fill(.25),4,2);}
});
