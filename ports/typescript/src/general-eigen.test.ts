import test from 'node:test';
import assert from 'node:assert/strict';
import {Matrix} from './matrix.js';
import {benchmark} from './runner.js';

function verify(a:Matrix) {
  const before=a.values, n=a.rows, r=a.eigenGeneral(), scale=Math.max(...before.map(Math.abs),1e-300);
  for(let col=0;col<n;col++) {
    let norm=0;
    for(let row=0;row<n;row++) {
      const re=r.vectorsReal.get(row,col),im=r.vectorsImag.get(row,col);
      norm+=re*re+im*im;
      let ar=0,ai=0;
      for(let k=0;k<n;k++) {ar+=a.get(row,k)/scale*r.vectorsReal.get(k,col);ai+=a.get(row,k)/scale*r.vectorsImag.get(k,col);}
      assert.ok(Math.hypot(ar-(r.valuesReal[col]/scale*re-r.valuesImag[col]/scale*im),
        ai-(r.valuesReal[col]/scale*im+r.valuesImag[col]/scale*re))<1e-10*n);
    }
    assert.ok(Math.abs(norm-1)<1e-12);
  }
  const sumReal=r.valuesReal.reduce((sum,x)=>sum+x,0),sumImag=r.valuesImag.reduce((sum,x)=>sum+x,0);
  assert.ok(Math.hypot((sumReal-a.trace())/scale,sumImag/scale)<1e-10*Math.max(n,1));
  assert.deepEqual(a.values,before);
  return r;
}
test('general real inputs return real, complex, and defective eigenpairs',()=>{
  const rotation=verify(new Matrix(2,2,[0,-1,1,0]));
  assert.deepEqual([...rotation.valuesImag],[-1,1]);
  assert.deepEqual([...verify(new Matrix(2,2,[2,1,0,2])).valuesReal],[2,2]);
  assert.deepEqual([...verify(new Matrix(3,3,[0,1,0,0,0,1,0,0,0])).valuesReal],[0,0,0]);
  verify(new Matrix(3,3,[1,-2,3,2,1,4,0,0,-3]));
  verify(new Matrix(3,3,[3,0,0,2,-1,0,1,4,2]));
  assert.equal(new Matrix().eigenGeneral().valuesReal.length,0);
});
test('balancing retains mixed scales and exact extreme diagonals',()=>{
  const r=verify(new Matrix(2,2,[0,1e300,-1e-300,0]));
  assert.ok(Math.abs(r.valuesImag[0]+1)<1e-14);
  assert.ok(Math.abs(r.valuesImag[1]-1)<1e-14);
  assert.deepEqual([...new Matrix(2,2,[1e300,0,0,1e-300]).eigenGeneral().valuesReal],[1e-300,1e300]);
  for(const scale of [1e-200,1e200]) verify(new Matrix(2,2,[scale,-2*scale,2*scale,scale]));
});
test('random nonsymmetric residuals, ownership, and explicit failure modes',()=>{
  let state=315;
  for(let n=1;n<10;n++) for(let j=0;j<8;j++) {
    const a=new Matrix(n,n,Array.from({length:n*n},()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return (state%101)-50;}));
    const r=verify(a);r.vectorsReal.set(0,0,900);assert.ok(Math.abs(a.eigenGeneral().vectorsReal.get(0,0))<=1);
  }
  const a=new Matrix(4,4,[4,1,2,3,0,5,1,2,2,1,6,1,3,0,1,7]);
  for(const max of [0,-1,1.5,NaN,Infinity,100001]) assert.throws(()=>a.eigenGeneral(max),RangeError);
  assert.throws(()=>a.eigenGeneral(1),/converge/);
  assert.throws(()=>new Matrix(2,3).eigenGeneral(),/square/);
  assert.throws(()=>new Matrix(2,2,[1e308,1e308,1e308,1e308]).eigenGeneral(),RangeError);
});
test('runner checksum consumes complex values and normalized eigenvectors',()=>{
  const result=benchmark(['eigen_general','3','2','17']) as {checksum:number};
  assert.ok(Math.abs(result.checksum-2*(1.5+3+3.375+3))<1e-12);
});
