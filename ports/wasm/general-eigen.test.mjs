import test from 'node:test';
import assert from 'node:assert/strict';
import { resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
import { createMatrixAPI } from './matrix.mjs';

const options = process.env.LINEAR_A_WASM_MODULE
  ? { moduleUrl: pathToFileURL(resolve(process.env.LINEAR_A_WASM_MODULE)) } : {};
const { Matrix } = await createMatrixAPI(options);
function dispose(result) { result?.vectorsReal.dispose(); result?.vectorsImag.dispose(); }
function invariants(a, result) {
  const n=a.rows, input=a.toArray(), real=result.vectorsReal.toArray(), imag=result.vectorsImag.toArray();
  assert.equal(result.valuesReal.length,n); assert.equal(result.valuesImag.length,n);
  for (const vectors of [result.vectorsReal,result.vectorsImag]) assert.deepEqual([vectors.rows,vectors.cols],[n,n]);
  const scale=Math.max(...input.map(Math.abs),0)||1;
  for (let col=0;col<n;col++) {
    const dr=result.valuesReal[col], di=result.valuesImag[col];
    assert.ok(Number.isFinite(dr)&&Number.isFinite(di));
    if(col) assert.ok(result.valuesReal[col-1]<dr || (result.valuesReal[col-1]===dr && result.valuesImag[col-1]<=di));
    let norm=0;
    for(let row=0;row<n;row++) {
      let ar=0,ai=0;
      for(let k=0;k<n;k++) {ar+=(input[row*n+k]/scale)*real[k*n+col];ai+=(input[row*n+k]/scale)*imag[k*n+col];}
      const index=row*n+col,br=real[index]*(dr/scale)-imag[index]*(di/scale),bi=imag[index]*(dr/scale)+real[index]*(di/scale);
      assert.ok(Math.hypot(ar-br,ai-bi)<5e-11*n,`complex eigenpair residual ${row},${col}`);
      norm=Math.hypot(norm,real[index],imag[index]);
    }
    assert.ok(Math.abs(norm-1)<1e-12*n);
  }
}
function decompose(n,input,verify) {
  const a=new Matrix(n,n,input);let result;
  try {result=a.eigenGeneral();invariants(a,result);verify?.(result,a);assert.deepEqual([...a.toArray()],input);}
  finally{dispose(result);a.dispose();}
}
test('general solver finds complex conjugates at ordinary, huge, tiny, and subnormal scales',()=>{
  for(const scale of [1,1e300,1e-300,Number.MIN_VALUE]) {
    decompose(2,[0,-scale,scale,0],result=>{
      assert.deepEqual([...result.valuesReal],[0,0]);
      assert.ok(Math.abs(result.valuesImag[0]/scale+1)<1e-14);
      assert.ok(Math.abs(result.valuesImag[1]/scale-1)<1e-14);
    });
  }
  decompose(2,[0,1e300,-1e-300,0],result=>{
    assert.ok(Math.abs(result.valuesImag[0]+1)<1e-14&&Math.abs(result.valuesImag[1]-1)<1e-14);
    const real=result.vectorsReal.toArray(),imag=result.vectorsImag.toArray();
    // This explicit unscaled check catches a zeroed tiny vector component,
    // which a global matrix-norm residual alone cannot detect.
    for(let col=0;col<2;col++) {
      const br=-imag[col]*result.valuesImag[col],bi=real[col]*result.valuesImag[col];
      assert.ok(Math.abs(1e300*real[2+col]-br)<1e-13);
      assert.ok(Math.abs(1e300*imag[2+col]-bi)<1e-13);
    }
  });
});
test('general solver preserves diagonal spectra and handles empty, zero, and defective matrices',()=>{
  decompose(0,[],result=>assert.equal(result.valuesReal.length,0));
  decompose(3,Array(9).fill(0),result=>assert.deepEqual([...result.valuesReal],[0,0,0]));
  decompose(3,[Number.MAX_VALUE,0,0,0,Number.MIN_VALUE,0,0,0,-Number.MAX_VALUE],result=>
    assert.deepEqual([...result.valuesReal],[-Number.MAX_VALUE,Number.MIN_VALUE,Number.MAX_VALUE]));
  for(const input of [[1e300,1,0,1e-300],[1e-300,0,1,1e300]])
    decompose(2,input,result=>assert.deepEqual([...result.valuesReal],[1e-300,1e300]));
  decompose(3,[1,1,0,0,1,1,0,0,1],result=>assert.deepEqual([...result.valuesReal],[1,1,1]));
  decompose(3,[0,1,0,0,0,1,0,0,0],result=>assert.deepEqual([...result.valuesReal],[0,0,0]));
});
test('general eigenpairs satisfy complex residuals for random non-symmetric matrices',()=>{
  let state=43;
  for(let trial=0;trial<40;trial++) {
    const input=Array.from({length:49},()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state%101-50;});
    decompose(7,input);
  }
});
test('general result arrays and both owned vector matrices are independent',()=>{
  const a=new Matrix(2,2,[0,-1,1,0]);let first,second;
  try {
    first=a.eigenGeneral();second=a.eigenGeneral();
    first.valuesReal[0]=99; first.valuesImag[0]=99; first.vectorsReal.set(0,0,99); first.vectorsImag.set(0,0,99);
    assert.deepEqual([...second.valuesReal],[0,0]); assert.deepEqual([...second.valuesImag],[-1,1]);
    assert.notEqual(second.vectorsReal.get(0,0),99); assert.notEqual(second.vectorsImag.get(0,0),99);
    a.dispose(); const original=new Matrix(2,2,[0,-1,1,0]);
    try {invariants(original,second);} finally {original.dispose();}
    assert.throws(()=>a.eigenGeneral(),/disposed/);
  } finally {dispose(first);dispose(second);a.dispose();}
});
test('general options, shape, nonconvergence, and nonfinite output fail and recover cleanly',()=>{
  const dense=new Matrix(3,3,[1,2,3,4,5,6,7,8,10]),rectangle=new Matrix(2,3),overflow=new Matrix(2,2,Array(4).fill(Number.MAX_VALUE));
  try {
    for(const maxIterations of [0,-1,1.5,100001,NaN,Infinity,'1000',null,true]) assert.throws(()=>dense.eigenGeneral({maxIterations}),/maxIterations/);
    assert.throws(()=>dense.eigenGeneral({maxIterations:1}),/converge/);
    assert.throws(()=>rectangle.eigenGeneral(),/dimensions/);
    assert.throws(()=>overflow.eigenGeneral(),/nonfinite/);
    for(let i=0;i<10;i++) {assert.throws(()=>dense.eigenGeneral({maxIterations:1}),/converge/);const result=dense.eigenGeneral();invariants(dense,result);dispose(result);}
  } finally {dense.dispose();rectangle.dispose();overflow.dispose();}
});
