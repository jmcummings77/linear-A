/** Educational applications; all decompositions use the public WASM matrix API. */
function points(input) {
  if(!Array.isArray(input)||input.length<2||input.length>80||input.some(p=>!Array.isArray(p)||p.length!==2||p.some(v=>!Number.isFinite(v)||Math.abs(v)>5)))throw new RangeError('Use 2–80 points with coordinates in [−5, 5].');
  return input;
}
export function pca(Matrix,input){
  const p=points(input),mean=[0,1].map(j=>p.reduce((s,v)=>s+v[j],0)/p.length);
  const centered=p.map(v=>v.map((x,j)=>x-mean[j]));
  const covariance=[0,1].flatMap(i=>[0,1].map(j=>centered.reduce((s,v)=>s+v[i]*v[j],0)/(p.length-1)));
  let a,e;
  try{a=new Matrix(2,2,covariance);e=a.eigenSymmetric();const vectors=Array.from(e.vectors.toArray()),values=Array.from(e.values).reverse().map(v=>Math.max(0,v)),axes=[[vectors[1],vectors[3]],[vectors[0],vectors[2]]];
    const total=values[0]+values[1];return {mean,covariance,values,axes,explained:total?values[0]/total:null,ambiguous:total===0||Math.abs(values[0]-values[1])<=total*1e-10,
      projected:centered.map(v=>{const s=v[0]*axes[0][0]+v[1]*axes[0][1];return [mean[0]+s*axes[0][0],mean[1]+s*axes[0][1]];})};
  }finally{e?.vectors.dispose();a?.dispose();}
}
export function fit(Matrix,input,degree){
  const p=points(input);if(!Number.isInteger(degree)||degree<1||degree>3||p.length<degree+1)throw new RangeError('A degree d fit needs at least d+1 points (degrees 1–3).');
  // Normalize x by the fixed plot extent; coefficients refer to t=x/5.
  const design=p.flatMap(([x])=>Array.from({length:degree+1},(_,j)=>(x/5)**j));
  let a,b,x;try{a=new Matrix(p.length,degree+1,design);b=new Matrix(p.length,1,p.map(v=>v[1]));x=a.leastSquares(b);
    const coefficients=Array.from(x.toArray()),predict=v=>coefficients.reduce((s,c,j)=>s+c*(v/5)**j,0),predicted=p.map(v=>predict(v[0])),residuals=p.map((v,i)=>v[1]-predicted[i]),sse=residuals.reduce((s,v)=>s+v*v,0);
    const mean=p.reduce((s,v)=>s+v[1],0)/p.length,tss=p.reduce((s,v)=>s+(v[1]-mean)**2,0);
    return {coefficients,predicted,residuals,rmse:Math.sqrt(sse/p.length),r2:tss?sse/tss===0?1:1-sse/tss:null,
      curve:Array.from({length:201},(_,i)=>{const x=-5+i/20;return [x,predict(x)];})};
  }finally{x?.dispose();b?.dispose();a?.dispose();}
}
export function compress(Matrix,pixels,size,rank){
  if(!Number.isInteger(size)||size<2||size>32||!Array.isArray(pixels)||pixels.length!==size*size||pixels.some(v=>!Number.isFinite(v)||v<0||v>1)||!Number.isInteger(rank)||rank<0||rank>size)throw new RangeError('Use a 2–32 square grayscale image and a valid rank.');
  let a,e;try{
    a=new Matrix(size,size,pixels);e=a.svd();
    const u=e.u.toArray(),vt=e.vt.toArray(),values=Array.from(e.values),reconstructed=Array(size*size).fill(0);
    for(let column=0;column<rank;column++)for(let i=0;i<size;i++)for(let j=0;j<size;j++)
      reconstructed[i*size+j]+=u[i*size+column]*values[column]*vt[column*size+j];
    const energy=pixels.reduce((s,v)=>s+v*v,0),error=pixels.reduce((s,v,i)=>s+(v-reconstructed[i])**2,0);
    return {reconstructed,singularValues:values,relativeError:energy?Math.sqrt(error/energy):0,
      rmse:Math.sqrt(error/pixels.length),retained:energy?Math.max(0,Math.min(1,1-error/energy)):null,
      originalScalars:size*size,factorScalars:rank*(2*size+1)};
  }finally{e?.u.dispose();e?.vt.dispose();a?.dispose();}
}
