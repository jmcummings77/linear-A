/** Actual V-cycle callbacks and independently checked PCG residual histories. */
export function exploreMultigrid(api,{width=15,pattern='mixed'}={}){
 if(![3,7,15,31].includes(width)||!['mixed','smooth','rough','broad'].includes(pattern))throw new RangeError('invalid grid settings');
 const n=width*width,norm=v=>v.reduce((s,x)=>Math.hypot(s,x),0),truth=Array.from({length:n},(_,i)=>{
  const x=(i%width+1)/(width+1),y=(Math.floor(i/width)+1)/(width+1),smooth=Math.sin(Math.PI*x)*Math.sin(Math.PI*y),rough=.25*Math.sin((width-1)*Math.PI*x)*Math.sin((width-1)*Math.PI*y);
  return pattern==='broad'?Math.sin((i+1)*1.731)+.2*Math.cos((i+3)*.413):pattern==='smooth'?smooth:pattern==='rough'?rough:smooth+rough;
 });
 let mg,a,ic;
 try{
  mg=new api.GeometricMultigrid(width);a=mg.matrix;const b=a.matvec(truth),trace=mg.trace(b);
  for(const f of trace.frames){
   let level,q;
   try{level=new api.GeometricMultigrid(f.width);q=level.matrix;const ax=q.matvec(f.x),r=f.b.map((v,i)=>v-ax[i]);if(norm(r.map((v,i)=>v-f.residual[i]))>1e-11*Math.max(1,norm(f.b)))throw new Error('V-cycle frame residual check failed');}
   finally{q?.dispose?.();level?.dispose?.();}
  }
  // Independently follow the error equation through each descent. Coarse RHS is
  // P^T r; its exact solution need not equal restricted fine-grid truth.
  const frames=trace.frames.map(f=>({...f,error:f.level===0?truth.map((v,i)=>v-f.x[i]):null}));
  ic=new api.IC0(a);const runs={};
  for(const method of ['cg','jacobi','ic0','multigrid']){
   const start=performance.now(),r=a.conjugateGradient(b,{rtol:1e-9,maxIterations:1000,capture:true,jacobi:method==='jacobi',preconditioner:method==='ic0'?ic:method==='multigrid'?mg:undefined}),solveMs=performance.now()-start;
   if(!r.converged)throw new Error(method+' failed: '+r.reason);
   r.iterates.forEach((x,k)=>{const ax=a.matvec(x),actual=norm(b.map((v,i)=>v-ax[i]));if(Math.abs(actual-r.residuals[k])>1e-11*Math.max(norm(b),1))throw new Error('CG frame residual check failed');});
   const ax=a.matvec(r.x);if(norm(b.map((v,i)=>v-ax[i]))>1e-9*norm(b)*(1+1e-5)+1e-13)throw new Error('CG final residual check failed');
   runs[method]={iterations:r.iterations,residuals:r.residuals,solveMs};
  }
  return {width,n,pattern,truth,b,frames,runs,levels:mg.levels};
 }finally{ic?.dispose?.();a?.dispose?.();mg?.dispose?.();}
}
