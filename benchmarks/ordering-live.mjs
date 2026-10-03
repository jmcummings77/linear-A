/** Live comparison in original coordinates; chart motion is a display interpolation. */
export function compareOrdering(api,grid,{scrambled=false,preconditioner='ilu'}={}){
 if(!['ilu','jacobi'].includes(preconditioner))throw new RangeError('invalid preconditioner');
 const {CSRMatrix,ILU0}=api,n=grid.n;
 let a=new CSRMatrix(n,n,grid.offsets,grid.indices,grid.values),b=grid.b.slice();
 try{
  if(scrambled){const p=Array.from({length:n},(_,i)=>i);let state=2026;for(let i=n-1;i>0;i--){state=(Math.imul(state,1664525)+1013904223)>>>0;const j=state%(i+1);[p[i],p[j]]=[p[j],p[i]];}const q=a.permuteSymmetric(p);b=CSRMatrix.permuteVector(p,b);a.dispose?.();a=q;}
  const start=performance.now(),p=a.reverseCuthillMcKee(),orderingMilliseconds=performance.now()-start;
  const points=[];const offsets=a.rowOffsets,indices=a.columnIndices,values=a.values;
  for(let i=0;i<n;i++)for(let k=offsets[i];k<offsets[i+1];k++)points.push([i,indices[k],values[k]]);
  const inverse=Array(n);p.forEach((old,i)=>inverse[old]=i);
  const width=points.reduce((v,[i,j])=>Math.max(v,Math.abs(i-j)),0),newWidth=points.reduce((v,[i,j])=>Math.max(v,Math.abs(inverse[i]-inverse[j])),0);
  const solve=reorder=>{
   let q=a,f;const totalStart=performance.now();let rhs=b;
   try{
    if(reorder){q=a.permuteSymmetric(p);rhs=CSRMatrix.permuteVector(p,b);}
    const setupStart=performance.now();if(preconditioner==='ilu')f=new ILU0(q);const setupMilliseconds=performance.now()-setupStart;
    const solveStart=performance.now(),r=q.gmres(rhs,{restart:20,rtol:1e-8,atol:0,maxIterations:800,preconditioner:f,jacobi:preconditioner==='jacobi'}),solveMilliseconds=performance.now()-solveStart;
    const x=reorder?CSRMatrix.permuteVector(p,r.x,true):r.x;
    const totalMilliseconds=performance.now()-totalStart+(reorder?orderingMilliseconds:0);
    const ax=a.matvec(x),residual=b.reduce((s,v,i)=>Math.hypot(s,v-ax[i]),0),threshold=1e-8*b.reduce((s,v)=>Math.hypot(s,v),0);
    if(r.converged&&!(residual<=threshold*(1+1e-5)+1e-300))throw new Error('Restored solution failed original-system residual check');
    return {iterations:r.iterations,converged:r.converged,reason:r.reason,residuals:r.residuals,residual,threshold,setupMilliseconds,solveMilliseconds,totalMilliseconds};
   }catch(e){return {converged:false,reason:e.message||String(e),residuals:[]};}
   finally{f?.dispose?.();if(q!==a)q.dispose?.();}
  };
  return {n,points,p,inverse,width,newWidth,orderingMilliseconds,preconditioner,natural:solve(false),rcm:solve(true)};
 }finally{a.dispose?.();}
}
