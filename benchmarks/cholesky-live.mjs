/** Actual symbolic fill history and independently checked numerical factors. */
export function compareCholesky(api,grid,{scrambled=false}={}){
 const {CSRMatrix,SparseCholeskySymbolic}=api,n=grid.n;
 let a=new CSRMatrix(n,n,grid.offsets,grid.indices,grid.values),b=grid.b.slice();
 const norm=v=>v.reduce((s,x)=>Math.hypot(s,x),0);
 try{
  if(scrambled){const p=Array.from({length:n},(_,i)=>i);let state=2026;for(let i=n-1;i>0;i--){state=(Math.imul(state,1664525)+1013904223)>>>0;const j=state%(i+1);[p[i],p[j]]=[p[j],p[i]];}const q=a.permuteSymmetric(p);b=CSRMatrix.permuteVector(p,b);a.dispose?.();a=q;}
  const solve=reorder=>{
   let q=a,plan,f,l;const start=performance.now();let rhs=b,p;
   try{
    if(reorder){p=a.reverseCuthillMcKee();q=a.permuteSymmetric(p);rhs=CSRMatrix.permuteVector(p,b);}
    const orderingMilliseconds=performance.now()-start,t1=performance.now();plan=new SparseCholeskySymbolic(q);const symbolicMilliseconds=performance.now()-t1;
    const t2=performance.now();f=plan.factorize(q);const factorMilliseconds=performance.now()-t2;
    const t3=performance.now(),y=f.solve(rhs);const solveMilliseconds=performance.now()-t3;
    const x=reorder?CSRMatrix.permuteVector(p,y,true):y,totalMilliseconds=performance.now()-start;
    const ax=a.matvec(x),residual=norm(b.map((v,i)=>v-ax[i]))/Math.max(norm(b),1e-300);
    l=f.lower;const offsets=l.rowOffsets,indices=l.columnIndices,values=l.values,steps=plan.fillSteps;
    // Sparse row intersections reconstruct every entry, including structural zeros.
    const qr=q.rowOffsets,qc=q.columnIndices,qv=q.values;let error=0;
    for(let i=0;i<n;i++)for(let j=0;j<=i;j++){
      let u=offsets[i],w=offsets[j],sum=0;
      while(u<offsets[i+1]&&w<offsets[j+1]){if(indices[u]===indices[w]){sum+=values[u]*values[w];u++;w++;}else if(indices[u]<indices[w])u++;else w++;}
      let value=0;for(let k=qr[i];k<qr[i+1];k++)if(qc[k]===j){value=qv[k];break;}
      error=Math.hypot(error,(sum-value)*(i===j?1:Math.SQRT2));
    }
    const reconstructionError=error/Math.max(norm(qv),1e-300);
    if(residual>1e-10||reconstructionError>1e-10||!Number.isFinite(residual+reconstructionError))throw new Error('Cholesky reconstruction or residual check failed');
    const points=[];for(let i=0;i<n;i++)for(let k=offsets[i];k<offsets[i+1];k++)points.push([i,indices[k],steps[k],values[k]]);
    return {ok:true,points,nnz:f.nnz,fillCount:plan.fillCount,logicalFactorBytes:16*f.nnz+8*(n+1),residual,reconstructionError,orderingMilliseconds,symbolicMilliseconds,factorMilliseconds,solveMilliseconds,totalMilliseconds};
   }catch(e){return {ok:false,error:e.message||String(e),points:[]};}
   finally{l?.dispose?.();f?.dispose?.();plan?.dispose?.();if(q!==a)q.dispose?.();}
  };
  return {n,natural:solve(false),rcm:solve(true)};
 }finally{a.dispose?.();}
}
