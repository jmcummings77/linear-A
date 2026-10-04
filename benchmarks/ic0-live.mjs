/** Captured PCG iterates; playback is algorithm progress, not physical or CPU time. */
export function compareIC0(api,grid,{ordering='natural',breakdown=false}={}){
 if(!['natural','rcm','amd'].includes(ordering)||typeof breakdown!=='boolean')throw new RangeError('invalid IC0 settings');
 let a,q,ic,plan,chol;const now=()=>performance.now(),norm=v=>v.reduce((s,x)=>Math.hypot(s,x),0);
 try{
  const n=breakdown?4:grid.n,size=breakdown?2:grid.size;
  a=breakdown?new api.CSRMatrix(4,4,[0,3,6,9,12],[0,1,3,0,1,2,1,2,3,0,2,3],[1.5,1,1,1,1.5,1,1,1.5,-1,1,-1,1.5]):new api.CSRMatrix(n,n,grid.offsets,grid.indices,grid.values);
  const b=breakdown?a.matvec([1,2,3,4]):grid.b;
  const begin=now(),p=ordering==='natural'?Array.from({length:n},(_,i)=>i):ordering==='rcm'?a.reverseCuthillMcKee():a.approximateMinimumDegree();
  q=ordering==='natural'?a:a.permuteSymmetric(p);const rhs=api.CSRMatrix.permuteVector(p,b),orderingMs=now()-begin,restore=x=>api.CSRMatrix.permuteVector(p,x,true),runs={};
  let start=now(),icError;
  try{ic=new api.IC0(q);}catch(e){icError=e.message;}
  const icSetupMs=now()-start;
  start=now();plan=new api.SparseCholeskySymbolic(q);chol=plan.factorize(q);const cholSetupMs=now()-start;
  const residual=x=>{const ax=a.matvec(x);return b.map((v,i)=>v-ax[i]);},threshold=1e-9*norm(b);
  for(const method of ['cg','jacobi','ic0','cholesky']){
   if(method==='ic0'&&!ic){runs[method]={ok:false,error:icError,setupMs:icSetupMs};continue;}
   start=now();const r=method==='cholesky'?{x:chol.solve(rhs),iterates:[],reason:'converged',iterations:0}:q.conjugateGradient(rhs,{rtol:1e-9,maxIterations:1500,jacobi:method==='jacobi',preconditioner:method==='ic0'?ic:undefined,capture:true});const solveMs=now()-start;
   const x=restore(r.x),frames=method==='cholesky'?[Array(n).fill(0),x]:r.iterates.map(restore),residualFields=frames.map(residual),residuals=residualFields.map(norm);
   if(r.reason==='converged'&&norm(residual(x))>threshold*(1+1e-5)+1e-12)throw new Error('independent residual check failed');
   if(method!=='cholesky'&&r.residuals.some((v,i)=>Math.abs(v-residuals[i])>1e-10*Math.max(norm(b),1e-300)))throw new Error('captured residual check failed');
   runs[method]={ok:true,reason:r.reason,iterations:r.iterations,setupMs:method==='ic0'?icSetupMs:method==='cholesky'?cholSetupMs:0,solveMs,nnz:method==='ic0'?ic.nnz:method==='cholesky'?chol.nnz:0,frames,residualFields,residuals};
  }
  return {n,size,ordering,orderingMs,breakdown,runs};
 }finally{ic?.dispose?.();chol?.dispose?.();plan?.dispose?.();if(q!==a)q?.dispose?.();a?.dispose?.();}
}
