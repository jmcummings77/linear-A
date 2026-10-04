/** Actual symbolic fill history and independently checked numerical factors. */
export function compareAMD(api,grid,{scrambled=false}={}){
 const {CSRMatrix,SparseCholeskySymbolic}=api,n=grid.n;
 let a=new CSRMatrix(n,n,grid.offsets,grid.indices,grid.values),b=grid.b.slice();
 const norm=v=>v.reduce((s,x)=>Math.hypot(s,x),0);
 try{
  if(scrambled){const p=Array.from({length:n},(_,i)=>i);let state=2026;for(let i=n-1;i>0;i--){state=(Math.imul(state,1664525)+1013904223)>>>0;const j=state%(i+1);[p[i],p[j]]=[p[j],p[i]];}const q=a.permuteSymmetric(p);b=CSRMatrix.permuteVector(p,b);a.dispose?.();a=q;}
  const solve=ordering=>{
   const reorder=ordering!=='natural';let q=a,plan,f,l;const start=performance.now();let rhs=b,p;
   try{
    if(reorder){p=ordering==='amd'?a.approximateMinimumDegree():a.reverseCuthillMcKee();q=a.permuteSymmetric(p);rhs=CSRMatrix.permuteVector(p,b);}
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
    const order=p||Array.from({length:n},(_,i)=>i),trace=eliminationTrace(a,order,points);
    return {ok:true,points,order,trace,nnz:f.nnz,fillCount:plan.fillCount,logicalFactorBytes:16*f.nnz+8*(n+1),residual,reconstructionError,orderingMilliseconds,symbolicMilliseconds,factorMilliseconds,solveMilliseconds,totalMilliseconds};
   }catch(e){return {ok:false,error:e.message||String(e),points:[]};}
   finally{l?.dispose?.();f?.dispose?.();plan?.dispose?.();if(q!==a)q.dispose?.();}
  };
  return {n,natural:solve('natural'),rcm:solve('rcm'),amd:solve('amd')};
 }finally{a.dispose?.();}
}

/** Explicit graph replay verifies the fill history independently of the WASM analysis. */
export function eliminationTrace(a,order,points){
 const n=a.rows,rp=a.rowOffsets,ci=a.columnIndices,g=Array.from({length:n},()=>new Set()),edges=[];
 for(let i=0;i<n;i++)for(let k=rp[i];k<rp[i+1];k++){const j=ci[k];if(i!==j){g[i].add(j);g[j].add(i);}}
 for(let i=0;i<n;i++)for(const j of g[i])if(i<j)edges.push([i,j,-1]);
 const active=new Set(order),frames=[];
 for(let step=0;step<n;step++){
  const pivot=order[step],neighbors=[...g[pivot]].filter(i=>active.has(i)).sort((a,b)=>a-b),fill=[];
  for(let u=0;u<neighbors.length;u++)for(let v=u+1;v<neighbors.length;v++){
   const i=neighbors[u],j=neighbors[v];if(!g[i].has(j)){g[i].add(j);g[j].add(i);fill.push([i,j]);edges.push([i,j,step]);}
  }
  const expected=points.filter(p=>p[2]===step).map(([i,j])=>[order[i],order[j]].sort((a,b)=>a-b).join(',')).sort();
  if(JSON.stringify(expected)!==JSON.stringify(fill.map(e=>e.join(',')).sort()))throw new Error('Elimination graph disagrees with symbolic fill');
  frames.push({pivot,neighbors,fill});active.delete(pivot);
 }
 return {edges,frames};
}

/** Deterministic small SPD examples: browser seed is intentionally independent of recorded data. */
export function graphProblem(size,kind){
 if(!Number.isInteger(size)||size<2||size>16||!['tree','irregular'].includes(kind))throw new RangeError('invalid graph problem');
 const n=size*size,g=Array.from({length:n},()=>new Set());let state=2026;
 const random=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state;};
 const edge=(i,j)=>{if(i!==j){g[i].add(j);g[j].add(i);}};
 for(let i=1;i<n;i++)edge(i,kind==='tree'?Math.floor((i-1)/2):random()%i);
 if(kind==='irregular')for(let k=0;k<2*n;k++){const i=random()%n,j=Math.floor(random()/65536)%n;edge(i,j);}
 const offsets=[0],indices=[],values=[],x=Array.from({length:n},(_,i)=>1+(i%7)/10),b=[];
 for(let i=0;i<n;i++){
  for(const j of [...g[i],i].sort((a,b)=>a-b)){indices.push(j);values.push(j===i?g[i].size+1:-1);}
  offsets.push(indices.length);b.push((g[i].size+1)*x[i]-[...g[i]].reduce((sum,j)=>sum+x[j],0));
 }
 return {n,offsets,indices,values,b};
}
