import {readFileSync} from 'node:fs';
export function sparseRun(args,CSRMatrix,Matrix,ILU0){
 if(args.length!==(args[0]==='gmres'?11:10))throw new Error('expected sparse OP ROWS COLS NNZ ITERATIONS RTOL ATOL LIMIT JACOBI CAPTURE');
 const [op]=args,[rows,cols,nnz,iters]=args.slice(1,5).map(Number),[rtol,atol]=args.slice(5,7).map(Number),[maxIterations,jacobi,capture]=args.slice(7,10).map(Number);
 const restart=op==='gmres'?Number(args[10]):30;
 if(!['spmv','dense','cg','gmres','ilu_setup','ilu_apply','rcm','permute','permutation_check','rcm_solve','ilu_solve'].includes(op)||[rows,cols,nnz,iters].some(x=>!Number.isSafeInteger(x)||x<0)||!(op==='gmres'?[0,1,2,3]:[0,1]).includes(jacobi)||![0,1].includes(capture))throw new Error('invalid sparse protocol');
 const raw=readFileSync(0,'utf8').trim(),data=raw?raw.split(/\s+/).map(Number):[],count=['cg','gmres'].includes(op)?rows:cols;
 if(data.length!==rows+1+2*nnz+count)throw new Error('incorrect sparse input count');
 const rp=data.slice(0,rows+1),ci=data.slice(rows+1,rows+1+nnz),v=data.slice(rows+1+nnz,rows+1+2*nnz),b=data.slice(rows+1+2*nnz);
 let a,dense,right,factor;
 try{
  a=new CSRMatrix(rows,cols,rp,ci,v);
  if(op==='dense'){const values=Array(rows*cols).fill(0);for(let i=0;i<rows;i++)for(let p=rp[i];p<rp[i+1];p++)values[i*cols+ci[p]]=v[p];dense=new Matrix(rows,cols,values);right=new Matrix(cols,1,b);}
  factor=op==='ilu_apply'||(op==='gmres'&&jacobi===2)?new ILU0(a):undefined;
  const compute=()=>{
   if(op==='rcm_solve'||op==='ilu_solve'){
    const p=op==='rcm_solve'?a.reverseCuthillMcKee():Array.from({length:rows},(_,i)=>i),q=op==='rcm_solve'?a.permuteSymmetric(p):a;let f;
    try{f=new ILU0(q);const rhs=op==='rcm_solve'?CSRMatrix.permuteVector(p,b):b,r=q.gmres(rhs,{restart:20,rtol,atol,maxIterations,preconditioner:f});if(!r.converged)throw new Error('ordering solve failed: '+r.reason);return [r.iterations,...(op==='rcm_solve'?CSRMatrix.permuteVector(p,r.x,true):r.x)];}
    finally{f?.dispose();if(q!==a)q.dispose();}
   }
   if(op==='rcm')return a.reverseCuthillMcKee();
   if(op==='permute'||op==='permutation_check'){
    const q=a.permuteSymmetric(b);
    try{if(op==='permute')return [...q.rowOffsets,...q.columnIndices,...q.values];
     const x=Array.from({length:rows},(_,i)=>i+1),y=CSRMatrix.permuteVector(b,x);
     return [...y,...CSRMatrix.permuteVector(b,y,true),...CSRMatrix.permuteVector(b,q.matvec(y),true)];
    }finally{q.dispose();}
   }
   if(op==='ilu_setup'){const f=new ILU0(a);try{return [f.size,f.nnz];}finally{f.dispose();}}
   if(op==='ilu_apply')return Array.from(factor.apply(b));
   if(op==='dense'){const out=dense.multiply(right);try{return Array.from(out.toArray?out.toArray():out.values);}finally{out.dispose?.();}}
   if(op==='spmv')return Array.from(a.matvec(b));
   const current=op==='gmres'&&jacobi===3?new ILU0(a):factor;
   let r;try{r=op==='gmres'?a.gmres(b,{restart,rtol,atol,maxIterations,jacobi:jacobi===1,capture:!!capture,preconditioner:current}):a.conjugateGradient(b,{rtol,atol,maxIterations,jacobi:!!jacobi,capture:!!capture});}finally{if(current!==factor)current?.dispose();}
   if(iters){if(!r.converged)throw new Error('benchmark CG did not converge: '+r.reason);return r.x;}
   return [['converged','iteration_limit','breakdown','nonfinite','stagnation'].indexOf(r.reason),r.iterations,r.residuals.length,...r.x,...r.residuals,...('estimatedResiduals' in r?[...r.estimatedResiduals,r.restarts.length,...r.restarts]:[]),...r.iterates.flat()];
  };
  if(!iters){const values=compute();return {rows:1,cols:values.length,values};}
  for(let i=0;i<3;i++)compute();let checksum=0;const start=process.hrtime.bigint();
  for(let i=0;i<iters;i++)for(const x of compute())checksum+=x;
  const elapsed_ns=Number(process.hrtime.bigint()-start);if(!Number.isFinite(checksum))throw new Error('nonfinite benchmark checksum');
  return {elapsed_ns,iterations:iters,checksum};
 }finally{factor?.dispose();a?.dispose?.();dense?.dispose?.();right?.dispose?.();}
}
