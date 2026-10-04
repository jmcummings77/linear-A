import {SparseCholeskySymbolic,IC0} from './cholesky.js';
import {GeometricMultigrid} from './multigrid.js';
import {ILU0} from "./ilu.js";
import type {GMRESResult} from "./gmres.js";
import {CSRMatrix} from "./sparse.js";
import {Matrix} from "./matrix.js";
import {readFileSync} from 'node:fs';
export function sparseRun(args:string[]){
 if(args.length!==(args[0]==='gmres'?11:10))throw new Error('expected sparse OP ROWS COLS NNZ ITERATIONS RTOL ATOL LIMIT JACOBI CAPTURE');
 const [op]=args,[rows,cols,nnz,iters]=args.slice(1,5).map(Number),[rtol,atol]=args.slice(5,7).map(Number),[maxIterations,jacobi,capture]=args.slice(7,10).map(Number);
 const restart=op==='gmres'?Number(args[10]):30;
 if(!['mg_setup','mg_matrix','mg_apply','spmv','dense','cg','gmres','ic0_factor','ic0_apply','ilu_setup','ilu_apply','rcm','amd','permute','permutation_check','rcm_solve','ilu_solve','chol_symbolic','chol_factor','chol_solve','chol_total','chol_rcm_total','chol_amd_total'].includes(op)||[rows,cols,nnz,iters].some(x=>!Number.isSafeInteger(x)||x<0)||!(op==='cg'?[0,1,2,3,4,5]:op==='gmres'?[0,1,2,3]:[0,1]).includes(jacobi)||![0,1].includes(capture))throw new Error('invalid sparse protocol');
 const raw=readFileSync(0,'utf8').trim(),data=raw?raw.split(/\s+/).map(Number):[],count=['cg','gmres'].includes(op)?rows:cols;
 if(data.length!==rows+1+2*nnz+count)throw new Error('incorrect sparse input count');
 const rp=data.slice(0,rows+1),ci=data.slice(rows+1,rows+1+nnz),v=data.slice(rows+1+nnz,rows+1+2*nnz),b=data.slice(rows+1+2*nnz);
 let a:CSRMatrix,dense:Matrix|undefined,right:Matrix|undefined;
 try{
  a=new CSRMatrix(rows,cols,rp,ci,v);
  if(op==='dense'){const values=Array(rows*cols).fill(0);for(let i=0;i<rows;i++)for(let p=rp[i];p<rp[i+1];p++)values[i*cols+ci[p]]=v[p];dense=new Matrix(rows,cols,values);right=new Matrix(cols,1,b);}
  const newMG=()=>{const w=Math.sqrt(rows);if(rows!==cols||!Number.isInteger(w))throw new Error('multigrid requires a square grid');return new GeometricMultigrid(w);};
  const mg=['mg_apply','mg_matrix'].includes(op)||(op==='cg'&&jacobi===4)?newMG():undefined;
  const factor=op==='ilu_apply'||(op==='gmres'&&jacobi===2)?new ILU0(a):undefined;
  const plan=['chol_factor','chol_solve'].includes(op)?new SparseCholeskySymbolic(a):undefined;
  const chol=op==='chol_solve'?plan!.factorize(a):undefined;
  const ic=op==='ic0_apply'||(op==='cg'&&jacobi===2)?new IC0(a):undefined;
  const compute=()=>{
   if(op==='mg_setup'){const m=newMG();return [m.size,m.levels];}
   if(op==='mg_matrix'){const m=mg!.matrix;return [...m.rowOffsets,...m.columnIndices,...m.values];}
   if(op==='mg_apply')return mg!.apply(b);
   if(op==='ic0_factor'){const f=new IC0(a);let l;try{l=f.lower;return [...l.rowOffsets,...l.columnIndices,...l.values];}finally{}}
   if(op==='ic0_apply')return Array.from(ic!.apply(b));
   if(op==='chol_symbolic'){const s=new SparseCholeskySymbolic(a);try{return [...s.rowOffsets,...s.columnIndices,...s.fillSteps];}finally{}}
   if(op==='chol_factor'){const f=plan!.factorize(a);let l;try{l=f.lower;return [...l.rowOffsets,...l.columnIndices,...l.values];}finally{}}
   if(op==='chol_solve')return chol!.solve(b);
   if(op==='chol_total'||op==='chol_rcm_total'||op==='chol_amd_total'){
    const p=op==='chol_amd_total'?a.approximateMinimumDegree():op==='chol_rcm_total'?a.reverseCuthillMcKee():undefined,q=p?a.permuteSymmetric(p):a;let s,f;
    try{s=new SparseCholeskySymbolic(q);f=s.factorize(q);const x=f.solve(p?CSRMatrix.permuteVector(p,b):b);return p?CSRMatrix.permuteVector(p,x,true):x;}finally{}
   }

   if(op==='rcm_solve'||op==='ilu_solve'){
    const p=op==='rcm_solve'?a.reverseCuthillMcKee():Array.from({length:rows},(_,i)=>i),q=op==='rcm_solve'?a.permuteSymmetric(p):a,f=new ILU0(q);
    try{const rhs=op==='rcm_solve'?CSRMatrix.permuteVector(p,b):b,r=q.gmres(rhs,{restart:20,rtol,atol,maxIterations,preconditioner:f});if(!r.converged)throw new Error('ordering solve failed: '+r.reason);return [r.iterations,...(op==='rcm_solve'?CSRMatrix.permuteVector(p,r.x,true):r.x)];}
    finally{/* Managed factors. */}
   }
   if(op==='amd')return a.approximateMinimumDegree();
   if(op==='rcm')return a.reverseCuthillMcKee();
   if(op==='permute'||op==='permutation_check'){
    const q=a.permuteSymmetric(b);
    try{if(op==='permute')return [...q.rowOffsets,...q.columnIndices,...q.values];
     const x=Array.from({length:rows},(_,i)=>i+1),y=CSRMatrix.permuteVector(b,x);
     return [...y,...CSRMatrix.permuteVector(b,y,true),...CSRMatrix.permuteVector(b,q.matvec(y),true)];
    }finally{/* Managed CSR. */}
   }
   if(op==='ilu_setup'){const f=new ILU0(a);return [f.size,f.nnz];}
   if(op==='ilu_apply')return factor!.apply(b);
   if(op==='dense'){const out=dense!.multiply(right!);return Array.from(out.values);}
   if(op==='spmv')return Array.from(a.matvec(b));
   const currentMG=op==='cg'?(jacobi===5?newMG():mg):undefined;
   const currentIC=op==='cg'&&jacobi===3?new IC0(a):ic;
   const current=op==='gmres'&&jacobi===3?new ILU0(a):factor;
   const r=op==='gmres'?a.gmres(b,{restart,rtol,atol,maxIterations,jacobi:jacobi===1,capture:!!capture,preconditioner:current}):a.conjugateGradient(b,{rtol,atol,maxIterations,jacobi:jacobi===1,capture:!!capture,preconditioner:currentMG??currentIC});
   if(iters){if(!r.converged)throw new Error('benchmark CG did not converge: '+r.reason);return r.x;}
   const extra=op==='gmres'?r as GMRESResult:undefined;
   return [['converged','iteration_limit','breakdown','nonfinite','stagnation'].indexOf(r.reason),r.iterations,r.residuals.length,...r.x,...r.residuals,...(extra?[...extra.estimatedResiduals,extra.restarts.length,...extra.restarts]:[]),...r.iterates.flat()];
  };
  if(!iters){const values=compute();return {rows:1,cols:values.length,values};}
  for(let i=0;i<3;i++)compute();let checksum=0;const start=process.hrtime.bigint();
  for(let i=0;i<iters;i++)for(const x of compute())checksum+=x;
  const elapsed_ns=Number(process.hrtime.bigint()-start);if(!Number.isFinite(checksum))throw new Error('nonfinite benchmark checksum');
  return {elapsed_ns,iterations:iters,checksum};
 }finally{/* Managed arrays require no explicit disposal. */}
}
