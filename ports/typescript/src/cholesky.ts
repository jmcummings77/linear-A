import {CSRMatrix} from './sparse.js';
/** Symbolic lower structure, including stored zeros. Fill steps are zero-based;
 * -1 marks original edges and diagonals. Apply an ordering before analysis. */
export class SparseCholeskySymbolic {
 readonly size:number;
 private sourceRP:number[];private sourceCI:number[];
 private rp:number[]=[0];private ci:number[]=[];private steps:number[]=[];
 constructor(a:CSRMatrix){
  if(a.rows!==a.cols)throw new RangeError('Cholesky requires square matrix');
  this.size=a.rows;this.sourceRP=a.rowOffsets;this.sourceCI=a.columnIndices;
  const g=Array.from({length:this.size},()=>new Map<number,number>());
  for(let i=0;i<this.size;i++)for(let p=this.sourceRP[i];p<this.sourceRP[i+1];p++){
   const j=this.sourceCI[p];if(i!==j){g[i].set(j,-1);g[j].set(i,-1);}
  }
  for(let k=0;k<this.size;k++){
   const ns=[...g[k].keys()].filter(j=>j>k).sort((a,b)=>a-b);
   for(let u=0;u<ns.length;u++)for(let w=0;w<u;w++){
    const i=ns[u],j=ns[w];if(!g[i].has(j)){g[i].set(j,k);g[j].set(i,k);}
   }
  }
  for(let i=0;i<this.size;i++){
   for(const j of [...g[i].keys()].filter(j=>j<i).sort((a,b)=>a-b)){this.ci.push(j);this.steps.push(g[i].get(j)!);}
   this.ci.push(i);this.steps.push(-1);this.rp.push(this.ci.length);
  }
 }
 get nnz(){return this.ci.length;}get fillCount(){return this.steps.filter(k=>k>=0).length;}
 get rowOffsets(){return this.rp.slice();}get columnIndices(){return this.ci.slice();}get fillSteps(){return this.steps.slice();}
 factorize(a:CSRMatrix):SparseCholesky{
  const ar=a.rowOffsets,ac=a.columnIndices,av=a.values;
  if(a.rows!==this.size||a.cols!==this.size||ar.length!==this.sourceRP.length||ac.length!==this.sourceCI.length||ar.some((x,i)=>x!==this.sourceRP[i])||ac.some((x,i)=>x!==this.sourceCI[i]))throw new RangeError('Cholesky symbolic pattern mismatch');
  const rows=Array.from({length:this.size},(_,i)=>new Map(ac.slice(ar[i],ar[i+1]).map((j,k)=>[j,av[ar[i]+k]])));
  for(let i=0;i<this.size;i++)for(const [j,v] of rows[i])if(v!==(rows[j].get(i)??0))throw new RangeError('Cholesky requires symmetric values');
  const {rp,ci}=this,v=new Array<number>(ci.length).fill(0);
  for(let i=0;i<this.size;i++)for(let p=rp[i];p<rp[i+1];p++){
   const j=ci[p];let s=rows[i].get(j)??0,u=rp[i],w=rp[j];
   while(u<p&&w<rp[j+1]-1){if(ci[u]===ci[w]){s-=v[u]*v[w];u++;w++;}else if(ci[u]<ci[w])u++;else w++;}
   if(!Number.isFinite(s))throw new RangeError('nonfinite Cholesky factor');
   if(i===j){if(s<=0)throw new RangeError('nonpositive Cholesky pivot');v[p]=Math.sqrt(s);}else v[p]=s/v[rp[j+1]-1];
   if(!Number.isFinite(v[p]))throw new RangeError('nonfinite Cholesky factor');
  }
  return new SparseCholesky(new CSRMatrix(this.size,this.size,rp,ci,v));
 }
}
/** Owned lower factor; normally constructed by a symbolic plan. */
export class SparseCholesky {
 private a:CSRMatrix;
 constructor(lower:CSRMatrix){
  const rp=lower.rowOffsets,ci=lower.columnIndices,v=lower.values;
  if(lower.rows!==lower.cols)throw new RangeError('invalid lower factor');
  for(let i=0;i<lower.rows;i++)if(rp[i]===rp[i+1]||ci[rp[i+1]-1]!==i||v[rp[i+1]-1]<=0)throw new RangeError('invalid lower factor');
  this.a=new CSRMatrix(lower.rows,lower.cols,rp,ci,v);
 }
 get size(){return this.a.rows;}get nnz(){return this.a.nnz;}
 get lower(){return new CSRMatrix(this.size,this.size,this.a.rowOffsets,this.a.columnIndices,this.a.values);}
 solve(b:ArrayLike<number>):number[]{
  const x=Array.from(b),rp=this.a.rowOffsets,ci=this.a.columnIndices,v=this.a.values;
  if(x.length!==this.size||x.some(z=>!Number.isFinite(z)))throw new RangeError('invalid Cholesky right-hand side');
  for(let i=0;i<this.size;i++){for(let p=rp[i];p<rp[i+1]-1;p++)x[i]-=v[p]*x[ci[p]];x[i]/=v[rp[i+1]-1];if(!Number.isFinite(x[i]))throw new RangeError('nonfinite Cholesky solve');}
  for(let i=this.size-1;i>=0;i--){x[i]/=v[rp[i+1]-1];if(!Number.isFinite(x[i]))throw new RangeError('nonfinite Cholesky solve');for(let p=rp[i];p<rp[i+1]-1;p++){x[ci[p]]-=v[p]*x[i];if(!Number.isFinite(x[ci[p]]))throw new RangeError('nonfinite Cholesky solve');}}
  return x;
 }
}
