import {GeometricMultigrid} from './multigrid.js';
import {IC0} from "./cholesky.js";
import {gmres, type GMRESOptions} from "./gmres.js";
/** Canonical CSR storage; zero-based, sorted unique column indices in each row. */
export interface CGResult { x:number[]; converged:boolean; iterations:number; reason:string; residuals:number[]; iterates:number[][]; }
export interface CGOptions { rtol?:number; atol?:number; maxIterations?:number; jacobi?:boolean; capture?:boolean; preconditioner?:IC0|GeometricMultigrid; }
export class CSRMatrix {
  readonly rows:number; readonly cols:number;
  private rp:number[];private ci:number[];private v:number[];
  constructor(rows:number,cols:number,offsets:ArrayLike<number>,indices:ArrayLike<number>,values:ArrayLike<number>){
    const rp=Array.from(offsets),ci=Array.from(indices),v=Array.from(values);
    if(!Number.isSafeInteger(rows)||!Number.isSafeInteger(cols)||rows<0||cols<0||rp.length!==rows+1||ci.length!==v.length||rp[0]!==0||rp[rows]!==v.length)throw new RangeError('invalid CSR dimensions or arrays');
    if(rp.some(x=>!Number.isSafeInteger(x)||x<0||x>v.length)||v.some(x=>!Number.isFinite(x)))throw new RangeError('invalid CSR offsets or values');
    for(let i=0;i<rows;i++){
      if(rp[i]>rp[i+1])throw new RangeError('CSR offsets must be monotone');
      let previous=-1;
      for(let p=rp[i];p<rp[i+1];p++){if(!Number.isSafeInteger(ci[p])||ci[p]<=previous||ci[p]>=cols)throw new RangeError('CSR columns must be sorted, unique and in range');previous=ci[p];}
    }
    this.rows=rows;this.cols=cols;this.rp=rp;this.ci=ci;this.v=v;
  }
  get nnz(){return this.v.length;}get rowOffsets(){return this.rp.slice();}get columnIndices(){return this.ci.slice();}get values(){return this.v.slice();}
  /** Quotient-graph AMD; deterministic new-to-old indices. */
  approximateMinimumDegree():number[]{
    if(this.rows!==this.cols)throw new RangeError('AMD requires square matrix');
    const n=this.rows,direct=Array.from({length:n},()=>new Set<number>()),elements=Array.from({length:n},()=>new Set<number>());
    for(let i=0;i<n;i++)for(let k=this.rp[i];k<this.rp[i+1];k++){const j=this.ci[k];if(i!==j){direct[i].add(j);direct[j].add(i);}}
    const active=new Set(Array.from({length:n},(_,i)=>i)),degree=direct.map(g=>g.size),order:number[]=[];
    while(active.size){
      const pivot=[...active].sort((a,b)=>Math.min(degree[a],active.size-1)-Math.min(degree[b],active.size-1)||a-b)[0];
      const neighbors=new Set(direct[pivot]);
      for(const e of elements)if(e.has(pivot)){for(const i of e)neighbors.add(i);e.clear();}
      neighbors.delete(pivot);active.delete(pivot);order.push(pivot);direct[pivot].clear();
      for(const i of neighbors){direct[i].delete(pivot);for(const j of neighbors)direct[i].delete(j);}
      elements[pivot]=neighbors;
      for(const i of neighbors){let bound=neighbors.size-1+direct[i].size;for(let e=0;e<n;e++)if(e!==pivot&&elements[e].has(i))for(const j of elements[e])if(!neighbors.has(j))bound++;degree[i]=Math.min(active.size-1,bound);}
    }return order;
  }

  reverseCuthillMcKee():number[]{
    if(this.rows!==this.cols)throw new RangeError('RCM requires square matrix');
    const graph=Array.from({length:this.rows},()=>new Set<number>());
    for(let i=0;i<this.rows;i++)for(let k=this.rp[i];k<this.rp[i+1];k++){const j=this.ci[k];if(i!==j){graph[i].add(j);graph[j].add(i);}}
    const compare=(a:number,b:number)=>graph[a].size-graph[b].size||a-b,seen=new Set<number>(),order:number[]=[];
    for(const start of Array.from({length:this.rows},(_,i)=>i).sort(compare)){
      if(seen.has(start))continue;const queue=[start];seen.add(start);
      for(let h=0;h<queue.length;h++)for(const j of [...graph[queue[h]]].filter(j=>!seen.has(j)).sort(compare)){seen.add(j);queue.push(j);}
      order.push(...queue);
    }return order.reverse();
  }
  static permuteVector(order:ArrayLike<number>,input:ArrayLike<number>,inverse=false):number[]{
    const p=Array.from(order),x=Array.from(input),n=x.length;
    if(p.length!==n||p.some(i=>!Number.isSafeInteger(i)||i<0||i>=n)||new Set(p).size!==n||x.some(v=>!Number.isFinite(v)))throw new RangeError('invalid permutation or vector');
    if(!inverse)return p.map(i=>x[i]);const out=Array(n).fill(0);p.forEach((j,i)=>out[j]=x[i]);return out;
  }
  permuteSymmetric(order:ArrayLike<number>):CSRMatrix{
    if(this.rows!==this.cols)throw new RangeError('permutation requires square matrix');
    const p=Array.from(order),inv=CSRMatrix.permuteVector(p,Array.from({length:this.rows},(_,i)=>i),true),rp=[0],ci:number[]=[],v:number[]=[];
    for(const i of p){const entries: [number,number][]=[];for(let k=this.rp[i];k<this.rp[i+1];k++)entries.push([inv[this.ci[k]],this.v[k]]);entries.sort((a,b)=>a[0]-b[0]);for(const [j,x] of entries){ci.push(j);v.push(x);}rp.push(v.length);}
    return new CSRMatrix(this.rows,this.cols,rp,ci,v);
  }
  static fromDense(a:{rows:number;cols:number;get(i:number,j:number):number}){
    const rp=[0],ci:number[]=[],v:number[]=[];
    for(let i=0;i<a.rows;i++){for(let j=0;j<a.cols;j++){const x=a.get(i,j);if(x!==0){ci.push(j);v.push(x);}}rp.push(v.length);}
    return new CSRMatrix(a.rows,a.cols,rp,ci,v);
  }
  matvec(input:ArrayLike<number>):number[]{
    const x=Array.from(input);if(x.length!==this.cols||x.some(v=>!Number.isFinite(v)))throw new RangeError('invalid vector');
    return Array.from({length:this.rows},(_,i)=>{let sum=0;for(let p=this.rp[i];p<this.rp[i+1];p++)sum+=this.v[p]*x[this.ci[p]];if(!Number.isFinite(sum))throw new RangeError('sparse multiplication outside float64 range');return sum;});
  }
  private find(row:number,col:number):number{let lo=this.rp[row],hi=this.rp[row+1];while(lo<hi){const mid=lo+Math.floor((hi-lo)/2);if(this.ci[mid]<col)lo=mid+1;else hi=mid;}return lo;}
  conjugateGradient(input:ArrayLike<number>,{rtol=1e-10,atol=0,maxIterations=1000,jacobi=false,capture=false,preconditioner}:CGOptions={}):CGResult{
    const b=Array.from(input),n=this.rows;
    if(this.cols!==n||b.length!==n||b.some(x=>!Number.isFinite(x)))throw new RangeError('CG requires square matrix and finite matching vector');
    if(!Number.isFinite(rtol)||rtol<0||rtol>=1||!Number.isFinite(atol)||atol<0||!Number.isInteger(maxIterations)||maxIterations<0||maxIterations>100000||typeof jacobi!=='boolean'||typeof capture!=='boolean')throw new RangeError('invalid CG options');
    if(preconditioner!==undefined&&(!(preconditioner instanceof IC0)&&!(preconditioner instanceof GeometricMultigrid)||preconditioner.size!==n||jacobi))throw new RangeError('invalid or conflicting CG preconditioner');
    const diagonal=Array(n).fill(1);
    for(let i=0;i<n;i++){
      for(let p=this.rp[i];p<this.rp[i+1];p++){const j=this.ci[p],q=this.find(j,i),other=q<this.rp[j+1]&&this.ci[q]===i?this.v[q]:0;if(this.v[p]!==other)throw new RangeError('CG requires exact symmetry');}
      if(jacobi){const q=this.find(i,i);if(q===this.rp[i+1]||this.ci[q]!==i||this.v[q]<=0)throw new RangeError('Jacobi requires positive diagonal');diagonal[i]=this.v[q];}
    }
    const norm=(v:number[])=>v.reduce((s,x)=>Math.hypot(s,x),0),dot=(a:number[],b:number[])=>a.reduce((s,x,i)=>s+x*b[i],0);
    let x:number[]=Array(n).fill(0),r=b.slice();const residuals=[norm(r)],iterates=capture?[x.slice()]:[],threshold=Math.max(atol,rtol*residuals[0]);
    const result=(reason:string):CGResult=>({x:x.slice(),residuals:residuals.slice(),iterates,iterations:residuals.length-1,reason,converged:reason==='converged'});
    if(!Number.isFinite(residuals[0]))return result('nonfinite');if(residuals[0]<=threshold)return result('converged');
    const apply=(r:number[])=>preconditioner?preconditioner.apply(r):r.map((v,i)=>v/diagonal[i]);
    let z:number[];try{z=apply(r);}catch{return result('nonfinite');}let p=z.slice(),rho=dot(r,z);
    for(let step=0;step<maxIterations;step++){
      let q:number[];try{q=this.matvec(p);}catch{return result('nonfinite');}
      const curvature=dot(p,q);if(!Number.isFinite(rho)||!Number.isFinite(curvature))return result('nonfinite');if(rho<=0||curvature<=0)return result('breakdown');
      const alpha=rho/curvature,candidate=x.map((v,i)=>v+alpha*p[i]);let ax:number[];
      try{ax=this.matvec(candidate);}catch{return result('nonfinite');}
      const residual=b.map((v,i)=>v-ax[i]),length=norm(residual);if(!Number.isFinite(length))return result('nonfinite');
      x=candidate;r=residual;residuals.push(length);if(capture)iterates.push(x.slice());if(length<=threshold)return result('converged');
      try{z=apply(r);}catch{return result('nonfinite');}const next=dot(r,z);if(!Number.isFinite(next))return result('nonfinite');if(next<=0)return result('breakdown');
      const beta=next/rho;p=z.map((v,i)=>v+beta*p[i]);rho=next;
    }
    return result('iteration_limit');
  }
  gmres(b:ArrayLike<number>,options:GMRESOptions={}){return gmres(this,b,options);}
}
