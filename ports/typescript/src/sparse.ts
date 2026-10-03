/** Canonical CSR storage; zero-based, sorted unique column indices in each row. */
export interface CGResult { x:number[]; converged:boolean; iterations:number; reason:string; residuals:number[]; iterates:number[][]; }
export interface CGOptions { rtol?:number; atol?:number; maxIterations?:number; jacobi?:boolean; capture?:boolean; }
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
  conjugateGradient(input:ArrayLike<number>,{rtol=1e-10,atol=0,maxIterations=1000,jacobi=false,capture=false}:CGOptions={}):CGResult{
    const b=Array.from(input),n=this.rows;
    if(this.cols!==n||b.length!==n||b.some(x=>!Number.isFinite(x)))throw new RangeError('CG requires square matrix and finite matching vector');
    if(!Number.isFinite(rtol)||rtol<0||rtol>=1||!Number.isFinite(atol)||atol<0||!Number.isInteger(maxIterations)||maxIterations<0||maxIterations>100000||typeof jacobi!=='boolean'||typeof capture!=='boolean')throw new RangeError('invalid CG options');
    const diagonal=Array(n).fill(1);
    for(let i=0;i<n;i++){
      for(let p=this.rp[i];p<this.rp[i+1];p++){const j=this.ci[p],q=this.find(j,i),other=q<this.rp[j+1]&&this.ci[q]===i?this.v[q]:0;if(this.v[p]!==other)throw new RangeError('CG requires exact symmetry');}
      if(jacobi){const q=this.find(i,i);if(q===this.rp[i+1]||this.ci[q]!==i||this.v[q]<=0)throw new RangeError('Jacobi requires positive diagonal');diagonal[i]=this.v[q];}
    }
    const norm=(v:number[])=>v.reduce((s,x)=>Math.hypot(s,x),0),dot=(a:number[],b:number[])=>a.reduce((s,x,i)=>s+x*b[i],0);
    let x:number[]=Array(n).fill(0),r=b.slice();const residuals=[norm(r)],iterates=capture?[x.slice()]:[],threshold=Math.max(atol,rtol*residuals[0]);
    const result=(reason:string):CGResult=>({x:x.slice(),residuals:residuals.slice(),iterates,iterations:residuals.length-1,reason,converged:reason==='converged'});
    if(!Number.isFinite(residuals[0]))return result('nonfinite');if(residuals[0]<=threshold)return result('converged');
    let z=r.map((v,i)=>v/diagonal[i]),p=z.slice(),rho=dot(r,z);
    for(let step=0;step<maxIterations;step++){
      let q:number[];try{q=this.matvec(p);}catch{return result('nonfinite');}
      const curvature=dot(p,q);if(!Number.isFinite(rho)||!Number.isFinite(curvature))return result('nonfinite');if(rho<=0||curvature<=0)return result('breakdown');
      const alpha=rho/curvature,candidate=x.map((v,i)=>v+alpha*p[i]);let ax:number[];
      try{ax=this.matvec(candidate);}catch{return result('nonfinite');}
      const residual=b.map((v,i)=>v-ax[i]),length=norm(residual);if(!Number.isFinite(length))return result('nonfinite');
      x=candidate;r=residual;residuals.push(length);if(capture)iterates.push(x.slice());if(length<=threshold)return result('converged');
      z=r.map((v,i)=>v/diagonal[i]);const next=dot(r,z);if(!Number.isFinite(next))return result('nonfinite');if(next<=0)return result('breakdown');
      const beta=next/rho;p=z.map((v,i)=>v+beta*p[i]);rho=next;
    }
    return result('iteration_limit');
  }
}
