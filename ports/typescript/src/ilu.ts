import type {CSRMatrix} from './sparse.js';
/** Owned zero-fill incomplete LU. No pivoting or diagonal shifts. */
export class ILU0 {
 readonly size:number;private rp:number[];private ci:number[];private v:number[];private d:number[];
 constructor(a:CSRMatrix){
  if(a.rows!==a.cols)throw new RangeError('ILU0 requires square matrix');
  this.size=a.rows;this.rp=a.rowOffsets;this.ci=a.columnIndices;this.v=a.values;this.d=[];
  const {rp,ci,v,d}=this;
  for(let i=0;i<this.size;i++){const p=this.find(i,i);if(p===rp[i+1]||ci[p]!==i)throw new RangeError('ILU0 requires stored diagonal');d.push(p);}
  for(let i=0;i<this.size;i++){
   for(let p=rp[i];p<d[i];p++){const j=ci[p];v[p]/=v[d[j]];if(!Number.isFinite(v[p]))throw new RangeError('nonfinite ILU0 factor');
    for(let q=d[j]+1;q<rp[j+1];q++){const k=this.find(i,ci[q]);if(k<rp[i+1]&&ci[k]===ci[q]){v[k]-=v[p]*v[q];if(!Number.isFinite(v[k]))throw new RangeError('nonfinite ILU0 factor');}}
   }
   if(v[d[i]]===0)throw new RangeError('zero ILU0 pivot');
  }
 }
 private find(i:number,j:number){let lo=this.rp[i],hi=this.rp[i+1];while(lo<hi){const m=lo+Math.floor((hi-lo)/2);if(this.ci[m]<j)lo=m+1;else hi=m;}return lo;}
 get nnz(){return this.v.length;}
 apply(b:ArrayLike<number>):number[]{
  const x=Array.from(b),{rp,ci,v,d}=this;if(x.length!==this.size||x.some(z=>!Number.isFinite(z)))throw new RangeError('invalid ILU0 vector');
  for(let i=0;i<this.size;i++){for(let p=rp[i];p<d[i];p++)x[i]-=v[p]*x[ci[p]];if(!Number.isFinite(x[i]))throw new RangeError('nonfinite ILU0 solve');}
  for(let i=this.size-1;i>=0;i--){for(let p=d[i]+1;p<rp[i+1];p++)x[i]-=v[p]*x[ci[p]];x[i]/=v[d[i]];if(!Number.isFinite(x[i]))throw new RangeError('nonfinite ILU0 solve');}return x;
 }
}
