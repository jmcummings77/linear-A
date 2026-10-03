import { Matrix } from './matrix.js';
/** Economy A = U diag(values) Vt; see ports/SVD.md for precision and ownership. */
export interface SingularValueDecomposition { u: Matrix; values: Float64Array; vt: Matrix; }
export function decompose(source: Matrix, tolerance=1e-12, maxSweeps=100): SingularValueDecomposition {
  if (!Number.isFinite(tolerance) || tolerance<=0 || tolerance>=1 || !Number.isInteger(maxSweeps) || maxSweeps<1 || maxSweeps>10000) throw new RangeError('invalid SVD options');
  const m=source.rows,n=source.cols;
  if(m<n){const r=decompose(source.transpose(),tolerance,maxSweeps);return {u:r.vt.transpose(),values:r.values,vt:r.u.transpose()};}
  const raw=source.values;let scale=0;
  for(const x of raw){if(!Number.isFinite(x))throw new RangeError('SVD requires finite input');scale=Math.max(scale,Math.abs(x));}
  scale=scale||1;
  const b=raw.map(x=>x/scale),v=new Float64Array(n*n);
  if(raw.some((x,i)=>x!==0&&b[i]===0))throw new RangeError('SVD scaling discards an entry');
  for(let j=0;j<n;j++)v[j*n+j]=1;
  const norm=(j:number)=>{let s=0;for(let i=0;i<m;i++)s=Math.hypot(s,b[i*n+j]);return s;};
  let converged=false;
  for(let sweep=0;sweep<=maxSweeps;sweep++){
    let changed=false;
    for(let p=0;p<n;p++)for(let q=p+1;q<n;q++){
      const np=norm(p),nq=norm(q);if(!np||!nq)continue;
      let corr=0;for(let i=0;i<m;i++)corr+=(b[i*n+p]/np)*(b[i*n+q]/nq);
      if(Math.abs(corr)<=tolerance)continue;
      changed=true;if(sweep===maxSweeps)continue;
      const pair=Math.max(np,nq),ap=np/pair,aq=nq/pair,delta=aq*aq-ap*ap,g=2*ap*aq*corr;
      const t=delta===0?(g<0?-1:1):g/(delta+(delta<0?-1:1)*Math.hypot(delta,g));
      if(Math.abs(t)<2.2250738585072014e-308){const small=np<nq?p:q;for(let i=0;i<m;i++)b[i*n+small]=0;continue;}
      const c=1/Math.hypot(1,t),s=c*t;
      for(const [data,count] of [[b,m],[v,n]] as const)for(let i=0;i<count;i++){const x=data[i*n+p],y=data[i*n+q];data[i*n+p]=c*x-s*y;data[i*n+q]=s*x+c*y;}
    }
    if(!changed){converged=true;break;}
  }
  if(!converged)throw new RangeError('SVD did not converge');
  const norms=Array.from({length:n},(_,j)=>norm(j)),order=Array.from({length:n},(_,j)=>j).sort((a,b)=>norms[b]-norms[a]);
  const u=new Float64Array(m*n),vt=new Float64Array(n*n),values=new Float64Array(n);
  for(let j=0;j<n;j++){
    const k=order[j];values[j]=norms[k]*scale;
    if(!Number.isFinite(values[j])||(norms[k]!==0&&values[j]===0))throw new RangeError('singular value outside float64 range');
    for(let i=0;i<n;i++)vt[j*n+i]=v[i*n+k];
    if(norms[k]){for(let i=0;i<m;i++)u[i*n+j]=b[i*n+k]/norms[k];}
    else{
      let found=false;
      for(let axis=0;axis<m;axis++){
        const candidate=Float64Array.from({length:m},(_,i)=>+(i===axis));
        for(let pass=0;pass<2;pass++)for(let col=0;col<j;col++){
          let dot=0;for(let i=0;i<m;i++)dot+=candidate[i]*u[i*n+col];
          for(let i=0;i<m;i++)candidate[i]-=dot*u[i*n+col];
        }
        let length=0;for(const x of candidate)length=Math.hypot(length,x);
        if(length>0.5/Math.sqrt(m)){for(let i=0;i<m;i++)u[i*n+j]=candidate[i]/length;found=true;break;}
      }
      if(!found)throw new RangeError('cannot complete SVD null basis');
    }
  }
  return {u:new Matrix(m,n,u),values,vt:new Matrix(n,n,vt)};
}
