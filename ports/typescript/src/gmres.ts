import type {ILU0} from './ilu.js';
import type {CSRMatrix} from './sparse.js';
export interface GMRESOptions {restart?:number;rtol?:number;atol?:number;maxIterations?:number;jacobi?:boolean;capture?:boolean;preconditioner?:ILU0;}
export interface GMRESResult {x:number[];converged:boolean;iterations:number;reason:string;residuals:number[];estimatedResiduals:number[];iterates:number[][];restarts:number[];}
export function gmres(a:CSRMatrix,input:ArrayLike<number>,{restart=30,rtol=1e-10,atol=0,maxIterations=1000,jacobi=false,capture=false,preconditioner}:GMRESOptions={}):GMRESResult {
 const b=Array.from(input),n=a.rows;
 if(a.cols!==n||b.length!==n||b.some(x=>!Number.isFinite(x)))throw new RangeError('GMRES requires square matrix and finite matching vector');
 if(!Number.isInteger(restart)||restart<1||restart>1024||!Number.isInteger(maxIterations)||maxIterations<0||maxIterations>100000||!Number.isFinite(rtol)||rtol<0||rtol>=1||!Number.isFinite(atol)||atol<0||typeof jacobi!=='boolean'||typeof capture!=='boolean')throw new RangeError('invalid GMRES options');
 if(preconditioner&&(jacobi||preconditioner.size!==n))throw new RangeError('incompatible preconditioner');
 const apply=(v:number[])=>preconditioner?preconditioner.apply(v):v.map((x,i)=>x/diagonal[i]);
 const diagonal=Array(n).fill(1);
 if(jacobi){const rp=a.rowOffsets,ci=a.columnIndices,v=a.values;for(let i=0;i<n;i++){let found=false;for(let p=rp[i];p<rp[i+1];p++)if(ci[p]===i){diagonal[i]=v[p];found=true;break;}if(!found||diagonal[i]===0)throw new RangeError('Jacobi requires a nonzero diagonal');}}
 const norm=(v:number[])=>v.reduce((s,x)=>Math.hypot(s,x),0),finite=(v:number[])=>v.every(Number.isFinite);
 let x:number[]=Array(n).fill(0),r=b.slice();const residuals=[norm(r)],estimatedResiduals=residuals.slice(),iterates=capture?[x.slice()]:[],restarts:number[]=[];
 const result=(reason:string):GMRESResult=>({x:x.slice(),converged:reason==='converged',iterations:residuals.length-1,reason,residuals,estimatedResiduals,iterates,restarts});
 const threshold=Math.max(atol,rtol*residuals[0]),m=Math.min(restart,n,maxIterations);
 if(!Number.isFinite(residuals[0]))return result('nonfinite');if(residuals[0]<=threshold)return result('converged');
 while(residuals.length-1<maxIterations){
  if(residuals.length>1)restarts.push(residuals.length-1);
  const base=x.slice(),beta=norm(r),basis=[r.map(v=>v/beta)],h=Array.from({length:m+1},()=>Array(m).fill(0)),cs=Array(m).fill(0),sn=Array(m).fill(0),g=[beta,...Array(m).fill(0)],steps=Math.min(m,maxIterations-(residuals.length-1));
  for(let j=0;j<steps;j++){
   let w:number[];try{w=a.matvec(apply(basis[j]));}catch{return result('nonfinite');}
   const original=norm(w);
   for(let pass=0;pass<2;pass++)for(let k=0;k<=j;k++){const dot=basis[k].reduce((s,v,i)=>s+v*w[i],0);h[k][j]+=dot;w=w.map((v,i)=>v-dot*basis[k][i]);}
   const tail=norm(w);if(!Number.isFinite(original)||!Number.isFinite(tail)||h.slice(0,j+1).some(row=>!Number.isFinite(row[j])))return result('nonfinite');
   const happy=tail<=8*Number.EPSILON*original;h[j+1][j]=happy?0:tail;if(!happy)basis.push(w.map(v=>v/tail));
   for(let k=0;k<j;k++){const top=cs[k]*h[k][j]+sn[k]*h[k+1][j];h[k+1][j]=-sn[k]*h[k][j]+cs[k]*h[k+1][j];h[k][j]=top;}
   const pivot=Math.hypot(h[j][j],h[j+1][j]);if(!Number.isFinite(pivot))return result('nonfinite');if(pivot===0)return result('breakdown');
   cs[j]=h[j][j]/pivot;sn[j]=h[j+1][j]/pivot;h[j][j]=pivot;h[j+1][j]=0;g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];const y=g.slice(0,j+1);
   for(let k=j;k>=0;k--){if(h[k][k]===0)return result('breakdown');let sum=0;for(let q=k+1;q<=j;q++)sum+=h[k][q]*y[q];y[k]=(y[k]-sum)/h[k][k];}
   let correction:number[];try{correction=apply(base.map((_,i)=>{let sum=0;for(let k=0;k<=j;k++)sum+=basis[k][i]*y[k];return sum;}));}catch{return result('nonfinite');}
   const candidate=base.map((v,i)=>v+correction[i]);
   if(!finite(y)||!finite(candidate)||!Number.isFinite(g[j+1]))return result('nonfinite');
   let ax:number[];try{ax=a.matvec(candidate);}catch{return result('nonfinite');}const residual=b.map((v,i)=>v-ax[i]),length=norm(residual);
   if(!Number.isFinite(length))return result('nonfinite');x=candidate;r=residual;residuals.push(length);estimatedResiduals.push(Math.abs(g[j+1]));if(capture)iterates.push(x.slice());
   if(length<=threshold)return result('converged');if(happy)return result('breakdown');
  }
  if(x.every((v,i)=>v===base[i]))return result('stagnation');
 }
 return result('iteration_limit');
}
