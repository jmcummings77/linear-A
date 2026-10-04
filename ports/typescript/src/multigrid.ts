import {CSRMatrix} from './sparse.js';
export interface MultigridFrame {level:number;width:number;phase:string;x:number[];b:number[];residual:number[];}
/** Symmetric V-cycle for the unit five-point Dirichlet Laplacian. */
export class GeometricMultigrid {
 readonly #width:number;
 get width(){return this.#width;}
 constructor(width:number){if(!Number.isInteger(width)||width<1||width>255||((width+1)&width)!==0)throw new RangeError('grid width must be 2**k-1 in 1..255');this.#width=width;}
 get size(){return this.width*this.width;}get levels(){return Math.log2(this.width+1);}
 private mv(w:number,x:number[]){return x.map((v,i)=>{const y=Math.floor(i/w),j=i%w;return 4*v-(j?x[i-1]:0)-(j+1<w?x[i+1]:0)-(y?x[i-w]:0)-(y+1<w?x[i+w]:0);});}
 get matrix(){const w=this.width,rp=[0],ci:number[]=[],v:number[]=[];for(let i=0;i<this.size;i++){const y=Math.floor(i/w),j=i%w,cols=[i];if(j)cols.push(i-1);if(j+1<w)cols.push(i+1);if(y)cols.push(i-w);if(y+1<w)cols.push(i+w);for(const k of cols.sort((a,b)=>a-b)){ci.push(k);v.push(k===i?4:-1);}rp.push(v.length);}return new CSRMatrix(this.size,this.size,rp,ci,v);}
 apply(b:ArrayLike<number>){return this.run(b,false).x;}trace(b:ArrayLike<number>){return this.run(b,true);}
 private run(input:ArrayLike<number>,capture:boolean){const b=Array.from(input);if(b.length!==this.size||b.some(v=>!Number.isFinite(v)))throw new RangeError('invalid multigrid RHS');const frames:MultigridFrame[]=[];
  const record=(level:number,width:number,phase:string,x:number[],b:number[])=>{if(x.some(v=>!Number.isFinite(v)))throw new RangeError('nonfinite multigrid cycle');if(capture)frames.push({level,width,phase,x:x.slice(),b:b.slice(),residual:this.mv(width,x).map((v,i)=>b[i]-v)});};
  const cycle=(w:number,b:number[],level:number):number[]=>{let x=Array(b.length).fill(0);record(level,w,'enter',x,b);if(w===1){x=[b[0]/4];record(level,w,'coarse_solve',x,b);return x;}
   const smooth=(v:number[])=>{for(let k=0;k<2;k++){const ax=this.mv(w,v);v=v.map((x,i)=>x+(b[i]-ax[i])/6);}return v;};x=smooth(x);record(level,w,'pre_smooth',x,b);const r=this.mv(w,x).map((v,i)=>b[i]-v),c=Math.floor(w/2),bc=Array(c*c).fill(0);
   for(let y=0;y<c;y++)for(let j=0;j<c;j++)for(let dy=-1;dy<=1;dy++)for(let dx=-1;dx<=1;dx++)bc[y*c+j]+=(dy===0?1:.5)*(dx===0?1:.5)*r[(2*y+1+dy)*w+2*j+1+dx];
   const ec=cycle(c,bc,level+1);for(let y=0;y<c;y++)for(let j=0;j<c;j++)for(let dy=-1;dy<=1;dy++)for(let dx=-1;dx<=1;dx++)x[(2*y+1+dy)*w+2*j+1+dx]+=(dy===0?1:.5)*(dx===0?1:.5)*ec[y*c+j];
   record(level,w,'correct',x,b);x=smooth(x);record(level,w,'post_smooth',x,b);return x;};return {x:cycle(this.width,b,0),frames};
 }
}
