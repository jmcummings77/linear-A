import { Matrix, finite } from "./matrix.js";
export type FactorAlgorithm = "lu" | "cholesky" | "qr";
const scaled = (value:number, scale:number):number => {
  const result=finite(value/scale);
  if(value!==0 && result===0)throw new RangeError("solver scaling would discard a nonzero value");
  return result;
};
/** Independent reusable factors. QR is column-pivoted and requires full column rank. */
export class Factorization {
  private readonly rows:number;
  private readonly cols:number;
  private readonly algorithm:FactorAlgorithm;
  private readonly data:Float64Array;
  private readonly permutation:number[];
  private readonly tau:Float64Array;
  private readonly scale:number;
  private readonly norm:number;
  constructor(matrix:Matrix,algorithm:FactorAlgorithm){
    if(!["lu","cholesky","qr"].includes(algorithm))throw new RangeError("unknown factorization");
    const m=matrix.rows,n=matrix.cols;
    if((algorithm!=="qr"&&m!==n)||m<n)throw new RangeError("factorization shape is invalid");
    this.rows=m;this.cols=n;this.algorithm=algorithm;
    const original=matrix.values;let maximum=0;
    for(const value of original)maximum=Math.max(maximum,Math.abs(finite(value)));
    this.scale=maximum?2**Math.min(1023,Math.floor(Math.log2(maximum))):1;
    const a=Float64Array.from(original,value=>scaled(value,this.scale));this.data=a;
    let matrixNorm=0;
    for(let i=0;i<m;i++){let sum=0;for(let j=0;j<n;j++)sum+=Math.abs(a[i*n+j]);matrixNorm=Math.max(matrixNorm,sum);}
    this.norm=matrixNorm;
    const p=Array.from({length:algorithm==="lu"?m:n},(_,i)=>i);this.permutation=p;
    const tau=new Float64Array(n);this.tau=tau;
    if(algorithm==="lu"){
      for(let k=0;k<n;k++){
        let pivot=k;for(let i=k+1;i<n;i++)if(Math.abs(a[i*n+k])>Math.abs(a[pivot*n+k]))pivot=i;
        if(a[pivot*n+k]===0)throw new RangeError("singular matrix: zero computed LU pivot");
        for(let j=0;j<n;j++){const t=a[k*n+j];a[k*n+j]=a[pivot*n+j];a[pivot*n+j]=t;}
        [p[k],p[pivot]]=[p[pivot],p[k]];
        for(let i=k+1;i<n;i++){
          a[i*n+k]=finite(a[i*n+k]/a[k*n+k]);
          for(let j=k+1;j<n;j++)a[i*n+j]=finite(a[i*n+j]-a[i*n+k]*a[k*n+j]);
        }
      }
    }else if(algorithm==="cholesky"){
      for(let i=0;i<n;i++)for(let j=0;j<i;j++)if(original[i*n+j]!==original[j*n+i])throw new RangeError("Cholesky requires exact symmetry");
      for(let i=0;i<n;i++)for(let j=0;j<=i;j++){
        let value=a[i*n+j];for(let k=0;k<j;k++)value=finite(value-a[i*n+k]*a[j*n+k]);
        if(i===j){if(value<=0)throw new RangeError("Cholesky requires positive computed pivots");a[i*n+j]=Math.sqrt(value);}
        else a[i*n+j]=finite(value/a[j*n+j]);
      }
    }else{
      let largest=0;for(let j=0;j<n;j++){let norm=0;for(let i=0;i<m;i++)norm=Math.hypot(norm,a[i*n+j]);largest=Math.max(largest,norm);}
      const threshold=Number.EPSILON*Math.max(m,n)*largest;
      for(let k=0;k<n;k++){
        let pivot=k,norm=-1;
        for(let j=k;j<n;j++){let candidate=0;for(let i=k;i<m;i++)candidate=Math.hypot(candidate,a[i*n+j]);if(candidate>norm){pivot=j;norm=candidate;}}
        if(norm<=threshold)throw new RangeError("QR input is numerically rank deficient");
        for(let i=0;i<m;i++){const t=a[i*n+k];a[i*n+k]=a[i*n+pivot];a[i*n+pivot]=t;}
        [p[k],p[pivot]]=[p[pivot],p[k]];
        const old=a[k*n+k],alpha=old<0||Object.is(old,-0)?norm:-norm,divisor=old-alpha;
        tau[k]=(alpha-old)/alpha;for(let i=k+1;i<m;i++)a[i*n+k]/=divisor;a[k*n+k]=alpha;
        for(let j=k+1;j<n;j++){
          let dot=a[k*n+j];for(let i=k+1;i<m;i++)dot+=a[i*n+k]*a[i*n+j];dot*=tau[k];
          a[k*n+j]=finite(a[k*n+j]-dot);for(let i=k+1;i<m;i++)a[i*n+j]=finite(a[i*n+j]-a[i*n+k]*dot);
        }
      }
    }
  }
  private solveInternal(rhs:Matrix,rescale:boolean):Matrix{
    const m=this.rows,n=this.cols,p=rhs.cols,a=this.data,b=rhs.values;
    if(rhs.rows!==m)throw new RangeError("right-hand side row count must match factorization");
    const work=new Float64Array(m*p);
    for(let i=0;i<m;i++)for(let j=0;j<p;j++){const v=b[(this.algorithm==="lu"?this.permutation[i]:i)*p+j];work[i*p+j]=rescale?scaled(v,this.scale):finite(v);}
    if(this.algorithm==="qr"){
      for(let k=0;k<n;k++)for(let j=0;j<p;j++){
        let dot=work[k*p+j];for(let i=k+1;i<m;i++)dot=finite(dot+a[i*n+k]*work[i*p+j]);dot=finite(dot*this.tau[k]);
        work[k*p+j]=finite(work[k*p+j]-dot);for(let i=k+1;i<m;i++)work[i*p+j]=finite(work[i*p+j]-a[i*n+k]*dot);
      }
    }else for(let i=0;i<n;i++)for(let j=0;j<p;j++){
      let value=work[i*p+j];for(let k=0;k<i;k++)value=finite(value-a[i*n+k]*work[k*p+j]);
      work[i*p+j]=this.algorithm==="cholesky"?finite(value/a[i*n+i]):value;
    }
    for(let i=n-1;i>=0;i--)for(let j=0;j<p;j++){
      let value=work[i*p+j];for(let k=i+1;k<n;k++)value=finite(value-(this.algorithm==="cholesky"?a[k*n+i]:a[i*n+k])*work[k*p+j]);
      work[i*p+j]=finite(value/a[i*n+i]);
    }
    const result=new Float64Array(n*p);
    for(let i=0;i<n;i++)for(let j=0;j<p;j++)result[(this.algorithm==="qr"?this.permutation[i]:i)*p+j]=work[i*p+j];
    return new Matrix(n,p,result);
  }
  solve(rhs:Matrix):Matrix{return this.solveInternal(rhs,true);}
  /** Computed infinity-norm rcond from n solves; O(n^3), not a certified bound. */
  reciprocalCondition():number{
    const n=this.cols;if(this.rows!==n)throw new RangeError("condition diagnostic requires square input");if(!n)return 1;
    try{
      const inverse=this.solveInternal(Matrix.identity(n),false).values;let norm=0;
      for(let i=0;i<n;i++){let sum=0;for(let j=0;j<n;j++)sum+=Math.abs(inverse[i*n+j]);norm=Math.max(norm,sum);}
      return Math.min(1,(1/this.norm)/norm);
    }catch(error){if(error instanceof RangeError && error.message === "matrix arithmetic requires finite float64 values")return 0;throw error;}
  }
}
