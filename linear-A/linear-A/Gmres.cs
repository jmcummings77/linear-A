#if NET10_0_OR_GREATER
using System;
using System.Collections.Generic;
using System.Linq;
namespace linear_A;
public sealed record GMRESResult(double[] X, bool Converged, int Iterations, string Reason, double[] Residuals, double[] EstimatedResiduals, double[][] Iterates, int[] Restarts);
public sealed partial class CSRMatrix
{
    public GMRESResult Gmres(double[] b, int restart=30, double relativeTolerance=1e-10, double absoluteTolerance=0, int maxIterations=1000, bool jacobi=false, bool capture=false, ILU0? preconditioner=null)
    {
        ArgumentNullException.ThrowIfNull(b);int n=Rows;
        if(Cols!=n||b.Length!=n||b.Any(v=>!double.IsFinite(v)))throw new ArgumentException("GMRES requires square matrix and finite matching vector.");
        if(restart<1||restart>1024||maxIterations<0||maxIterations>100000||!double.IsFinite(relativeTolerance)||relativeTolerance<0||relativeTolerance>=1||!double.IsFinite(absoluteTolerance)||absoluteTolerance<0)throw new ArgumentOutOfRangeException(nameof(restart),"Invalid GMRES options.");
        if(preconditioner is not null && (jacobi || preconditioner.Size!=n))throw new ArgumentException("Incompatible preconditioner.");
        var diagonal=Enumerable.Repeat(1.0,n).ToArray();
        if(jacobi)for(int i=0;i<n;i++){bool found=false;for(int p=offsets[i];p<offsets[i+1];p++)if(indices[p]==i){diagonal[i]=values[p];found=true;break;}if(!found||diagonal[i]==0)throw new ArgumentException("Jacobi requires a nonzero diagonal.");}
        double[] Apply(double[] v)=>preconditioner is null?v.Select((z,i)=>z/diagonal[i]).ToArray():preconditioner.Apply(v);
        static double Norm(double[] v){double s=0;foreach(var z in v)s=double.Hypot(s,z);return s;}
        var x=new double[n];var r=(double[])b.Clone();var history=new List<double>{Norm(r)};var estimates=new List<double>(history);var frames=new List<double[]>();var restarts=new List<int>();if(capture)frames.Add((double[])x.Clone());
        GMRESResult Result(string reason)=>new((double[])x.Clone(),reason=="converged",history.Count-1,reason,history.ToArray(),estimates.ToArray(),frames.ToArray(),restarts.ToArray());
        double threshold=Math.Max(absoluteTolerance,relativeTolerance*history[0]);int m=Math.Min(restart,Math.Min(n,maxIterations));
        if(!double.IsFinite(history[0]))return Result("nonfinite");if(history[0]<=threshold)return Result("converged");
        while(history.Count-1<maxIterations){
            if(history.Count>1)restarts.Add(history.Count-1);var start=(double[])x.Clone();double beta=Norm(r);var basis=new List<double[]>{r.Select(v=>v/beta).ToArray()};var h=new double[m+1,m];var cs=new double[m];var sn=new double[m];var g=new double[m+1];g[0]=beta;int steps=Math.Min(m,maxIterations-(history.Count-1));
            for(int j=0;j<steps;j++){
                double[] w;try{w=Matvec(Apply(basis[j]));}catch(ArithmeticException){return Result("nonfinite");}catch(ArgumentException){return Result("nonfinite");}double original=Norm(w);
                for(int pass=0;pass<2;pass++)for(int k=0;k<=j;k++){double dot=0;for(int i=0;i<n;i++)dot+=basis[k][i]*w[i];h[k,j]+=dot;for(int i=0;i<n;i++)w[i]-=dot*basis[k][i];}
                double tail=Norm(w);if(!double.IsFinite(original)||!double.IsFinite(tail))return Result("nonfinite");for(int k=0;k<=j;k++)if(!double.IsFinite(h[k,j]))return Result("nonfinite");
                bool happy=tail<=8*2.220446049250313e-16*original;h[j+1,j]=happy?0:tail;if(!happy)basis.Add(w.Select(v=>v/tail).ToArray());
                for(int k=0;k<j;k++){double top=cs[k]*h[k,j]+sn[k]*h[k+1,j];h[k+1,j]=-sn[k]*h[k,j]+cs[k]*h[k+1,j];h[k,j]=top;}
                double pivot=double.Hypot(h[j,j],h[j+1,j]);if(!double.IsFinite(pivot))return Result("nonfinite");if(pivot==0)return Result("breakdown");
                cs[j]=h[j,j]/pivot;sn[j]=h[j+1,j]/pivot;h[j,j]=pivot;h[j+1,j]=0;g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];var y=g[..(j+1)];
                for(int k=j;k>=0;k--){if(h[k,k]==0)return Result("breakdown");double sum=0;for(int q=k+1;q<=j;q++)sum+=h[k,q]*y[q];y[k]=(y[k]-sum)/h[k,k];}
                var candidate=new double[n];for(int i=0;i<n;i++){double sum=0;for(int k=0;k<=j;k++)sum+=basis[k][i]*y[k];candidate[i]=sum;}
                try{candidate=Apply(candidate);}catch(ArithmeticException){return Result("nonfinite");}for(int i=0;i<n;i++)candidate[i]+=start[i];
                if(y.Any(v=>!double.IsFinite(v))||candidate.Any(v=>!double.IsFinite(v))||!double.IsFinite(g[j+1]))return Result("nonfinite");
                double[] ax;try{ax=Matvec(candidate);}catch(ArithmeticException){return Result("nonfinite");}catch(ArgumentException){return Result("nonfinite");}var residual=b.Select((v,i)=>v-ax[i]).ToArray();double length=Norm(residual);
                if(!double.IsFinite(length))return Result("nonfinite");x=candidate;r=residual;history.Add(length);estimates.Add(Math.Abs(g[j+1]));if(capture)frames.Add((double[])x.Clone());
                if(length<=threshold)return Result("converged");if(happy)return Result("breakdown");
            }
            if(x.SequenceEqual(start))return Result("stagnation");
        }
        return Result("iteration_limit");
    }
}
#endif
