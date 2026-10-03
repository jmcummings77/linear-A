#if NET10_0_OR_GREATER
using System;
using System.Linq;
namespace linear_A;
/// <summary>Owned ILU(0), without pivoting, fill or diagonal shifts.</summary>
public sealed class ILU0
{
    private readonly int[] rp,ci,d;
    private readonly double[] v;
    public int Size { get; }
    public int NNZ => v.Length;
    public ILU0(CSRMatrix a)
    {
        ArgumentNullException.ThrowIfNull(a);
        if(a.Rows!=a.Cols)throw new ArgumentException("ILU0 requires square matrix.");
        Size=a.Rows;rp=a.RowOffsets;ci=a.ColumnIndices;v=a.Values;d=new int[Size];
        for(int i=0;i<Size;i++){d[i]=Find(i,i);if(d[i]==rp[i+1]||ci[d[i]]!=i)throw new ArgumentException("ILU0 requires stored diagonal.");}
        for(int i=0;i<Size;i++){
            for(int p=rp[i];p<d[i];p++){int j=ci[p];v[p]/=v[d[j]];Finite(v[p]);
                for(int q=d[j]+1;q<rp[j+1];q++){int k=Find(i,ci[q]);if(k<rp[i+1]&&ci[k]==ci[q]){v[k]-=v[p]*v[q];Finite(v[k]);}}
            }
            if(v[d[i]]==0)throw new ArithmeticException("Zero ILU0 pivot.");
        }
    }
    private static void Finite(double value){if(!double.IsFinite(value))throw new ArithmeticException("Nonfinite ILU0 arithmetic.");}
    private int Find(int i,int j){int lo=rp[i],hi=rp[i+1];while(lo<hi){int m=lo+(hi-lo)/2;if(ci[m]<j)lo=m+1;else hi=m;}return lo;}
    public double[] Apply(double[] b)
    {
        ArgumentNullException.ThrowIfNull(b);
        if(b.Length!=Size||b.Any(x=>!double.IsFinite(x)))throw new ArgumentException("Invalid ILU0 vector.");
        var x=(double[])b.Clone();
        for(int i=0;i<Size;i++){for(int p=rp[i];p<d[i];p++)x[i]-=v[p]*x[ci[p]];Finite(x[i]);}
        for(int i=Size-1;i>=0;i--){for(int p=d[i]+1;p<rp[i+1];p++)x[i]-=v[p]*x[ci[p]];x[i]/=v[d[i]];Finite(x[i]);}
        return x;
    }
}
#endif
