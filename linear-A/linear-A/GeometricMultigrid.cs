#if NET10_0_OR_GREATER
using System;
using System.Collections.Generic;
using System.Linq;
namespace linear_A;

/// <summary>A fixed symmetric positive-definite preconditioner with matching vector size.</summary>
public interface ISymmetricPreconditioner
{
    int Size { get; }
    double[] Apply(double[] b);
}
/// <summary>Symmetric V-cycle for the unit five-point Dirichlet Laplacian.</summary>
public sealed class GeometricMultigrid : ISymmetricPreconditioner
{
    public int Width { get; }
    public int Size => Width * Width;
    public int Levels => int.Log2(Width + 1);
    public GeometricMultigrid(int width)
    {
        if (width < 1 || width > 255 || (width & (width + 1)) != 0) throw new ArgumentException("Grid width must be 2^k-1 in 1..255.");
        Width = width;
    }
    public CSRMatrix Matrix
    {
        get
        {
            var rp = new List<int> { 0 }; var ci = new List<int>(); var v = new List<double>();
            for (int i = 0; i < Size; i++)
            {
                void Add(int j) { ci.Add(j); v.Add(j == i ? 4 : -1); }
                if (i >= Width) Add(i - Width); if (i % Width > 0) Add(i - 1); Add(i); if (i % Width + 1 < Width) Add(i + 1); if (i + Width < Size) Add(i + Width); rp.Add(v.Count);
            }
            return new CSRMatrix(Size, Size, rp.ToArray(), ci.ToArray(), v.ToArray());
        }
    }
    private static double[] Multiply(int w, double[] x)
    {
        var r = new double[x.Length];
        for (int i = 0; i < x.Length; i++) r[i] = 4 * x[i] - (i % w > 0 ? x[i - 1] : 0) - (i % w + 1 < w ? x[i + 1] : 0) - (i >= w ? x[i - w] : 0) - (i + w < x.Length ? x[i + w] : 0);
        return r;
    }
    private static double[] Cycle(int w, double[] b)
    {
        if (w == 1) return [b[0] / 4];
        var x = new double[b.Length];
        void Smooth() { for (int k = 0; k < 2; k++) { var ax = Multiply(w, x); for (int i = 0; i < x.Length; i++) x[i] += (b[i] - ax[i]) / 6; } }
        Smooth(); var r = Multiply(w, x); for (int i = 0; i < r.Length; i++) r[i] = b[i] - r[i];
        int c = w / 2; var bc = new double[c * c];
        for (int y = 0; y < c; y++) for (int j = 0; j < c; j++) for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) bc[y * c + j] += (dy == 0 ? 1 : .5) * (dx == 0 ? 1 : .5) * r[(2 * y + 1 + dy) * w + 2 * j + 1 + dx];
        var ec = Cycle(c, bc);
        for (int y = 0; y < c; y++) for (int j = 0; j < c; j++) for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) x[(2 * y + 1 + dy) * w + 2 * j + 1 + dx] += (dy == 0 ? 1 : .5) * (dx == 0 ? 1 : .5) * ec[y * c + j];
        Smooth(); if (x.Any(v => !double.IsFinite(v))) throw new ArithmeticException("Nonfinite multigrid cycle."); return x;
    }
    public double[] Apply(double[] b)
    {
        ArgumentNullException.ThrowIfNull(b); if (b.Length != Size || b.Any(v => !double.IsFinite(v))) throw new ArgumentException("Invalid multigrid RHS.");
        return Cycle(Width, b);
    }
}
#endif
