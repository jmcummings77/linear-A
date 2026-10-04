#if NET10_0_OR_GREATER
using System;
using System.Collections.Generic;
using System.Linq;

namespace linear_A;

/// <summary>Reusable symbolic lower pattern. FillSteps uses -1 for original entries.</summary>
public sealed class SparseCholeskySymbolic
{
    private readonly int[] sourceRP, sourceCI, rp, ci, steps;
    public int Size { get; }
    public int NNZ => ci.Length;
    public int FillCount => steps.Count(k => k >= 0);
    public int[] RowOffsets => (int[])rp.Clone();
    public int[] ColumnIndices => (int[])ci.Clone();
    public int[] FillSteps => (int[])steps.Clone();
    public SparseCholeskySymbolic(CSRMatrix a) : this(a, false) { }
    internal SparseCholeskySymbolic(CSRMatrix a, bool incomplete)
    {
        ArgumentNullException.ThrowIfNull(a);
        if (a.Rows != a.Cols) throw new ArgumentException("Cholesky requires square matrix.");
        Size = a.Rows; sourceRP = a.RowOffsets; sourceCI = a.ColumnIndices;
        var g = Enumerable.Range(0, Size).Select(_ => new SortedDictionary<int, int>()).ToArray();
        for (int i = 0; i < Size; i++) for (int p = sourceRP[i]; p < sourceRP[i + 1]; p++)
            { int j = sourceCI[p]; if (i != j) { g[i][j] = -1; g[j][i] = -1; } }
        for (int k = 0; !incomplete && k < Size; k++)
        {
            var ns = g[k].Keys.Where(j => j > k).ToArray();
            for (int u = 0; u < ns.Length; u++) for (int w = 0; w < u; w++)
                { int i = ns[u], j = ns[w]; if (!g[i].ContainsKey(j)) { g[i][j] = k; g[j][i] = k; } }
        }
        rp = new int[Size + 1]; var columns = new List<int>(); var births = new List<int>();
        for (int i = 0; i < Size; i++)
        {
            foreach (var pair in g[i]) if (pair.Key < i) { columns.Add(pair.Key); births.Add(pair.Value); }
            columns.Add(i); births.Add(-1); rp[i + 1] = columns.Count;
        }
        ci = columns.ToArray(); steps = births.ToArray();
    }
    public SparseCholesky Factorize(CSRMatrix a)
    {
        ArgumentNullException.ThrowIfNull(a);
        if (a.Rows != Size || a.Cols != Size || !a.RowOffsets.SequenceEqual(sourceRP) || !a.ColumnIndices.SequenceEqual(sourceCI)) throw new ArgumentException("Cholesky symbolic pattern mismatch.");
        var av = a.Values;
        var rows = Enumerable.Range(0, Size).Select(_ => new Dictionary<int, double>()).ToArray();
        for (int i = 0; i < Size; i++) for (int p = sourceRP[i]; p < sourceRP[i + 1]; p++) rows[i][sourceCI[p]] = av[p];
        for (int i = 0; i < Size; i++) foreach (var pair in rows[i]) if (pair.Value != rows[pair.Key].GetValueOrDefault(i)) throw new ArgumentException("Cholesky requires symmetric values.");
        var v = new double[ci.Length];
        for (int i = 0; i < Size; i++) for (int p = rp[i]; p < rp[i + 1]; p++)
        {
            int j = ci[p], u = rp[i], w = rp[j]; double s = rows[i].GetValueOrDefault(j);
            while (u < p && w < rp[j + 1] - 1) { if (ci[u] == ci[w]) { s -= v[u] * v[w]; u++; w++; } else if (ci[u] < ci[w]) u++; else w++; }
            if (!double.IsFinite(s)) throw new ArithmeticException("Nonfinite Cholesky factor.");
            if (i == j) { if (s <= 0) throw new ArithmeticException("Nonpositive Cholesky pivot."); v[p] = Math.Sqrt(s); }
            else v[p] = s / v[rp[j + 1] - 1];
            if (!double.IsFinite(v[p])) throw new ArithmeticException("Nonfinite Cholesky factor.");
        }
        return new SparseCholesky(new CSRMatrix(Size, Size, rp, ci, v));
    }
}
/// <summary>Owned lower factor supporting repeated right-hand sides.</summary>
public sealed class SparseCholesky
{
    private readonly CSRMatrix lower;
    internal SparseCholesky(CSRMatrix lower) { this.lower = lower; }
    public int Size => lower.Rows;
    public int NNZ => lower.NNZ;
    public CSRMatrix Lower => new(Size, Size, lower.RowOffsets, lower.ColumnIndices, lower.Values);
    public double[] Solve(double[] b)
    {
        ArgumentNullException.ThrowIfNull(b);
        if (b.Length != Size || b.Any(z => !double.IsFinite(z))) throw new ArgumentException("Invalid Cholesky right-hand side.");
        var x = (double[])b.Clone(); var rp = lower.RowOffsets; var ci = lower.ColumnIndices; var v = lower.Values;
        for (int i = 0; i < Size; i++)
        {
            for (int p = rp[i]; p < rp[i + 1] - 1; p++) x[i] -= v[p] * x[ci[p]];
            x[i] /= v[rp[i + 1] - 1]; if (!double.IsFinite(x[i])) throw new ArithmeticException("Nonfinite Cholesky solve.");
        }
        for (int i = Size - 1; i >= 0; i--)
        {
            x[i] /= v[rp[i + 1] - 1]; if (!double.IsFinite(x[i])) throw new ArithmeticException("Nonfinite Cholesky solve.");
            for (int p = rp[i]; p < rp[i + 1] - 1; p++) { x[ci[p]] -= v[p] * x[i]; if (!double.IsFinite(x[ci[p]])) throw new ArithmeticException("Nonfinite Cholesky solve."); }
        }
        return x;
    }
}
/// <summary>Owned zero-fill incomplete Cholesky; no shifts or pivoting.</summary>
public sealed class IC0 : ISymmetricPreconditioner
{
    private readonly SparseCholesky factor;
    public IC0(CSRMatrix a) { factor = new SparseCholeskySymbolic(a, true).Factorize(a); }
    public int Size => factor.Size;
    public int NNZ => factor.NNZ;
    public CSRMatrix Lower => factor.Lower;
    public double[] Apply(double[] b) => factor.Solve(b);
}
#endif
