#if NET10_0_OR_GREATER
using System;
using System.Collections.Generic;
using System.Linq;
namespace linear_A;

/// <summary>Owned CG solution, true residual history, and optional solution snapshots.</summary>
public sealed record CGResult(double[] X, bool Converged, int Iterations, string Reason, double[] Residuals, double[][] Iterates);

/// <summary>Canonical compressed sparse rows; sorted unique zero-based column indices.</summary>
public sealed partial class CSRMatrix
{
    private readonly int[] offsets, indices;
    private readonly double[] values;
    public int Rows { get; }
    public int Cols { get; }
    public int NNZ => values.Length;
    public int[] RowOffsets => (int[])offsets.Clone();
    public int[] ColumnIndices => (int[])indices.Clone();
    public double[] Values => (double[])values.Clone();
    public CSRMatrix(int rows, int cols, int[] rowOffsets, int[] columnIndices, double[] entries)
    {
        ArgumentNullException.ThrowIfNull(rowOffsets); ArgumentNullException.ThrowIfNull(columnIndices); ArgumentNullException.ThrowIfNull(entries);
        offsets = (int[])rowOffsets.Clone(); indices = (int[])columnIndices.Clone(); values = (double[])entries.Clone();
        if (rows < 0 || rows == int.MaxValue || cols < 0 || offsets.Length != rows + 1 || indices.Length != values.Length || offsets[0] != 0 || offsets[rows] != values.Length)
            throw new ArgumentException("Invalid CSR dimensions or array lengths.");
        if (offsets.Any(x => x < 0 || x > values.Length) || values.Any(x => !double.IsFinite(x))) throw new ArgumentException("Invalid CSR offsets or values.");
        for (var i = 0; i < rows; i++)
        {
            if (offsets[i] > offsets[i + 1]) throw new ArgumentException("CSR offsets must be monotone.");
            var previous = -1;
            for (var p = offsets[i]; p < offsets[i + 1]; p++)
            {
                if (indices[p] <= previous || indices[p] >= cols) throw new ArgumentException("CSR columns must be sorted, unique and in range.");
                previous = indices[p];
            }
        }
        Rows = rows; Cols = cols;
    }
    public static CSRMatrix FromDense(Matrix<double> a)
    {
        ArgumentNullException.ThrowIfNull(a);
        var rp = new List<int> { 0 }; var ci = new List<int>(); var v = new List<double>();
        for (var i = 0; i < a.RowCount; i++)
        {
            for (var j = 0; j < a.ColumnCount; j++) if (a[i, j] != 0) { ci.Add(j); v.Add(a[i, j]); }
            rp.Add(v.Count);
        }
        return new(a.RowCount, a.ColumnCount, rp.ToArray(), ci.ToArray(), v.ToArray());
    }
    public double[] Matvec(double[] x)
    {
        ArgumentNullException.ThrowIfNull(x);
        if (x.Length != Cols || x.Any(v => !double.IsFinite(v))) throw new ArgumentException("Invalid vector.");
        var result = new double[Rows];
        for (var i = 0; i < Rows; i++)
        {
            for (var p = offsets[i]; p < offsets[i + 1]; p++) result[i] += values[p] * x[indices[p]];
            if (!double.IsFinite(result[i])) throw new ArithmeticException("Sparse multiplication outside float64 range.");
        }
        return result;
    }
    private int Find(int row, int col)
    {
        var lo = offsets[row]; var hi = offsets[row + 1];
        while (lo < hi) { var mid = lo + (hi - lo) / 2; if (indices[mid] < col) lo = mid + 1; else hi = mid; }
        return lo;
    }
    public CGResult ConjugateGradient(double[] b, double relativeTolerance = 1e-10, double absoluteTolerance = 0, int maxIterations = 1000, bool jacobi = false, bool capture = false)
    {
        ArgumentNullException.ThrowIfNull(b); var n = Rows;
        if (Cols != n || b.Length != n || b.Any(x => !double.IsFinite(x))) throw new ArgumentException("CG requires square matrix and finite matching vector.");
        if (!double.IsFinite(relativeTolerance) || relativeTolerance < 0 || relativeTolerance >= 1 || !double.IsFinite(absoluteTolerance) || absoluteTolerance < 0 || maxIterations < 0 || maxIterations > 100000) throw new ArgumentException("Invalid CG options.");
        var diagonal = Enumerable.Repeat(1.0, n).ToArray();
        for (var i = 0; i < n; i++)
        {
            for (var p = offsets[i]; p < offsets[i + 1]; p++)
            {
                var j = indices[p]; var q = Find(j, i); var other = q < offsets[j + 1] && indices[q] == i ? values[q] : 0;
                if (values[p] != other) throw new ArgumentException("CG requires exact symmetry.");
            }
            if (jacobi) { var q = Find(i, i); if (q == offsets[i + 1] || indices[q] != i || values[q] <= 0) throw new ArgumentException("Jacobi requires a positive diagonal."); diagonal[i] = values[q]; }
        }
        static double Norm(double[] v) { var s = 0.0; foreach (var x in v) s = double.Hypot(s, x); return s; }
        static double Dot(double[] a, double[] b) { var s = 0.0; for (var i = 0; i < a.Length; i++) s += a[i] * b[i]; return s; }
        var x = new double[n]; var r = (double[])b.Clone(); var history = new List<double> { Norm(r) }; var frames = new List<double[]>();
        if (capture) frames.Add((double[])x.Clone()); var threshold = Math.Max(absoluteTolerance, relativeTolerance * history[0]);
        CGResult Result(string reason) => new((double[])x.Clone(), reason == "converged", history.Count - 1, reason, history.ToArray(), frames.ToArray());
        if (!double.IsFinite(history[0])) return Result("nonfinite"); if (history[0] <= threshold) return Result("converged");
        var z = r.Select((v, i) => v / diagonal[i]).ToArray(); var direction = (double[])z.Clone(); var rho = Dot(r, z);
        for (var step = 0; step < maxIterations; step++)
        {
            double[] q; try { q = Matvec(direction); } catch (Exception e) when (e is ArithmeticException or ArgumentException) { return Result("nonfinite"); }
            var curvature = Dot(direction, q); if (!double.IsFinite(rho) || !double.IsFinite(curvature)) return Result("nonfinite"); if (rho <= 0 || curvature <= 0) return Result("breakdown");
            var alpha = rho / curvature; var candidate = x.Select((v, i) => v + alpha * direction[i]).ToArray(); double[] ax;
            try { ax = Matvec(candidate); } catch (Exception e) when (e is ArithmeticException or ArgumentException) { return Result("nonfinite"); }
            var residual = b.Select((v, i) => v - ax[i]).ToArray(); var length = Norm(residual); if (!double.IsFinite(length)) return Result("nonfinite");
            x = candidate; r = residual; history.Add(length); if (capture) frames.Add((double[])x.Clone()); if (length <= threshold) return Result("converged");
            z = r.Select((v, i) => v / diagonal[i]).ToArray(); var next = Dot(r, z); if (!double.IsFinite(next)) return Result("nonfinite"); if (next <= 0) return Result("breakdown");
            var beta = next / rho; direction = z.Select((v, i) => v + beta * direction[i]).ToArray(); rho = next;
        }
        return Result("iteration_limit");
    }
}
#endif
