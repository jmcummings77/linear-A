#if NET10_0_OR_GREATER
using System;
using System.Linq;

namespace linear_A;

/// <summary>Economy factors A = U diag(Values) Vt with descending nonnegative values.</summary>
public sealed record SingularValueDecomposition(Matrix<double> U, double[] Values, Matrix<double> Vt);

/// <summary>One-sided Jacobi SVD for finite real matrices; see ports/SVD.md.</summary>
public static class MatrixSvd
{
    /// <summary>Returns independent economy factors without modifying the input.</summary>
    public static SingularValueDecomposition Svd(this Matrix<double> source, double tolerance = 1e-12, int maxSweeps = 100)
    {
        ArgumentNullException.ThrowIfNull(source);
        if (!double.IsFinite(tolerance) || tolerance <= 0 || tolerance >= 1 || maxSweeps < 1 || maxSweeps > 10000)
            throw new ArgumentOutOfRangeException(nameof(tolerance), "Invalid SVD options.");
        var m = source.RowCount; var n = source.ColumnCount;
        if (m < n)
        {
            var r = source.Transpose().Svd(tolerance, maxSweeps);
            return new(r.Vt.Transpose(), r.Values, r.U.Transpose());
        }
        var b = new double[checked(m * n)]; var v = new double[checked(n * n)];
        var scale = 0.0;
        for (var i = 0; i < m; i++)
            for (var j = 0; j < n; j++)
            {
                var x = source[i, j];
                if (!double.IsFinite(x)) throw new ArgumentException("SVD requires finite input.", nameof(source));
                scale = Math.Max(scale, Math.Abs(x));
            }
        if (scale == 0) scale = 1;
        for (var i = 0; i < m; i++)
            for (var j = 0; j < n; j++)
            {
                b[i * n + j] = source[i, j] / scale;
                if (source[i, j] != 0 && b[i * n + j] == 0) throw new ArithmeticException("SVD scaling discards an entry.");
            }
        for (var j = 0; j < n; j++) v[j * n + j] = 1;
        double Norm(int j)
        {
            var sum = 0.0;
            for (var i = 0; i < m; i++) sum = double.Hypot(sum, b[i * n + j]);
            return sum;
        }
        var converged = false;
        for (var sweep = 0; sweep <= maxSweeps; sweep++)
        {
            var changed = false;
            for (var p = 0; p < n; p++)
                for (var q = p + 1; q < n; q++)
                {
                    var np = Norm(p); var nq = Norm(q);
                    if (np == 0 || nq == 0) continue;
                    var corr = 0.0;
                    for (var i = 0; i < m; i++) corr += (b[i * n + p] / np) * (b[i * n + q] / nq);
                    if (Math.Abs(corr) <= tolerance) continue;
                    changed = true;
                    if (sweep == maxSweeps) continue;
                    var pair = Math.Max(np, nq); var ap = np / pair; var aq = nq / pair;
                    var delta = aq * aq - ap * ap; var g = 2 * ap * aq * corr;
                    var t = delta == 0 ? Math.CopySign(1, g) : g / (delta + Math.CopySign(double.Hypot(delta, g), delta));
                    if (Math.Abs(t) < 2.2250738585072014e-308)
                    {
                        var small = np < nq ? p : q;
                        for (var i = 0; i < m; i++) b[i * n + small] = 0;
                        continue;
                    }
                    var c = 1 / double.Hypot(1, t); var s = c * t;
                    void Rotate(double[] data, int count)
                    {
                        for (var i = 0; i < count; i++)
                        {
                            var x = data[i * n + p]; var y = data[i * n + q];
                            data[i * n + p] = c * x - s * y; data[i * n + q] = s * x + c * y;
                        }
                    }
                    Rotate(b, m); Rotate(v, n);
                }
            if (!changed) { converged = true; break; }
        }
        if (!converged) throw new ArithmeticException("SVD did not converge.");
        var norms = Enumerable.Range(0, n).Select(Norm).ToArray();
        var order = Enumerable.Range(0, n).OrderByDescending(j => norms[j]).ToArray();
        var u = new Matrix<double>(m, n); var vt = new Matrix<double>(n, n); var values = new double[n];
        for (var j = 0; j < n; j++)
        {
            var k = order[j]; values[j] = norms[k] * scale;
            if (!double.IsFinite(values[j]) || norms[k] != 0 && values[j] == 0)
                throw new ArithmeticException("Singular value outside double range.");
            for (var i = 0; i < n; i++) vt[j, i] = v[i * n + k];
            if (norms[k] != 0)
                for (var i = 0; i < m; i++) u[i, j] = b[i * n + k] / norms[k];
            else
            {
                var found = false;
                for (var axis = 0; axis < m; axis++)
                {
                    var candidate = new double[m]; candidate[axis] = 1;
                    for (var pass = 0; pass < 2; pass++)
                        for (var col = 0; col < j; col++)
                        {
                            var dot = 0.0;
                            for (var i = 0; i < m; i++) dot += candidate[i] * u[i, col];
                            for (var i = 0; i < m; i++) candidate[i] -= dot * u[i, col];
                        }
                    var length = 0.0;
                    foreach (var x in candidate) length = double.Hypot(length, x);
                    if (length > 0.5 / Math.Sqrt(m))
                    {
                        for (var i = 0; i < m; i++) u[i, j] = candidate[i] / length;
                        found = true; break;
                    }
                }
                if (!found) throw new ArithmeticException("Cannot complete SVD null basis.");
            }
        }
        return new(u, values, vt);
    }
}
#endif
