#if NET10_0_OR_GREATER
using System;

namespace linear_A;

/// <summary>Reusable double-precision factorization algorithms.</summary>
public enum FactorizationAlgorithm { Lu, Cholesky, Qr }

/// <summary>An owned snapshot for solving one or more right-hand sides.</summary>
public sealed class MatrixFactorization
{
    private readonly int rows, cols;
    private readonly FactorizationAlgorithm algorithm;
    private readonly double scale, norm;
    private readonly double[] data, tau;
    private readonly int[] permutation;

    private static double Finite(double value) => double.IsFinite(value) ? value : throw new ArithmeticException("Solver arithmetic exceeds the finite double range.");
    private static double Scaled(double value, double scale)
    {
        var result = Finite(value / scale);
        if (value != 0 && result == 0) throw new ArithmeticException("Solver scaling would discard a nonzero value.");
        return result;
    }

    /// <summary>QR requires rows >= columns and numerical full column rank; other algorithms require square input.</summary>
    public MatrixFactorization(Matrix<double> source, FactorizationAlgorithm algorithm)
    {
        ArgumentNullException.ThrowIfNull(source);
        if (!Enum.IsDefined(algorithm)) throw new ArgumentOutOfRangeException(nameof(algorithm));
        rows = source.RowCount; cols = source.ColumnCount; this.algorithm = algorithm;
        if (rows < cols || algorithm != FactorizationAlgorithm.Qr && rows != cols)
            throw new ArgumentException("Factorization requires square input, or rows >= columns for QR.", nameof(source));
        var m = rows; var n = cols;
        data = new double[checked(m * n)]; tau = new double[n]; permutation = new int[m];
        var maximum = 0.0;
        for (var i = 0; i < m; i++)
        {
            permutation[i] = i;
            for (var j = 0; j < n; j++) maximum = Math.Max(maximum, Math.Abs(Finite(source[i, j])));
        }
        scale = maximum == 0 ? 1 : Math.ScaleB(1, Math.ILogB(maximum));
        for (var i = 0; i < m; i++)
        {
            var sum = 0.0;
            for (var j = 0; j < n; j++) { data[i * n + j] = Scaled(source[i, j], scale); sum += Math.Abs(data[i * n + j]); }
            norm = Math.Max(norm, sum);
        }
        var a = data;
        if (algorithm == FactorizationAlgorithm.Lu)
        {
            for (var k = 0; k < n; k++)
            {
                var pivot = k;
                for (var i = k + 1; i < n; i++) if (Math.Abs(a[i * n + k]) > Math.Abs(a[pivot * n + k])) pivot = i;
                if (a[pivot * n + k] == 0) throw new ArgumentException("Singular matrix: zero computed LU pivot.");
                for (var j = 0; j < n; j++) (a[k * n + j], a[pivot * n + j]) = (a[pivot * n + j], a[k * n + j]);
                (permutation[k], permutation[pivot]) = (permutation[pivot], permutation[k]);
                for (var i = k + 1; i < n; i++)
                {
                    a[i * n + k] = Finite(a[i * n + k] / a[k * n + k]);
                    for (var j = k + 1; j < n; j++) a[i * n + j] = Finite(a[i * n + j] - a[i * n + k] * a[k * n + j]);
                }
            }
        }
        else if (algorithm == FactorizationAlgorithm.Cholesky)
        {
            for (var i = 0; i < n; i++)
                for (var j = 0; j < i; j++)
                    if (source[i, j] != source[j, i]) throw new ArgumentException("Cholesky requires exact symmetry.");
            for (var i = 0; i < n; i++)
                for (var j = 0; j <= i; j++)
                {
                    var value = a[i * n + j];
                    for (var k = 0; k < j; k++) value = Finite(value - a[i * n + k] * a[j * n + k]);
                    if (i == j)
                    {
                        if (value <= 0) throw new ArgumentException("Cholesky requires positive computed pivots.");
                        a[i * n + j] = Math.Sqrt(value);
                    }
                    else a[i * n + j] = Finite(value / a[j * n + j]);
                }
        }
        else
        {
            var largest = 0.0;
            for (var j = 0; j < n; j++)
            {
                var length = 0.0;
                for (var i = 0; i < m; i++) length = double.Hypot(length, a[i * n + j]);
                largest = Math.Max(largest, length);
            }
            var threshold = 2.220446049250313e-16 * m * largest;
            for (var k = 0; k < n; k++)
            {
                var pivot = k; var length = 0.0;
                for (var j = k; j < n; j++)
                {
                    var candidate = 0.0;
                    for (var i = k; i < m; i++) candidate = double.Hypot(candidate, a[i * n + j]);
                    if (candidate > length) { pivot = j; length = candidate; }
                }
                if (length <= threshold) throw new ArgumentException("QR input is numerically rank deficient.");
                for (var i = 0; i < m; i++) (a[i * n + k], a[i * n + pivot]) = (a[i * n + pivot], a[i * n + k]);
                (permutation[k], permutation[pivot]) = (permutation[pivot], permutation[k]);
                var old = a[k * n + k]; var alpha = -Math.CopySign(length, old); var divisor = old - alpha;
                tau[k] = (alpha - old) / alpha;
                for (var i = k + 1; i < m; i++) a[i * n + k] /= divisor;
                a[k * n + k] = alpha;
                for (var j = k + 1; j < n; j++)
                {
                    var dot = a[k * n + j];
                    for (var i = k + 1; i < m; i++) dot += a[i * n + k] * a[i * n + j];
                    dot *= tau[k]; a[k * n + j] = Finite(a[k * n + j] - dot);
                    for (var i = k + 1; i < m; i++) a[i * n + j] = Finite(a[i * n + j] - a[i * n + k] * dot);
                }
            }
        }
    }

    /// <summary>Solves using these factors without mutating either input or factor storage.</summary>
    public Matrix<double> Solve(Matrix<double> rhs) => SolveInternal(rhs, true);
    private Matrix<double> SolveInternal(Matrix<double> rhs, bool rescale)
    {
        ArgumentNullException.ThrowIfNull(rhs);
        if (rhs.RowCount != rows) throw new ArgumentException("Right-hand side row count must match factorization.", nameof(rhs));
        var m = rows; var n = cols; var p = rhs.ColumnCount; var a = data;
        var work = new double[checked(m * p)];
        for (var i = 0; i < m; i++)
            for (var j = 0; j < p; j++)
            {
                var value = Finite(rhs[algorithm == FactorizationAlgorithm.Lu ? permutation[i] : i, j]);
                work[i * p + j] = rescale ? Scaled(value, scale) : value;
            }
        if (algorithm == FactorizationAlgorithm.Qr)
        {
            for (var k = 0; k < n; k++)
                for (var j = 0; j < p; j++)
                {
                    var dot = work[k * p + j];
                    for (var i = k + 1; i < m; i++) dot = Finite(dot + a[i * n + k] * work[i * p + j]);
                    dot = Finite(dot * tau[k]); work[k * p + j] = Finite(work[k * p + j] - dot);
                    for (var i = k + 1; i < m; i++) work[i * p + j] = Finite(work[i * p + j] - a[i * n + k] * dot);
                }
        }
        else
        {
            for (var i = 0; i < n; i++)
                for (var j = 0; j < p; j++)
                {
                    var value = work[i * p + j];
                    for (var k = 0; k < i; k++) value = Finite(value - a[i * n + k] * work[k * p + j]);
                    work[i * p + j] = algorithm == FactorizationAlgorithm.Cholesky ? Finite(value / a[i * n + i]) : value;
                }
        }
        for (var i = n - 1; i >= 0; i--)
            for (var j = 0; j < p; j++)
            {
                var value = work[i * p + j];
                for (var k = i + 1; k < n; k++) value = Finite(value - (algorithm == FactorizationAlgorithm.Cholesky ? a[k * n + i] : a[i * n + k]) * work[k * p + j]);
                work[i * p + j] = Finite(value / a[i * n + i]);
            }
        var result = new Matrix<double>(n, p);
        for (var i = 0; i < n; i++)
            for (var j = 0; j < p; j++) result[algorithm == FactorizationAlgorithm.Qr ? permutation[i] : i, j] = work[i * p + j];
        return result;
    }

    /// <summary>Computed infinity-norm reciprocal condition using an inverse, O(n^3), not a certified bound.</summary>
    public double ReciprocalCondition()
    {
        if (rows != cols) throw new InvalidOperationException("Condition diagnostic requires a square matrix.");
        if (cols == 0) return 1;
        var identity = new Matrix<double>(cols, cols);
        for (var i = 0; i < cols; i++) identity[i, i] = 1;
        Matrix<double> inverse;
        try { inverse = SolveInternal(identity, false); }
        catch (ArithmeticException) { return 0; }
        var inverseNorm = 0.0;
        for (var i = 0; i < cols; i++)
        {
            var sum = 0.0;
            for (var j = 0; j < cols; j++) sum += Math.Abs(inverse[i, j]);
            inverseNorm = Math.Max(inverseNorm, sum);
        }
        return Math.Min(1, (1 / norm) / inverseNorm);
    }
}

/// <summary>Double-precision system solving; integer and decimal matrices are not implicitly converted.</summary>
public static class MatrixSolvers
{
    public static MatrixFactorization FactorLu(this Matrix<double> matrix) => new(matrix, FactorizationAlgorithm.Lu);
    public static MatrixFactorization FactorCholesky(this Matrix<double> matrix) => new(matrix, FactorizationAlgorithm.Cholesky);
    public static MatrixFactorization FactorQr(this Matrix<double> matrix) => new(matrix, FactorizationAlgorithm.Qr);
    public static Matrix<double> Solve(this Matrix<double> matrix, Matrix<double> rhs) => matrix.FactorLu().Solve(rhs);
    public static Matrix<double> LeastSquares(this Matrix<double> matrix, Matrix<double> rhs) => matrix.FactorQr().Solve(rhs);
}
#endif
