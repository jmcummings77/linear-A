using System;

namespace linear_A;

/// <summary>The real eigendecomposition of a symmetric matrix.</summary>
/// <remarks>
/// Eigenvalues are in ascending order. Column j of <see cref="EigenVectors"/>
/// is a unit eigenvector for EigenValues[j]. Repeated eigenvalues may have any
/// orthonormal basis of their eigenspace. The result owns its mutable storage.
/// </remarks>
public sealed class SymmetricEigenDecomposition
{
    internal SymmetricEigenDecomposition(double[] values, DoubleMatrix vectors)
    {
        EigenValues = values;
        EigenVectors = vectors;
    }

    /// <summary>Gets the ascending eigenvalues.</summary>
    public double[] EigenValues { get; }

    /// <summary>Gets the corresponding orthonormal eigenvectors as columns.</summary>
    public DoubleMatrix EigenVectors { get; }
}

internal static class SymmetricEigenSolver
{
    internal static SymmetricEigenDecomposition Solve(double[,] source, double tolerance, int maxSweeps)
    {
        if (double.IsNaN(tolerance) || double.IsInfinity(tolerance) || tolerance <= 0 || tolerance >= 1)
            throw new ArgumentOutOfRangeException(nameof(tolerance), "Tolerance must be finite and between zero and one.");
        if (maxSweeps <= 0) throw new ArgumentOutOfRangeException(nameof(maxSweeps));
        var n = source.GetLength(0);
        if (n != source.GetLength(1)) throw new NotSquareMatrixException();
        if (n == 0) return new SymmetricEigenDecomposition(Array.Empty<double>(), new DoubleMatrix());
        var scale = 0d;
        var diagonal = true;
        for (var row = 0; row < n; row++)
        {
            for (var column = 0; column < n; column++)
            {
                var value = source[row, column];
                if (double.IsNaN(value) || double.IsInfinity(value))
                    throw new ArgumentException("Eigendecomposition requires finite entries.", nameof(source));
                if (value != source[column, row])
                    throw new ArgumentException("Eigendecomposition requires an exactly symmetric matrix.", nameof(source));
                scale = Math.Max(scale, Math.Abs(value));
                if (row != column && value != 0) diagonal = false;
            }
        }

        var vectors = new DoubleMatrix(n, true);
        var work = (double[,])source.Clone();
        if (diagonal) return Finish(work, vectors, 1);

        // A power of two avoids extra rounding when normal entries are scaled.
        // Subnormal scales use their finite maximum directly.
        var exponentBits = (BitConverter.DoubleToInt64Bits(scale) >> 52) & 0x7ff;
        if (exponentBits != 0) scale = BitConverter.Int64BitsToDouble(exponentBits << 52);
        var norm = 0d;
        for (var row = 0; row < n; row++)
            for (var column = 0; column < n; column++)
            {
                work[row, column] /= scale;
                norm = Hypot(norm, work[row, column]);
            }
        var threshold = tolerance * norm;
        var pairThreshold = threshold / (2d * n);

        for (var sweep = 0; ; sweep++)
        {
            var offNorm = 0d;
            for (var row = 0; row < n; row++)
                for (var column = row + 1; column < n; column++)
                    offNorm = Hypot(offNorm, Math.Sqrt(2) * work[row, column]);
            if (offNorm <= threshold) return Finish(work, vectors, scale);
            if (sweep == maxSweeps)
                throw new InvalidOperationException("Symmetric eigendecomposition did not converge within the sweep limit.");

            for (var p = 0; p < n - 1; p++)
                for (var q = p + 1; q < n; q++)
                {
                    var b = work[p, q];
                    if (Math.Abs(b) <= pairThreshold) continue;
                    var delta = (work[q, q] - work[p, p]) / 2;
                    var hypotenuse = Hypot(delta, b);
                    var t = b / (delta + (delta < 0 ? -hypotenuse : hypotenuse));
                    var c = 1 / Math.Sqrt(1 + t * t);
                    var s = t * c;
                    work[p, p] -= t * b;
                    work[q, q] += t * b;
                    work[p, q] = work[q, p] = 0;
                    for (var k = 0; k < n; k++)
                    {
                        if (k != p && k != q)
                        {
                            var x = work[k, p];
                            var y = work[k, q];
                            work[k, p] = work[p, k] = c * x - s * y;
                            work[k, q] = work[q, k] = s * x + c * y;
                        }
                        var vp = vectors[k, p];
                        var vq = vectors[k, q];
                        vectors[k, p] = c * vp - s * vq;
                        vectors[k, q] = s * vp + c * vq;
                    }
                }
        }
    }

    private static double Hypot(double a, double b)
    {
        a = Math.Abs(a);
        b = Math.Abs(b);
        var largest = Math.Max(a, b);
        if (largest == 0) return 0;
        var ratio = Math.Min(a, b) / largest;
        return largest * Math.Sqrt(1 + ratio * ratio);
    }

    private static SymmetricEigenDecomposition Finish(double[,] work, DoubleMatrix vectors, double scale)
    {
        var n = work.GetLength(0);
        var values = new double[n];
        var order = new int[n];
        for (var i = 0; i < n; i++)
        {
            values[i] = work[i, i] * scale;
            if (double.IsNaN(values[i]) || double.IsInfinity(values[i]))
                throw new ArithmeticException("An eigenvalue is outside the finite double range.");
            order[i] = i;
        }
        Array.Sort(order, (left, right) =>
        {
            var comparison = values[left].CompareTo(values[right]);
            return comparison != 0 ? comparison : left.CompareTo(right);
        });
        var sortedValues = new double[n];
        var sortedVectors = new DoubleMatrix(n, false);
        for (var column = 0; column < n; column++)
        {
            var originalColumn = order[column];
            sortedValues[column] = values[originalColumn];
            var lengthSquared = 0d;
            var largestRow = 0;
            for (var row = 0; row < n; row++)
            {
                var entry = vectors[row, originalColumn];
                lengthSquared += entry * entry;
                if (Math.Abs(entry) > Math.Abs(vectors[largestRow, originalColumn])) largestRow = row;
            }
            var length = Math.Sqrt(lengthSquared);
            var divisor = vectors[largestRow, originalColumn] < 0 ? -length : length;
            for (var row = 0; row < n; row++)
                sortedVectors[row, column] = vectors[row, originalColumn] / divisor;
        }
        return new SymmetricEigenDecomposition(sortedValues, sortedVectors);
    }
}
