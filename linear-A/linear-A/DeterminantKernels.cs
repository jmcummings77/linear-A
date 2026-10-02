using System;
using System.Numerics;

namespace linear_A;

/// <summary>Determinant arithmetic shared by generic and compatibility matrix types.</summary>
/// <remarks>Public callers validate square dimensions before entering these kernels.</remarks>
internal static class DeterminantKernels
{
    /// <summary>Copies the input and performs scaled partial-pivot LU elimination.</summary>
    internal static double Lu(double[,] matrix, out bool isInvertible)
    {
        var size = matrix.GetLength(0);
        var work = (double[,])matrix.Clone();
        isInvertible = false;

        for (var row = 0; row < size; row++)
            for (var column = 0; column < size; column++)
                if (double.IsNaN(work[row, column]) || double.IsInfinity(work[row, column]))
                    return double.NaN;

        long rowScaleExponent = 0;
        for (var row = 0; row < size; row++)
        {
            var largest = 0d;
            for (var column = 0; column < size; column++)
                largest = Math.Max(largest, Math.Abs(work[row, column]));
            if (largest == 0d) return 0d;

            var encodedExponent = (int)((BitConverter.DoubleToInt64Bits(largest) >> 52) & 0x7ff);
            if (encodedExponent == 0)
            {
                // Normalize subnormal rows before constructing their reciprocal scale.
                for (var column = 0; column < size; column++) work[row, column] *= 18014398509481984d;
                largest *= 18014398509481984d;
                encodedExponent = (int)((BitConverter.DoubleToInt64Bits(largest) >> 52) & 0x7ff);
                rowScaleExponent -= 54;
            }

            var exponent = encodedExponent - 1023;
            var scale = Math.Pow(2d, -exponent);
            var preservesEntries = true;
            for (var column = 0; column < size && preservesEntries; column++)
                preservesEntries = work[row, column] * scale / scale == work[row, column];
            if (!preservesEntries) continue;

            for (var column = 0; column < size; column++) work[row, column] *= scale;
            rowScaleExponent += exponent;
        }

        var sign = 1;
        for (var pivotColumn = 0; pivotColumn < size; pivotColumn++)
        {
            var pivotRow = pivotColumn;
            for (var row = pivotColumn + 1; row < size; row++)
                if (Math.Abs(work[row, pivotColumn]) > Math.Abs(work[pivotRow, pivotColumn]))
                    pivotRow = row;

            var pivot = work[pivotRow, pivotColumn];
            if (double.IsNaN(pivot) || double.IsInfinity(pivot)) return double.NaN;
            if (pivot == 0d) return 0d;

            if (pivotRow != pivotColumn)
            {
                for (var column = pivotColumn; column < size; column++)
                {
                    var temporary = work[pivotColumn, column];
                    work[pivotColumn, column] = work[pivotRow, column];
                    work[pivotRow, column] = temporary;
                }
                sign = -sign;
            }

            for (var row = pivotColumn + 1; row < size; row++)
            {
                var numerator = work[row, pivotColumn];
                var factor = numerator / pivot;
                var scaledUpdate = numerator != 0d && Math.Abs(factor) < 2.2250738585072014e-308;
                work[row, pivotColumn] = 0d;
                for (var column = pivotColumn + 1; column < size; column++)
                {
                    var product = scaledUpdate
                        ? MultiplyQuotient(numerator, pivot, work[pivotColumn, column])
                        : factor * work[pivotColumn, column];
                    work[row, column] -= product;
                    if (double.IsNaN(work[row, column]) || double.IsInfinity(work[row, column]))
                        return double.NaN;
                }
            }
        }

        isInvertible = true;
        return MultiplyDiagonal(work, sign, rowScaleExponent);
    }

    private static double MultiplyQuotient(double numerator, double denominator, double value)
    {
        if (value == 0d) return 0d;
        // A tiny quotient can still have a representable product. Keep the
        // binary exponents separate instead of rounding the factor to zero.
        var a = NormalizedMantissa(numerator, out var aExponent);
        var b = NormalizedMantissa(denominator, out var bExponent);
        var c = NormalizedMantissa(value, out var cExponent);
        var mantissa = a / b * c;
        var exponent = aExponent - bExponent + cExponent;
        if (Math.Abs(mantissa) >= 2d) { mantissa *= 0.5d; exponent++; }
        else if (Math.Abs(mantissa) < 1d) { mantissa *= 2d; exponent--; }
        if (exponent > 1023) return mantissa > 0d ? double.PositiveInfinity : double.NegativeInfinity;
        if (exponent < -1075) return mantissa * 0d;
        if (exponent < -1022)
            return mantissa * Math.Pow(2d, exponent + 1022) * Math.Pow(2d, -1022);
        return mantissa * BitConverter.Int64BitsToDouble((long)(exponent + 1023) << 52);
    }

    private static double NormalizedMantissa(double value, out int exponent)
    {
        var bits = BitConverter.DoubleToInt64Bits(value);
        var encodedExponent = (int)((bits >> 52) & 0x7ff);
        exponent = 0;
        if (encodedExponent == 0)
        {
            value *= 18014398509481984d;
            bits = BitConverter.DoubleToInt64Bits(value);
            encodedExponent = (int)((bits >> 52) & 0x7ff);
            exponent -= 54;
        }
        exponent += encodedExponent - 1023;
        return BitConverter.Int64BitsToDouble((bits & unchecked((long)0x800fffffffffffffUL)) | (1023L << 52));
    }

    private static double MultiplyDiagonal(double[,] matrix, int sign, long exponent)
    {
        var mantissa = (double)sign;
        for (var i = 0; i < matrix.GetLength(0); i++)
        {
            var value = matrix[i, i];
            var bits = BitConverter.DoubleToInt64Bits(value);
            var valueExponent = (int)((bits >> 52) & 0x7ff);
            if (valueExponent == 0)
            {
                // Scaling by 2^54 normalizes any nonzero subnormal double exactly.
                value *= 18014398509481984d;
                bits = BitConverter.DoubleToInt64Bits(value);
                valueExponent = (int)((bits >> 52) & 0x7ff);
                exponent -= 54;
            }

            exponent += valueExponent - 1023;
            bits = (bits & unchecked((long)0x800fffffffffffffUL)) | (1023L << 52);
            mantissa *= BitConverter.Int64BitsToDouble(bits);
            if (Math.Abs(mantissa) >= 2d)
            {
                mantissa *= 0.5d;
                exponent++;
            }
        }

        if (exponent > 1023) return mantissa > 0d ? double.PositiveInfinity : double.NegativeInfinity;
        if (exponent < -1075) return mantissa * 0d;
        if (exponent < -1022)
            return mantissa * Math.Pow(2d, exponent + 1022) * Math.Pow(2d, -1022);

        return mantissa * BitConverter.Int64BitsToDouble((exponent + 1023) << 52);
    }

    /// <summary>Copies integer input and computes its exact determinant.</summary>
    internal static BigInteger Bareiss(int[,] matrix)
    {
        var size = matrix.GetLength(0);
        if (size == 0) return BigInteger.One;
        if (size == 1) return matrix[0, 0];
        var work = new BigInteger[size, size];
        for (var row = 0; row < size; row++)
            for (var column = 0; column < size; column++)
                work[row, column] = matrix[row, column];
        return BareissInPlace(work);
    }

    /// <summary>Computes an exact determinant, consuming a caller-owned square workspace.</summary>
    internal static BigInteger BareissInPlace(BigInteger[,] work)
    {
        var size = work.GetLength(0);
        if (size == 0) return BigInteger.One;
        var previousPivot = BigInteger.One;
        var sign = 1;
        for (var pivotColumn = 0; pivotColumn < size - 1; pivotColumn++)
        {
            var pivotRow = pivotColumn;
            while (pivotRow < size && work[pivotRow, pivotColumn].IsZero) pivotRow++;
            if (pivotRow == size) return BigInteger.Zero;

            if (pivotRow != pivotColumn)
            {
                for (var column = pivotColumn; column < size; column++)
                {
                    var temporary = work[pivotColumn, column];
                    work[pivotColumn, column] = work[pivotRow, column];
                    work[pivotRow, column] = temporary;
                }
                sign = -sign;
            }

            var pivot = work[pivotColumn, pivotColumn];
            for (var row = pivotColumn + 1; row < size; row++)
            {
                for (var column = pivotColumn + 1; column < size; column++)
                    work[row, column] =
                        (pivot * work[row, column] - work[row, pivotColumn] * work[pivotColumn, column]) /
                        previousPivot;
                work[row, pivotColumn] = BigInteger.Zero;
            }
            previousPivot = pivot;
        }

        return sign * work[size - 1, size - 1];
    }
}
