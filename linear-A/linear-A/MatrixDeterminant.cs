using System;
using System.Numerics;

namespace linear_A;

internal static class MatrixDeterminant<T> where T : INumber<T>
{
    private static readonly bool IsInteger = typeof(T) == typeof(byte) || typeof(T) == typeof(sbyte)
        || typeof(T) == typeof(short) || typeof(T) == typeof(ushort)
        || typeof(T) == typeof(int) || typeof(T) == typeof(uint)
        || typeof(T) == typeof(long) || typeof(T) == typeof(ulong)
        || typeof(T) == typeof(nint) || typeof(T) == typeof(nuint)
        || typeof(T) == typeof(Int128) || typeof(T) == typeof(UInt128) || typeof(T) == typeof(BigInteger);
    private static readonly bool IsDecimal = typeof(T) == typeof(decimal);
    private static readonly bool IsBinaryFloat = typeof(T) == typeof(double)
        || typeof(T) == typeof(float) || typeof(T) == typeof(Half);

    internal static DeterminantAlgorithm Resolve(DeterminantAlgorithm algorithm)
    {
        if (algorithm == DeterminantAlgorithm.Auto)
            return IsInteger || IsDecimal ? DeterminantAlgorithm.Bareiss
                : IsBinaryFloat ? DeterminantAlgorithm.Lu : DeterminantAlgorithm.Cofactor;
        if (algorithm == DeterminantAlgorithm.Cofactor) return algorithm;
        if (algorithm == DeterminantAlgorithm.Lu || algorithm == DeterminantAlgorithm.Cholesky)
        {
            if (!IsBinaryFloat) throw new NotSupportedException("LU and Cholesky determinants support only double, float, and Half.");
            return algorithm;
        }
        if (algorithm == DeterminantAlgorithm.Bareiss)
        {
            if (!IsInteger && !IsDecimal) throw new NotSupportedException("Bareiss determinants support only built-in integers and decimal.");
            return algorithm;
        }
        throw new ArgumentOutOfRangeException(nameof(algorithm));
    }

    internal static T Lu(T[,] values)
    {
        var size = values.GetLength(0);
        if (size == 0) return T.One;
        // Reuse the established power-of-two scaling and separately accumulated
        // pivot exponents, including the protections for mixed huge/tiny rows.
        var work = new double[size, size];
        for (var row = 0; row < size; row++)
            for (var column = 0; column < size; column++)
                work[row, column] = double.CreateChecked(values[row, column]);
        return T.CreateChecked(DeterminantKernels.Lu(work, out _));
    }

    internal static T Bareiss(T[,] values)
    {
        var size = values.GetLength(0);
        var work = new BigInteger[size, size];
        var decimalScale = 0;
        for (var row = 0; row < size; row++)
        {
            var rowScale = 0;
            if (IsDecimal)
                for (var column = 0; column < size; column++)
                    rowScale = Math.Max(rowScale, decimal.GetBits(decimal.CreateChecked(values[row, column]))[3] >> 16 & 0xff);
            decimalScale += rowScale;
            for (var column = 0; column < size; column++)
            {
                if (!IsDecimal)
                {
                    work[row, column] = BigInteger.CreateChecked(values[row, column]);
                    continue;
                }
                var bits = decimal.GetBits(decimal.CreateChecked(values[row, column]));
                var coefficient = (BigInteger)(uint)bits[0] | ((BigInteger)(uint)bits[1] << 32) | ((BigInteger)(uint)bits[2] << 64);
                if (bits[3] < 0) coefficient = -coefficient;
                var scale = bits[3] >> 16 & 0xff;
                work[row, column] = coefficient * BigInteger.Pow(10, rowScale - scale);
            }
        }
        var determinant = DeterminantKernels.BareissInPlace(work);
        return IsDecimal ? T.CreateChecked(ToDecimal(determinant, decimalScale)) : T.CreateChecked(determinant);
    }

    internal static T Cholesky(T[,] values)
    {
        var size = values.GetLength(0);
        var work = new double[size, size];
        for (var row = 0; row < size; row++)
            for (var column = 0; column < size; column++)
            {
                var value = double.CreateChecked(values[row, column]);
                if (!double.IsFinite(value) || values[row, column] != values[column, row])
                    throw new InvalidOperationException("Cholesky requires a finite, exactly symmetric positive-definite matrix.");
                work[row, column] = value;
            }
        var mantissa = 1d;
        long exponent = 0;
        for (var row = 0; row < size; row++)
        {
            for (var column = 0; column <= row; column++)
            {
                var sum = work[row, column];
                for (var k = 0; k < column; k++) sum -= work[row, k] * work[column, k];
                if (!double.IsFinite(sum) || row == column && sum <= 0d)
                    throw new InvalidOperationException("Cholesky encountered a nonpositive or nonfinite pivot; the matrix is not numerically positive definite.");
                work[row, column] = row == column ? Math.Sqrt(sum) : sum / work[column, column];
                if (!double.IsFinite(work[row, column]))
                    throw new InvalidOperationException("Cholesky encountered a nonfinite factor.");
            }
            // det(A) = product(L[i,i]) squared. Feed each factor separately so
            // even representable mixed-scale results avoid an intermediate square.
            var factor = work[row, row];
            var factorExponent = Math.ILogB(factor);
            var fraction = Math.ScaleB(factor, -factorExponent);
            for (var repeat = 0; repeat < 2; repeat++)
            {
                mantissa *= fraction;
                exponent += factorExponent;
                if (mantissa >= 2d) { mantissa *= 0.5d; exponent++; }
            }
        }
        var determinant = exponent > 1023 ? double.PositiveInfinity
            : exponent < -1075 ? 0d : Math.ScaleB(mantissa, (int)exponent);
        return T.CreateChecked(determinant);
    }

    private static decimal ToDecimal(BigInteger coefficient, int scale)
    {
        var negative = coefficient.Sign < 0;
        coefficient = BigInteger.Abs(coefficient);
        var maximum = (BigInteger.One << 96) - 1;
        var removedDigits = Math.Max(0, scale - 28);
        var divisor = BigInteger.Pow(10, removedDigits);
        while (true)
        {
            var rounded = BigInteger.DivRem(coefficient, divisor, out var remainder);
            var midpoint = (remainder * 2).CompareTo(divisor);
            if (midpoint > 0 || midpoint == 0 && !rounded.IsEven) rounded++;
            if (rounded <= maximum)
                return new decimal(unchecked((int)(uint)(rounded & uint.MaxValue)),
                    unchecked((int)(uint)((rounded >> 32) & uint.MaxValue)),
                    unchecked((int)(uint)((rounded >> 64) & uint.MaxValue)), negative, (byte)(scale - removedDigits));
            if (removedDigits == scale) throw new OverflowException("The determinant exceeds the decimal range.");
            removedDigits++;
            divisor *= 10;
        }
    }
}
