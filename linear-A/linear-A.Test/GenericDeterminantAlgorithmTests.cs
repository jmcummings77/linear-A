using System;
using System.Numerics;
using NUnit.Framework;

namespace linear_A.Test;

public class GenericDeterminantAlgorithmTests
{
    private static Matrix<T> Convert<T>(int[,] values) where T : INumber<T>
    {
        var matrix = new Matrix<T>(values.GetLength(0), values.GetLength(1));
        for (var row = 0; row < matrix.RowCount; row++)
            for (var col = 0; col < matrix.ColumnCount; col++)
                matrix[row, col] = T.CreateChecked(values[row, col]);
        return matrix;
    }

    // Independent Leibniz formula: enumerate permutations and count inversions.
    private static BigInteger PermutationDeterminant(int[,] values)
    {
        var size = values.GetLength(0);
        var permutation = new int[size];
        var used = new bool[size];
        var result = BigInteger.Zero;
        void Visit(int row)
        {
            if (row == size)
            {
                var product = BigInteger.One;
                var inversions = 0;
                for (var i = 0; i < size; i++)
                {
                    product *= values[i, permutation[i]];
                    for (var j = i + 1; j < size; j++)
                        if (permutation[i] > permutation[j]) inversions++;
                }
                result += inversions % 2 == 0 ? product : -product;
                return;
            }
            for (var col = 0; col < size; col++)
            {
                if (used[col]) continue;
                used[col] = true;
                permutation[row] = col;
                Visit(row + 1);
                used[col] = false;
            }
        }
        Visit(0);
        return result;
    }

    [Test]
    public void AllAlgorithmsMatchIndependentExactOracleForSeededSizesZeroThroughSix()
    {
        var random = new Random(8201);
        for (var size = 0; size <= 6; size++)
            for (var sample = 0; sample < 8; sample++)
            {
                var values = new int[size, size];
                for (var row = 0; row < size; row++)
                    for (var col = 0; col < size; col++)
                        values[row, col] = random.Next(-3, 4);
                var expected = PermutationDeterminant(values);
                var integers = Convert<int>(values);
                var big = Convert<BigInteger>(values);
                var decimals = Convert<decimal>(values);
                var floats = Convert<double>(values);
                foreach (var algorithm in new[] { DeterminantAlgorithm.Auto, DeterminantAlgorithm.Bareiss, DeterminantAlgorithm.Cofactor })
                {
                    Assert.That(integers.GetDeterminant(algorithm), Is.EqualTo((int)expected), $"int {size}/{sample}/{algorithm}");
                    Assert.That(big.GetDeterminant(algorithm), Is.EqualTo(expected), $"BigInteger {size}/{sample}/{algorithm}");
                    Assert.That(decimals.GetDeterminant(algorithm), Is.EqualTo((decimal)expected), $"decimal {size}/{sample}/{algorithm}");
                }
                foreach (var algorithm in new[] { DeterminantAlgorithm.Auto, DeterminantAlgorithm.Lu, DeterminantAlgorithm.Cofactor })
                    Assert.That(floats.GetDeterminant(algorithm), Is.EqualTo((double)expected).Within(1e-9 * Math.Max(1, (double)BigInteger.Abs(expected))), $"double {size}/{sample}/{algorithm}");
                for (var row = 0; row < size; row++)
                    for (var col = 0; col < size; col++)
                    {
                        Assert.That(big[row, col], Is.EqualTo((BigInteger)values[row, col]));
                        Assert.That(floats[row, col], Is.EqualTo((double)values[row, col]));
                    }
            }
    }

    [Test]
    public void BareissAvoidsBoundedIntermediateOverflowButChecksItsFinalIntegerResult()
    {
        var n = int.MaxValue;
        var matrix = new Matrix<int>(new[,] { { n, n - 1 }, { n - 1, n - 2 } });
        Assert.That(matrix.GetDeterminant(), Is.EqualTo(-1));
        Assert.Throws<OverflowException>(() => matrix.GetDeterminant(DeterminantAlgorithm.Cofactor));
        Assert.Throws<OverflowException>(() => new Matrix<int>(new[,] { { n, 0 }, { 0, 2 } }).GetDeterminant());
        Assert.Throws<OverflowException>(() => new Matrix<uint>(new uint[,] { { 0, 1 }, { 1, 0 } }).GetDeterminant());
        Assert.That(new Matrix<uint>(new uint[,] { { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 0 } }).GetDeterminant(), Is.EqualTo(1u));
    }

    [Test]
    public void AllBuiltInIntegerDomainsUseExactElimination()
    {
        void Check<T>() where T : INumber<T>
        {
            var matrix = Convert<T>(new[,] { { 2, 3 }, { 1, 2 } });
            Assert.That(matrix.GetDeterminant(), Is.EqualTo(T.One));
            Assert.That(matrix.GetDeterminant(DeterminantAlgorithm.Bareiss), Is.EqualTo(T.One));
        }
        Check<byte>(); Check<sbyte>(); Check<short>(); Check<ushort>(); Check<int>(); Check<uint>();
        Check<long>(); Check<ulong>(); Check<nint>(); Check<nuint>(); Check<Int128>(); Check<UInt128>(); Check<BigInteger>();
    }

    [Test]
    public void BigIntegerEliminationRetainsHugeCancellationAndPivotSigns()
    {
        var n = BigInteger.Pow(10, 150);
        var matrix = new Matrix<BigInteger>(new[,] { { n + 1, n, BigInteger.Zero }, { n, n - 1, BigInteger.Zero }, { BigInteger.Zero, BigInteger.Zero, BigInteger.One } });
        Assert.That(matrix.GetDeterminant(), Is.EqualTo(BigInteger.MinusOne));
        var swapped = new Matrix<BigInteger>(new[,] { { BigInteger.Zero, BigInteger.Zero, BigInteger.One }, { n, n - 1, BigInteger.Zero }, { n + 1, n, BigInteger.Zero } });
        Assert.That(swapped.GetDeterminant(), Is.EqualTo(BigInteger.One));
    }

    [Test]
    public void DecimalEliminationRetainsFractionsAndCancellationBeyondItsStoragePrecision()
    {
        var n = decimal.MaxValue;
        var matrix = new Matrix<decimal>(new[,] { { n, n - 1 }, { n - 1, n - 2 } });
        Assert.That(matrix.GetDeterminant(), Is.EqualTo(-1m));
        Assert.Throws<OverflowException>(() => matrix.GetDeterminant(DeterminantAlgorithm.Cofactor));
        var scaled = 1000000000000000000000000000.1m;
        Assert.That(new Matrix<decimal>(new[,] { { scaled, scaled - 0.1m }, { scaled - 0.1m, scaled - 0.2m } }).GetDeterminant(), Is.EqualTo(-0.01m));
        var fractions = new Matrix<decimal>(new[,] { { 0.1m, 0.2m, 0.3m }, { -0.25m, 0.5m, 0.75m }, { 1.5m, 0.25m, -0.5m } });
        Assert.That(fractions.GetDeterminant(), Is.EqualTo(fractions.GetDeterminant(DeterminantAlgorithm.Cofactor)));
        Assert.That(new Matrix<decimal>(new[,] { { n, 0m }, { 0m, 1m } }).GetDeterminant(), Is.EqualTo(n));
        Assert.Throws<OverflowException>(() => new Matrix<decimal>(new[,] { { n, 0m }, { 0m, 2m } }).GetDeterminant());
    }

    [Test]
    public void DecimalFinalRoundingUsesNearestTiesToEvenIncludingNegativeUnderflow()
    {
        var unit = 0.0000000000000000000000000001m;
        foreach (var sign in new[] { -1m, 1m })
            foreach (var (coefficient, expected) in new[] { (1m, 0m), (3m, 2m), (5m, 2m), (7m, 4m) })
            {
                var matrix = new Matrix<decimal>(new[,] { { sign * coefficient * unit, 0m }, { 0m, 0.5m } });
                Assert.That(matrix.GetDeterminant(), Is.EqualTo(sign * expected * unit));
            }
        // Precision is reduced from the original exact quotient, not by rounding
        // repeatedly: .49 must not become .5 and then round the odd maximum up.
        var max = decimal.MaxValue;
        var boundary = new Matrix<decimal>(new[,] { { max / 10m, 0m }, { 0m, 10m } });
        Assert.That(boundary.GetDeterminant(), Is.EqualTo(max));
        foreach (var sign in new[] { -1m, 1m })
        {
            var belowHalf = new Matrix<decimal>(new[,] { { sign * max, -sign * 0.49m }, { 1m, 1m } });
            var aboveHalf = new Matrix<decimal>(new[,] { { sign * max, -sign * 0.51m }, { 1m, 1m } });
            Assert.That(belowHalf.GetDeterminant(), Is.EqualTo(sign * max));
            Assert.Throws<OverflowException>(() => aboveHalf.GetDeterminant());
        }
    }

    [Test]
    public void LuPreservesMixedScalesAndUsesDoubleWorkingPrecisionForSmallerFloats()
    {
        var mixed = new Matrix<double>(new[,] { { 1e308, 1e308 }, { 1e-308, 2e-308 } });
        Assert.That(mixed.GetDeterminant(), Is.EqualTo(1d).Within(1e-12));
        Assert.That(mixed.Transpose().GetDeterminant(), Is.EqualTo(1d).Within(1e-12));
        var diagonal = new Matrix<double>(new[,] { { 1e-200, 0d, 0d, 0d }, { 0d, 1e-200, 0d, 0d }, { 0d, 0d, 1e200, 0d }, { 0d, 0d, 0d, 1e200 } });
        Assert.That(diagonal.GetDeterminant(), Is.EqualTo(1d).Within(1e-12));
        Assert.That(new Matrix<float>(new[,] { { 1e30f, 0f }, { 0f, 1e-30f } }).GetDeterminant(), Is.EqualTo(1f).Within(1e-6f));
        Assert.That(new Matrix<Half>(new[,] { { (Half)2, (Half)3 }, { (Half)1, (Half)2 } }).GetDeterminant(), Is.EqualTo((Half)1));
        Assert.That(new Matrix<float>(new[,] { { float.MaxValue, 0f }, { 0f, 2f } }).GetDeterminant(), Is.EqualTo(float.PositiveInfinity));
        Assert.That(new Matrix<float>(new[,] { { float.Epsilon, 0f }, { 0f, 0.5f } }).GetDeterminant(), Is.Zero);
        Assert.That(new Matrix<Half>(new[,] { { Half.MaxValue, (Half)0 }, { (Half)0, (Half)2 } }).GetDeterminant(), Is.EqualTo(Half.PositiveInfinity));
        Assert.That(new Matrix<Half>(new[,] { { Half.Epsilon, (Half)0 }, { (Half)0, (Half)0.5 } }).GetDeterminant(), Is.EqualTo((Half)0));
        Assert.That(new Matrix<double>(new[,] { { 1e-200, 1e-200 }, { 3e-124, 6e-124 } }).GetDeterminant(), Is.EqualTo(double.Epsilon));
        var unscalable = new Matrix<double>(new[,] { { 1e308, 1e308, 1e-308 }, { 1e-308, 2e-308, 1e308 }, { 0d, 0d, 1e-308 } });
        Assert.That(unscalable.GetDeterminant() / 1e-308, Is.EqualTo(1d).Within(1e-12));
    }

    [Test]
    public void CholeskyMatchesSpdDeterminantsAndPreservesInputsAndExtremeScales()
    {
        var matrix = new Matrix<double>(new[,] { { 4d, 2d }, { 2d, 3d } });
        Assert.That(matrix.GetDeterminant(DeterminantAlgorithm.Cholesky), Is.EqualTo(8d).Within(1e-12));
        Assert.That(matrix.GetRow(0), Is.EqualTo(new[] { 4d, 2d }));
        Assert.That(matrix.GetRow(1), Is.EqualTo(new[] { 2d, 3d }));
        Assert.That(new Matrix<double>().GetDeterminant(DeterminantAlgorithm.Cholesky), Is.EqualTo(1d));
        var mixed = new Matrix<double>(new[,] { { 1e300, 0d, 0d, 0d }, { 0d, 1e300, 0d, 0d }, { 0d, 0d, 1e-300, 0d }, { 0d, 0d, 0d, 1e-300 } });
        Assert.That(mixed.GetDeterminant(DeterminantAlgorithm.Cholesky), Is.EqualTo(1d).Within(1e-12));
        Assert.That(new Matrix<double>(new[,] { { double.Epsilon } }).GetDeterminant(DeterminantAlgorithm.Cholesky), Is.EqualTo(double.Epsilon));
        Assert.That(new Matrix<float>(new[,] { { 4f, 2f }, { 2f, 3f } }).GetDeterminant(DeterminantAlgorithm.Cholesky), Is.EqualTo(8f).Within(1e-6f));
        Assert.That(new Matrix<Half>(new[,] { { (Half)4, (Half)2 }, { (Half)2, (Half)3 } }).GetDeterminant(DeterminantAlgorithm.Cholesky), Is.EqualTo((Half)8));
    }

    [Test]
    public void CholeskyRejectsNonsymmetricIndefiniteSingularAndNonfiniteInputs()
    {
        foreach (var values in new[] {
            new[,] { { 4d, 2d }, { 1d, 3d } }, new[,] { { 1d, 2d }, { 2d, 1d } },
            new[,] { { 1d, 1d }, { 1d, 1d } }, new[,] { { -1d, 0d }, { 0d, -1d } },
            new[,] { { 1d, double.NaN }, { double.NaN, 1d } }, new[,] { { double.PositiveInfinity, 0d }, { 0d, 1d } },
        })
            Assert.Throws<InvalidOperationException>(() => new Matrix<double>(values).GetDeterminant(DeterminantAlgorithm.Cholesky));
        Assert.Throws<NotSupportedException>(() => new Matrix<int>().GetDeterminant(DeterminantAlgorithm.Cholesky));
        Assert.Throws<NotSupportedException>(() => new Matrix<decimal>().GetDeterminant(DeterminantAlgorithm.Cholesky));
    }

    [Test]
    public void NonfiniteLuInputIsNaNEvenAwayFromTheDiagonal()
    {
        foreach (var value in new[] { double.NaN, double.PositiveInfinity, double.NegativeInfinity })
        {
            var matrix = new Matrix<double>(new[,] { { 1d, value }, { 0d, 2d } });
            Assert.That(matrix.GetDeterminant(), Is.NaN);
            Assert.That(matrix.TryGetDeterminant(DeterminantAlgorithm.Lu, out var result), Is.True);
            Assert.That(result, Is.NaN);
        }
    }

    [Test]
    public void SelectionRejectsUnknownAlgorithmsAndUnsupportedDomainsBeforeSizeShortcuts()
    {
        Assert.Throws<NotSupportedException>(() => new Matrix<int>().GetDeterminant(DeterminantAlgorithm.Lu));
        Assert.Throws<NotSupportedException>(() => new Matrix<decimal>().GetDeterminant(DeterminantAlgorithm.Lu));
        Assert.Throws<NotSupportedException>(() => new Matrix<double>().GetDeterminant(DeterminantAlgorithm.Bareiss));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Matrix<int>().GetDeterminant((DeterminantAlgorithm)999));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Matrix<int>(2, 3).TryGetDeterminant((DeterminantAlgorithm)999, out _));
        var rectangle = new Matrix<int>(2, 3);
        Assert.That(rectangle.TryGetDeterminant(DeterminantAlgorithm.Bareiss, out var determinant), Is.False);
        Assert.That(determinant, Is.Zero);
        Assert.Throws<NotSquareMatrixException>(() => rectangle.GetDeterminant(DeterminantAlgorithm.Cofactor));
    }

    [Test]
    public void OtherGenericMathProvidersKeepTheCofactorFallback()
    {
        // NFloat supplies INumber but is not one of the specialized T domains.
        var matrix = new Matrix<System.Runtime.InteropServices.NFloat>(new[,]
        {
            { new System.Runtime.InteropServices.NFloat(1.5), new System.Runtime.InteropServices.NFloat(2) },
            { new System.Runtime.InteropServices.NFloat(0.25), new System.Runtime.InteropServices.NFloat(3) },
        });
        Assert.That(matrix.GetDeterminant(), Is.EqualTo(new System.Runtime.InteropServices.NFloat(4)));
        Assert.That(matrix.GetDeterminant(), Is.EqualTo(matrix.GetDeterminant(DeterminantAlgorithm.Cofactor)));
        Assert.Throws<NotSupportedException>(() => matrix.GetDeterminant(DeterminantAlgorithm.Lu));
        Assert.Throws<NotSupportedException>(() => matrix.GetDeterminant(DeterminantAlgorithm.Bareiss));
    }
}
