using System;
using System.Linq;
using System.Numerics;
using NUnit.Framework;

namespace linear_A.Test;

public class SymmetricEigenTests
{
    private static DoubleMatrix Legacy(double[,] values)
    {
        var result = new DoubleMatrix(values.GetLength(0), values.GetLength(1));
        for (var i = 0; i < result.RowCount; i++)
            for (var j = 0; j < result.ColumnCount; j++) result[i, j] = values[i, j];
        return result;
    }

    private static void Verify(double[,] a, SymmetricEigenDecomposition result, double tolerance = 2e-10)
    {
        var n = a.GetLength(0);
        var values = result.EigenValues;
        var q = result.EigenVectors;
        Assert.That(values.Length, Is.EqualTo(n));
        Assert.That(q.RowCount, Is.EqualTo(n));
        Assert.That(q.ColumnCount, Is.EqualTo(n));
        Assert.That(values, Is.Ordered);
        var scale = 0d;
        foreach (var value in a) scale = Math.Max(scale, Math.Abs(value));
        if (scale == 0) scale = 1;
        var residualSquared = 0d;
        var reconstructionSquared = 0d;
        var inputSquared = 0d;
        for (var row = 0; row < n; row++)
            for (var column = 0; column < n; column++)
            {
                var aq = 0d;
                var qtq = 0d;
                var reconstructed = 0d;
                for (var k = 0; k < n; k++)
                {
                    aq += (a[row, k] / scale) * q[k, column];
                    qtq += q[k, row] * q[k, column];
                    reconstructed += q[row, k] * (values[k] / scale) * q[column, k];
                }
                var residual = aq - (values[column] / scale) * q[row, column];
                var difference = reconstructed - a[row, column] / scale;
                residualSquared += residual * residual;
                reconstructionSquared += difference * difference;
                inputSquared += (a[row, column] / scale) * (a[row, column] / scale);
                Assert.That(qtq, Is.EqualTo(row == column ? 1d : 0d).Within(tolerance), $"Q^TQ[{row},{column}]");
            }
        Assert.That(Math.Sqrt(residualSquared), Is.LessThanOrEqualTo(tolerance * Math.Max(1, Math.Sqrt(inputSquared))), "AQ = QΛ");
        Assert.That(Math.Sqrt(reconstructionSquared), Is.LessThanOrEqualTo(tolerance * Math.Max(1, Math.Sqrt(inputSquared))), "A = QΛQ^T");
        for (var column = 0; column < n; column++)
        {
            var largest = 0;
            for (var row = 0; row < n; row++)
                if (Math.Abs(q[row, column]) > Math.Abs(q[largest, column])) largest = row;
            Assert.That(q[largest, column], Is.GreaterThan(0));
        }
    }

    // A Householder reflection gives an independent, prescribed-spectrum oracle.
    private static double[,] WithSpectrum(double[] eigenvalues)
    {
        var n = eigenvalues.Length;
        var q = new double[n, n];
        var lengthSquared = Enumerable.Range(1, n).Sum(value => (double)value * value);
        for (var i = 0; i < n; i++)
            for (var j = 0; j < n; j++)
                q[i, j] = (i == j ? 1 : 0) - 2d * (i + 1) * (j + 1) / lengthSquared;
        var a = new double[n, n];
        for (var i = 0; i < n; i++)
            for (var j = i; j < n; j++)
            {
                for (var k = 0; k < n; k++) a[i, j] += q[i, k] * eigenvalues[k] * q[j, k];
                a[j, i] = a[i, j];
            }
        return a;
    }

    [Test]
    public void AnalyticPairsHaveSortedEigenvaluesAndOrthonormalColumns()
    {
        foreach (var a in new[] { new[,] { { 2d, 1d }, { 1d, 2d } }, new[,] { { 2d, -1d }, { -1d, 2d } } })
        {
            var generic = new Matrix<double>(a);
            var result = generic.GetSymmetricEigenDecomposition();
            Assert.That(result.EigenValues[0], Is.EqualTo(1d).Within(1e-14));
            Assert.That(result.EigenValues[1], Is.EqualTo(3d).Within(1e-14));
            Verify(a, result);
            Verify(a, Legacy(a).GetSymmetricEigenDecomposition());
            for (var row = 0; row < 2; row++)
                for (var column = 0; column < 2; column++) Assert.That(generic[row, column], Is.EqualTo(a[row, column]));
            result.EigenVectors[0, 0] = 999;
            Assert.That(generic[0, 0], Is.EqualTo(a[0, 0]));
        }
    }

    [Test]
    public void RepeatedNegativeZeroAndDensePrescribedSpectraReconstruct()
    {
        for (var n = 1; n <= 12; n++)
        {
            var expected = Enumerable.Range(0, n).Select(i => (double)(i / 2 - n / 4)).ToArray();
            var a = WithSpectrum(expected);
            var result = new Matrix<double>(a).GetSymmetricEigenDecomposition();
            Verify(a, result);
            for (var i = 0; i < n; i++) Assert.That(result.EigenValues[i], Is.EqualTo(expected[i]).Within(1e-10));
        }
    }

    [Test]
    public void EmptyZeroAndMixedScaleDiagonalSpectraRemainExact()
    {
        Verify(new double[0, 0], new Matrix<double>().GetSymmetricEigenDecomposition());
        Verify(new double[4, 4], new Matrix<double>(4).GetSymmetricEigenDecomposition());
        var a = new[,] { { 1e300, 0d, 0d }, { 0d, -1e-300, 0d }, { 0d, 0d, double.Epsilon } };
        var result = new Matrix<double>(a).GetSymmetricEigenDecomposition();
        Assert.That(result.EigenValues, Is.EqualTo(new[] { -1e-300, double.Epsilon, 1e300 }));
        Verify(a, result);
    }

    [TestCase(1e300)]
    [TestCase(1e-300)]
    [TestCase(1e-310)]
    public void ScalingSupportsHugeTinyAndSubnormalDenseMatrices(double scale)
    {
        var a = new[,] { { 2 * scale, scale }, { scale, 2 * scale } };
        var result = new Matrix<double>(a).GetSymmetricEigenDecomposition();
        Assert.That(result.EigenValues[0] / scale, Is.EqualTo(1d).Within(1e-12));
        Assert.That(result.EigenValues[1] / scale, Is.EqualTo(3d).Within(1e-12));
        Verify(a, result);
    }

    [Test]
    public void IntegerDecimalAndBigIntegerInputsReturnRealSpectra()
    {
        var a = new[,] { { 1d, 1d }, { 1d, 0d } };
        var integer = new IntegerMatrix(2, false);
        integer[0, 0] = integer[0, 1] = integer[1, 0] = 1;
        var results = new[] {
            integer.GetSymmetricEigenDecomposition(),
            new Matrix<int>(new[,] { { 1, 1 }, { 1, 0 } }).GetSymmetricEigenDecomposition(),
            new Matrix<decimal>(new[,] { { 1m, 1m }, { 1m, 0m } }).GetSymmetricEigenDecomposition(),
            new Matrix<BigInteger>(new[,] { { BigInteger.One, BigInteger.One }, { BigInteger.One, BigInteger.Zero } }).GetSymmetricEigenDecomposition(),
        };
        foreach (var result in results)
        {
            Verify(a, result);
            Assert.That(result.EigenValues[0], Is.EqualTo((1 - Math.Sqrt(5)) / 2).Within(1e-14));
            Assert.That(result.EigenValues[1], Is.EqualTo((1 + Math.Sqrt(5)) / 2).Within(1e-14));
        }
    }

    [Test]
    public void RejectsUnsupportedShapesNonfiniteAndAsymmetricInput()
    {
        Assert.Throws<NotSquareMatrixException>(() => new Matrix<double>(2, 3).GetSymmetricEigenDecomposition());
        Assert.Throws<ArgumentException>(() => new Matrix<double>(new[,] { { 1d, 2d }, { 3d, 4d } }).GetSymmetricEigenDecomposition());
        foreach (var value in new[] { double.NaN, double.PositiveInfinity, double.NegativeInfinity })
            Assert.Throws<ArgumentException>(() => new Matrix<double>(new[,] { { value } }).GetSymmetricEigenDecomposition());
        var large = BigInteger.One << 60;
        Assert.Throws<ArgumentException>(() => new Matrix<BigInteger>(new[,] { { large, large }, { large + 1, large } }).GetSymmetricEigenDecomposition());
        Assert.Throws<ArgumentException>(() => new Matrix<BigInteger>(new[,] { { BigInteger.One << 2000 } }).GetSymmetricEigenDecomposition());
        var overflowing = new Matrix<double>(new[,] { { double.MaxValue, double.MaxValue }, { double.MaxValue, double.MaxValue } });
        Assert.Throws<ArithmeticException>(() => overflowing.GetSymmetricEigenDecomposition());
    }

    [Test]
    public void ValidatesToleranceAndReportsNonconvergence()
    {
        var matrix = new Matrix<double>(WithSpectrum(new[] { -4d, 1d, 2d, 3d, 8d }));
        foreach (var tolerance in new[] { 0d, -1d, 1d, double.NaN, double.PositiveInfinity })
            Assert.Throws<ArgumentOutOfRangeException>(() => matrix.GetSymmetricEigenDecomposition(tolerance));
        Assert.Throws<ArgumentOutOfRangeException>(() => matrix.GetSymmetricEigenDecomposition(maxSweeps: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => matrix.GetSymmetricEigenDecomposition(maxSweeps: -1));
        Assert.Throws<InvalidOperationException>(() => matrix.GetSymmetricEigenDecomposition(maxSweeps: 1));
        Verify(WithSpectrum(new[] { -4d, 1d, 2d, 3d, 8d }), matrix.GetSymmetricEigenDecomposition());
    }

    [Test]
    public void TightToleranceDoesNotLoseOffDiagonalNormToUnderflow()
    {
        var a = new[,] { { 1d, 1e-200 }, { 1e-200, 2d } };
        var result = new Matrix<double>(a).GetSymmetricEigenDecomposition(1e-250, 1);
        Assert.That(result.EigenVectors[0, 1] / 1e-200, Is.EqualTo(1d).Within(1e-14));
        Assert.That(result.EigenVectors[1, 0] / 1e-200, Is.EqualTo(-1d).Within(1e-14));
        Verify(a, result);
    }

    [Test]
    public void LegacyDoubleHelpersHaveDocumentedLayoutAndTrySemantics()
    {
        var matrix = Legacy(new[,] { { 2d, 1d }, { 1d, 2d } });
        var decomposition = matrix.GetSymmetricEigenDecomposition();
        Assert.That(matrix.GetEigenValues(), Is.EqualTo(decomposition.EigenValues));
        var flattened = matrix.GetEigenVectors();
        for (var row = 0; row < 2; row++)
            for (var column = 0; column < 2; column++)
                Assert.That(flattened[row * 2 + column], Is.EqualTo(decomposition.EigenVectors[row, column]));
        Assert.That(matrix.TryGetEigenValues(out var values), Is.True);
        Assert.That(values, Is.EqualTo(decomposition.EigenValues));
        foreach (var invalid in new[] { new DoubleMatrix(2, 3), Legacy(new[,] { { 1d, 2d }, { 3d, 4d } }), Legacy(new[,] { { double.NaN } }) })
        {
            Assert.That(invalid.TryGetEigenValues(out values), Is.False);
            Assert.That(values, Is.Empty);
        }
    }

    [Test]
    public void LegacyIntegerSignaturesNeverSilentlyRoundRealResults()
    {
        var matrix = new IntegerMatrix(2, true);
#pragma warning disable CS0618 // Exercise the retained legacy signatures.
        Assert.Throws<NotSupportedException>(() => matrix.GetEigenValues());
        Assert.Throws<NotSupportedException>(() => matrix.GetEigenVectors());
        Assert.That(matrix.TryGetEigenValues(out var values), Is.False);
#pragma warning restore CS0618
        Assert.That(values, Is.Empty);
    }
}
