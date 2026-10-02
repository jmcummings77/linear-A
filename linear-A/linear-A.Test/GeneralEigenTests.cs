using System;
using System.Numerics;
using NUnit.Framework;

namespace linear_A.Test;

public class GeneralEigenTests
{
    private static void Verify(double[,] source, EigenDecomposition result)
    {
        var n = source.GetLength(0);
        Assert.That(result.EigenValues.Length, Is.EqualTo(n));
        Assert.That(result.EigenVectors.GetLength(0), Is.EqualTo(n));
        Assert.That(result.EigenVectors.GetLength(1), Is.EqualTo(n));
        var scale = 0d;
        foreach (var value in source) scale = Math.Max(scale, Math.Abs(value));
        if (scale == 0) scale = 1;
        for (var col = 0; col < n; col++)
        {
            var lambda = result.EigenValues[col];
            Assert.That(double.IsFinite(lambda.Real) && double.IsFinite(lambda.Imaginary), Is.True);
            if (col > 0)
            {
                var previous = result.EigenValues[col - 1];
                Assert.That(previous.Real < lambda.Real || previous.Real == lambda.Real && previous.Imaginary <= lambda.Imaginary, Is.True);
            }
            var norm = 0d;
            var residual = 0d;
            for (var row = 0; row < n; row++)
            {
                var vector = result.EigenVectors[row, col];
                norm += vector.Magnitude * vector.Magnitude;
                var product = Complex.Zero;
                for (var k = 0; k < n; k++) product += (source[row, k] / scale) * result.EigenVectors[k, col];
                residual = double.Hypot(residual, (product - (lambda / scale) * vector).Magnitude);
            }
            Assert.That(norm, Is.EqualTo(1).Within(1e-12));
            Assert.That(residual, Is.LessThan(1e-10));
        }
    }

    [Test]
    public void RealInputCanReturnConjugateComplexEigenpairs()
    {
        var source = new[,] { { 0d, -1d }, { 1d, 0d } };
        var matrix = new Matrix<double>(source);
        var result = matrix.GetEigenDecomposition();
        Assert.That(result.EigenValues, Is.EqualTo(new[] { -Complex.ImaginaryOne, Complex.ImaginaryOne }));
        Verify(source, result);
        Verify(source, DoubleMatrix.FromArray(source).GetEigenDecomposition());
        Verify(source, IntegerMatrix.FromArray(new[,] { { 0, -1 }, { 1, 0 } }).GetEigenDecomposition());
        Verify(source, new Matrix<int>(new[,] { { 0, -1 }, { 1, 0 } }).GetEigenDecomposition());
        result.EigenVectors[0, 0] = 123;
        Assert.That(matrix[0, 0], Is.Zero);
        Assert.That(matrix[0, 1], Is.EqualTo(-1));
    }

    [Test]
    public void DefectiveMatricesReturnUnitEigenvectorsWithoutPromisingAnIndependentBasis()
    {
        foreach (var source in new[]
        {
            new[,] { { 2d, 1d }, { 0d, 2d } },
            new[,] { { 0d, 1d, 0d }, { 0d, 0d, 1d }, { 0d, 0d, 0d } },
            new[,] { { 2d, 0d }, { 1d, 2d } }
        }) Verify(source, new Matrix<double>(source).GetEigenDecomposition());
    }

    [Test]
    public void DenseSimilarityAndMixedSpectrumHaveSmallResiduals()
    {
        var source = new[,] { { 1d, -2d, 3d }, { 2d, 1d, 4d }, { 0d, 0d, -3d } };
        var result = new Matrix<double>(source).GetEigenDecomposition();
        Assert.That(result.EigenValues, Is.EqualTo(new[] { new Complex(-3, 0), new Complex(1, -2), new Complex(1, 2) }));
        Verify(source, result);
        var random = new Random(279);
        for (var n = 2; n < 12; n++)
        {
            var dense = new double[n, n];
            for (var i = 0; i < n; i++) for (var j = 0; j < n; j++) dense[i, j] = random.NextDouble() * 2 - 1;
            Verify(dense, new Matrix<double>(dense).GetEigenDecomposition());
        }
    }

    [Test]
    public void BalancingPreservesExtremeOpposingOffDiagonalEntries()
    {
        var source = new[,] { { 0d, 1e300 }, { -1e-300, 0d } };
        var result = new Matrix<double>(source).GetEigenDecomposition();
        Assert.That(result.EigenValues[0].Imaginary, Is.EqualTo(-1).Within(1e-14));
        Assert.That(result.EigenValues[1].Imaginary, Is.EqualTo(1).Within(1e-14));
        Verify(source, result);
        Assert.That(Math.Abs(result.EigenVectors[1, 0].Imaginary), Is.EqualTo(1e-300).Within(1e-314));
    }

    [Test]
    public void ScalingAndTriangularShortcutsPreserveExtremeSpectra()
    {
        foreach (var scale in new[] { 1e-200, 1e200 })
        {
            var source = new[,] { { scale, -2 * scale }, { 2 * scale, scale } };
            var result = new Matrix<double>(source).GetEigenDecomposition();
            Assert.That(result.EigenValues[0].Real / scale, Is.EqualTo(1).Within(1e-14));
            Assert.That(result.EigenValues[0].Imaginary / scale, Is.EqualTo(-2).Within(1e-14));
            Verify(source, result);
        }
        foreach (var source in new[]
        {
            new[,] { { 1e300, 0d }, { 0d, 1e-300 } },
            new[,] { { 1e300, 1d }, { 0d, 1e-300 } },
            new[,] { { 1e300, 0d }, { 1d, 1e-300 } }
        })
        {
            var result = new Matrix<double>(source).GetEigenDecomposition();
            Assert.That(result.EigenValues, Is.EqualTo(new[] { new Complex(1e-300, 0), new Complex(1e300, 0) }));
            Verify(source, result);
        }
        var tiny = new Matrix<double>(new[,] { { double.Epsilon } }).GetEigenDecomposition();
        Assert.That(tiny.EigenValues[0].Real, Is.EqualTo(double.Epsilon));
    }

    [Test]
    public void EmptyAndZeroInputsAreValid()
    {
        Verify(new double[0, 0], new Matrix<double>().GetEigenDecomposition());
        Verify(new double[0, 0], new DoubleMatrix().GetEigenDecomposition());
        Verify(new double[3, 3], new Matrix<double>(3, 3).GetEigenDecomposition());
    }

    [Test]
    public void InvalidInputsAndIterationBoundsFailExplicitly()
    {
        Assert.Throws<NotSquareMatrixException>(() => new Matrix<double>(2, 3).GetEigenDecomposition());
        foreach (var value in new[] { double.NaN, double.PositiveInfinity, double.NegativeInfinity })
            Assert.Throws<ArgumentException>(() => new Matrix<double>(new[,] { { value } }).GetEigenDecomposition());
        foreach (var limit in new[] { -1, 0, 100001 })
            Assert.Throws<ArgumentOutOfRangeException>(() => new Matrix<double>(1, 1).GetEigenDecomposition(limit));
        var difficult = new Matrix<double>(new[,] { { 1d, 2d, 3d }, { 4d, 5d, 6d }, { 7d, 8d, 10d } });
        Assert.Throws<InvalidOperationException>(() => difficult.GetEigenDecomposition(1));
        Assert.That(difficult[2, 2], Is.EqualTo(10));
    }
}
