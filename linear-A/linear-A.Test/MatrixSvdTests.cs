using System;
using NUnit.Framework;

namespace linear_A.Test;

public class MatrixSvdTests
{
    [Test]
    public void EconomyFactorsOwnStorageAndRejectInvalidOptions()
    {
        var a = new Matrix<double>(new double[,] { { 3, 0, 0 }, { 0, 4, 0 } });
        var result = a.Svd();
        Assert.That(result.Values, Is.EqualTo(new[] { 4.0, 3.0 }));
        Assert.That(result.U.RowCount, Is.EqualTo(2));
        Assert.That(result.Vt.ColumnCount, Is.EqualTo(3));
        result.U[0, 0] = 99; result.Values[0] = 99;
        Assert.That(a[0, 0], Is.EqualTo(3));
        Assert.Throws<ArgumentOutOfRangeException>(() => a.Svd(tolerance: double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => a.Svd(maxSweeps: 0));
        a[0, 0] = double.PositiveInfinity;
        Assert.Throws<ArgumentException>(() => a.Svd());
    }

    [Test]
    public void ExhaustedIterationsFailExplicitly()
    {
        var a = new Matrix<double>(new double[,] { { 1, 2, 3 }, { 4, 5, 7 }, { 6, 8, 9 } });
        Assert.Throws<ArithmeticException>(() => a.Svd(maxSweeps: 1));
        Assert.That(a.Svd().Values.Length, Is.EqualTo(3));
    }
}
