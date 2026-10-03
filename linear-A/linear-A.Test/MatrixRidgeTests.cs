using System;
using NUnit.Framework;

namespace linear_A.Test;

public class MatrixRidgeTests
{
    [Test]
    public void PenaltyShrinksSolutionAndDoesNotModifyInputs()
    {
        var a = new Matrix<double>(new double[,] { { 2 } });
        var b = new Matrix<double>(new double[,] { { 3 } });
        Assert.That(a.SolveRidge(b, 1)[0, 0], Is.EqualTo(1.2).Within(1e-14));
        Assert.That(a.SolveRidge(b, 0)[0, 0], Is.EqualTo(a.SolveMinimumNorm(b)[0, 0]));
        Assert.That(a[0, 0], Is.EqualTo(2));
        Assert.That(b[0, 0], Is.EqualTo(3));
        foreach (var lambda in new[] { -1.0, double.NaN, double.PositiveInfinity })
            Assert.Throws<ArgumentOutOfRangeException>(() => a.SolveRidge(b, lambda));
        Assert.Throws<ArgumentNullException>(() => a.SolveRidge(null!, 1));
        b[0, 0] = double.NaN;
        Assert.Throws<ArgumentException>(() => a.SolveRidge(b, 1));
    }

    [Test]
    public void FilterDoesNotUnderflowBeforeMultiplyingLargeRightHandSide()
    {
        var a = new Matrix<double>(new double[,] { { 1e-300 } });
        var b = new Matrix<double>(new double[,] { { 1e300 } });
        Assert.That(a.SolveRidge(b, 1e300)[0, 0] / 1e-300, Is.EqualTo(1).Within(1e-14));
    }
}
