using System;
using NUnit.Framework;

namespace linear_A.Test;

public class MatrixPseudoinverseTests
{
    [Test]
    public void StrictCutoffControlsRankAndMultipleRightHandSides()
    {
        var a = new Matrix<double>(new double[,] { { 4, 0 }, { 0, 1 } });
        var b = new Matrix<double>(new double[,] { { 4, 8 }, { 1, 2 } });
        var p = a.Pseudoinverse(.25);
        var x = a.SolveMinimumNorm(b, .25);
        Assert.That(p[0, 0], Is.EqualTo(.25));
        Assert.That(p[1, 1], Is.Zero);
        Assert.That(x[0, 1], Is.EqualTo(2));
        Assert.That(x[1, 1], Is.Zero);
        Assert.That(a.SpectralDiagnostics(.25), Is.EqualTo(new SpectralDiagnostics(1, .25, 1)));
        a[0, 0] = 8;
        Assert.That(p[0, 0], Is.EqualTo(.25));
        Assert.That(b[0, 0], Is.EqualTo(4));
        Assert.Throws<ArgumentOutOfRangeException>(() => a.Pseudoinverse(double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => a.SpectralDiagnostics(1.1));
        b[0, 0] = double.NaN;
        Assert.Throws<ArgumentException>(() => a.SolveMinimumNorm(b, 1));
    }

    [Test]
    public void DirectSolveDoesNotRequireARepresentableInverse()
    {
        var a = new Matrix<double>(new double[,] { { 1e-310 } });
        Assert.Throws<ArithmeticException>(() => a.Pseudoinverse(0));
        Assert.That(a.SolveMinimumNorm(a, 0)[0, 0], Is.EqualTo(1).Within(1e-14));
    }
}
