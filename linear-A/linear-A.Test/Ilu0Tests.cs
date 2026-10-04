using System;
using NUnit.Framework;

namespace linear_A.Test;

public class Ilu0Tests
{
    [Test]
    public void FactorOwnsStorageAndSupportsMultipleRightHandSides()
    {
        var values = new[] { 4.0, 1, 1, 3 };
        var a = new CSRMatrix(2, 2, new[] { 0, 2, 4 }, new[] { 0, 1, 0, 1 }, values);
        var f = new ILU0(a); values[0] = 99;
        var x = f.Apply(new[] { 6.0, 7 }); Assert.That(x[0], Is.EqualTo(1).Within(1e-12));
        x[0] = 99; var y = f.Apply(new[] { 11.0, 13 }); Assert.That(y[0], Is.EqualTo(20.0 / 11).Within(1e-12));
        Assert.That(a.Gmres(new[] { 6.0, 7 }, restart: 1, preconditioner: f).Iterations, Is.EqualTo(1));
        Assert.Throws<ArgumentException>(() => a.Gmres(new[] { 6.0, 7 }, jacobi: true, preconditioner: f));
    }
    [Test]
    public void ZeroPivotIsExplicitEvenForNonsingularMatrix()
    {
        var a = new CSRMatrix(2, 2, new[] { 0, 2, 4 }, new[] { 0, 1, 0, 1 }, new[] { 0.0, 1, 1, 2 });
        Assert.Throws<ArithmeticException>(() => new ILU0(a));
    }

    [Test]
    public void IncompleteFactorDropsUpdatesOutsideTheStoredPattern()
    {
        // A = [4 0 2; 2 3 0; 0 1 5]. Eliminating column 0 would add
        // U[1,2] = -1, which ILU(0) drops. The retained factors are
        // L = [1 0 0; 1/2 1 0; 0 1/3 1], U = [4 0 2; 0 3 0; 0 0 5].
        // Their product applied to [1,2,3] is [10,11,17].
        var a = new CSRMatrix(3, 3, [0, 2, 4, 6], [0, 2, 0, 1, 1, 2], [4, 2, 2, 3, 1, 5]);
        var factor = new ILU0(a);
        double[] rhs = [10, 11, 17];

        Assert.That(factor.NNZ, Is.EqualTo(6));
        Assert.That(factor.Apply(rhs), Is.EqualTo(new double[] { 1, 2, 3 }).Within(1e-12));
        Assert.That(rhs, Is.EqualTo(new double[] { 10, 11, 17 }));
        Assert.That(a.Values, Is.EqualTo(new double[] { 4, 2, 2, 3, 1, 5 }));
    }

    [Test]
    public void TinyNonzeroPivotIsAcceptedButSolveOverflowIsExplicit()
    {
        var factor = new ILU0(new CSRMatrix(1, 1, [0, 1], [0], [1e-320]));

        Assert.That(factor.Apply([1e-320]), Is.EqualTo(new double[] { 1 }));
        Assert.Throws<ArithmeticException>(() => factor.Apply([1]));
    }

    [TestCase(1e-320, 1, 1)]
    [TestCase(1, double.MaxValue, 2)]
    public void FactorizationRejectsDivisionOrUpdateOverflow(double pivot, double upper, double lower)
    {
        var a = new CSRMatrix(2, 2, [0, 2, 4], [0, 1, 0, 1], [pivot, upper, lower, 1]);

        Assert.Throws<ArithmeticException>(() => new ILU0(a));
    }
}
