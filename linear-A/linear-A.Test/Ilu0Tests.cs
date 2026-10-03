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
}
