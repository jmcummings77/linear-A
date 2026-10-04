using System;
using System.Linq;
using NUnit.Framework;
namespace linear_A.Test;

[TestFixture]
public class MultigridTests
{
    [TestCase(1)]
    [TestCase(3)]
    [TestCase(7)]
    public void ReusableCycleSolvesTheGeneratedOperator(int width)
    {
        var m = new GeometricMultigrid(width);
        var a = m.Matrix;
        var x = Enumerable.Range(0, m.Size).Select(i => Math.Sin(i + 1.0)).ToArray();
        var b = a.Matvec(x);
        var saved = (double[])b.Clone();
        for (var k = 0; k < 2; k++)
        {
            var r = a.ConjugateGradient(b, preconditioner: m, capture: true);
            Assert.That(r.Converged, Is.True);
            Assert.That(r.X, Is.EqualTo(x).Within(1e-8));
            Assert.That(b, Is.EqualTo(saved));
        }
        Assert.Throws<ArgumentException>(() => a.ConjugateGradient(b, jacobi: true, preconditioner: m));
        Assert.Throws<ArgumentException>(() => m.Apply([]));
    }
    [Test]
    public void InvalidGridsAndNonfiniteInputsAreRejected()
    {
        foreach (var width in new[] { -1, 0, 2, 4, 256 }) Assert.Throws<ArgumentException>(() => new GeometricMultigrid(width));
        Assert.Throws<ArgumentException>(() => new GeometricMultigrid(1).Apply([double.NaN]));
        Assert.That(new GeometricMultigrid(1).Apply([8]), Is.EqualTo(new double[] { 2 }));
    }
}
