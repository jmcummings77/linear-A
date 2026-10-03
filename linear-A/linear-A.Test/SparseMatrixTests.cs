using System;
using NUnit.Framework;
namespace linear_A.Test;

public class SparseMatrixTests
{
    [Test]
    public void StorageAndFramesOwnTheirValues()
    {
        var values = new[] { 4.0, 1, 1, 3 };
        var a = new CSRMatrix(2, 2, new[] { 0, 2, 4 }, new[] { 0, 1, 0, 1 }, values);
        values[0] = 99;
        Assert.That(a.Matvec(new[] { 1.0, 2 }), Is.EqualTo(new[] { 6.0, 7 }));
        var result = a.ConjugateGradient(new[] { 6.0, 7 }, jacobi: true, capture: true);
        Assert.That(result.Converged, Is.True);
        Assert.That(result.X[0], Is.EqualTo(1).Within(1e-12));
        Assert.That(result.Iterates.Length, Is.EqualTo(result.Iterations + 1));
        result.X[0] = 99;
        Assert.That(result.Iterates[^1][0], Is.EqualTo(1).Within(1e-12));
    }
    [Test]
    public void LimitAndCurvatureAreExplicitResults()
    {
        var a = new CSRMatrix(1, 1, new[] { 0, 1 }, new[] { 0 }, new[] { -1.0 });
        Assert.That(a.ConjugateGradient(new[] { 1.0 }, maxIterations: 0).Reason, Is.EqualTo("iteration_limit"));
        Assert.That(a.ConjugateGradient(new[] { 1.0 }).Reason, Is.EqualTo("breakdown"));
        Assert.Throws<ArgumentException>(() => new CSRMatrix(1, 1, new[] { 0, 2 }, new[] { 0, 0 }, new[] { 1.0, 2 }));
    }
}
