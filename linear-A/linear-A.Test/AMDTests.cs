using System;
using NUnit.Framework;

namespace linear_A.Test;

[TestFixture]
public class AMDTests
{
    [Test]
    public void OrderingIsOwnedDeterministicAndSolvesOriginalSystem()
    {
        var a = new CSRMatrix(3, 3, [0, 3, 5, 7], [0, 1, 2, 0, 1, 0, 2], [4, -1, -1, -1, 2, -1, 2]);
        var p = a.ApproximateMinimumDegree();
        Assert.That(p, Is.EqualTo(new[] { 1, 0, 2 }));
        var q = a.PermuteSymmetric(p);
        var b = a.Matvec([1, 2, 3]);
        var x = new SparseCholeskySymbolic(q).Factorize(q).Solve(CSRMatrix.PermuteVector(p, b));
        Assert.That(CSRMatrix.PermuteVector(p, x, true), Is.EqualTo(new double[] { 1, 2, 3 }).Within(1e-12));
        p[0] = 999;
        Assert.That(a.ApproximateMinimumDegree(), Is.EqualTo(new[] { 1, 0, 2 }));
    }
    [Test]
    public void PatternIncludesStoredZerosAndRejectsRectangularMatrices()
    {
        Assert.That(new CSRMatrix(3, 3, [0, 2, 2, 2], [1, 2], [0, 0]).ApproximateMinimumDegree(), Is.EqualTo(new[] { 1, 0, 2 }));
        Assert.That(new CSRMatrix(0, 0, [0], [], []).ApproximateMinimumDegree(), Is.Empty);
        Assert.Throws<ArgumentException>(() => new CSRMatrix(1, 2, [0, 0], [], []).ApproximateMinimumDegree());
    }
}
