using System;
using NUnit.Framework;
namespace linear_A.Test;

public class OrderingTests
{
    [Test]
    public void PermutationPreservesStoredZerosAndRestoresOriginalCoordinates()
    {
        var a = new CSRMatrix(3, 3, new[] { 0, 2, 3, 5 }, new[] { 0, 2, 1, 0, 2 }, new[] { 4.0, 0, 5, 2, 6 });
        var p = new[] { 2, 0, 1 }; var q = a.PermuteSymmetric(p); var x = new[] { 1.0, 2, 3 }; var y = CSRMatrix.PermuteVector(p, x);
        Assert.That(q.NNZ, Is.EqualTo(5)); Assert.That(CSRMatrix.PermuteVector(p, y, true), Is.EqualTo(x));
        Assert.That(CSRMatrix.PermuteVector(p, q.Matvec(y), true), Is.EqualTo(a.Matvec(x)));
        Assert.That(a.ReverseCuthillMcKee(), Is.EqualTo(new[] { 2, 0, 1 }));
        p[0] = 0; Assert.That(q.Matvec(y), Is.EqualTo(new[] { 20.0, 4, 10 }));
        Assert.Throws<ArgumentException>(() => a.PermuteSymmetric(p));
        Assert.Throws<ArgumentException>(() => CSRMatrix.PermuteVector(new[] { 0 }, new[] { double.NaN }));
        Assert.Throws<ArgumentException>(() => new CSRMatrix(1, 2, new[] { 0, 0 }, Array.Empty<int>(), Array.Empty<double>()).ReverseCuthillMcKee());
    }
}
