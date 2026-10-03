using System;
using NUnit.Framework;

namespace linear_A.Test;

[TestFixture]
public class SparseCholeskyTests
{
    [Test]
    public void PlansAndFactorsOwnStorageAndCanBeReused()
    {
        var a = new CSRMatrix(2, 2, [0, 2, 4], [0, 1, 0, 1], [4, 2, 2, 3]);
        var plan = new SparseCholeskySymbolic(a);
        var factor = plan.Factorize(a);
        plan.RowOffsets[1] = 99;
        plan.ColumnIndices[0] = 99;
        var scaled = new CSRMatrix(2, 2, [0, 2, 4], [0, 1, 0, 1], [8, 4, 4, 6]);
        Assert.That(factor.Solve([8, 8]), Is.EqualTo(new double[] { 1, 2 }).Within(1e-12));
        Assert.That(plan.Factorize(scaled).Solve([16, 16]), Is.EqualTo(new double[] { 1, 2 }).Within(1e-12));
        Assert.Throws<ArgumentException>(() => plan.Factorize(new CSRMatrix(2, 2, [0, 1, 2], [0, 1], [4, 3])));
        Assert.Throws<ArgumentException>(() => factor.Solve([double.NaN, 1]));
        Assert.Throws<ArithmeticException>(() => factor.Solve([double.MaxValue, -double.MaxValue]));
        Assert.That(factor.Solve([0, 0]), Is.EqualTo(new double[] { 0, 0 }));
    }
    [Test]
    public void StarCreatesThreeEntriesAtFirstPivot()
    {
        var a = new CSRMatrix(4, 4, [0, 4, 6, 8, 10], [0, 1, 2, 3, 0, 1, 0, 2, 0, 3], [5, 1, 1, 1, 1, 3, 1, 3, 1, 3]);
        var s = new SparseCholeskySymbolic(a);
        Assert.That(s.FillCount, Is.EqualTo(3));
        Assert.That(s.FillSteps, Is.EqualTo(new[] { -1, -1, -1, -1, 0, -1, -1, 0, 0, -1 }));
        Assert.That(s.Factorize(a).Solve(a.Matvec([1, 2, 3, 4])), Is.EqualTo(new double[] { 1, 2, 3, 4 }).Within(1e-12));
    }
    [Test]
    public void NumericalPhaseRejectsInvalidPivotsAndSymmetry()
    {
        foreach (var values in new[] { new double[] { 1, 2, 2, 1 }, new double[] { 1, 1, 1, 1 } })
        {
            var a = new CSRMatrix(2, 2, [0, 2, 4], [0, 1, 0, 1], values);
            Assert.Throws<ArithmeticException>(() => new SparseCholeskySymbolic(a).Factorize(a));
        }
        var asymmetric = new CSRMatrix(2, 2, [0, 2, 4], [0, 1, 0, 1], [4, 2, 1, 3]);
        Assert.Throws<ArgumentException>(() => new SparseCholeskySymbolic(asymmetric).Factorize(asymmetric));
    }
}
