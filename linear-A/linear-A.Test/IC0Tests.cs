using System;
using NUnit.Framework;

namespace linear_A.Test;

[TestFixture]
public class IC0Tests
{
    [Test]
    public void FactorReuseAndConflictingOptions()
    {
        var a = new CSRMatrix(2, 2, [0, 2, 4], [0, 1, 0, 1], [4, 1, 1, 3]);
        var f = new IC0(a);
        f.Lower.Values[0] = 99;
        foreach (var rhs in new double[][] { [6, 7], [11, 13] })
        {
            var r = a.ConjugateGradient(rhs, preconditioner: f, capture: true);
            Assert.That(r.Converged, Is.True);
            Assert.That(r.Iterations, Is.EqualTo(1));
            Assert.That(a.Matvec(r.X), Is.EqualTo(rhs).Within(1e-12));
        }
        Assert.Throws<ArgumentException>(() => a.ConjugateGradient([0, 0], jacobi: true, preconditioner: f));
        var small = new IC0(new CSRMatrix(1, 1, [0, 1], [0], [1]));
        Assert.Throws<ArgumentException>(() => a.ConjugateGradient([0, 0], preconditioner: small));
    }
    [Test]
    public void PositiveDefiniteInputCanHaveAnIncompletePivotFailure()
    {
        var a = new CSRMatrix(4, 4, [0, 3, 6, 9, 12], [0, 1, 3, 0, 1, 2, 1, 2, 3, 0, 2, 3], [1.5, 1, 1, 1, 1.5, 1, 1, 1.5, -1, 1, -1, 1.5]);
        Assert.DoesNotThrow(() => new SparseCholeskySymbolic(a).Factorize(a));
        Assert.Throws<ArithmeticException>(() => new IC0(a));
    }
    [Test]
    public void ApplyOverflowKeepsTheInitialIterate()
    {
        var a = new CSRMatrix(1, 1, [0, 1], [0], [1e-320]);
        var r = a.ConjugateGradient([1], preconditioner: new IC0(a), capture: true);
        Assert.That(r.Reason, Is.EqualTo("nonfinite"));
        Assert.That(r.Iterations, Is.Zero);
        Assert.That(r.X, Is.EqualTo(new double[] { 0 }));
        Assert.That(r.Residuals, Is.EqualTo(new double[] { 1 }));
    }

    [Test]
    public void IncompleteFactorOmitsTheFillCreatedByTheFirstPivot()
    {
        // L = [2 0 0; 1/2 sqrt(15/4) 0; 1/2 0 sqrt(15/4)].
        // Dropping L[2,1] makes (L L^T)[1,2] = 1/4, so applying
        // L L^T to [1,2,3] gives [9,39/4,27/2].
        var a = new CSRMatrix(3, 3, [0, 3, 5, 7], [0, 1, 2, 0, 1, 0, 2], [4, 1, 1, 1, 4, 1, 4]);
        var factor = new IC0(a);

        Assert.That(new SparseCholeskySymbolic(a).FillCount, Is.EqualTo(1));
        Assert.That(factor.NNZ, Is.EqualTo(5));
        Assert.That(factor.Lower.RowOffsets, Is.EqualTo(new[] { 0, 1, 3, 5 }));
        Assert.That(factor.Lower.ColumnIndices, Is.EqualTo(new[] { 0, 0, 1, 0, 2 }));
        Assert.That(factor.Apply([9, 9.75, 13.5]), Is.EqualTo(new double[] { 1, 2, 3 }).Within(1e-12));
    }
}
