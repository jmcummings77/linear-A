using System;
using NUnit.Framework;

namespace linear_A.Test;

public class MatrixSolverTests
{
    [TestCase(FactorizationAlgorithm.Lu)]
    [TestCase(FactorizationAlgorithm.Cholesky)]
    [TestCase(FactorizationAlgorithm.Qr)]
    public void FactorsOwnSnapshotsAndSupportRepeatedRightHandSides(FactorizationAlgorithm algorithm)
    {
        var a = new Matrix<double>(new[,] { { 4.0, 1 }, { 1, 3 } });
        var b = new Matrix<double>(new[,] { { 6.0, 5 }, { 7, 4 } });
        var factor = new MatrixFactorization(a, algorithm);
        a[0, 0] = 99;
        for (var run = 0; run < 3; run++)
        {
            var x = factor.Solve(b);
            Assert.That(x[0, 0], Is.EqualTo(1).Within(1e-12));
            Assert.That(x[1, 0], Is.EqualTo(2).Within(1e-12));
            Assert.That(x[0, 1], Is.EqualTo(1).Within(1e-12));
            Assert.That(x[1, 1], Is.EqualTo(1).Within(1e-12));
        }
        Assert.That(factor.ReciprocalCondition(), Is.EqualTo(.44).Within(1e-12));
        Assert.That(b[0, 0], Is.EqualTo(6));
        Assert.Throws<ArgumentException>(() => factor.Solve(new Matrix<double>(1, 1)));
        b[0, 0] = double.NaN;
        Assert.Throws<ArithmeticException>(() => factor.Solve(b));
    }
}
