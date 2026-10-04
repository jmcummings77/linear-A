using System;
using NUnit.Framework;

namespace linear_A.Test;

public class GmresTests
{
    [Test]
    public void RestartedNonsymmetricSolveOwnsItsFrames()
    {
        var a = new CSRMatrix(2, 2, new[] { 0, 2, 3 }, new[] { 0, 1, 1 }, new[] { 3.0, 1, 2 });
        var b = new[] { 5.0, 4 };
        var r = a.Gmres(b, restart: 1, capture: true);

        Assert.That(r.Converged, Is.True);
        Assert.That(r.X, Is.EqualTo(new[] { 1.0, 2 }).Within(1e-8));
        Assert.That(r.Restarts.Length, Is.EqualTo(r.Iterations - 1));
        for (var i = 0; i < r.Restarts.Length; i++)
        {
            Assert.That(r.Restarts[i], Is.EqualTo(i + 1));
        }
        AssertCapturedResiduals(a, b, r);

        var firstStep = r.Iterates[1][0];
        r.Iterates[0][0] = 98;
        r.Iterates[^1][0] = 99;
        Assert.That(r.Iterates[1][0], Is.EqualTo(firstStep));
        Assert.That(r.X[0], Is.EqualTo(1).Within(1e-8));
        Assert.That(b, Is.EqualTo(new[] { 5.0, 4 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => a.Gmres(new[] { 5.0, 4 }, restart: 0));
    }

    [Test]
    public void FullKrylovCycleSolvesNonsymmetricSystemWithoutCapturingFrames()
    {
        var a = new CSRMatrix(2, 2, new[] { 0, 2, 3 }, new[] { 0, 1, 1 }, new[] { 3.0, 1, 2 });
        var r = a.Gmres(new[] { 5.0, 4 }, restart: 2);

        Assert.That(r.Converged, Is.True);
        Assert.That(r.Iterations, Is.EqualTo(2));
        Assert.That(r.X, Is.EqualTo(new[] { 1.0, 2 }).Within(1e-12));
        Assert.That(r.Restarts, Is.Empty);
        Assert.That(r.Iterates, Is.Empty);
        Assert.That(r.Residuals.Length, Is.EqualTo(3));
        Assert.That(r.EstimatedResiduals.Length, Is.EqualTo(3));
    }

    [TestCase(0.0, 0.0, true)]
    [TestCase(2.0, 2.0, true)]
    [TestCase(2.0, 1.0, false)]
    public void InitialResidualIsCheckedBeforeTheIterationLimit(double rhs, double tolerance, bool converged)
    {
        var a = new CSRMatrix(1, 1, new[] { 0, 1 }, new[] { 0 }, new[] { 2.0 });
        var r = a.Gmres(new[] { rhs }, absoluteTolerance: tolerance, maxIterations: 0, capture: true);

        Assert.That(r.Converged, Is.EqualTo(converged));
        Assert.That(r.Reason, Is.EqualTo(converged ? "converged" : "iteration_limit"));
        Assert.That(r.Iterations, Is.Zero);
        Assert.That(r.X, Is.EqualTo(new[] { 0.0 }));
        Assert.That(r.Residuals, Is.EqualTo(new[] { rhs }));
        Assert.That(r.EstimatedResiduals, Is.EqualTo(r.Residuals));
        Assert.That(r.Restarts, Is.Empty);
        AssertCapturedResiduals(a, new[] { rhs }, r);
    }

    [TestCase(0.3, 0.0, true)]
    [TestCase(0.0, 2.0, true)]
    [TestCase(0.0, 1.0, false)]
    public void TrueResidualControlsToleranceAtTheLastAllowedIteration(double relative, double absolute, bool converged)
    {
        var a = new CSRMatrix(2, 2, new[] { 0, 2, 3 }, new[] { 0, 1, 1 }, new[] { 3.0, 1, 2 });
        var b = new[] { 5.0, 4 };
        var r = a.Gmres(b, relativeTolerance: relative, absoluteTolerance: absolute, maxIterations: 1, capture: true);

        // The first minimizer is x = b * ((Ab)^T b / ||Ab||^2), with Ab = (19, 8).
        Assert.That(r.X, Is.EqualTo(new[] { 635.0 / 425, 508.0 / 425 }).Within(1e-12));
        Assert.That(r.Converged, Is.EqualTo(converged));
        Assert.That(r.Reason, Is.EqualTo(converged ? "converged" : "iteration_limit"));
        Assert.That(r.Iterations, Is.EqualTo(1));
        AssertCapturedResiduals(a, b, r);
    }

    [Test]
    public void IterationLimitCanEndPartwayThroughARestartCycle()
    {
        var a = new CSRMatrix(3, 3, new[] { 0, 2, 4, 5 }, new[] { 0, 1, 1, 2, 2 }, new[] { 3.0, 1, 2, 1, 1 });
        var b = new[] { 5.0, 7, 3 };
        var r = a.Gmres(b, restart: 2, relativeTolerance: 0, maxIterations: 3, capture: true);

        Assert.That(r.Converged, Is.False);
        Assert.That(r.Reason, Is.EqualTo("iteration_limit"));
        Assert.That(r.Iterations, Is.EqualTo(3));
        Assert.That(r.Restarts, Is.EqualTo(new[] { 2 }));
        AssertCapturedResiduals(a, b, r);
    }

    [Test]
    public void JacobiAppliesTheRightPreconditionerToTheSolutionCorrection()
    {
        var a = new CSRMatrix(2, 2, new[] { 0, 1, 2 }, new[] { 0, 1 }, new[] { 2.0, 8 });
        var b = new[] { 2.0, 16 };
        var r = a.Gmres(b, restart: 1, maxIterations: 1, jacobi: true, capture: true);

        Assert.That(r.Converged, Is.True);
        Assert.That(r.Iterations, Is.EqualTo(1));
        Assert.That(r.X, Is.EqualTo(new[] { 1.0, 2 }).Within(1e-12));
        Assert.That(r.Residuals[0], Is.EqualTo(double.Hypot(2, 16)));
        AssertCapturedResiduals(a, b, r);
    }

    [Test]
    public void IluAppliesTheRightPreconditionerToANonsymmetricSystem()
    {
        var a = new CSRMatrix(2, 2, new[] { 0, 2, 3 }, new[] { 0, 1, 1 }, new[] { 3.0, 1, 2 });
        var b = new[] { 5.0, 4 };
        var r = a.Gmres(b, restart: 1, maxIterations: 1, preconditioner: new ILU0(a), capture: true);

        Assert.That(r.Converged, Is.True);
        Assert.That(r.Iterations, Is.EqualTo(1));
        Assert.That(r.X, Is.EqualTo(new[] { 1.0, 2 }).Within(1e-12));
        AssertCapturedResiduals(a, b, r);
    }

    [Test]
    public void InvalidPreconditionersAreRejectedBeforeEarlyConvergence()
    {
        var missingDiagonal = new CSRMatrix(2, 2, new[] { 0, 1, 1 }, new[] { 1 }, new[] { 1.0 });
        var zeroDiagonal = new CSRMatrix(2, 2, new[] { 0, 1, 2 }, new[] { 0, 1 }, new[] { 1.0, 0 });
        var smallFactor = new ILU0(new CSRMatrix(1, 1, new[] { 0, 1 }, new[] { 0 }, new[] { 1.0 }));

        Assert.Throws<ArgumentException>(() => missingDiagonal.Gmres(new double[2], jacobi: true));
        Assert.Throws<ArgumentException>(() => zeroDiagonal.Gmres(new double[2], jacobi: true));
        Assert.Throws<ArgumentException>(() => zeroDiagonal.Gmres(new double[2], preconditioner: smallFactor));
    }

    [Test]
    public void SingularSystemCannotClaimHappyConvergence()
    {
        var a = new CSRMatrix(2, 2, new[] { 0, 1, 1 }, new[] { 0 }, new[] { 1.0 });
        var b = new[] { 1.0, 1 };
        var r = a.Gmres(b, capture: true);

        Assert.That(r.Converged, Is.False);
        Assert.That(r.Reason, Is.EqualTo("breakdown"));
        Assert.That(r.Residuals[^1], Is.GreaterThanOrEqualTo(1));
        AssertCapturedResiduals(a, b, r);
    }

    [Test]
    public void ZeroOperatorBreaksDownWithoutAcceptingAnIterate()
    {
        var a = new CSRMatrix(1, 1, new[] { 0, 0 }, Array.Empty<int>(), Array.Empty<double>());
        var r = a.Gmres(new[] { 1.0 }, capture: true);

        AssertUnacceptedStep(r, "breakdown", 1);
        Assert.That(r.Residuals, Is.EqualTo(new[] { 1.0 }));
    }

    [Test]
    public void OverflowingInitialNormIsReportedWithoutStartingACycle()
    {
        var a = new CSRMatrix(2, 2, new[] { 0, 1, 2 }, new[] { 0, 1 }, new[] { 1.0, 1 });
        var r = a.Gmres(new[] { double.MaxValue, double.MaxValue }, capture: true);

        AssertUnacceptedStep(r, "nonfinite", 2);
        Assert.That(r.Residuals[0], Is.EqualTo(double.PositiveInfinity));
    }

    [Test]
    public void OverflowingArnoldiProductDoesNotAppendAHistoryEntry()
    {
        var a = new CSRMatrix(2, 2, new[] { 0, 2, 3 }, new[] { 0, 1, 1 }, new[] { double.MaxValue, double.MaxValue, 1.0 });
        var r = a.Gmres(new[] { 1.0, 1 }, capture: true);

        AssertUnacceptedStep(r, "nonfinite", 2);
        Assert.That(r.Residuals[0], Is.EqualTo(double.Hypot(1, 1)));
    }

    [TestCase(false)]
    [TestCase(true)]
    public void OverflowingPreconditionerDoesNotAppendAHistoryEntry(bool ilu)
    {
        var a = new CSRMatrix(1, 1, new[] { 0, 1 }, new[] { 0 }, new[] { 1e-320 });
        var r = a.Gmres(new[] { 1.0 }, jacobi: !ilu, preconditioner: ilu ? new ILU0(a) : null, capture: true);

        AssertUnacceptedStep(r, "nonfinite", 1);
        Assert.That(r.Residuals, Is.EqualTo(new[] { 1.0 }));
    }

    private static void AssertUnacceptedStep(GMRESResult result, string reason, int size)
    {
        Assert.That(result.Converged, Is.False);
        Assert.That(result.Reason, Is.EqualTo(reason));
        Assert.That(result.Iterations, Is.Zero);
        Assert.That(result.X, Is.EqualTo(new double[size]));
        Assert.That(result.Residuals.Length, Is.EqualTo(1));
        Assert.That(result.EstimatedResiduals, Is.EqualTo(result.Residuals));
        Assert.That(result.Iterates.Length, Is.EqualTo(1));
        Assert.That(result.Iterates[0], Is.EqualTo(result.X));
        Assert.That(result.Restarts, Is.Empty);
    }

    private static void AssertCapturedResiduals(CSRMatrix matrix, double[] rhs, GMRESResult result)
    {
        Assert.That(result.Residuals.Length, Is.EqualTo(result.Iterations + 1));
        Assert.That(result.EstimatedResiduals.Length, Is.EqualTo(result.Residuals.Length));
        Assert.That(result.Iterates.Length, Is.EqualTo(result.Residuals.Length));
        Assert.That(result.Iterates[0], Is.EqualTo(new double[rhs.Length]));
        Assert.That(result.Iterates[^1], Is.EqualTo(result.X));

        for (var step = 0; step < result.Iterates.Length; step++)
        {
            var product = matrix.Matvec(result.Iterates[step]);
            var squaredResidual = 0.0;
            for (var i = 0; i < rhs.Length; i++)
            {
                var residual = rhs[i] - product[i];
                squaredResidual += residual * residual;
            }
            Assert.That(result.Residuals[step], Is.EqualTo(Math.Sqrt(squaredResidual)).Within(1e-12));
            Assert.That(double.IsFinite(result.EstimatedResiduals[step]), Is.True);
            Assert.That(result.EstimatedResiduals[step], Is.GreaterThanOrEqualTo(0));
        }
    }
}
