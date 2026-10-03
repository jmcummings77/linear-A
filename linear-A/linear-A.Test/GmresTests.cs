using System;
using NUnit.Framework;
namespace linear_A.Test;

public class GmresTests
{
    [Test]
    public void RestartedNonsymmetricSolveOwnsItsFrames()
    {
        var a = new CSRMatrix(2, 2, new[] { 0, 2, 3 }, new[] { 0, 1, 1 }, new[] { 3.0, 1, 2 });
        var r = a.Gmres(new[] { 5.0, 4 }, restart: 1, capture: true);
        Assert.That(r.Converged, Is.True);
        Assert.That(r.Restarts.Length, Is.GreaterThan(0));
        Assert.That(r.X[0], Is.EqualTo(1).Within(1e-8));
        Assert.That(r.EstimatedResiduals.Length, Is.EqualTo(r.Iterations + 1));
        r.Iterates[^1][0] = 99;
        Assert.That(r.X[0], Is.EqualTo(1).Within(1e-8));
        Assert.Throws<ArgumentOutOfRangeException>(() => a.Gmres(new[] { 5.0, 4 }, restart: 0));
    }
    [Test]
    public void SingularSystemCannotClaimHappyConvergence()
    {
        var a = new CSRMatrix(2, 2, new[] { 0, 1, 1 }, new[] { 0 }, new[] { 1.0 });
        var r = a.Gmres(new[] { 1.0, 1 });
        Assert.That(r.Converged, Is.False);
        Assert.That(r.Reason, Is.EqualTo("breakdown"));
    }
}
