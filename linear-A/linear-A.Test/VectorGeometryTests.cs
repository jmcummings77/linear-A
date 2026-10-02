using System;
using System.Numerics;
using NUnit.Framework;

namespace linear_A.Test;

public class VectorGeometryTests
{
    private static Matrix<T> Vector<T>(T x, T y, T z, bool row = false) where T : INumber<T>
        => new(row ? new[,] { { x, y, z } } : new[,] { { x }, { y }, { z } });

    private static DoubleMatrix Legacy(double x, double y, double z, bool row = false)
    {
        var result = new DoubleMatrix(row ? 1 : 3, row ? 3 : 1);
        var values = new[] { x, y, z };
        for (var i = 0; i < 3; i++) result[row ? 0 : i, row ? i : 0] = values[i];
        return result;
    }

    private static void AssertMatrix(DoubleMatrix actual, DoubleMatrix expected, double tolerance = 2e-14)
    {
        Assert.That(actual.RowCount, Is.EqualTo(expected.RowCount));
        Assert.That(actual.ColumnCount, Is.EqualTo(expected.ColumnCount));
        for (var i = 0; i < actual.RowCount; i++)
            for (var j = 0; j < actual.ColumnCount; j++)
                Assert.That(actual[i, j], Is.EqualTo(expected[i, j]).Within(tolerance));
    }

    private static void AssertRotation(DoubleMatrix rotation)
    {
        AssertMatrix(rotation.Transpose() * rotation, new DoubleMatrix(rotation.RowCount, true));
        Assert.That(rotation.GetDeterminant(), Is.EqualTo(1d).Within(2e-14));
    }

    [Test]
    public void CrossIsRightHandedAndPreservesLeftOrientationAndInputs()
    {
        foreach (var leftRow in new[] { false, true })
            foreach (var rightRow in new[] { false, true })
            {
                var a = Vector(1d, 2d, 3d, leftRow);
                var b = Vector(4d, 5d, 6d, rightRow);
                var cross = a.CrossProduct(b);
                Assert.That(cross.RowCount, Is.EqualTo(a.RowCount));
                Assert.That(cross.ColumnCount, Is.EqualTo(a.ColumnCount));
                Assert.That(leftRow ? cross.GetRow(0) : cross.GetColumn(0), Is.EqualTo(new[] { -3d, 6d, -3d }));
                Assert.That(leftRow ? a.GetRow(0) : a.GetColumn(0), Is.EqualTo(new[] { 1d, 2d, 3d }));
                Assert.That(rightRow ? b.GetRow(0) : b.GetColumn(0), Is.EqualTo(new[] { 4d, 5d, 6d }));
                cross[0, 0] = 99;
                Assert.That(a[0, 0], Is.EqualTo(1));
            }
        Assert.That(Vector(1, 0, 0).CrossProduct(Vector(0, 1, 0)).GetColumn(0), Is.EqualTo(new[] { 0, 0, 1 }));
        Assert.That(Vector(0, 1, 0).CrossProduct(Vector(1, 0, 0)).GetColumn(0), Is.EqualTo(new[] { 0, 0, -1 }));
        Assert.That(Vector(1, 2, 3).CrossProduct(Vector(2, 4, 6)).GetColumn(0), Is.EqualTo(new[] { 0, 0, 0 }));
        var big = BigInteger.One << 100;
        Assert.That(Vector(big, BigInteger.Zero, BigInteger.Zero).CrossProduct(Vector(BigInteger.Zero, big, BigInteger.Zero))[2, 0], Is.EqualTo(big * big));
        Assert.That(Vector(0.5m, 0m, 0m).CrossProduct(Vector(0m, 0.25m, 0m))[2, 0], Is.EqualTo(0.125m));
    }

    [Test]
    public void CrossRejectsInvalidAndNonfiniteArithmeticWithoutMutation()
    {
        var a = Vector(1d, 2d, 3d);
        Assert.Throws<ArgumentNullException>(() => a.CrossProduct(null!));
        foreach (var shape in new[] { new Matrix<double>(), new Matrix<double>(3), new Matrix<double>(2, 1), new Matrix<double>(1, 4) })
            Assert.Throws<ArgumentException>(() => a.CrossProduct(shape));
        foreach (var value in new[] { double.NaN, double.PositiveInfinity })
            Assert.Throws<ArgumentException>(() => a.CrossProduct(Vector(value, 0d, 0d)));
        Assert.Throws<ArithmeticException>(() => Vector(double.MaxValue, 0d, 0d).CrossProduct(Vector(0d, 2d, 0d)));
        var integers = Vector(int.MaxValue, 0, 0);
        Assert.Throws<OverflowException>(() => integers.CrossProduct(Vector(0, 2, 0)));
        Assert.That(integers[0, 0], Is.EqualTo(int.MaxValue));
        Assert.That(a.GetColumn(0), Is.EqualTo(new[] { 1d, 2d, 3d }));
    }

    [Test]
    public void LegacyCrossHasTransactionalInPlaceAndIndependentResultForms()
    {
        var a = Legacy(1, 2, 3, true);
        var b = Legacy(4, 5, 6);
        Assert.That(a.GetCrossProduct(b).GetRow(0), Is.EqualTo(new[] { -3d, 6d, -3d }));
        Assert.That(a.GetRow(0), Is.EqualTo(new[] { 1d, 2d, 3d }));
        Assert.That(a.CrossProduct(b), Is.True);
        Assert.That(a.GetRow(0), Is.EqualTo(new[] { -3d, 6d, -3d }));
        Assert.That(a.CrossProduct(null), Is.False);
        Assert.That(a.CrossProduct(new DoubleMatrix(3, true)), Is.False);
        Assert.That(a.CrossProduct(Legacy(double.NaN, 0, 0)), Is.False);
        Assert.That(a.GetRow(0), Is.EqualTo(new[] { -3d, 6d, -3d }));
        var huge = Legacy(double.MaxValue, 0, 0);
        Assert.That(huge.CrossProduct(Legacy(0, 2, 0)), Is.False);
        Assert.That(huge[0, 0], Is.EqualTo(double.MaxValue));

        var integer = new IntegerMatrix(3, 1);
        var right = new IntegerMatrix(1, 3);
        integer[0, 0] = 1;
        right[0, 1] = 1;
        Assert.That(integer.GetCrossProduct(right).GetColumn(0), Is.EqualTo(new[] { 0, 0, 1 }));
        Assert.That(integer.CrossProduct(right), Is.True);
        Assert.That(integer.GetColumn(0), Is.EqualTo(new[] { 0, 0, 1 }));
        integer[0, 0] = int.MaxValue;
        right[0, 1] = 2;
        Assert.That(integer.CrossProduct(right), Is.False);
        Assert.That(integer[0, 0], Is.EqualTo(int.MaxValue));
        Assert.That(integer.CrossProduct(null), Is.False);
    }

    [Test]
    public void PrincipalRotationsHaveCorrectHandedness()
    {
        var quarter = Math.PI / 2;
        var plane = MatrixRotation.Create2D(quarter);
        Assert.That(plane[0, 0], Is.EqualTo(0).Within(1e-15));
        Assert.That(plane[1, 0], Is.EqualTo(1).Within(1e-15));
        AssertMatrix(MatrixRotation.CreateX(quarter) * Legacy(0, 1, 0), Legacy(0, 0, 1));
        AssertMatrix(MatrixRotation.CreateY(quarter) * Legacy(0, 0, 1), Legacy(1, 0, 0));
        AssertMatrix(MatrixRotation.CreateZ(quarter) * Legacy(1, 0, 0), Legacy(0, 1, 0));
        foreach (var angle in new[] { 0d, 1e-12, -0.5, Math.PI, 2 * Math.PI, double.MaxValue })
        {
            foreach (var create in new Func<double, DoubleMatrix>[] { MatrixRotation.Create2D, MatrixRotation.CreateX, MatrixRotation.CreateY, MatrixRotation.CreateZ })
            {
                var rotation = create(angle);
                AssertRotation(rotation);
                AssertMatrix(create(-angle), rotation.Transpose());
            }
        }
    }

    [Test]
    public void AxisAnglePreservesAxisAndHandlesExtremeScales()
    {
        var expected = MatrixRotation.CreateAxisAngle(Legacy(1, 2, 3), 0.5);
        foreach (var scale in new[] { 1d, 1e300, 1e-300, double.Epsilon })
            foreach (var row in new[] { false, true })
            {
                var axis = Legacy(scale, 2 * scale, 3 * scale, row);
                var rotation = MatrixRotation.CreateAxisAngle(axis, 0.5);
                AssertRotation(rotation);
                AssertMatrix(rotation, expected);
                AssertMatrix(rotation * Legacy(1, 2, 3), Legacy(1, 2, 3));
                AssertMatrix(MatrixRotation.CreateAxisAngle(axis, -0.5) * rotation, new DoubleMatrix(3, true));
                AssertMatrix(MatrixRotation.CreateAxisAngle(axis, 0), new DoubleMatrix(3, true), 0);
                Assert.That(axis[0, 0], Is.EqualTo(scale));
            }
        AssertMatrix(MatrixRotation.CreateAxisAngle(Vector(1, 2, 3), 0.5), expected);
        AssertMatrix(MatrixRotation.CreateAxisAngle(Legacy(0, 0, 1), 0.5), MatrixRotation.CreateZ(0.5));
        AssertMatrix(MatrixRotation.CreateAxisAngle(Legacy(-1, -2, -3), -0.5), expected);
    }

    [Test]
    public void RotationRejectsInvalidAxisAndAngle()
    {
        foreach (var value in new[] { double.NaN, double.PositiveInfinity, double.NegativeInfinity })
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => MatrixRotation.Create2D(value));
            Assert.Throws<ArgumentOutOfRangeException>(() => MatrixRotation.CreateX(value));
            Assert.Throws<ArgumentOutOfRangeException>(() => MatrixRotation.CreateAxisAngle(Legacy(1, 2, 3), value));
            Assert.Throws<ArgumentException>(() => MatrixRotation.CreateAxisAngle(Legacy(value, 0, 0), 0.5));
        }
        Assert.Throws<ArgumentException>(() => MatrixRotation.CreateAxisAngle(Legacy(0, 0, 0), 0));
        Assert.Throws<ArgumentException>(() => MatrixRotation.CreateAxisAngle(new DoubleMatrix(3, true), 0.5));
        Assert.Throws<ArgumentNullException>(() => MatrixRotation.CreateAxisAngle((DoubleMatrix)null!, 0.5));
    }
}
