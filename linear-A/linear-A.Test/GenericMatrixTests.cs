using System;
using System.Collections;
using System.Numerics;
using NUnit.Framework;

namespace linear_A.Test;

[TestFixture(typeof(int))]
[TestFixture(typeof(double))]
[TestFixture(typeof(decimal))]
[TestFixture(typeof(BigInteger))]
public class GenericMatrixTests<T> where T : INumber<T>
{
    private static T Number(int value) => T.CreateChecked(value);

    private static void AssertDeterminant(T actual, int expected, string? message = null)
    {
        if (typeof(T) == typeof(double))
            Assert.That(double.CreateChecked(actual), Is.EqualTo((double)expected).Within(1e-10 * Math.Max(1, Math.Abs(expected))), message);
        else
            Assert.That(actual, Is.EqualTo(Number(expected)), message);
    }

    private static Matrix<T> From(int[,] values)
    {
        var result = new Matrix<T>(values.GetLength(0), values.GetLength(1));
        for (var row = 0; row < result.RowCount; row++)
            for (var column = 0; column < result.ColumnCount; column++)
                result[row, column] = Number(values[row, column]);
        return result;
    }

    private static void AssertMatrix(Matrix<T> actual, int[,] expected)
    {
        Assert.That(actual.RowCount, Is.EqualTo(expected.GetLength(0)));
        Assert.That(actual.ColumnCount, Is.EqualTo(expected.GetLength(1)));
        for (var row = 0; row < actual.RowCount; row++)
            for (var column = 0; column < actual.ColumnCount; column++)
                Assert.That(actual[row, column], Is.EqualTo(Number(expected[row, column])), $"Element [{row}, {column}]");
    }

    [Test]
    public void ConstructorsInitializeZerosAndIdentity()
    {
        AssertMatrix(new Matrix<T>(2, 3), new int[,] { { 0, 0, 0 }, { 0, 0, 0 } });
        AssertMatrix(new Matrix<T>(3, true), new int[,] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } });
        AssertMatrix(new Matrix<T>(2), new int[,] { { 0, 0 }, { 0, 0 } });
    }

    [Test]
    public void CopiesDoNotShareArrayStorage()
    {
        var values = new[,] { { Number(1), Number(2) }, { Number(3), Number(4) } };
        var matrix = new Matrix<T>(values);
        var copy = new Matrix<T>(matrix);
        values[0, 0] = Number(9);
        matrix[1, 1] = Number(8);

        AssertMatrix(copy, new int[,] { { 1, 2 }, { 3, 4 } });
        Assert.That(matrix[0, 0], Is.EqualTo(Number(1)));
    }

    [Test]
    public void ConstructorsRejectNegativeDimensionsAndNullSources()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new Matrix<T>(-1, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Matrix<T>(2, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Matrix<T>(-1, true));
        Assert.Throws<ArgumentNullException>(() => new Matrix<T>((Matrix<T>)null!));
        Assert.Throws<ArgumentNullException>(() => new Matrix<T>((T[,])null!));
    }

    [Test]
    public void ArrayConstructorRejectsNonZeroLowerBounds()
    {
        var values = (T[,])Array.CreateInstance(typeof(T), new[] { 2, 2 }, new[] { 1, 1 });
        Assert.Throws<ArgumentException>(() => new Matrix<T>(values));
    }

    [Test]
    public void AdditionAndSubtractionOperatorsReturnIndependentResults()
    {
        var left = From(new int[,] { { 1, 2, 3 }, { -4, 5, 6 } });
        var right = From(new int[,] { { 6, 5, -4 }, { 3, 2, 1 } });
        var sum = left + right;
        var difference = left - right;

        AssertMatrix(sum, new int[,] { { 7, 7, -1 }, { -1, 7, 7 } });
        AssertMatrix(difference, new int[,] { { -5, -3, 7 }, { -7, 3, 5 } });
        sum[0, 0] = Number(99);
        difference[1, 1] = Number(99);
        AssertMatrix(left, new int[,] { { 1, 2, 3 }, { -4, 5, 6 } });
        AssertMatrix(right, new int[,] { { 6, 5, -4 }, { 3, 2, 1 } });
    }

    [Test]
    public void InPlaceOperationsHandleSelfAndShapeMismatch()
    {
        var matrix = From(new int[,] { { 1, 2 }, { 3, 4 } });
        Assert.That(matrix.TryAddMatrix(new Matrix<T>(1, 4)), Is.False);
        Assert.That(matrix.TrySubtractMatrix(new Matrix<T>(4, 1)), Is.False);
        Assert.That(matrix.TryAddMatrix(null), Is.False);
        Assert.That(matrix.TrySubtractMatrix(null), Is.False);
        AssertMatrix(matrix, new int[,] { { 1, 2 }, { 3, 4 } });

        Assert.That(matrix.TryAddMatrix(matrix), Is.True);
        AssertMatrix(matrix, new int[,] { { 2, 4 }, { 6, 8 } });
        Assert.That(matrix.TrySubtractMatrix(matrix), Is.True);
        AssertMatrix(matrix, new int[,] { { 0, 0 }, { 0, 0 } });
    }

    [Test]
    public void MatrixProductUsesRowsOfLeftAndColumnsOfRight()
    {
        var left = From(new int[,] { { 1, 2, 3 }, { 4, 5, 6 } });
        var right = From(new int[,] { { 7, 8 }, { 9, 10 }, { 11, 12 } });
        var expected = new int[,] { { 58, 64 }, { 139, 154 } };

        AssertMatrix(left * right, expected);
        AssertMatrix(left.DotProduct(right), expected);
        AssertMatrix(left, new int[,] { { 1, 2, 3 }, { 4, 5, 6 } });
        AssertMatrix(right, new int[,] { { 7, 8 }, { 9, 10 }, { 11, 12 } });
    }

    [Test]
    public void ProductIsNotElementwiseOrCommutative()
    {
        var left = From(new int[,] { { 1, 2 }, { 0, 1 } });
        var right = From(new int[,] { { 1, 0 }, { 3, 1 } });
        AssertMatrix(left * right, new int[,] { { 7, 2 }, { 3, 1 } });
        AssertMatrix(right * left, new int[,] { { 1, 2 }, { 3, 7 } });
    }

    [Test]
    public void OperatorsRejectNullAndIncompatibleDimensions()
    {
        var matrix = new Matrix<T>(2, 3);
        Assert.Throws<ArgumentException>(() => { _ = matrix + new Matrix<T>(3, 2); });
        Assert.Throws<ArgumentException>(() => { _ = matrix - new Matrix<T>(2, 2); });
        Assert.Throws<ArgumentException>(() => { _ = matrix * new Matrix<T>(2, 3); });
        Assert.Throws<ArgumentNullException>(() => { _ = matrix + null!; });
        Assert.Throws<ArgumentNullException>(() => { _ = (Matrix<T>)null! - matrix; });
        Assert.Throws<ArgumentNullException>(() => { _ = (Matrix<T>)null! * matrix; });
        Assert.Throws<ArgumentNullException>(() => matrix.DotProduct(null!));
    }

    [Test]
    public void RectangularTransposeSwapsDimensionsAndCopiesValues()
    {
        var matrix = From(new int[,] { { 1, 2, 3 }, { 4, 5, 6 } });
        var transpose = matrix.Transpose();
        AssertMatrix(transpose, new int[,] { { 1, 4 }, { 2, 5 }, { 3, 6 } });
        AssertMatrix(transpose.Transpose(), new int[,] { { 1, 2, 3 }, { 4, 5, 6 } });
        transpose[0, 0] = Number(99);
        Assert.That(matrix[0, 0], Is.EqualTo(Number(1)));
    }

    [Test]
    public void ScaleChangesEveryEntry()
    {
        var matrix = From(new int[,] { { 1, 2, -3 }, { 0, 5, 6 } });
        matrix.Scale(Number(-2));
        AssertMatrix(matrix, new int[,] { { -2, -4, 6 }, { 0, -10, -12 } });
    }

    [Test]
    public void RowsAndColumnsAreCopiesAndEnumerationIsInRowOrder()
    {
        var matrix = From(new int[,] { { 1, 2, 3 }, { 4, 5, 6 } });
        var row = matrix.GetRow(1);
        var column = matrix.GetColumn(1);
        Assert.That(row, Is.EqualTo(new[] { Number(4), Number(5), Number(6) }));
        Assert.That(column, Is.EqualTo(new[] { Number(2), Number(5) }));
        row[0] = Number(99);
        column[0] = Number(99);
        Assert.That((IEnumerable)matrix, Is.EqualTo(new[] { Number(1), Number(2), Number(3), Number(4), Number(5), Number(6) }));
        Assert.Throws<ArgumentOutOfRangeException>(() => matrix.GetRow(-1));
        Assert.Throws<ArgumentOutOfRangeException>(() => matrix.GetRow(2));
        Assert.Throws<ArgumentOutOfRangeException>(() => matrix.GetColumn(-1));
        Assert.Throws<ArgumentOutOfRangeException>(() => matrix.GetColumn(3));
    }

    [Test]
    public void EmptyDimensionsObeyMatrixIdentities()
    {
        var empty = new Matrix<T>();
        Assert.That(empty.GetTrace(), Is.EqualTo(T.Zero));
        Assert.That(empty.GetDeterminant(), Is.EqualTo(T.One));
        Assert.That(empty.IsTriangular(), Is.True);
        Assert.That(empty.TryGetTrace(out var trace), Is.True);
        Assert.That(trace, Is.EqualTo(T.Zero));
        Assert.That(empty.TryGetDeterminant(out var determinant), Is.True);
        Assert.That(determinant, Is.EqualTo(T.One));
        Assert.That((empty + empty).RowCount, Is.Zero);
        Assert.That((empty * empty).ColumnCount, Is.Zero);

        AssertMatrix(new Matrix<T>(2, 0) * new Matrix<T>(0, 3), new int[2, 3]);
        var transpose = new Matrix<T>(0, 3).Transpose();
        Assert.That(transpose.RowCount, Is.EqualTo(3));
        Assert.That(transpose.ColumnCount, Is.Zero);
        Assert.That(transpose.GetRow(2), Is.Empty);
        Assert.Throws<ArgumentOutOfRangeException>(() => empty.GetRow(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => empty.GetColumn(0));
    }

    [Test]
    public void TraceAndDeterminantHaveKnownValuesAtEachRecursionSize()
    {
        var singleton = From(new int[,] { { -7 } });
        Assert.That(singleton.GetTrace(), Is.EqualTo(Number(-7)));
        AssertDeterminant(singleton.GetDeterminant(), -7);

        var two = From(new int[,] { { 1, 2 }, { 3, 4 } });
        Assert.That(two.GetTrace(), Is.EqualTo(Number(5)));
        AssertDeterminant(two.GetDeterminant(), -2);

        var three = From(new int[,] { { 6, 1, 1 }, { 4, -2, 5 }, { 2, 8, 7 } });
        Assert.That(three.GetTrace(), Is.EqualTo(Number(11)));
        Assert.That(three.TryGetDeterminant(out var determinant), Is.True);
        AssertDeterminant(determinant, -306);
        Assert.That(three.TryGetTrace(out var trace), Is.True);
        Assert.That(trace, Is.EqualTo(Number(11)));

        var four = From(new int[,] { { 3, 2, 0, 1 }, { 4, 0, 1, 2 }, { 3, 0, 2, 1 }, { 9, 2, 3, 1 } });
        AssertDeterminant(four.GetDeterminant(), 24);
    }

    [Test]
    public void DeterminantHandlesSingularityAndRowSwaps()
    {
        AssertDeterminant(From(new int[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 1, 2, 3 } }).GetDeterminant(), 0);
        Assert.That(From(new int[,] { { 0, 1, 0 }, { 1, 0, 0 }, { 0, 0, 1 } }).GetDeterminant(), Is.EqualTo(Number(-1)));
        Assert.That(From(new int[,] { { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 0 } }).GetDeterminant(), Is.EqualTo(T.One));
        Assert.That(new Matrix<T>(4, true).GetDeterminant(), Is.EqualTo(T.One));
        Assert.That(new Matrix<T>(4, false).GetDeterminant(), Is.EqualTo(T.Zero));
    }

    [Test]
    public void DeterminantMatchesIndependentLeibnizExpansionForSmallMatrices()
    {
        var random = new Random(1729);
        for (var sample = 0; sample < 30; sample++)
        {
            var values = new int[4, 4];
            for (var row = 0; row < 4; row++)
                for (var column = 0; column < 4; column++)
                    values[row, column] = random.Next(-3, 4);

            // Independent permutation formula: sign is the parity of inversions.
            var expected = 0;
            for (var a = 0; a < 4; a++)
                for (var b = 0; b < 4; b++)
                    for (var c = 0; c < 4; c++)
                        for (var d = 0; d < 4; d++)
                        {
                            if (a == b || a == c || a == d || b == c || b == d || c == d) continue;
                            var permutation = new[] { a, b, c, d };
                            var inversions = 0;
                            for (var i = 0; i < 4; i++)
                                for (var j = i + 1; j < 4; j++)
                                    if (permutation[i] > permutation[j]) inversions++;
                            var product = values[0, a] * values[1, b] * values[2, c] * values[3, d];
                            expected += inversions % 2 == 0 ? product : -product;
                        }

            var matrix = From(values);
            AssertDeterminant(matrix.GetDeterminant(), expected, $"Sample {sample}");
            AssertDeterminant(matrix.Transpose().GetDeterminant(), expected, $"Transpose {sample}");
        }
    }

    [Test]
    public void NonSquareMatricesRejectTraceAndDeterminant()
    {
        var matrix = new Matrix<T>(2, 3);
        Assert.Throws<NotSquareMatrixException>(() => matrix.GetTrace());
        Assert.Throws<NotSquareMatrixException>(() => matrix.GetDeterminant());
        Assert.That(matrix.TryGetTrace(out var trace), Is.False);
        Assert.That(trace, Is.EqualTo(T.Zero));
        Assert.That(matrix.TryGetDeterminant(out var determinant), Is.False);
        Assert.That(determinant, Is.EqualTo(T.Zero));
        Assert.That(matrix.IsSquare(), Is.False);
        Assert.That(matrix.IsTriangular(), Is.False);
        Assert.That(matrix.IsUpperTriangular(), Is.False);
        Assert.That(matrix.IsLowerTriangular(), Is.False);
    }

    [Test]
    public void TriangularPredicatesInspectTheCorrectSideOfTheDiagonal()
    {
        var upper = From(new int[,] { { 2, 3, 4 }, { 0, -3, 5 }, { 0, 0, 7 } });
        Assert.That(upper.IsUpperTriangular(), Is.True);
        Assert.That(upper.IsLowerTriangular(), Is.False);
        Assert.That(upper.IsTriangular(), Is.True);
        var lower = upper.Transpose();
        Assert.That(lower.IsUpperTriangular(), Is.False);
        Assert.That(lower.IsLowerTriangular(), Is.True);
        AssertDeterminant(upper.GetDeterminant(), -42);
        AssertDeterminant(lower.GetDeterminant(), -42);
        upper[2, 0] = T.One;
        Assert.That(upper.IsTriangular(), Is.False);
    }
}

public class GenericMatrixNumericBehaviorTests
{
    [Test]
    public void DecimalArithmeticRetainsFractionalValues()
    {
        var matrix = new Matrix<decimal>(new[,] { { 1.5m, 2m }, { 0.25m, 3m } });
        Assert.That(matrix.GetTrace(), Is.EqualTo(4.5m));
        Assert.That(matrix.GetDeterminant(), Is.EqualTo(4m));
        matrix.Scale(0.5m);
        Assert.That(matrix.GetRow(0), Is.EqualTo(new[] { 0.75m, 1m }));
        Assert.That(matrix.GetDeterminant(), Is.EqualTo(1m));
    }

    [Test]
    public void DoubleArithmeticRetainsFractionalValues()
    {
        var left = new Matrix<double>(new[,] { { 0.5, 0.25 }, { -0.5, 1.5 } });
        var right = new Matrix<double>(new[,] { { 2.0 }, { 4.0 } });
        Assert.That((left * right).GetColumn(0), Is.EqualTo(new[] { 2.0, 5.0 }));
        Assert.That(left.GetDeterminant(), Is.EqualTo(0.875));
    }

    [Test]
    public void BigIntegerDeterminantRemainsExactThroughHugeCancellation()
    {
        var huge = BigInteger.Pow(10, 50);
        var matrix = new Matrix<BigInteger>(new[,] { { huge + 1, huge }, { huge, huge - 1 } });
        Assert.That(matrix.GetDeterminant(), Is.EqualTo(BigInteger.MinusOne));
        Assert.That(matrix.GetTrace(), Is.EqualTo(2 * huge));
    }

    [Test]
    public void IntegerArithmeticThrowsInsteadOfWrapping()
    {
        var maximum = new Matrix<int>(new[,] { { int.MaxValue } });
        var one = new Matrix<int>(1, true);
        Assert.Throws<OverflowException>(() => { _ = maximum + one; });
        Assert.Throws<OverflowException>(() => { _ = maximum * new Matrix<int>(new[,] { { 2 } }); });
        Assert.Throws<OverflowException>(() => { _ = new Matrix<int>(new[,] { { int.MinValue } }) - one; });
        Assert.Throws<OverflowException>(() => new Matrix<int>(new[,] { { int.MaxValue, 0 }, { 0, 1 } }).GetTrace());
        Assert.Throws<OverflowException>(() => new Matrix<int>(new[,] { { int.MaxValue, 0 }, { 0, 2 } }).GetDeterminant());
        Assert.Throws<OverflowException>(() => new Matrix<int>(new[,] { { int.MaxValue, 1 } })
            .DotProduct(new Matrix<int>(new[,] { { 1 }, { 1 } })));
    }

    [Test]
    public void OverflowLeavesInPlaceOperationsUnchanged()
    {
        var matrix = new Matrix<int>(new[,] { { 1, int.MaxValue } });
        Assert.Throws<OverflowException>(() => matrix.TryAddMatrix(new Matrix<int>(new[,] { { 2, 1 } })));
        Assert.That(matrix.GetRow(0), Is.EqualTo(new[] { 1, int.MaxValue }));
        Assert.Throws<OverflowException>(() => matrix.Scale(2));
        Assert.That(matrix.GetRow(0), Is.EqualTo(new[] { 1, int.MaxValue }));

        matrix = new Matrix<int>(new[,] { { 1, int.MinValue } });
        Assert.Throws<OverflowException>(() => matrix.TrySubtractMatrix(new Matrix<int>(new[,] { { 2, 1 } })));
        Assert.That(matrix.GetRow(0), Is.EqualTo(new[] { 1, int.MinValue }));
    }

    [Test]
    public void UnsignedNegativeDeterminantsReportOverflow()
    {
        var matrix = new Matrix<uint>(new uint[,] { { 0, 1 }, { 1, 0 } });
        Assert.Throws<OverflowException>(() => matrix.GetDeterminant());
    }

    [Test]
    public void FloatingPointSpecialValuesFollowTheNumericType()
    {
        var matrix = new Matrix<double>(new[,] { { double.NaN, 1.0 }, { 0.0, 1.0 } });
        Assert.That(double.IsNaN(matrix.GetDeterminant()), Is.True);
        Assert.That(matrix.TryGetDeterminant(out var determinant), Is.True);
        Assert.That(double.IsNaN(determinant), Is.True);
        matrix = new Matrix<double>(new[,] { { double.MaxValue } });
        matrix.Scale(2.0);
        Assert.That(matrix[0, 0], Is.EqualTo(double.PositiveInfinity));
    }

    [Test]
    public void TriangularChecksUseExactZerosAndRejectNaNInZeroPositions()
    {
        var matrix = new Matrix<double>(new[,] { { 1.0, 1.0 }, { double.Epsilon, 1.0 } });
        Assert.That(matrix.IsUpperTriangular(), Is.False);
        matrix[1, 0] = double.NaN;
        Assert.That(matrix.IsUpperTriangular(), Is.False);
        matrix[1, 0] = -0.0;
        Assert.That(matrix.IsUpperTriangular(), Is.True);
    }
}
