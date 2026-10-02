using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace linear_A.Test;

public class LegacyMatrixCompatibilityTests
{
    [Test]
    public void NullCopyAndDefaultConstructorsRetainEmptyMatrixBehavior()
    {
        foreach (var matrix in new[] { new IntegerMatrix(), new IntegerMatrix((IntegerMatrix?)null) })
        {
            Assert.That(matrix.RowCount, Is.Zero);
            Assert.That(matrix.ColumnCount, Is.Zero);
            Assert.That(matrix.GetTrace(), Is.Zero);
            Assert.That(matrix.GetDeterminant(), Is.EqualTo(1));
            Assert.That(matrix.IsInvertible(), Is.True);
        }
        foreach (var matrix in new[] { new DoubleMatrix(), new DoubleMatrix((DoubleMatrix?)null) })
        {
            Assert.That(matrix.RowCount, Is.Zero);
            Assert.That(matrix.ColumnCount, Is.Zero);
            Assert.That(matrix.GetTrace(), Is.Zero);
            Assert.That(matrix.GetDeterminant(), Is.EqualTo(1d));
            Assert.That(matrix.IsInvertible(), Is.True);
        }
    }

    [TestCase(0, 0)]
    [TestCase(0, 3)]
    [TestCase(3, 0)]
    [TestCase(-1, 2)]
    [TestCase(2, -1)]
    public void DimensionConstructorsContinueToRequirePositiveDimensions(int rows, int columns)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new IntegerMatrix(rows, columns));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DoubleMatrix(rows, columns));
        Assert.Throws<ArgumentOutOfRangeException>(() => new IntegerMatrix(0, false));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DoubleMatrix(0, true));
    }

    [Test]
    public void BasicIntegerArithmeticRemainsUncheckedWhileGenericMathRemainsChecked()
    {
        var maximum = IntegerMatrix.FromArray(new[,] { { int.MaxValue } });
        var one = IntegerMatrix.FromArray(new[,] { { 1 } });
        var two = IntegerMatrix.FromArray(new[,] { { 2 } });
        Assert.That((maximum + one)[0, 0], Is.EqualTo(int.MinValue));
        Assert.That((IntegerMatrix.FromArray(new[,] { { int.MinValue } }) - one)[0, 0], Is.EqualTo(int.MaxValue));
        Assert.That((maximum * two)[0, 0], Is.EqualTo(-2));
        Assert.That(IntegerMatrix.FromArray(new[,] { { int.MaxValue, 0 }, { 0, 1 } }).GetTrace(), Is.EqualTo(int.MinValue));
        var scaled = new IntegerMatrix(maximum);
        scaled.Scale(2);
        Assert.That(scaled[0, 0], Is.EqualTo(-2));
        var added = new IntegerMatrix(maximum);
        Assert.That(added.TryAddMatrix(one), Is.True);
        Assert.That(added[0, 0], Is.EqualTo(int.MinValue));
        Assert.That(added.TrySubtractMatrix(one), Is.True);
        Assert.That(added[0, 0], Is.EqualTo(int.MaxValue));
        Assert.That(maximum[0, 0], Is.EqualTo(int.MaxValue));

        var generic = maximum.ToMatrix();
        Assert.Throws<OverflowException>(() => { _ = generic + one.ToMatrix(); });
        Assert.Throws<OverflowException>(() => { _ = generic * two.ToMatrix(); });
        Assert.Throws<OverflowException>(() => generic.Scale(2));
        Assert.Throws<OverflowException>(() => generic.TryAddMatrix(one.ToMatrix()));
        Assert.That(generic[0, 0], Is.EqualTo(int.MaxValue));
        Assert.Throws<OverflowException>(() => IntegerMatrix.FromArray(new[,] { { int.MaxValue }, { 0 }, { 0 } })
            .GetCrossProduct(IntegerMatrix.FromArray(new[,] { { 0 }, { 2 }, { 0 } })));
    }

    [Test]
    public void PublicAndTypedEnumerationPreserveRowMajorOrderAndOriginalSignature()
    {
        var integers = IntegerMatrix.FromArray(new[,] { { 1, 2, 3 }, { 4, 5, 6 } });
        var doubles = DoubleMatrix.FromArray(new[,] { { 0.5, -1d, 2d }, { 3d, 4d, 5d } });
        Assert.That(typeof(IntegerMatrix).GetMethod(nameof(IntegerMatrix.GetEnumerator))!.ReturnType, Is.EqualTo(typeof(IEnumerator)));
        Assert.That(typeof(DoubleMatrix).GetMethod(nameof(DoubleMatrix.GetEnumerator))!.ReturnType, Is.EqualTo(typeof(IEnumerator)));
        Assert.That(Read(integers.GetEnumerator()), Is.EqualTo(new object[] { 1, 2, 3, 4, 5, 6 }));
        Assert.That(Read(doubles.GetEnumerator()), Is.EqualTo(new object[] { 0.5, -1d, 2d, 3d, 4d, 5d }));
        Assert.That(((IEnumerable<int>)integers).ToArray(), Is.EqualTo(new[] { 1, 2, 3, 4, 5, 6 }));
        Assert.That(((IEnumerable<double>)doubles).ToArray(), Is.EqualTo(new[] { 0.5, -1d, 2d, 3d, 4d, 5d }));
        Assert.That(((IEnumerable)integers).Cast<int>().ToArray(), Is.EqualTo(new[] { 1, 2, 3, 4, 5, 6 }));
        Assert.That(((IEnumerable<int>)IntegerMatrix.FromArray(new int[0, 3])).Any(), Is.False);
        Assert.That(((IEnumerable<double>)DoubleMatrix.FromArray(new double[3, 0])).Any(), Is.False);
    }

    private static object[] Read(IEnumerator enumerator)
    {
        var values = new List<object>();
        while (enumerator.MoveNext()) values.Add(enumerator.Current!);
        return values.ToArray();
    }

    [Test]
    public void ArrayFactoriesCopiesAndConversionsOwnIndependentStorage()
    {
        var source = new[,] { { 1, 2 }, { 3, 4 } };
        var legacy = IntegerMatrix.FromArray(source);
        var array = legacy.ToArray();
        var copy = new IntegerMatrix(legacy);
        var generic = legacy.ToMatrix();
        var returned = IntegerMatrix.FromMatrix(generic);
        source[0, 0] = 11;
        array[0, 0] = 12;
        copy[0, 0] = 13;
        generic[0, 0] = 14;
        returned[0, 0] = 15;
        Assert.That(legacy[0, 0], Is.EqualTo(1));
        Assert.That(source[0, 0], Is.EqualTo(11));
        Assert.That(array[0, 0], Is.EqualTo(12));
        Assert.That(copy[0, 0], Is.EqualTo(13));
        Assert.That(generic[0, 0], Is.EqualTo(14));
        Assert.That(returned[0, 0], Is.EqualTo(15));

        var doubles = new[,] { { 0.5, 2d }, { 3d, 4d } };
        var floating = DoubleMatrix.FromArray(doubles);
        var floatingArray = floating.ToArray();
        var floatingCopy = new DoubleMatrix(floating);
        var floatingGeneric = floating.ToMatrix();
        var floatingReturned = DoubleMatrix.FromMatrix(floatingGeneric);
        doubles[0, 0] = 11;
        floatingArray[0, 0] = 12;
        floatingCopy[0, 0] = 13;
        floatingGeneric[0, 0] = 14;
        floatingReturned[0, 0] = 15;
        Assert.That(floating[0, 0], Is.EqualTo(0.5));
        Assert.That(doubles[0, 0], Is.EqualTo(11));
        Assert.That(floatingArray[0, 0], Is.EqualTo(12));
        Assert.That(floatingCopy[0, 0], Is.EqualTo(13));
        Assert.That(floatingGeneric[0, 0], Is.EqualTo(14));
        Assert.That(floatingReturned[0, 0], Is.EqualTo(15));
    }

    [TestCase(0, 0)]
    [TestCase(0, 3)]
    [TestCase(3, 0)]
    public void FactoriesAndConversionsPreserveEmptyRectangularShapes(int rows, int columns)
    {
        var integers = IntegerMatrix.FromArray(new int[rows, columns]);
        var doubles = DoubleMatrix.FromArray(new double[rows, columns]);
        foreach (var copy in new[] { integers, new IntegerMatrix(integers), IntegerMatrix.FromMatrix(integers.ToMatrix()) })
        {
            Assert.That((copy.RowCount, copy.ColumnCount), Is.EqualTo((rows, columns)));
            Assert.That((copy.ToArray().GetLength(0), copy.ToArray().GetLength(1)), Is.EqualTo((rows, columns)));
            Assert.That((copy.Transpose().RowCount, copy.Transpose().ColumnCount), Is.EqualTo((columns, rows)));
            Assert.That((copy + integers).ToArray(), Is.EqualTo(new int[rows, columns]));
            if (rows > 0) Assert.That(copy.GetRow(0), Is.Empty);
            if (columns > 0) Assert.That(copy.GetColumn(0), Is.Empty);
        }
        foreach (var copy in new[] { doubles, new DoubleMatrix(doubles), DoubleMatrix.FromMatrix(doubles.ToMatrix()) })
        {
            Assert.That((copy.RowCount, copy.ColumnCount), Is.EqualTo((rows, columns)));
            Assert.That((copy.ToArray().GetLength(0), copy.ToArray().GetLength(1)), Is.EqualTo((rows, columns)));
            Assert.That((copy.Transpose().RowCount, copy.Transpose().ColumnCount), Is.EqualTo((columns, rows)));
            Assert.That((copy + doubles).ToArray(), Is.EqualTo(new double[rows, columns]));
            if (rows > 0) Assert.That(copy.GetRow(0), Is.Empty);
            if (columns > 0) Assert.That(copy.GetColumn(0), Is.Empty);
        }
        Assert.That(IntegerMatrix.FromMatrix(new Matrix<int>(rows, columns)).ToArray(), Is.EqualTo(new int[rows, columns]));
        Assert.That(DoubleMatrix.FromMatrix(new Matrix<double>(rows, columns)).ToArray(), Is.EqualTo(new double[rows, columns]));
    }

    [Test]
    public void EmptyInnerDimensionProductsHaveTheCorrectNonemptyOutputShape()
    {
        var integers = IntegerMatrix.FromArray(new int[2, 0]) * IntegerMatrix.FromArray(new int[0, 3]);
        var doubles = DoubleMatrix.FromArray(new double[2, 0]) * DoubleMatrix.FromArray(new double[0, 3]);
        Assert.That((integers.RowCount, integers.ColumnCount), Is.EqualTo((2, 3)));
        Assert.That((doubles.RowCount, doubles.ColumnCount), Is.EqualTo((2, 3)));
        Assert.That(integers.ToArray(), Is.EqualTo(new int[2, 3]));
        Assert.That(doubles.ToArray(), Is.EqualTo(new double[2, 3]));
        var empty = DoubleMatrix.FromArray(new double[0, 2]) * new DoubleMatrix(2, 3);
        Assert.That((empty.RowCount, empty.ColumnCount), Is.EqualTo((0, 3)));
    }

    [Test]
    public void DoubleConversionsPreserveSpecialValuesAndTheirBits()
    {
        var bits = new[] { 0L, long.MinValue, 1L, 0x7ff0000000000000L,
            unchecked((long)0xfff0000000000000UL), 0x7ff8000000000042L };
        var source = new double[2, 3];
        for (var i = 0; i < bits.Length; i++) source[i / 3, i % 3] = BitConverter.Int64BitsToDouble(bits[i]);
        var legacy = DoubleMatrix.FromArray(source);
        var generic = legacy.ToMatrix();
        foreach (var array in new[] { legacy.ToArray(), new DoubleMatrix(legacy).ToArray(), DoubleMatrix.FromMatrix(generic).ToArray() })
            for (var i = 0; i < bits.Length; i++)
                Assert.That(BitConverter.DoubleToInt64Bits(array[i / 3, i % 3]), Is.EqualTo(bits[i]));
        for (var i = 0; i < bits.Length; i++)
            Assert.That(BitConverter.DoubleToInt64Bits(generic[i / 3, i % 3]), Is.EqualTo(bits[i]));
    }

    [Test]
    public void NewFactoriesRejectNullAndNonZeroArrayLowerBounds()
    {
        Assert.Throws<ArgumentNullException>(() => IntegerMatrix.FromArray(null!));
        Assert.Throws<ArgumentNullException>(() => DoubleMatrix.FromArray(null!));
        Assert.Throws<ArgumentNullException>(() => IntegerMatrix.FromMatrix(null!));
        Assert.Throws<ArgumentNullException>(() => DoubleMatrix.FromMatrix(null!));
        var integers = (int[,])Array.CreateInstance(typeof(int), new[] { 2, 2 }, new[] { 1, 1 });
        var doubles = (double[,])Array.CreateInstance(typeof(double), new[] { 2, 2 }, new[] { 1, 1 });
        Assert.Throws<ArgumentException>(() => IntegerMatrix.FromArray(integers));
        Assert.Throws<ArgumentException>(() => DoubleMatrix.FromArray(doubles));
    }
}
