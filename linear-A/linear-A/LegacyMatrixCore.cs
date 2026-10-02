using System;
using System.Collections.Generic;

namespace linear_A;

// Shared implementation for the two compatibility APIs. Arithmetic policies keep
// their historical overflow rules available on .NET Standard without generic math.
internal sealed class LegacyMatrixCore<T, TArithmetic>
    where T : struct
    where TArithmetic : struct, ILegacyArithmetic<T>
{
    private static readonly TArithmetic Arithmetic = default;

    internal T[,] Items { get; }
    internal int RowCount => Items.GetLength(0);
    internal int ColumnCount => Items.GetLength(1);

    internal LegacyMatrixCore(int rows, int columns) => Items = new T[rows, columns];

    private LegacyMatrixCore(T[,] ownedItems) => Items = ownedItems;

    internal static LegacyMatrixCore<T, TArithmetic> Square(int squareDimensions, bool initializeAsIdentity)
    {
        if (squareDimensions < 1)
            throw new ArgumentOutOfRangeException(nameof(squareDimensions), "Dimensions must be greater than zero");
        var result = new LegacyMatrixCore<T, TArithmetic>(squareDimensions, squareDimensions);
        if (initializeAsIdentity)
            for (var i = 0; i < squareDimensions; i++) result.Items[i, i] = Arithmetic.One;
        return result;
    }

    internal static LegacyMatrixCore<T, TArithmetic> Rectangle(int rowCount, int columnCount)
    {
        if (rowCount < 1 || columnCount < 1)
        {
            var parameter = rowCount < 1 && columnCount < 1 ? "rowCount & columnCount"
                : rowCount < 1 ? nameof(rowCount) : nameof(columnCount);
            throw new ArgumentOutOfRangeException(parameter, "Dimensions must be greater than zero");
        }
        return new LegacyMatrixCore<T, TArithmetic>(rowCount, columnCount);
    }

    internal static LegacyMatrixCore<T, TArithmetic> FromArray(T[,] values)
    {
        if (values == null) throw new ArgumentNullException(nameof(values));
        if (values.GetLowerBound(0) != 0 || values.GetLowerBound(1) != 0)
            throw new ArgumentException("Array indices must start at zero.", nameof(values));
        return new LegacyMatrixCore<T, TArithmetic>((T[,])values.Clone());
    }

    internal LegacyMatrixCore<T, TArithmetic> Copy() => FromArray(Items);

    internal bool TryCombine(LegacyMatrixCore<T, TArithmetic>? other, bool subtract)
    {
        if (other == null || RowCount != other.RowCount || ColumnCount != other.ColumnCount) return false;
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < ColumnCount; column++)
                Items[row, column] = subtract
                    ? Arithmetic.Subtract(Items[row, column], other.Items[row, column])
                    : Arithmetic.Add(Items[row, column], other.Items[row, column]);
        return true;
    }

    internal void Scale(T scalar)
    {
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < ColumnCount; column++)
                Items[row, column] = Arithmetic.Multiply(Items[row, column], scalar);
    }

    internal LegacyMatrixCore<T, TArithmetic> Transpose()
    {
        var result = new LegacyMatrixCore<T, TArithmetic>(ColumnCount, RowCount);
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < ColumnCount; column++)
                result.Items[column, row] = Items[row, column];
        return result;
    }

    internal LegacyMatrixCore<T, TArithmetic> DotProduct(LegacyMatrixCore<T, TArithmetic>? multiplicand)
    {
        if (multiplicand == null) throw new ArgumentNullException(nameof(multiplicand));
        if (ColumnCount != multiplicand.RowCount)
            throw new ArgumentOutOfRangeException(nameof(multiplicand), "The left column count must equal the right row count.");
        var result = new LegacyMatrixCore<T, TArithmetic>(RowCount, multiplicand.ColumnCount);
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < multiplicand.ColumnCount; column++)
                for (var k = 0; k < ColumnCount; k++)
                    result.Items[row, column] = Arithmetic.Add(result.Items[row, column],
                        Arithmetic.Multiply(Items[row, k], multiplicand.Items[k, column]));
        return result;
    }

    internal LegacyMatrixCore<T, TArithmetic> GetCrossProduct(LegacyMatrixCore<T, TArithmetic>? multiplicand)
    {
        if (multiplicand == null) throw new ArgumentNullException(nameof(multiplicand));
        if (!MatrixRotation.IsVector3(RowCount, ColumnCount)
            || !MatrixRotation.IsVector3(multiplicand.RowCount, multiplicand.ColumnCount))
            throw new ArgumentException("Cross product requires two 3-by-1 or 1-by-3 vectors.", nameof(multiplicand));
        var a = RowCount == 1 ? GetRow(0) : GetColumn(0);
        var b = multiplicand.RowCount == 1 ? multiplicand.GetRow(0) : multiplicand.GetColumn(0);
        for (var i = 0; i < 3; i++)
            if (!Arithmetic.IsFinite(a[i]) || !Arithmetic.IsFinite(b[i]))
                throw new ArgumentException("Cross product requires finite components.");

        var result = new LegacyMatrixCore<T, TArithmetic>(RowCount, ColumnCount);
        for (var i = 0; i < 3; i++)
        {
            var j = (i + 1) % 3;
            var k = (i + 2) % 3;
            var value = Arithmetic.CrossComponent(a[j], b[k], a[k], b[j]);
            if (!Arithmetic.IsFinite(value)) throw new ArithmeticException("Cross product produced a nonfinite component.");
            result.Items[RowCount == 1 ? 0 : i, RowCount == 1 ? i : 0] = value;
        }
        return result;
    }

    internal bool CrossProduct(LegacyMatrixCore<T, TArithmetic>? multiplicand)
    {
        try
        {
            var result = GetCrossProduct(multiplicand);
            Array.Copy(result.Items, Items, Items.Length);
            return true;
        }
        catch (ArgumentException) { return false; }
        catch (ArithmeticException) { return false; }
    }

    internal T[] GetRow(int rowNumber)
    {
        if (rowNumber < 0 || rowNumber >= RowCount) throw new ArgumentOutOfRangeException(nameof(rowNumber));
        var result = new T[ColumnCount];
        for (var column = 0; column < ColumnCount; column++) result[column] = Items[rowNumber, column];
        return result;
    }

    internal T[] GetColumn(int columnNumber)
    {
        if (columnNumber < 0 || columnNumber >= ColumnCount) throw new ArgumentOutOfRangeException(nameof(columnNumber));
        var result = new T[RowCount];
        for (var row = 0; row < RowCount; row++) result[row] = Items[row, columnNumber];
        return result;
    }

    internal bool IsSquare() => RowCount == ColumnCount;

    internal bool IsTriangular(bool upper)
    {
        if (!IsSquare()) return false;
        for (var row = 0; row < RowCount; row++)
            for (var column = row + 1; column < ColumnCount; column++)
                if (!Arithmetic.IsZero(upper ? Items[column, row] : Items[row, column])) return false;
        return true;
    }

    internal T GetTrace()
    {
        if (!IsSquare()) throw new NotSquareMatrixException();
        var result = default(T);
        for (var i = 0; i < RowCount; i++) result = Arithmetic.Add(result, Items[i, i]);
        return result;
    }

    internal bool TryGetTrace(out T trace)
    {
        trace = default;
        if (!IsSquare()) return false;
        trace = GetTrace();
        return true;
    }

    internal IEnumerator<T> GetTypedEnumerator()
    {
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < ColumnCount; column++)
                yield return Items[row, column];
    }
}

internal interface ILegacyArithmetic<T>
{
    T One { get; }
    T Add(T left, T right);
    T Subtract(T left, T right);
    T Multiply(T left, T right);
    T CrossComponent(T a, T b, T c, T d);
    bool IsZero(T value);
    bool IsFinite(T value);
}

internal readonly struct LegacyIntegerArithmetic : ILegacyArithmetic<int>
{
    public int One => 1;
    public int Add(int left, int right) => unchecked(left + right);
    public int Subtract(int left, int right) => unchecked(left - right);
    public int Multiply(int left, int right) => unchecked(left * right);
    public int CrossComponent(int a, int b, int c, int d) => checked(a * b - c * d);
    public bool IsZero(int value) => value == 0;
    public bool IsFinite(int value) => true;
}

internal readonly struct LegacyDoubleArithmetic : ILegacyArithmetic<double>
{
    public double One => 1d;
    public double Add(double left, double right) => left + right;
    public double Subtract(double left, double right) => left - right;
    public double Multiply(double left, double right) => left * right;
    public double CrossComponent(double a, double b, double c, double d) => a * b - c * d;
    public bool IsZero(double value) => value == 0d;
    public bool IsFinite(double value) => MatrixRotation.IsFinite(value);
}
