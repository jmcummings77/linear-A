using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;

namespace linear_A;

/// <summary>A mutable, fixed-size matrix using .NET generic math.</summary>
/// <remarks>
/// Arithmetic uses the checked operators provided by <typeparamref name="T"/>.
/// Built-in bounded integers and decimals throw on overflow; floating-point types
/// retain their rounding, infinity, and NaN behavior. Custom types control their
/// own checked arithmetic. Zero-size dimensions are supported. This type is not
/// thread-safe.
/// </remarks>
public class Matrix<T> : IEnumerable<T> where T : INumber<T>
{
    private readonly T[,] _items;

    /// <summary>Creates the empty 0-by-0 matrix.</summary>
    public Matrix() : this(0, 0) { }

    /// <summary>Creates a zero or identity square matrix.</summary>
    public Matrix(int squareDimensions, bool initializeAsIdentity = false)
        : this(squareDimensions, squareDimensions)
    {
        if (initializeAsIdentity)
            for (var i = 0; i < squareDimensions; i++)
                _items[i, i] = T.One;
    }

    /// <summary>Creates a matrix filled with <see cref="INumberBase{TSelf}.Zero"/>.</summary>
    public Matrix(int rowCount, int columnCount)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(rowCount);
        ArgumentOutOfRangeException.ThrowIfNegative(columnCount);
        RowCount = rowCount;
        ColumnCount = columnCount;
        _items = new T[rowCount, columnCount];
        // Generic numeric types need not have default(T) equal to their zero.
        for (var row = 0; row < rowCount; row++)
            for (var column = 0; column < columnCount; column++)
                _items[row, column] = T.Zero;
    }

    /// <summary>Copies the values of another matrix into independent storage.</summary>
    public Matrix(Matrix<T> matrixToCopy)
    {
        ArgumentNullException.ThrowIfNull(matrixToCopy);
        RowCount = matrixToCopy.RowCount;
        ColumnCount = matrixToCopy.ColumnCount;
        _items = (T[,])matrixToCopy._items.Clone();
    }

    /// <summary>Copies a zero-based, two-dimensional array into independent storage.</summary>
    public Matrix(T[,] values)
    {
        ArgumentNullException.ThrowIfNull(values);
        if (values.GetLowerBound(0) != 0 || values.GetLowerBound(1) != 0)
            throw new ArgumentException("Array indices must start at zero.", nameof(values));
        RowCount = values.GetLength(0);
        ColumnCount = values.GetLength(1);
        _items = (T[,])values.Clone();
    }

    /// <summary>Gets the number of rows.</summary>
    public int RowCount { get; }

    /// <summary>Gets the number of columns.</summary>
    public int ColumnCount { get; }

    /// <summary>Gets whether the matrix has a fixed size.</summary>
    public bool IsFixedSize => true;

    /// <summary>Gets whether the matrix is read-only.</summary>
    public bool IsReadOnly => false;

    /// <summary>Gets or sets a value using zero-based indices.</summary>
    public T this[int rowIndex, int columnIndex]
    {
        get => _items[rowIndex, columnIndex];
        set => _items[rowIndex, columnIndex] = value;
    }

    /// <summary>Adds a matrix in place, returning false for null or unequal dimensions.</summary>
    /// <remarks>Arithmetic exceptions propagate and leave this matrix unchanged.</remarks>
    public bool TryAddMatrix(Matrix<T>? matrixToAdd)
    {
        if (matrixToAdd is null || !HasSameDimensions(matrixToAdd)) return false;
        var result = this + matrixToAdd;
        Array.Copy(result._items, _items, _items.Length);
        return true;
    }

    /// <summary>Subtracts a matrix in place, returning false for null or unequal dimensions.</summary>
    /// <remarks>Arithmetic exceptions propagate and leave this matrix unchanged.</remarks>
    public bool TrySubtractMatrix(Matrix<T>? matrixToSubtract)
    {
        if (matrixToSubtract is null || !HasSameDimensions(matrixToSubtract)) return false;
        var result = this - matrixToSubtract;
        Array.Copy(result._items, _items, _items.Length);
        return true;
    }

    /// <summary>Multiplies every element by a scalar in place.</summary>
    /// <remarks>Arithmetic exceptions leave this matrix unchanged.</remarks>
    public void Scale(T scalar)
    {
        var result = new T[RowCount, ColumnCount];
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < ColumnCount; column++)
                result[row, column] = checked(_items[row, column] * scalar);
        Array.Copy(result, _items, _items.Length);
    }

    /// <summary>Returns the transpose in independent storage.</summary>
    public Matrix<T> Transpose()
    {
        var result = new Matrix<T>(ColumnCount, RowCount);
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < ColumnCount; column++)
                result[column, row] = _items[row, column];
        return result;
    }

    /// <summary>Returns the matrix product of this matrix and <paramref name="multiplicand"/>.</summary>
    public Matrix<T> DotProduct(Matrix<T> multiplicand)
    {
        ArgumentNullException.ThrowIfNull(multiplicand);
        if (ColumnCount != multiplicand.RowCount)
            throw new ArgumentException("The left column count must equal the right row count.", nameof(multiplicand));

        var result = new Matrix<T>(RowCount, multiplicand.ColumnCount);
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < multiplicand.ColumnCount; column++)
            {
                var sum = T.Zero;
                for (var k = 0; k < ColumnCount; k++)
                    sum = checked(sum + checked(_items[row, k] * multiplicand._items[k, column]));
                result._items[row, column] = sum;
            }
        return result;
    }

    /// <summary>Returns the right-handed cross product of two three-component vectors.</summary>
    /// <remarks>
    /// Each operand may be 3-by-1 or 1-by-3; the result retains this matrix's shape.
    /// Inputs are unchanged. Arithmetic is checked, and nonfinite inputs/results fail.
    /// </remarks>
    public Matrix<T> CrossProduct(Matrix<T> multiplicand)
    {
        ArgumentNullException.ThrowIfNull(multiplicand);
        if (!MatrixRotation.IsVector3(RowCount, ColumnCount) || !MatrixRotation.IsVector3(multiplicand.RowCount, multiplicand.ColumnCount))
            throw new ArgumentException("Cross product requires two 3-by-1 or 1-by-3 vectors.", nameof(multiplicand));
        var a = RowCount == 1 ? GetRow(0) : GetColumn(0);
        var b = multiplicand.RowCount == 1 ? multiplicand.GetRow(0) : multiplicand.GetColumn(0);
        for (var i = 0; i < 3; i++)
            if (!T.IsFinite(a[i]) || !T.IsFinite(b[i])) throw new ArgumentException("Cross product requires finite components.");
        var result = new Matrix<T>(RowCount, ColumnCount);
        for (var i = 0; i < 3; i++)
        {
            var j = (i + 1) % 3;
            var k = (i + 2) % 3;
            var value = checked(checked(a[j] * b[k]) - checked(a[k] * b[j]));
            if (!T.IsFinite(value)) throw new ArithmeticException("Cross product produced a nonfinite component.");
            result[RowCount == 1 ? 0 : i, RowCount == 1 ? i : 0] = value;
        }
        return result;
    }

    /// <summary>Returns a copy of one row.</summary>
    public T[] GetRow(int rowNumber)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(rowNumber);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(rowNumber, RowCount);
        var result = new T[ColumnCount];
        for (var column = 0; column < ColumnCount; column++)
            result[column] = _items[rowNumber, column];
        return result;
    }

    /// <summary>Returns a copy of one column.</summary>
    public T[] GetColumn(int columnNumber)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(columnNumber);
        ArgumentOutOfRangeException.ThrowIfGreaterThanOrEqual(columnNumber, ColumnCount);
        var result = new T[RowCount];
        for (var row = 0; row < RowCount; row++)
            result[row] = _items[row, columnNumber];
        return result;
    }

    /// <summary>Returns whether the row and column counts match.</summary>
    public bool IsSquare() => RowCount == ColumnCount;

    /// <summary>Returns whether this square matrix has exact zeros below its diagonal.</summary>
    public bool IsUpperTriangular()
    {
        if (!IsSquare()) return false;
        for (var row = 1; row < RowCount; row++)
            for (var column = 0; column < row; column++)
                if (!T.IsZero(_items[row, column])) return false;
        return true;
    }

    /// <summary>Returns whether this square matrix has exact zeros above its diagonal.</summary>
    public bool IsLowerTriangular()
    {
        if (!IsSquare()) return false;
        for (var row = 0; row < RowCount; row++)
            for (var column = row + 1; column < ColumnCount; column++)
                if (!T.IsZero(_items[row, column])) return false;
        return true;
    }

    /// <summary>Returns whether this matrix is upper or lower triangular using exact zero checks.</summary>
    public bool IsTriangular() => IsUpperTriangular() || IsLowerTriangular();

    /// <summary>Returns the sum of the diagonal; the empty matrix has trace zero.</summary>
    /// <exception cref="NotSquareMatrixException">The matrix is not square.</exception>
    public T GetTrace()
    {
        if (!IsSquare()) throw new NotSquareMatrixException();
        var result = T.Zero;
        for (var i = 0; i < RowCount; i++)
            result = checked(result + _items[i, i]);
        return result;
    }

    /// <summary>Returns false with trace zero for a non-square matrix.</summary>
    /// <remarks>Arithmetic exceptions propagate; floating-point NaN and infinity follow the numeric type.</remarks>
    public bool TryGetTrace(out T trace)
    {
        trace = T.Zero;
        if (!IsSquare()) return false;
        trace = GetTrace();
        return true;
    }

    /// <summary>Returns the determinant; the empty matrix has determinant one.</summary>
    /// <remarks>
    /// Auto uses cubic-time elimination for built-in numeric types: exact BigInteger
    /// Bareiss for integers and decimal, and partial-pivoting LU with double working
    /// precision for double, float, and Half. Integers convert back with overflow checks;
    /// decimal rounds only the final result to nearest, ties to even. Custom types use
    /// checked cofactor expansion. Floating results remain subject to rounding and
    /// conditioning; nonfinite LU inputs or elimination values produce NaN. A nonzero
    /// determinant is not a guarantee that an inverse can be represented by T.
    /// </remarks>
    /// <exception cref="NotSquareMatrixException">The matrix is not square.</exception>
    public T GetDeterminant() => GetDeterminant(DeterminantAlgorithm.Auto);

    /// <summary>Returns the determinant using the selected algorithm without changing this matrix.</summary>
    /// <remarks>
    /// Cofactor uses checked arithmetic with factorial cost and can overflow in intermediate
    /// steps. Bareiss uses exact BigInteger intermediates for built-in integers and decimal;
    /// its cubic arithmetic-operation count excludes the cost of growing integer operands.
    /// LU and Cholesky support only double, float, and Half. Cholesky additionally requires
    /// finite, exactly symmetric input and finite positive pivots; floating-point rounding
    /// can cause it to reject an ill-conditioned positive-definite matrix. Auto chooses
    /// according to the numeric type and does not automatically select Cholesky.
    /// </remarks>
    /// <exception cref="NotSupportedException">The algorithm does not support T.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The algorithm is not a defined value.</exception>
    /// <exception cref="NotSquareMatrixException">The matrix is not square.</exception>
    /// <exception cref="InvalidOperationException">Cholesky's symmetry, finiteness, or positive-pivot requirements failed.</exception>
    public T GetDeterminant(DeterminantAlgorithm algorithm)
    {
        algorithm = MatrixDeterminant<T>.Resolve(algorithm);
        if (!IsSquare()) throw new NotSquareMatrixException();
        return algorithm switch
        {
            DeterminantAlgorithm.Lu => MatrixDeterminant<T>.Lu(_items),
            DeterminantAlgorithm.Bareiss => MatrixDeterminant<T>.Bareiss(_items),
            DeterminantAlgorithm.Cholesky => MatrixDeterminant<T>.Cholesky(_items),
            _ => CofactorDeterminant(_items),
        };
    }

    /// <summary>Returns false with determinant zero for a non-square matrix.</summary>
    /// <remarks>Uses Auto; arithmetic exceptions propagate.</remarks>
    public bool TryGetDeterminant(out T determinant) => TryGetDeterminant(DeterminantAlgorithm.Auto, out determinant);

    /// <summary>Returns false with determinant zero for a non-square matrix, using the selected algorithm.</summary>
    /// <remarks>Invalid algorithms, unsupported numeric types, and arithmetic exceptions propagate.</remarks>
    public bool TryGetDeterminant(DeterminantAlgorithm algorithm, out T determinant)
    {
        determinant = T.Zero;
        MatrixDeterminant<T>.Resolve(algorithm);
        if (!IsSquare()) return false;
        determinant = GetDeterminant(algorithm);
        return true;
    }

    /// <summary>Computes real eigenvalues and orthonormal eigenvector columns for a symmetric matrix.</summary>
    /// <remarks>
    /// Entries must be finite and exactly symmetric before conversion to double.
    /// Numeric values are approximated as doubles; conversion can lose precision.
    /// The Jacobi iteration uses a relative Frobenius-norm tolerance, so tiny
    /// eigenvalues in a matrix with much larger entries need not have small relative error.
    /// The input is unchanged. Nonconvergence or unrepresentable eigenvalues throw.
    /// </remarks>
    public SymmetricEigenDecomposition GetSymmetricEigenDecomposition(double tolerance = 1e-12, int maxSweeps = 50)
    {
        if (!IsSquare()) throw new NotSquareMatrixException();
        var values = new double[RowCount, ColumnCount];
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < ColumnCount; column++)
            {
                if (!T.IsFinite(_items[row, column]))
                    throw new ArgumentException("Eigendecomposition requires finite entries.");
                if (_items[row, column] != _items[column, row])
                    throw new ArgumentException("Eigendecomposition requires an exactly symmetric matrix.");
                values[row, column] = double.CreateChecked(_items[row, column]);
            }
        return SymmetricEigenSolver.Solve(values, tolerance, maxSweeps);
    }

    /// <summary>Computes approximate complex right eigenpairs of a real square matrix.</summary>
    /// <remarks>Entries are converted to double and can lose precision. Values sort
    /// by real part then imaginary part; normalized right eigenvectors are columns.
    /// No orthogonality or independent basis is promised for defective matrices.
    /// Nonconvergence and nonfinite input/results throw. The bounded QR iteration
    /// limit applies per unresolved root and must be in 1..100000.</remarks>
    public EigenDecomposition GetEigenDecomposition(int maxIterations = 1000)
    {
        if (!IsSquare()) throw new NotSquareMatrixException();
        var values = new double[RowCount, ColumnCount];
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < ColumnCount; column++) values[row, column] = double.CreateChecked(_items[row, column]);
        return GeneralEigenSolver.Solve(values, maxIterations);
    }

    private static T CofactorDeterminant(T[,] values)
    {
        var size = values.GetLength(0);
        if (size == 0) return T.One;
        if (size == 1) return values[0, 0];
        if (size == 2)
            return checked(checked(values[0, 0] * values[1, 1]) - checked(values[0, 1] * values[1, 0]));

        var result = T.Zero;
        var minor = new T[size - 1, size - 1];
        for (var excludedColumn = 0; excludedColumn < size; excludedColumn++)
        {
            for (var row = 1; row < size; row++)
            {
                var minorColumn = 0;
                for (var column = 0; column < size; column++)
                {
                    if (column == excludedColumn) continue;
                    minor[row - 1, minorColumn++] = values[row, column];
                }
            }
            var term = checked(values[0, excludedColumn] * CofactorDeterminant(minor));
            result = excludedColumn % 2 == 0 ? checked(result + term) : checked(result - term);
        }
        return result;
    }

    /// <summary>Returns the sum without modifying either operand.</summary>
    public static Matrix<T> operator +(Matrix<T> left, Matrix<T> right) => Combine(left, right, false);

    /// <summary>Returns the difference without modifying either operand.</summary>
    public static Matrix<T> operator -(Matrix<T> left, Matrix<T> right) => Combine(left, right, true);

    /// <summary>Returns the matrix product without modifying either operand.</summary>
    public static Matrix<T> operator *(Matrix<T> left, Matrix<T> right)
    {
        ArgumentNullException.ThrowIfNull(left);
        return left.DotProduct(right);
    }

    private static Matrix<T> Combine(Matrix<T> left, Matrix<T> right, bool subtract)
    {
        ArgumentNullException.ThrowIfNull(left);
        ArgumentNullException.ThrowIfNull(right);
        if (!left.HasSameDimensions(right))
            throw new ArgumentException("Matrix dimensions must match.", nameof(right));

        var result = new Matrix<T>(left.RowCount, left.ColumnCount);
        for (var row = 0; row < left.RowCount; row++)
            for (var column = 0; column < left.ColumnCount; column++)
                result._items[row, column] = subtract
                    ? checked(left._items[row, column] - right._items[row, column])
                    : checked(left._items[row, column] + right._items[row, column]);
        return result;
    }

    private bool HasSameDimensions(Matrix<T> other) =>
        RowCount == other.RowCount && ColumnCount == other.ColumnCount;

    /// <summary>Enumerates elements one row at a time.</summary>
    public IEnumerator<T> GetEnumerator()
    {
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < ColumnCount; column++)
                yield return _items[row, column];
    }

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
}
