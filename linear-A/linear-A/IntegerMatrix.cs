// Copyright John-Michael Cummings 2021.

using System;
using System.Collections;
using System.Collections.Generic;

namespace linear_A;

/// <summary>A mutable, fixed-size matrix of 32-bit integers.</summary>
/// <remarks>
/// Available on .NET Standard 2.1 and .NET 10. Prefer Matrix&lt;int&gt; for new .NET 10 code.
/// Addition, subtraction, scaling, multiplication, and trace use unchecked int arithmetic.
/// Determinants use exact intermediates; cross products use checked arithmetic.
/// Instances are not thread-safe.
/// </remarks>
public class IntegerMatrix : IEnumerable<int>
{
    private readonly LegacyMatrixCore<int, LegacyIntegerArithmetic> _core;

    /// <summary>Creates the empty 0-by-0 matrix.</summary>
    public IntegerMatrix() : this(new LegacyMatrixCore<int, LegacyIntegerArithmetic>(0, 0)) { }

    /// <summary>Creates a zero or identity square matrix with strictly positive dimensions.</summary>
    public IntegerMatrix(int squareDimensions, bool initializeAsIdentity)
        : this(LegacyMatrixCore<int, LegacyIntegerArithmetic>.Square(squareDimensions, initializeAsIdentity)) { }

    /// <summary>Creates a zero matrix with strictly positive row and column counts.</summary>
    public IntegerMatrix(int rowCount, int columnCount)
        : this(LegacyMatrixCore<int, LegacyIntegerArithmetic>.Rectangle(rowCount, columnCount)) { }

    /// <summary>Copies a matrix into independent storage; null creates an empty 0-by-0 matrix.</summary>
    public IntegerMatrix(IntegerMatrix? matrixToCopy)
        : this(matrixToCopy?._core.Copy() ?? new LegacyMatrixCore<int, LegacyIntegerArithmetic>(0, 0)) { }

    private IntegerMatrix(LegacyMatrixCore<int, LegacyIntegerArithmetic> core) => _core = core;

    /// <summary>Copies a zero-based array, preserving its shape, including zero-size dimensions.</summary>
    public static IntegerMatrix FromArray(int[,] values)
        => new(LegacyMatrixCore<int, LegacyIntegerArithmetic>.FromArray(values));

    /// <summary>Returns an independent two-dimensional copy of the entries.</summary>
    public int[,] ToArray() => (int[,])_core.Items.Clone();

#if NET10_0_OR_GREATER
    /// <summary>Copies the entries and shape into a generic matrix.</summary>
    public Matrix<int> ToMatrix() => new(_core.Items);

    /// <summary>Copies a generic matrix, preserving its shape, including zero-size dimensions.</summary>
    public static IntegerMatrix FromMatrix(Matrix<int> matrix)
    {
        if (matrix == null) throw new ArgumentNullException(nameof(matrix));
        var core = new LegacyMatrixCore<int, LegacyIntegerArithmetic>(matrix.RowCount, matrix.ColumnCount);
        for (var row = 0; row < matrix.RowCount; row++)
            for (var column = 0; column < matrix.ColumnCount; column++)
                core.Items[row, column] = matrix[row, column];
        return new IntegerMatrix(core);
    }
#endif

    /// <summary>Gets the number of rows.</summary>
    public int RowCount => _core.RowCount;

    /// <summary>Gets the number of columns.</summary>
    public int ColumnCount => _core.ColumnCount;

    /// <summary>Gets whether the matrix has a fixed size.</summary>
    public bool IsFixedSize => true;

    /// <summary>Gets whether the matrix is read-only.</summary>
    public bool IsReadOnly => false;

    /// <summary>Gets or sets a value using zero-based indices.</summary>
    public int this[int rowIndex, int columnIndex]
    {
        get => _core.Items[rowIndex, columnIndex];
        set => _core.Items[rowIndex, columnIndex] = value;
    }

    /// <summary>Adds in place; null or unequal shapes return false without changing this matrix.</summary>
    public bool TryAddMatrix(IntegerMatrix? matrixToAdd) => _core.TryCombine(matrixToAdd?._core, false);

    /// <summary>Subtracts in place; null or unequal shapes return false without changing this matrix.</summary>
    public bool TrySubtractMatrix(IntegerMatrix? matrixToSubtract) => _core.TryCombine(matrixToSubtract?._core, true);

    /// <summary>Returns an independent transpose with swapped dimensions.</summary>
    public IntegerMatrix Transpose() => new(_core.Transpose());

    /// <summary>Multiplies each entry by the scalar in place.</summary>
    public void Scale(int scalar) => _core.Scale(scalar);

    /// <summary>Returns the matrix product without modifying either operand.</summary>
    /// <remarks>The left column count must equal the right row count.</remarks>
    public IntegerMatrix DotProduct(IntegerMatrix multiplicand) => new(_core.DotProduct(multiplicand?._core));

    /// <summary>Replaces this vector with its right-handed cross product on success.</summary>
    /// <remarks>Null, invalid shapes, nonfinite entries/results, or checked overflow return false without mutation.</remarks>
    public bool CrossProduct(IntegerMatrix? multiplicand) => _core.CrossProduct(multiplicand?._core);

    /// <summary>Returns an independent right-handed cross product, retaining this vector's shape.</summary>
    /// <remarks>Each operand must be a finite 3-by-1 or 1-by-3 vector. Invalid input or checked overflow throws.</remarks>
    public IntegerMatrix GetCrossProduct(IntegerMatrix multiplicand) => new(_core.GetCrossProduct(multiplicand?._core));

    /// <summary>Returns an independent copy of the requested row.</summary>
    public int[] GetRow(int rowNumber) => _core.GetRow(rowNumber);

    /// <summary>Returns an independent copy of the requested column.</summary>
    public int[] GetColumn(int columnNumber) => _core.GetColumn(columnNumber);

    /// <summary>Returns whether row and column counts match.</summary>
    public bool IsSquare() => _core.IsSquare();

    /// <summary>Checks whether a square matrix is upper or lower triangular using exact zero comparisons.</summary>
    public bool IsTriangular() => IsUpperTriangular() || IsLowerTriangular();

    /// <summary>Checks whether every entry below the diagonal of a square matrix is exactly zero.</summary>
    public bool IsUpperTriangular() => _core.IsTriangular(true);

    /// <summary>Checks whether every entry above the diagonal of a square matrix is exactly zero.</summary>
    public bool IsLowerTriangular() => _core.IsTriangular(false);

    /// <summary>Returns the diagonal sum; nonsquare matrices throw NotSquareMatrixException.</summary>
    public int GetTrace() => _core.GetTrace();

    /// <summary>Returns false and zero for a nonsquare matrix; otherwise returns the diagonal sum.</summary>
    public bool TryGetTrace(out int trace) => _core.TryGetTrace(out trace);

    /// <summary>Computes real eigenvalues and orthonormal eigenvector columns as doubles.</summary>
    /// <remarks>
    /// Requires exact symmetry. Integer matrices can have noninteger eigenvalues
    /// and eigenvectors. The input is unchanged; this method uses approximate
    /// double arithmetic with a relative Frobenius-norm tolerance.
    /// </remarks>
    public SymmetricEigenDecomposition GetSymmetricEigenDecomposition(double tolerance = 1e-12, int maxSweeps = 50)
    {
        if (!IsSquare()) throw new NotSquareMatrixException();
        var values = new double[RowCount, ColumnCount];
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < ColumnCount; column++)
                values[row, column] = _core.Items[row, column];
        return SymmetricEigenSolver.Solve(values, tolerance, maxSweeps);
    }

    /// <summary>Computes approximate complex right eigenpairs of a real square matrix.</summary>
    /// <remarks>Values sort by real part then imaginary part; unit vectors are columns.
    /// Defective matrices need not have an independent eigenbasis. Nonconvergence
    /// and nonfinite results throw. The iteration limit must be in 1..100000.</remarks>
    public EigenDecomposition GetEigenDecomposition(int maxIterations = 1000)
    {
        var values = new double[RowCount, ColumnCount];
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < ColumnCount; column++) values[row, column] = _core.Items[row, column];
        return GeneralEigenSolver.Solve(values, maxIterations);
    }

    /// <summary>This legacy integer result cannot represent general eigenvalues.</summary>
    [Obsolete("Use GetSymmetricEigenDecomposition(), which returns real-valued eigenpairs.")]
    public int[] GetEigenValues() => throw new NotSupportedException("Integer eigenvalue arrays cannot represent general spectra. Use GetSymmetricEigenDecomposition().");

    /// <summary>This legacy integer result cannot represent normalized eigenvectors.</summary>
    [Obsolete("Use GetSymmetricEigenDecomposition(), which returns real-valued eigenpairs.")]
    public int[] GetEigenVectors() => throw new NotSupportedException("Integer arrays cannot represent normalized eigenvectors. Use GetSymmetricEigenDecomposition().");

    /// <summary>
    ///     The TryGetEigenValues method.
    /// </summary>
    /// <param name="eigenValues">
    ///     An out variable for the eigen values.
    /// </param>
    /// <returns>
    ///     Always false because this legacy signature cannot represent general real spectra.
    /// </returns>
    [Obsolete("Use GetSymmetricEigenDecomposition(), which returns real-valued eigenpairs.")]
    public bool TryGetEigenValues(out int[] eigenValues)
    {
        eigenValues = Array.Empty<int>();
        return false;
    }

    /// <summary>Tests for a nonzero exact determinant over the rational or real numbers.</summary>
    /// <remarks>The inverse need not contain integers. A nonsquare matrix returns false.</remarks>
    public bool IsInvertible() => IsSquare() && !DeterminantKernels.Bareiss(_core.Items).IsZero;

    /// <summary>Returns the exact determinant using Bareiss elimination; the empty determinant is one.</summary>
    /// <exception cref="OverflowException">The final determinant is outside the int range.</exception>
    public int GetDeterminant()
    {
        if (!IsSquare()) throw new NotSquareMatrixException();
        return (int)DeterminantKernels.Bareiss(_core.Items);
    }

    /// <summary>Returns false and zero for a nonsquare matrix or a determinant outside the int range.</summary>
    public bool TryGetDeterminant(out int determinant)
    {
        determinant = 0;
        if (!IsSquare()) return false;
        var exact = DeterminantKernels.Bareiss(_core.Items);
        if (exact < int.MinValue || exact > int.MaxValue) return false;
        determinant = (int)exact;
        return true;
    }

    /// <summary>Returns the product without modifying either operand.</summary>
    public static IntegerMatrix operator *(IntegerMatrix a, IntegerMatrix b)
    {
        if (a == null) throw new ArgumentNullException(nameof(a));
        if (b == null) throw new ArgumentNullException(nameof(b));
        return a.DotProduct(b);
    }

    /// <summary>Returns the sum without modifying either operand.</summary>
    public static IntegerMatrix operator +(IntegerMatrix a, IntegerMatrix b) => Combine(a, b, false);

    /// <summary>Returns the difference without modifying either operand.</summary>
    public static IntegerMatrix operator -(IntegerMatrix a, IntegerMatrix b) => Combine(a, b, true);

    private static IntegerMatrix Combine(IntegerMatrix a, IntegerMatrix b, bool subtract)
    {
        if (a == null) throw new ArgumentNullException(nameof(a));
        if (b == null) throw new ArgumentNullException(nameof(b));
        var result = a._core.Copy();
        if (!result.TryCombine(b._core, subtract)) throw new ArgumentOutOfRangeException();
        return new IntegerMatrix(result);
    }

    /// <summary>Enumerates entries in row-major order, retaining the original public return type.</summary>
    public IEnumerator GetEnumerator() => _core.Items.GetEnumerator();

    IEnumerator<int> IEnumerable<int>.GetEnumerator() => _core.GetTypedEnumerator();
}
