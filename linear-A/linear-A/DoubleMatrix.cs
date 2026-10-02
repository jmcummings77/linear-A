// Copyright John-Michael Cummings 2021.

using System;
using System.Collections;
using System.Collections.Generic;

namespace linear_A;

/// <summary>A mutable, fixed-size matrix of doubles.</summary>
/// <remarks>
/// Available on .NET Standard 2.1 and .NET 10. Prefer Matrix&lt;double&gt; for new .NET 10 code.
/// Basic arithmetic follows IEEE double rounding, infinity, and NaN rules.
/// Instances are not thread-safe.
/// </remarks>
public class DoubleMatrix : IEnumerable<double>
{
    private readonly LegacyMatrixCore<double, LegacyDoubleArithmetic> _core;

    /// <summary>Creates the empty 0-by-0 matrix.</summary>
    public DoubleMatrix() : this(new LegacyMatrixCore<double, LegacyDoubleArithmetic>(0, 0)) { }

    /// <summary>Creates a zero or identity square matrix with strictly positive dimensions.</summary>
    public DoubleMatrix(int squareDimensions, bool initializeAsIdentity)
        : this(LegacyMatrixCore<double, LegacyDoubleArithmetic>.Square(squareDimensions, initializeAsIdentity)) { }

    /// <summary>Creates a zero matrix with strictly positive row and column counts.</summary>
    public DoubleMatrix(int rowCount, int columnCount)
        : this(LegacyMatrixCore<double, LegacyDoubleArithmetic>.Rectangle(rowCount, columnCount)) { }

    /// <summary>Copies a matrix into independent storage; null creates an empty 0-by-0 matrix.</summary>
    public DoubleMatrix(DoubleMatrix? matrixToCopy)
        : this(matrixToCopy?._core.Copy() ?? new LegacyMatrixCore<double, LegacyDoubleArithmetic>(0, 0)) { }

    private DoubleMatrix(LegacyMatrixCore<double, LegacyDoubleArithmetic> core) => _core = core;

    /// <summary>Copies a zero-based array, preserving its shape, including zero-size dimensions.</summary>
    public static DoubleMatrix FromArray(double[,] values)
        => new(LegacyMatrixCore<double, LegacyDoubleArithmetic>.FromArray(values));

    /// <summary>Returns an independent two-dimensional copy of the entries.</summary>
    public double[,] ToArray() => (double[,])_core.Items.Clone();

#if NET10_0_OR_GREATER
    /// <summary>Copies the entries and shape into a generic matrix.</summary>
    public Matrix<double> ToMatrix() => new(_core.Items);

    /// <summary>Copies a generic matrix, preserving its shape, including zero-size dimensions.</summary>
    public static DoubleMatrix FromMatrix(Matrix<double> matrix)
    {
        if (matrix == null) throw new ArgumentNullException(nameof(matrix));
        var core = new LegacyMatrixCore<double, LegacyDoubleArithmetic>(matrix.RowCount, matrix.ColumnCount);
        for (var row = 0; row < matrix.RowCount; row++)
            for (var column = 0; column < matrix.ColumnCount; column++)
                core.Items[row, column] = matrix[row, column];
        return new DoubleMatrix(core);
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
    public double this[int rowIndex, int columnIndex]
    {
        get => _core.Items[rowIndex, columnIndex];
        set => _core.Items[rowIndex, columnIndex] = value;
    }

    /// <summary>Adds in place; null or unequal shapes return false without changing this matrix.</summary>
    public bool TryAddMatrix(DoubleMatrix? matrixToAdd) => _core.TryCombine(matrixToAdd?._core, false);

    /// <summary>Subtracts in place; null or unequal shapes return false without changing this matrix.</summary>
    public bool TrySubtractMatrix(DoubleMatrix? matrixToSubtract) => _core.TryCombine(matrixToSubtract?._core, true);

    /// <summary>Returns an independent transpose with swapped dimensions.</summary>
    public DoubleMatrix Transpose() => new(_core.Transpose());

    /// <summary>Multiplies each entry by the scalar in place.</summary>
    public void Scale(int scalar) => _core.Scale(scalar);

    /// <summary>Multiplies each entry by a double-precision scalar in place.</summary>
    public void Scale(double scalar) => _core.Scale(scalar);

    /// <summary>Returns the matrix product without modifying either operand.</summary>
    /// <remarks>The left column count must equal the right row count.</remarks>
    public DoubleMatrix DotProduct(DoubleMatrix multiplicand) => new(_core.DotProduct(multiplicand?._core));

    /// <summary>Replaces this vector with its right-handed cross product on success.</summary>
    /// <remarks>Null, invalid shapes, nonfinite entries/results, or checked overflow return false without mutation.</remarks>
    public bool CrossProduct(DoubleMatrix? multiplicand) => _core.CrossProduct(multiplicand?._core);

    /// <summary>Returns an independent right-handed cross product, retaining this vector's shape.</summary>
    /// <remarks>Each operand must be a finite 3-by-1 or 1-by-3 vector. Invalid input or checked overflow throws.</remarks>
    public DoubleMatrix GetCrossProduct(DoubleMatrix multiplicand) => new(_core.GetCrossProduct(multiplicand?._core));

    /// <summary>Returns an independent copy of the requested row.</summary>
    public double[] GetRow(int rowNumber) => _core.GetRow(rowNumber);

    /// <summary>Returns an independent copy of the requested column.</summary>
    public double[] GetColumn(int columnNumber) => _core.GetColumn(columnNumber);

    /// <summary>Returns whether row and column counts match.</summary>
    public bool IsSquare() => _core.IsSquare();

    /// <summary>Checks whether a square matrix is upper or lower triangular using exact zero comparisons.</summary>
    public bool IsTriangular() => IsUpperTriangular() || IsLowerTriangular();

    /// <summary>Checks whether every entry below the diagonal of a square matrix is exactly zero.</summary>
    public bool IsUpperTriangular() => _core.IsTriangular(true);

    /// <summary>Checks whether every entry above the diagonal of a square matrix is exactly zero.</summary>
    public bool IsLowerTriangular() => _core.IsTriangular(false);

    /// <summary>Returns the diagonal sum; nonsquare matrices throw NotSquareMatrixException.</summary>
    public double GetTrace() => _core.GetTrace();

    /// <summary>Returns false and zero for a nonsquare matrix; otherwise returns the diagonal sum.</summary>
    public bool TryGetTrace(out double trace) => _core.TryGetTrace(out trace);

    /// <summary>Computes ascending real eigenvalues and orthonormal eigenvector columns.</summary>
    /// <remarks>
    /// Requires a finite, exactly symmetric square matrix. Jacobi iteration uses
    /// a relative Frobenius-norm tolerance; small eigenvalues in mixed-scale matrices
    /// need not have small relative error. The matrix is unchanged. Nonconvergence
    /// or eigenvalues outside the finite double range throw.
    /// </remarks>
    public SymmetricEigenDecomposition GetSymmetricEigenDecomposition(double tolerance = 1e-12, int maxSweeps = 50)
        => SymmetricEigenSolver.Solve(_core.Items, tolerance, maxSweeps);

    /// <summary>Computes normalized complex right eigenvector columns of a finite real square matrix.</summary>
    /// <remarks>Uses Hessenberg reduction and double-shift QR. The input is unchanged.
    /// Values sort by real part then imaginary part. Defective matrices need not
    /// have an independent eigenbasis. Nonconvergence and nonfinite results throw.
    /// The iteration limit applies per unresolved root and must be in 1..100000.</remarks>
    public EigenDecomposition GetEigenDecomposition(int maxIterations = 1000)
        => GeneralEigenSolver.Solve(_core.Items, maxIterations);

    /// <summary>Gets ascending real eigenvalues of a finite, exactly symmetric square matrix.</summary>
    public double[] GetEigenValues() => GetSymmetricEigenDecomposition().EigenValues;

    /// <summary>Gets a row-major flattened n-by-n matrix of unit eigenvector columns.</summary>
    /// <remarks>Column j corresponds to the jth ascending eigenvalue. Requires exact symmetry.</remarks>
    public double[] GetEigenVectors()
    {
        var vectors = GetSymmetricEigenDecomposition().EigenVectors;
        var result = new double[checked(RowCount * ColumnCount)];
        for (var row = 0; row < RowCount; row++)
            for (var column = 0; column < ColumnCount; column++)
                result[row * ColumnCount + column] = vectors[row, column];
        return result;
    }

    /// <summary>
    ///     The TryGetEigenValues method.
    /// </summary>
    /// <param name="eigenValues">
    ///     An out variable for the eigen values.
    /// </param>
    /// <returns>
    ///     True if finite real eigenvalues were obtained for a symmetric square matrix.
    /// </returns>
    public bool TryGetEigenValues(out double[] eigenValues)
    {
        eigenValues = Array.Empty<double>();
        if (!IsSquare()) return false;
        try
        {
            eigenValues = GetEigenValues();
            return true;
        }
        catch (ArgumentException) { return false; }
        catch (ArithmeticException) { return false; }
        catch (InvalidOperationException) { return false; }
    }

    /// <summary>Tests whether partial-pivoting elimination produces finite nonzero pivots.</summary>
    /// <remarks>
    /// Uses no tolerance or conditioning estimate. Rounding can affect nearly singular inputs.
    /// Nonfinite inputs or elimination results return false; determinant underflow alone does not.
    /// </remarks>
    public bool IsInvertible()
    {
        if (!IsSquare()) return false;
        DeterminantKernels.Lu(_core.Items, out var isInvertible);
        return isInvertible;
    }

    /// <summary>Returns the determinant using scaled partial-pivoting LU; the empty determinant is one.</summary>
    /// <remarks>
    /// The input is unchanged. Rounding, overflow, and underflow can affect the result,
    /// especially for ill-conditioned matrices. Nonfinite inputs or elimination results produce NaN.
    /// </remarks>
    public double GetDeterminant()
    {
        if (!IsSquare()) throw new NotSquareMatrixException();
        return DeterminantKernels.Lu(_core.Items, out _);
    }

    /// <summary>Returns false and zero for a nonsquare matrix; a square result may be NaN or infinite.</summary>
    public bool TryGetDeterminant(out double determinant)
    {
        determinant = 0d;
        if (!IsSquare()) return false;
        determinant = DeterminantKernels.Lu(_core.Items, out _);
        return true;
    }

    /// <summary>Returns the product without modifying either operand.</summary>
    public static DoubleMatrix operator *(DoubleMatrix a, DoubleMatrix b)
    {
        if (a == null) throw new ArgumentNullException(nameof(a));
        if (b == null) throw new ArgumentNullException(nameof(b));
        return a.DotProduct(b);
    }

    /// <summary>Returns the sum without modifying either operand.</summary>
    public static DoubleMatrix operator +(DoubleMatrix a, DoubleMatrix b) => Combine(a, b, false);

    /// <summary>Returns the difference without modifying either operand.</summary>
    public static DoubleMatrix operator -(DoubleMatrix a, DoubleMatrix b) => Combine(a, b, true);

    private static DoubleMatrix Combine(DoubleMatrix a, DoubleMatrix b, bool subtract)
    {
        if (a == null) throw new ArgumentNullException(nameof(a));
        if (b == null) throw new ArgumentNullException(nameof(b));
        var result = a._core.Copy();
        if (!result.TryCombine(b._core, subtract)) throw new ArgumentOutOfRangeException();
        return new DoubleMatrix(result);
    }

    /// <summary>Enumerates entries in row-major order, retaining the original public return type.</summary>
    public IEnumerator GetEnumerator() => _core.Items.GetEnumerator();

    IEnumerator<double> IEnumerable<double>.GetEnumerator() => _core.GetTypedEnumerator();
}
