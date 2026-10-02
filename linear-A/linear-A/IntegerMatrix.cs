#region Copyright

// --------------------------------------------------------------------------------------------------------------------
// <copyright file="IntegerMatrix.cs" company="John-Michael Cummings">
//   John-Michael Cummings 2021
// </copyright>
// <summary>
//   A matrix of 32-bit integers.
// </summary>
// --------------------------------------------------------------------------------------------------------------------

#endregion

namespace linear_A
{

    #region Using Statements

    using System;
    using System.Collections;
    using System.Numerics;
    // ReSharper disable NotAccessedField.Local
    // ReSharper disable UnusedMember.Global
    // ReSharper disable UnusedParameter.Global

    #endregion

    /// <summary>
    ///     A matrix of integers with basic matrix operations.
    /// </summary>
    /// <remarks>
    ///     Addition, subtraction, scaling, multiplication and trace retain unchecked 32-bit
    ///     arithmetic. Determinants use exact intermediates and report an out-of-range result.
    /// </remarks>
    public class IntegerMatrix
    {
        /// <summary>
        ///     The integer byte size. Used for fast array copying.
        /// </summary>
        private const int IntegerByteSize = 4;

        /// <summary>
        ///     The items stored in the matrix.
        /// </summary>
        private readonly int[,] _items;

        /// <summary>
        ///     Gets the row count.
        /// </summary>
        public int RowCount { get; }

        /// <summary>
        ///     Gets the column count.
        /// </summary>
        public int ColumnCount { get; }

        /// <summary>
        ///     The index accessor method.
        /// <param name="rowIndex">
        ///     The row index.
        /// </param>
        /// <param name="columnIndex">
        ///     The column index.
        /// </param>
        /// <returns>
        ///     The <see cref="int" />.
        /// </returns>
        /// </summary>
        public int this[int rowIndex, int columnIndex]
        {
            get => _items[rowIndex, columnIndex];
            set => _items[rowIndex, columnIndex] = value;
        }

        /// <summary>
        ///     The is fixed size.
        /// </summary>
        public bool IsFixedSize => true;

        /// <summary>
        ///     The is read only.
        /// </summary>
        public bool IsReadOnly => false;

        /// <summary>
        ///     Initializes a new instance of the <see cref="IntegerMatrix" /> class.
        /// </summary>
        /// <param name="squareDimensions">
        ///     The dimension for both rows and columns for a square matrix.
        /// </param>
        /// <param name="initializeAsIdentity">
        ///     The initialize as identity flag
        ///     Creates a matrix with default zero values if false
        ///     Or an identity matrix with dimensions specified by the squareDimensions param
        /// </param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// </exception>
        public IntegerMatrix(int squareDimensions, bool initializeAsIdentity)
        {
            if (squareDimensions < 1)
                throw new ArgumentOutOfRangeException(nameof(squareDimensions), "Dimensions must be greater than zero");

            RowCount = squareDimensions;
            ColumnCount = squareDimensions;
            _items = new int[RowCount, ColumnCount];

            if (!initializeAsIdentity) return;

            for (var i = 0; i < RowCount; i++)
                _items[i, i] = 1;
        }

        /// <summary>
        ///     Initializes a new instance of the <see cref="IntegerMatrix" /> class.
        /// </summary>
        /// <param name="rowCount">
        ///     The row count.
        /// </param>
        /// <param name="columnCount">
        ///     The column count.
        /// </param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// </exception>
        public IntegerMatrix(int rowCount, int columnCount)
        {
            if (rowCount < 1 || columnCount < 1)
            {
                if (rowCount < 1 && columnCount < 1)
                    throw new ArgumentOutOfRangeException(
                        $"{nameof(rowCount)} & {nameof(columnCount)}",
                        "Dimensions must be greater than zero");

                if (rowCount < 1)
                    throw new ArgumentOutOfRangeException(nameof(rowCount), "Dimensions must be greater than zero");

                throw new ArgumentOutOfRangeException(nameof(columnCount), "Dimensions must be greater than zero");
            }

            RowCount = rowCount;
            ColumnCount = columnCount;
            _items = new int[RowCount, ColumnCount];
        }

        /// <summary>
        ///     Initializes a new instance of the <see cref="IntegerMatrix" /> class as a duplicate of the provided matrix.
        /// </summary>
        /// <param name="matrixToCopy">Matrix to copy</param>
        public IntegerMatrix(IntegerMatrix? matrixToCopy)
        {
            if (matrixToCopy == null)
            {
                _items = new int[0, 0];
                return;
            }

            RowCount = matrixToCopy.RowCount;
            ColumnCount = matrixToCopy.ColumnCount;
            _items = (int[,])matrixToCopy._items.Clone();
        }

        /// <summary>
        ///     Initializes a default instance of the <see cref="IntegerMatrix" /> class.
        /// </summary>
        public IntegerMatrix()
        {
            _items = new int[0, 0];
        }

        /// <summary>
        ///     The TrySubtract method.
        ///     Returns false without changing this matrix if the argument is null or the dimensions do not match.
        ///     Otherwise subtracts the provided matrix from this matrix.
        /// </summary>
        /// <param name="matrixToSubtract">
        ///     The matrix to subtract.
        /// </param>
        /// <returns>
        ///     The <see cref="bool" />.
        /// </returns>
        public bool TrySubtractMatrix(IntegerMatrix? matrixToSubtract)
        {
            if (matrixToSubtract == null || matrixToSubtract.RowCount != RowCount || matrixToSubtract.ColumnCount != ColumnCount) return false;

            SubtractMatrix(matrixToSubtract);

            return true;
        }

        /// <summary>
        ///     The unsafe subtract method.
        ///     Only works if the matrices are of the correct dimensions.
        /// </summary>
        /// <param name="matrixToSubtract"></param>
        private void SubtractMatrix(IntegerMatrix matrixToSubtract)
        {
            for (var i = 0; i < RowCount; i++)
                for (var j = 0; j < ColumnCount; j++)
                    _items[i, j] -= matrixToSubtract[i, j];
        }

        /// <summary>
        ///     The TryAdd method.
        ///     Returns false and leaves this matrix unchanged if the argument is null or has different dimensions.
        ///     Otherwise returns true and adds the values from the provided matrix to this matrix.
        /// </summary>
        /// <param name="matrixToAdd">
        ///     The matrix to add to this matrix.
        /// </param>
        /// <returns>
        ///     The success status, which returns false when matrices are not the same size <see cref="bool" />.
        /// </returns>
        public bool TryAddMatrix(IntegerMatrix? matrixToAdd)
        {
            if (matrixToAdd == null || matrixToAdd.RowCount != RowCount || matrixToAdd.ColumnCount != ColumnCount) return false;

            AddMatrix(matrixToAdd);

            return true;
        }

        /// <summary>
        ///     The unsafe add method. Only works correctly if the matrices are of the same dimensions.
        /// </summary>
        /// <param name="matrixToAdd"></param>
        private void AddMatrix(IntegerMatrix matrixToAdd)
        {
            for (var i = 0; i < RowCount; i++)
                for (var j = 0; j < ColumnCount; j++)
                    _items[i, j] += matrixToAdd[i, j];
        }

        /// <summary>
        ///     The transpose method.
        /// </summary>
        /// <returns>
        ///     A newly instantiated matrix consisting of the transpose of this matrix <see cref="IntegerMatrix" />.
        /// </returns>
        public IntegerMatrix Transpose()
        {
            if (RowCount == 0) return new IntegerMatrix();

            var result = new IntegerMatrix(ColumnCount, RowCount);
            for (var i = 0; i < RowCount; i++)
                for (var j = 0; j < ColumnCount; j++)
                    result[j, i] = _items[i, j];

            return result;
        }

        /// <summary>
        ///     The scale method. Multiplies each element of the matrix by the provided scalar.
        /// </summary>
        /// <param name="scalar">
        ///     The scalar.
        /// </param>
        public void Scale(int scalar)
        {
            for (var i = 0; i < RowCount; i++)
                for (var j = 0; j < ColumnCount; j++)
                    _items[i, j] *= scalar;
        }

        /// <summary>
        ///     The dot product method.
        /// </summary>
        /// <param name="multiplicand">
        ///     The multiplicand.
        /// </param>
        /// <returns>
        ///     A newly instantiated matrix <see cref="IntegerMatrix" />. consisting of the dot product
        /// of this matrix and the provided multiplicand
        /// </returns>
        /// <exception cref="ArgumentOutOfRangeException">
        /// </exception>
        public IntegerMatrix DotProduct(IntegerMatrix multiplicand)
        {
            if (multiplicand == null) throw new ArgumentNullException(nameof(multiplicand));
            if (ColumnCount != multiplicand.RowCount)
                throw new ArgumentOutOfRangeException(nameof(multiplicand), "The left column count must equal the right row count.");

            if (RowCount == 0) return new IntegerMatrix();

            var result = new IntegerMatrix(RowCount, multiplicand.ColumnCount);

            for (var i = 0; i < RowCount; i++)
                for (var j = 0; j < multiplicand.ColumnCount; j++)
                    for (var k = 0; k < ColumnCount; k++)
                        result._items[i, j] += _items[i, k] * multiplicand._items[k, j];

            return result;
        }

        /// <summary>
        ///     The cross product.
        /// </summary>
        /// <param name="multiplicand">
        ///     The multiplicand.
        /// </param>
        /// <returns>
        ///     The <see cref="bool" />.
        /// </returns>
        public bool CrossProduct(IntegerMatrix multiplicand) => false;

        /// <summary>
        ///     The get enumerator method for implementing the IEnumerable interface.
        /// </summary>
        /// <returns>
        ///     The the enumerator for the underlying matrix <see cref="IEnumerator" />.
        /// </returns>
        public IEnumerator GetEnumerator() => _items.GetEnumerator();

        /// <summary>
        ///     The get row method.
        /// </summary>
        /// <param name="rowNumber">
        ///     The row number.
        /// </param>
        /// <returns>
        ///     Returns an array of <see cref="int" /> corresponding to the requested row in the matrix.
        /// </returns>
        public int[] GetRow(int rowNumber)
        {
            if (rowNumber < 0 || rowNumber >= RowCount) throw new ArgumentOutOfRangeException(nameof(rowNumber));

            var result = new int[ColumnCount];
            Buffer.BlockCopy(
                _items,
                IntegerByteSize * ColumnCount * rowNumber,
                result,
                0,
                IntegerByteSize * ColumnCount);
            return result;
        }

        /// <summary>
        ///     The get column method.
        /// </summary>
        /// <param name="columnNumber">
        ///     The column number.
        /// </param>
        /// <returns>
        ///     Returns an array of <see cref="int" /> corresponding to the requested column in the matrix.
        /// </returns>
        public int[] GetColumn(int columnNumber)
        {
            if (columnNumber < 0 || columnNumber >= ColumnCount) throw new ArgumentOutOfRangeException(nameof(columnNumber));

            var result = new int[RowCount];
            for (var i = 0; i < RowCount; i++) result[i] = _items[i, columnNumber];

            return result;
        }

        /// <summary>
        ///     The IsSquare property.
        /// </summary>
        /// <returns>
        ///     A <see cref="bool" /> indicating whether the matrix is square.
        /// </returns>
        public bool IsSquare() => RowCount == ColumnCount;

        /// <summary>
        ///     The IsTriangular property.
        /// </summary>
        /// <returns>
        ///     A <see cref="bool" />  indicating whether the matrix is triangular.
        /// </returns>
        public bool IsTriangular() => IsSquare() && (IsUpperTriangularUnsafe() || IsLowerTriangularUnsafe());

        /// <summary>
        ///     The IsUpperTriangular property.
        /// </summary>
        /// <returns>
        ///     A <see cref="bool" />  indicating whether the matrix is upper triangular.
        /// </returns>
        public bool IsUpperTriangular() => IsSquare() && IsUpperTriangularUnsafe();

        /// <summary>
        ///     The IsLowerTriangular property.
        /// </summary>
        /// <returns>
        ///     A <see cref="bool" />  indicating whether the matrix is lower triangular.
        /// </returns>
        public bool IsLowerTriangular() => IsSquare() && IsLowerTriangularUnsafe();

        /// <summary>
        ///     Determines whether this matrix has an inverse over the rational or real numbers.
        ///     The inverse need not have integer entries.
        /// </summary>
        /// <returns>True for a square matrix with a nonzero exact determinant.</returns>
        public bool IsInvertible()
        {
            if (!IsSquare()) return false;

            return !GetDeterminantUnsafe(_items).IsZero;
        }

        /// <summary>
        ///     Unsafe method for checking if the matrix is lower triangular.
        /// </summary>
        /// <returns>
        ///     A <see cref="bool" />  indicating whether the matrix is lower triangular.
        /// </returns>
        private bool IsLowerTriangularUnsafe()
        {
            for (var i = 0; i < RowCount; i++)
                for (var j = i + 1; j < ColumnCount; j++)
                    if (_items[i, j] != 0)
                        return false;

            return true;
        }

        /// <summary>
        ///     Unsafe method for checking if the matrix is upper triangular.
        /// </summary>
        /// <returns>
        ///     A <see cref="bool" />  indicating whether the matrix is upper triangular.
        /// </returns>
        private bool IsUpperTriangularUnsafe()
        {
            for (var i = 0; i < RowCount; i++)
                for (var j = i + 1; j < ColumnCount; j++)
                    if (_items[j, i] != 0)
                        return false;

            return true;
        }

        /// <summary>
        ///     The GetTrace method.
        /// </summary>
        /// <returns>
        ///     An <see cref="int" /> corresponding to the matrix's trace.
        /// </returns>
        public int GetTrace()
        {
            if (!IsSquare()) throw new NotSquareMatrixException();

            return GetTraceUnsafe();
        }

        /// <summary>
        ///     Unsafe GetTrace method.
        /// </summary>
        /// <returns>
        ///     An <see cref="int" /> corresponding to the matrix's trace.
        /// </returns>
        private int GetTraceUnsafe()
        {
            var result = 0;
            for (var i = 0; i < RowCount; i++)
                result += _items[i, i];
            return result;
        }

        /// <summary>
        ///     The TryGet trace method.
        /// </summary>
        /// <param name="trace">
        ///     The out variable for returning trace.
        /// </param>
        /// <returns>
        ///     Flag <see cref="bool" /> indicating whether the trace can be calculated.
        /// </returns>
        public bool TryGetTrace(out int trace)
        {
            trace = 0;

            if (!IsSquare()) return false;

            trace = GetTraceUnsafe();

            return true;
        }

        /// <summary>
        ///     The get eigen values method.
        /// </summary>
        /// <returns>
        ///     An array of <see cref="int" /> corresponding to the matrix eigenvalues.
        /// </returns>
        public int[] GetEigenValues()
        {
            if (!IsSquare()) throw new NotSquareMatrixException();

            return GetEigenValuesUnsafe();
        }

        /// <summary>
        ///     The unsafe get eigen values method.
        /// </summary>
        /// <returns>
        ///     An array of <see cref="int" /> corresponding to the matrix eigenvalues.
        /// </returns>
        private static int[] GetEigenValuesUnsafe() => throw new NotImplementedException();

        /// <summary>
        ///     The get eigen vectors method.
        /// </summary>
        /// <returns>
        ///     An array of <see cref="int" /> corresponding to the matrix eigen vectors.
        /// </returns>
        public int[] GetEigenVectors() => throw new NotImplementedException();

        /// <summary>
        ///     The TryGetEigenValues method.
        /// </summary>
        /// <param name="eigenValues">
        ///     An out variable for the eigen values.
        /// </param>
        /// <returns>
        ///     Flag indicating whether the values were retrieved <see cref="bool" />.
        /// </returns>
        public bool TryGetEigenValues(out int[] eigenValues)
        {
            eigenValues = new int[0];

            if (!IsSquare()) return false;

            eigenValues = GetEigenValuesUnsafe();

            return true;
        }

        /// <summary>
        ///     The TryGetDeterminant method.
        /// </summary>
        /// <param name="determinant">
        ///     An out variable for the determinant.
        /// </param>
        /// <returns>
        ///     True if this matrix is square and its exact determinant fits in an <see cref="int" />.
        ///     Otherwise returns false and sets <paramref name="determinant" /> to zero.
        /// </returns>
        public bool TryGetDeterminant(out int determinant)
        {
            determinant = 0;

            if (!IsSquare()) return false;

            var exactDeterminant = GetDeterminantUnsafe(_items);
            if (exactDeterminant < int.MinValue || exactDeterminant > int.MaxValue) return false;

            determinant = (int)exactDeterminant;
            return true;
        }

        /// <summary>
        ///     The GetDeterminant method.
        /// </summary>
        /// <returns>
        ///     The exact determinant, with the determinant of the empty matrix defined as one.
        /// </returns>
        /// <exception cref="OverflowException">The determinant is outside the range of <see cref="int" />.</exception>
        public int GetDeterminant()
        {
            if (!IsSquare()) throw new NotSquareMatrixException();

            return (int)GetDeterminantUnsafe(_items);
        }

        /// <summary>
        ///     Computes an exact determinant using fraction-free Bareiss elimination.
        ///     Arbitrary-precision intermediates prevent overflow even when large terms cancel.
        /// </summary>
        private static BigInteger GetDeterminantUnsafe(int[,] matrix)
        {
            var size = matrix.GetLength(0);
            if (size == 0) return BigInteger.One;
            if (size == 1) return matrix[0, 0];

            var work = new BigInteger[size, size];
            for (var row = 0; row < size; row++)
                for (var column = 0; column < size; column++)
                    work[row, column] = matrix[row, column];

            var previousPivot = BigInteger.One;
            var sign = 1;
            for (var pivotColumn = 0; pivotColumn < size - 1; pivotColumn++)
            {
                var pivotRow = pivotColumn;
                while (pivotRow < size && work[pivotRow, pivotColumn].IsZero) pivotRow++;
                if (pivotRow == size) return BigInteger.Zero;

                if (pivotRow != pivotColumn)
                {
                    for (var column = pivotColumn; column < size; column++)
                    {
                        var temporary = work[pivotColumn, column];
                        work[pivotColumn, column] = work[pivotRow, column];
                        work[pivotRow, column] = temporary;
                    }
                    sign = -sign;
                }

                var pivot = work[pivotColumn, pivotColumn];
                for (var row = pivotColumn + 1; row < size; row++)
                {
                    for (var column = pivotColumn + 1; column < size; column++)
                        work[row, column] =
                            (pivot * work[row, column] - work[row, pivotColumn] * work[pivotColumn, column]) /
                            previousPivot;
                    work[row, pivotColumn] = BigInteger.Zero;
                }
                previousPivot = pivot;
            }

            return sign * work[size - 1, size - 1];
        }

        /// <summary>
        ///     Overloads * operator for matrix class
        /// </summary>
        /// <param name="a">The left side matrix</param>
        /// <param name="b">The right matrix</param>
        /// <returns>The dot product</returns>
        public static IntegerMatrix operator *(IntegerMatrix a, IntegerMatrix b)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));

            return a.DotProduct(b);
        }

        /// <summary>
        ///     Overloads + operator for addition
        /// </summary>
        /// <param name="a">The left side matrix</param>
        /// <param name="b">The right side matrix</param>
        /// <returns></returns>
        public static IntegerMatrix operator +(IntegerMatrix a, IntegerMatrix b)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));

            var result = new IntegerMatrix(a);
            if (result.TryAddMatrix(b)) return result;

            throw new ArgumentOutOfRangeException();
        }

        /// <summary>
        ///     Overloads - operator for subtraction
        /// </summary>
        /// <param name="a">The left side matrix</param>
        /// <param name="b">The right side matrix</param>
        /// <returns></returns>
        public static IntegerMatrix operator -(IntegerMatrix a, IntegerMatrix b)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));

            var result = new IntegerMatrix(a);
            if (result.TrySubtractMatrix(b)) return result;

            throw new ArgumentOutOfRangeException();
        }
    }
}
