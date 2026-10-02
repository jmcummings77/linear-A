#region Copyright

// --------------------------------------------------------------------------------------------------------------------
// <copyright file="DoubleMatrix.cs" company="John-Michael Cummings">
//   John-Michael Cummings 2021
// </copyright>
// <summary>
//   The double matrix.
// </summary>
// --------------------------------------------------------------------------------------------------------------------

#endregion

// ReSharper disable NotAccessedField.Local
// ReSharper disable UnusedMember.Global
// ReSharper disable UnusedParameter.Global
namespace linear_A
{

    #region Using Statements

    using System;
    using System.Collections;

    #endregion

    /// <summary>
    ///     A matrix of doubles with basic matrix operations.
    /// </summary>
    public class DoubleMatrix
    {
        /// <summary>
        ///     The double byte size. Used for fast array copying.
        /// </summary>
        private const int DoubleByteSize = 8;

        /// <summary>
        ///     The items stored in the matrix.
        /// </summary>
        private readonly double[,] _items;

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
        ///     The <see cref="double" />.
        /// </returns>
        /// </summary>
        public double this[int rowIndex, int columnIndex]
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
        ///     Initializes a new instance of the <see cref="DoubleMatrix" /> class.
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
        public DoubleMatrix(int squareDimensions, bool initializeAsIdentity)
        {
            if (squareDimensions < 1)
                throw new ArgumentOutOfRangeException(nameof(squareDimensions), "Dimensions must be greater than zero");

            RowCount = squareDimensions;
            ColumnCount = squareDimensions;
            _items = new double[RowCount, ColumnCount];

            if (!initializeAsIdentity) return;

            for (var i = 0; i < RowCount; i++)
                _items[i, i] = 1;
        }

        /// <summary>
        ///     Initializes a new instance of the <see cref="DoubleMatrix" /> class.
        /// </summary>
        /// <param name="rowCount">
        ///     The row count.
        /// </param>
        /// <param name="columnCount">
        ///     The column count.
        /// </param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// </exception>
        public DoubleMatrix(int rowCount, int columnCount)
        {
            if (rowCount < 1 || columnCount < 1)
            {
                if (rowCount < 1 && columnCount < 1)
                    throw new ArgumentOutOfRangeException(
                        nameof(rowCount) + " & " + nameof(columnCount),
                        "Dimensions must be greater than zero");

                if (rowCount < 1)
                    throw new ArgumentOutOfRangeException(nameof(rowCount), "Dimensions must be greater than zero");

                throw new ArgumentOutOfRangeException(nameof(columnCount), "Dimensions must be greater than zero");
            }

            RowCount = rowCount;
            ColumnCount = columnCount;
            _items = new double[RowCount, ColumnCount];
        }

        /// <summary>
        ///     Initializes a new instance of the <see cref="DoubleMatrix" /> class as a duplicate of the provided matrix.
        /// </summary>
        /// <param name="matrixToCopy">Matrix to copy</param>
        public DoubleMatrix(DoubleMatrix? matrixToCopy)
        {
            if (matrixToCopy == null)
            {
                _items = new double[0, 0];
                return;
            }

            RowCount = matrixToCopy.RowCount;
            ColumnCount = matrixToCopy.ColumnCount;
            _items = (double[,])matrixToCopy._items.Clone();
        }

        /// <summary>
        ///     Initializes a default instance of the <see cref="DoubleMatrix" /> class.
        /// </summary>
        public DoubleMatrix()
        {
            _items = new double[0, 0];
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
        public bool TrySubtractMatrix(DoubleMatrix? matrixToSubtract)
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
        private void SubtractMatrix(DoubleMatrix matrixToSubtract)
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
        public bool TryAddMatrix(DoubleMatrix? matrixToAdd)
        {
            if (matrixToAdd == null || matrixToAdd.RowCount != RowCount || matrixToAdd.ColumnCount != ColumnCount) return false;

            AddMatrix(matrixToAdd);

            return true;
        }

        /// <summary>
        ///     The unsafe add method. Only works correctly if the matrices are of the same dimensions.
        /// </summary>
        /// <param name="matrixToAdd"></param>
        private void AddMatrix(DoubleMatrix matrixToAdd)
        {
            for (var i = 0; i < RowCount; i++)
                for (var j = 0; j < ColumnCount; j++)
                    _items[i, j] += matrixToAdd[i, j];
        }

        /// <summary>
        ///     The transpose method.
        /// </summary>
        /// <returns>
        ///     A newly instantiated matrix consisting of the transpose of this matrix <see cref="DoubleMatrix" />.
        /// </returns>
        public DoubleMatrix Transpose()
        {
            if (RowCount == 0) return new DoubleMatrix();

            var result = new DoubleMatrix(ColumnCount, RowCount);
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
        public void Scale(int scalar) => Scale((double)scalar);

        /// <summary>
        ///     Multiplies each element by a double-precision scalar.
        /// </summary>
        /// <param name="scalar">The scalar.</param>
        public void Scale(double scalar)
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
        ///     A newly instantiated matrix <see cref="DoubleMatrix" />. consisting of the dot product
        /// of this matrix and the provided multiplicand
        /// </returns>
        /// <exception cref="ArgumentOutOfRangeException">
        /// </exception>
        public DoubleMatrix DotProduct(DoubleMatrix multiplicand)
        {
            if (multiplicand == null) throw new ArgumentNullException(nameof(multiplicand));
            if (ColumnCount != multiplicand.RowCount)
                throw new ArgumentOutOfRangeException(nameof(multiplicand), "The left column count must equal the right row count.");

            if (RowCount == 0) return new DoubleMatrix();

            var result = new DoubleMatrix(RowCount, multiplicand.ColumnCount);

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
        public bool CrossProduct(DoubleMatrix multiplicand) => false;

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
        ///     Returns an array of <see cref="double" /> corresponding to the requested row in the matrix.
        /// </returns>
        public double[] GetRow(int rowNumber)
        {
            if (rowNumber < 0 || rowNumber >= RowCount) throw new ArgumentOutOfRangeException(nameof(rowNumber));

            var result = new double[ColumnCount];
            Buffer.BlockCopy(
                _items,
                DoubleByteSize * ColumnCount * rowNumber,
                result,
                0,
                DoubleByteSize * ColumnCount);
            return result;
        }

        /// <summary>
        ///     The get column method.
        /// </summary>
        /// <param name="columnNumber">
        ///     The column number.
        /// </param>
        /// <returns>
        ///     Returns an array of <see cref="double" /> corresponding to the requested column in the matrix.
        /// </returns>
        public double[] GetColumn(int columnNumber)
        {
            if (columnNumber < 0 || columnNumber >= ColumnCount) throw new ArgumentOutOfRangeException(nameof(columnNumber));

            var result = new double[RowCount];
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
        ///     Tests whether partial-pivoting elimination produces nonzero finite pivots.
        /// </summary>
        /// <remarks>
        ///     This floating-point test uses no tolerance and does not estimate conditioning.
        ///     Rounding can affect the result for nearly singular matrices. Nonfinite inputs or
        ///     elimination results return false. Determinant underflow does not imply singularity.
        /// </remarks>
        /// <returns>True when every computed pivot is finite and nonzero.</returns>
        public bool IsInvertible()
        {
            if (!IsSquare()) return false;

            GetDeterminantUnsafe(_items, out var isInvertible);
            return isInvertible;
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
                    if (_items[i, j] != 0d)
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
                    if (_items[j, i] != 0d)
                        return false;

            return true;
        }

        /// <summary>
        ///     The GetTrace method.
        /// </summary>
        /// <returns>
        ///     An <see cref="double" /> corresponding to the matrix's trace.
        /// </returns>
        public double GetTrace()
        {
            if (!IsSquare()) throw new NotSquareMatrixException();

            return GetTraceUnsafe();
        }

        /// <summary>
        ///     Unsafe GetTrace method.
        /// </summary>
        /// <returns>
        ///     An <see cref="double" /> corresponding to the matrix's trace.
        /// </returns>
        private double GetTraceUnsafe()
        {
            var result = 0d;
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
        public bool TryGetTrace(out double trace)
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
        ///     An array of <see cref="double" /> corresponding to the matrix eigenvalues.
        /// </returns>
        public double[] GetEigenValues()
        {
            if (!IsSquare()) throw new NotSquareMatrixException();

            return GetEigenValuesUnsafe();
        }

        /// <summary>
        ///     The unsafe get eigen values method.
        /// </summary>
        /// <returns>
        ///     An array of <see cref="double" /> corresponding to the matrix eigenvalues.
        /// </returns>
        private static double[] GetEigenValuesUnsafe() => throw new NotImplementedException();

        /// <summary>
        ///     The get eigen vectors method.
        /// </summary>
        /// <returns>
        ///     An array of <see cref="double" /> corresponding to the matrix eigen vectors.
        /// </returns>
        public double[] GetEigenVectors() => throw new NotImplementedException();

        /// <summary>
        ///     The TryGetEigenValues method.
        /// </summary>
        /// <param name="eigenValues">
        ///     An out variable for the eigen values.
        /// </param>
        /// <returns>
        ///     Flag indicating whether the values were retrieved <see cref="bool" />.
        /// </returns>
        public bool TryGetEigenValues(out double[] eigenValues)
        {
            eigenValues = new double[0];

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
        ///     True if this matrix is square. The result may be NaN or infinite; the same
        ///     floating-point limitations as <see cref="GetDeterminant" /> apply.
        /// </returns>
        public bool TryGetDeterminant(out double determinant)
        {
            determinant = 0;

            if (!IsSquare()) return false;

            determinant = GetDeterminantUnsafe(_items, out _);

            return true;
        }

        /// <summary>
        ///     The GetDeterminant method.
        /// </summary>
        /// <returns>
        ///     The determinant, with the determinant of the empty matrix defined as one.
        /// </returns>
        /// <remarks>
        ///     Uses partial-pivoting elimination with power-of-two row scaling where it preserves
        ///     the entries. Results remain subject to floating-point rounding, overflow and underflow,
        ///     especially for ill-conditioned matrices. Nonfinite input or elimination results produce NaN.
        /// </remarks>
        public double GetDeterminant()
        {
            if (!IsSquare()) throw new NotSquareMatrixException();

            return GetDeterminantUnsafe(_items, out _);
        }

        /// <summary>
        ///     Computes a determinant by Gaussian elimination with partial pivoting.
        ///     The input is copied so inspecting a matrix never changes its entries.
        /// </summary>
        private static double GetDeterminantUnsafe(double[,] matrix, out bool isInvertible)
        {
            var size = matrix.GetLength(0);
            var work = (double[,])matrix.Clone();
            isInvertible = false;

            for (var row = 0; row < size; row++)
                for (var column = 0; column < size; column++)
                    if (double.IsNaN(work[row, column]) || double.IsInfinity(work[row, column]))
                        return double.NaN;

            long rowScaleExponent = 0;
            for (var row = 0; row < size; row++)
            {
                var largest = 0d;
                for (var column = 0; column < size; column++)
                    largest = Math.Max(largest, Math.Abs(work[row, column]));
                if (largest == 0d) return 0d;

                var encodedExponent = (int)((BitConverter.DoubleToInt64Bits(largest) >> 52) & 0x7ff);
                if (encodedExponent == 0)
                {
                    // Normalize subnormal rows before constructing their reciprocal scale.
                    for (var column = 0; column < size; column++) work[row, column] *= 18014398509481984d;
                    largest *= 18014398509481984d;
                    encodedExponent = (int)((BitConverter.DoubleToInt64Bits(largest) >> 52) & 0x7ff);
                    rowScaleExponent -= 54;
                }

                var exponent = encodedExponent - 1023;
                var scale = Math.Pow(2d, -exponent);
                var preservesEntries = true;
                for (var column = 0; column < size && preservesEntries; column++)
                    preservesEntries = work[row, column] * scale / scale == work[row, column];
                if (!preservesEntries) continue;

                for (var column = 0; column < size; column++) work[row, column] *= scale;
                rowScaleExponent += exponent;
            }

            var sign = 1;
            for (var pivotColumn = 0; pivotColumn < size; pivotColumn++)
            {
                var pivotRow = pivotColumn;
                for (var row = pivotColumn + 1; row < size; row++)
                    if (Math.Abs(work[row, pivotColumn]) > Math.Abs(work[pivotRow, pivotColumn]))
                        pivotRow = row;

                var pivot = work[pivotRow, pivotColumn];
                if (double.IsNaN(pivot) || double.IsInfinity(pivot)) return double.NaN;
                if (pivot == 0d) return 0d;

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

                for (var row = pivotColumn + 1; row < size; row++)
                {
                    var factor = work[row, pivotColumn] / pivot;
                    work[row, pivotColumn] = 0d;
                    for (var column = pivotColumn + 1; column < size; column++)
                    {
                        work[row, column] -= factor * work[pivotColumn, column];
                        if (double.IsNaN(work[row, column]) || double.IsInfinity(work[row, column]))
                            return double.NaN;
                    }
                }
            }

            isInvertible = true;
            return MultiplyDiagonal(work, sign, rowScaleExponent);
        }

        /// <summary>
        ///     Multiplies the pivots with a separate binary exponent so mixed large and small
        ///     pivots do not cause avoidable intermediate overflow or underflow.
        /// </summary>
        private static double MultiplyDiagonal(double[,] matrix, int sign, long exponent)
        {
            var mantissa = (double)sign;
            for (var i = 0; i < matrix.GetLength(0); i++)
            {
                var value = matrix[i, i];
                var bits = BitConverter.DoubleToInt64Bits(value);
                var valueExponent = (int)((bits >> 52) & 0x7ff);
                if (valueExponent == 0)
                {
                    // Scaling by 2^54 normalizes any nonzero subnormal double exactly.
                    value *= 18014398509481984d;
                    bits = BitConverter.DoubleToInt64Bits(value);
                    valueExponent = (int)((bits >> 52) & 0x7ff);
                    exponent -= 54;
                }

                exponent += valueExponent - 1023;
                bits = (bits & unchecked((long)0x800fffffffffffffUL)) | (1023L << 52);
                mantissa *= BitConverter.Int64BitsToDouble(bits);
                if (Math.Abs(mantissa) >= 2d)
                {
                    mantissa *= 0.5d;
                    exponent++;
                }
            }

            if (exponent > 1023) return mantissa > 0d ? double.PositiveInfinity : double.NegativeInfinity;
            if (exponent < -1075) return mantissa * 0d;
            if (exponent < -1022)
                return mantissa * Math.Pow(2d, exponent + 1022) * Math.Pow(2d, -1022);

            return mantissa * BitConverter.Int64BitsToDouble((exponent + 1023) << 52);
        }

        /// <summary>
        ///     Overloads * operator for matrix class
        /// </summary>
        /// <param name="a">The left side matrix</param>
        /// <param name="b">The right matrix</param>
        /// <returns>The dot product</returns>
        public static DoubleMatrix operator *(DoubleMatrix a, DoubleMatrix b)
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
        public static DoubleMatrix operator +(DoubleMatrix a, DoubleMatrix b)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));

            var result = new DoubleMatrix(a);
            if (result.TryAddMatrix(b)) return result;

            throw new ArgumentOutOfRangeException();
        }

        /// <summary>
        ///     Overloads - operator for subtraction
        /// </summary>
        /// <param name="a">The left side matrix</param>
        /// <param name="b">The right side matrix</param>
        /// <returns></returns>
        public static DoubleMatrix operator -(DoubleMatrix a, DoubleMatrix b)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));

            var result = new DoubleMatrix(a);
            if (result.TrySubtractMatrix(b)) return result;

            throw new ArgumentOutOfRangeException();
        }
    }
}
