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
    using System.Collections.Generic;

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
        ///     Gets the count of items in the matrix.
        /// </summary>
        private int _count;

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
        public double this[int rowIndex, int columnIndex]
        {
            get => _items[rowIndex, columnIndex];
            set
            {
                _valuesHaveChanged = true;
                _items[rowIndex, columnIndex] = value;
            }
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
        ///     Flag indicating whether the matrix values have changed, requiring re-computation of cached characteristics values
        /// </summary>
#pragma warning disable 414
        private bool _valuesHaveChanged;
#pragma warning restore 414

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
            _count = RowCount * ColumnCount;
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
            _count = RowCount * ColumnCount;
            _items = new double[RowCount, ColumnCount];
        }

        /// <summary>
        ///     Initializes a new instance of the <see cref="DoubleMatrix" /> class as a duplicate of the provided matrix.
        /// </summary>
        /// <param name="matrixToCopy">Matrix to copy</param>
        public DoubleMatrix(DoubleMatrix matrixToCopy)
        {
            if (matrixToCopy == null)
            {
                _items = new double[0,0];
                return;
            }

            RowCount = matrixToCopy.RowCount;
            ColumnCount = matrixToCopy.ColumnCount;
            _count = RowCount * ColumnCount;
            _items = (double[,]) matrixToCopy._items.Clone();
        }

        /// <summary>
        ///     Initializes a default instance of the <see cref="DoubleMatrix" /> class.
        /// </summary>
        public DoubleMatrix()
        {
            _items = new double[0,0];
        }

        /// <summary>
        ///     The TrySubtract method.
        ///     Returns false and does not execute the subtraction if the dimensions do not match.
        ///     Otherwise subtracts the provided matrix from this matrix.
        /// </summary>
        /// <param name="matrixToSubtract">
        ///     The matrix to subtract.
        /// </param>
        /// <returns>
        ///     The <see cref="bool" />.
        /// </returns>
        public bool TrySubtractMatrix(DoubleMatrix matrixToSubtract)
        {
            if (matrixToSubtract.RowCount != RowCount || matrixToSubtract.ColumnCount != ColumnCount) return false;

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
        ///     Returns false and leaves this matrix unchanged if the matrix to add is not of the correct dimensions.
        ///     Otherwise returns true and adds the values from the provided matrix to this matrix.
        /// </summary>
        /// <param name="matrixToAdd">
        ///     The matrix to add to this matrix.
        /// </param>
        /// <returns>
        ///     The success status, which returns false when matrices are not the same size <see cref="bool" />.
        /// </returns>
        public bool TryAddMatrix(DoubleMatrix matrixToAdd)
        {
            if (matrixToAdd.RowCount != RowCount || matrixToAdd.ColumnCount != ColumnCount) return false;

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
            var result = new DoubleMatrix(ColumnCount, RowCount);
            for (var i = 0; i < RowCount; i++)
            for (var j = 0; j < ColumnCount; j++)
                result[i, j] += _items[j, i];

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
        ///     A newly instantiated matrix <see cref="DoubleMatrix" />. consisting of the dot product
        /// of this matrix and the provided multiplicand
        /// </returns>
        /// <exception cref="ArgumentOutOfRangeException">
        /// </exception>
        public DoubleMatrix DotProduct(DoubleMatrix multiplicand)
        {
            if (ColumnCount != multiplicand.RowCount) throw new ArgumentOutOfRangeException();

            var result = new DoubleMatrix(RowCount, multiplicand.ColumnCount);

            for (var i = 0; i < RowCount; i++)
            for (var j = 0; j < multiplicand.ColumnCount; j++)
                result[i, j] = VectorDotProduct(GetRow(i), multiplicand.GetColumn(j));

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
        ///     The vector dot product.
        /// </summary>
        /// <param name="vec1">
        ///     The first vector.
        /// </param>
        /// <param name="vec2">
        ///     The second vector.
        /// </param>
        /// <returns>
        ///     The vector dot product <see cref="int" />.
        /// </returns>
        /// <exception cref="ArgumentNullException">
        /// </exception>
        /// <exception cref="ArgumentOutOfRangeException">
        /// </exception>
        private static double VectorDotProduct(IReadOnlyList<double> vec1, IReadOnlyList<double> vec2)
        {
            if (vec1 == null || vec2 == null) throw new ArgumentNullException();

            var length = vec1.Count;
            if (vec2.Count != length) throw new ArgumentOutOfRangeException();

            var result = 0d;
            for (var i = 0; i < length; i++) result += vec1[i] * vec2[i];

            return result;
        }

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
        ///     The IsInvertible property.
        /// </summary>
        ///     A <see cref="bool" />  indicating whether the matrix is invertible.
        /// <returns></returns>
        public bool IsInvertible()
        {
            if (!IsSquare()) return false;

            return Math.Abs(GetDeterminantUnsafe(_items)) < double.Epsilon;
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
                if (Math.Abs(_items[i, j]) > double.Epsilon)
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
                if (Math.Abs(_items[j, i]) > double.Epsilon)
                    return false;

            return true;
        }

        /// <summary>
        ///     The GetTrace method.
        /// </summary>
        /// <returns>
        ///     An <see cref="int" /> corresponding to the matrix's trace.
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
        ///     An <see cref="int" /> corresponding to the matrix's trace.
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
        ///     An array of <see cref="int" /> corresponding to the matrix eigenvalues.
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
        ///     An array of <see cref="int" /> corresponding to the matrix eigenvalues.
        /// </returns>
        private static double[] GetEigenValuesUnsafe() => throw new NotImplementedException();

        /// <summary>
        ///     The get eigen vectors method.
        /// </summary>
        /// <returns>
        ///     An array of <see cref="int" /> corresponding to the matrix eigen vectors.
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
        ///     Flag indicating whether the determinant could be calculated <see cref="bool" />.
        /// </returns>
        public bool TryGetDeterminant(out double determinant)
        {
            determinant = 0;

            if (!IsSquare()) return false;

            determinant = GetDeterminantUnsafe(_items);

            return true;
        }

        /// <summary>
        ///     The GetDeterminant method.
        /// </summary>
        /// <returns>
        ///     The determinant for the matrix <see cref="int" />.
        /// </returns>
        public double GetDeterminant()
        {
            if (!IsSquare()) throw new NotSquareMatrixException();

            return GetDeterminantUnsafe(_items);
        }

        /// <summary>
        ///     The unsafe GetDeterminant method.
        ///     Recursively calculates the matrix determinant.
        /// </summary>
        /// <param name="subarray">
        ///     The matrix or sub-matrix for which to calculate the determinant.
        /// </param>
        /// <returns>
        ///     The <see cref="int" /> matrix determinant.
        /// </returns>
        private double GetDeterminantUnsafe(double[,] subarray)
        {
            if (!IsSquare()) throw new NotSquareMatrixException();

            if (subarray.GetLength(0) == 2) return subarray[0, 0] * subarray[1, 1] - subarray[0, 1] * subarray[1, 0];

            var rowLength = subarray.GetLength(0);
            var subsubarray = new double[rowLength - 1, rowLength - 1];
            var results = new double[rowLength];

            for (var i = 0; i < rowLength; i++)
            {
                for (var column = 0; column < rowLength; column++)
                {
                    if (column < i)
                        for (var row = 1; row < rowLength; row++)
                            subsubarray[row - 1, column] = subarray[row, column];

                    if (column <= i) continue;

                    for (var row = 1; row < rowLength; row++)
                        subsubarray[row - 1, column - 1] = subarray[row, column];
                }

                results[i] = subarray[0, i] * GetDeterminantUnsafe(subsubarray);
            }

            var result = 0d;
            var sign = 1;
            for (var i = 0; i < rowLength; i++)
            {
                result += sign * results[i];
                sign *= -1;
            }

            return result;
        }

        /// <summary>
        ///     Overloads * operator for matrix class
        /// </summary>
        /// <param name="a">The left side matrix</param>
        /// <param name="b">The right matrix</param>
        /// <returns>The dot product</returns>
        public static DoubleMatrix operator *(DoubleMatrix a, DoubleMatrix b) => a.DotProduct(b);

        /// <summary>
        ///     Overloads + operator for addition
        /// </summary>
        /// <param name="a">The left side matrix</param>
        /// <param name="b">The right side matrix</param>
        /// <returns></returns>
        public static DoubleMatrix operator +(DoubleMatrix a, DoubleMatrix b)
        {
            if (a == null || b == null) throw new ArgumentNullException();

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
            if (a == null || b == null) throw new ArgumentNullException();

            var result = new DoubleMatrix(a);
            if (result.TrySubtractMatrix(b)) return result;

            throw new ArgumentOutOfRangeException();
        }
    }
}