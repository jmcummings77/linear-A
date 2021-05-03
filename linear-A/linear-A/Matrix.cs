#region Copyright

// --------------------------------------------------------------------------------------------------------------------
// <copyright file="Matrix.cs" company="John-Michael Cummings">
//   John-Michael Cummings 2021
// </copyright>
// <summary>
//   Generic matrix and vector base classes
// </summary>
// --------------------------------------------------------------------------------------------------------------------

#endregion

namespace linear_A
{
    using System;
    using System.Collections;
    using System.Collections.Generic;
    using System.Linq;

    /// <summary>
    /// 
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class Vector<T> : IEquatable<Vector<T>>, IList<T> where T : struct, IEquatable<T>
    {
        /// <summary>
        ///     Initializes a vector with the provided length
        /// </summary>
        /// <param name="length"></param>
        public Vector(int length)
        {
            _items = new T[length];
            Count = length;
        }

        /// <summary>
        ///     The index accessor method.
        /// <param name="index">
        ///     The item index.
        /// </param>
        /// <returns>
        ///     The <see cref="int" />.
        /// </returns>
        /// </summary>
        public T this[int index]
        {
            get => _items[index];
            set => _items[index] = value;
        }
        
        /// <summary>
        ///     Gets the row count.
        /// </summary>
        public int Count { get; private set; }
        
        /// <summary>
        ///     Gets the raw matrix data storage.
        /// </summary>
        private T[] _items;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="item"></param>
        /// <returns></returns>
        public int IndexOf(T item)
        {
            return Array.IndexOf(_items, item);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="index"></param>
        /// <param name="item"></param>
        public void Insert(int index, T item)
        {
            if (index > Count + 1)
            {
                throw new IndexOutOfRangeException("Index must be less than the vector length plus one.");
            }
            
            Count++;
            var items = new T[Count];
            _items.CopyTo(items, index - 1);
            items[index] = item;

            for (var i = index + 1; i < Count; i++)
            {
                items[i] = _items[i - 1];
            }

            _items = items;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="index"></param>
        public void RemoveAt(int index)
        {
            if (index > Count - 1)
            {
                throw new IndexOutOfRangeException("Index is greater than vector length.");
            }
            
            if (index < 0)
            {
                throw new IndexOutOfRangeException("Index cannot be negative.");
            }
            Count--;
            var items = new T[Count];
            _items.CopyTo(items, index - 1);

            for (var i = index; i < Count - 1; i++)
            {
                _items[i] = _items[i + 1];
            }

            _items = items;
        }
        
        /// <summary>
        /// 
        /// </summary>
        /// <param name="item"></param>
        public void Add(T item)
        {
            Count++;
            var items = new T[Count];
            _items.CopyTo(items, 0);
            items[Count - 1] = item;
        }

        /// <summary>
        /// 
        /// </summary>
        public void Clear() => _items = new T[Count];

        /// <summary>
        /// 
        /// </summary>
        /// <param name="item"></param>
        /// <returns></returns>
        public bool Contains(T item) => _items.Contains(item);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="array"></param>
        /// <param name="arrayIndex"></param>
        public void CopyTo(T[] array, int arrayIndex) => _items.CopyTo(array, arrayIndex);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="item"></param>
        /// <returns></returns>
        public bool Remove(T item)
        {
            var list = _items.ToList();
            if (!list.Remove(item)) return false;
            _items = list.ToArray();
            return true;
        }

        /// <summary>
        /// 
        /// </summary>
        public bool IsReadOnly => false;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Equals(Vector<T> other)
        {
            if (other == null) return false;
            if (other.Count != Count) return false;
            for (var i = 0; i < Count; i++)
            {
                if (!other[i].Equals(_items[i])) return false;
            }

            return true;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public IEnumerator<T> GetEnumerator() => (IEnumerator<T>) _items.GetEnumerator();

        IEnumerator IEnumerable.GetEnumerator() => _items.GetEnumerator();
    }

    /// <summary>
    /// 
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class Matrix<T> : IEquatable<Matrix<T>> where T : struct, IEquatable<T>
    {
        /// <summary>
        /// Initializes a new instance of the Matrix class.
        /// </summary>
        public Matrix(int rowCount, int columnCount)
        {
            _items = new T[rowCount * columnCount];
            RowCount = rowCount;
            ColumnCount = ColumnCount;
        }

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
        public T this[int rowIndex, int columnIndex]
        {
            get => _items[columnIndex * RowCount + rowIndex];
            set => _items[columnIndex * RowCount + rowIndex] = value;
        }
        
        /// <summary>
        ///     Gets the raw matrix data storage.
        /// </summary>
        private readonly T[] _items;

        /// <summary>
        ///     Gets the total count of elements.
        /// </summary>
        public int Count => RowCount * ColumnCount;
        
        /// <summary>
        ///     Gets the row count.
        /// </summary>
        public int RowCount { get; }

        /// <summary>
        ///     Gets the column count.
        /// </summary>
        public int ColumnCount { get; }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Equals(Matrix<T> other)
        {
            if (other == null) return false;
            if (other.RowCount != RowCount) return false;
            if (other.ColumnCount != ColumnCount) return false;
            
            for (var i = 0; i < RowCount; i++)
            {
                for (var j = 0; j < ColumnCount; j++)
                {
                    if (!other[i,j].Equals(_items[j * RowCount + i])) return false;
                }
            }

            return true;
        }
    }
}