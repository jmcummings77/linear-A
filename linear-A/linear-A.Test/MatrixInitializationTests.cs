#region Copyright

// --------------------------------------------------------------------------------------------------------------------
// <copyright file="MatrixInitializationTests.cs" company="John-Michael Cummings">
//   John-Michael Cummings 2021
// </copyright>
// <summary>
//   The matrix initialization tests.
// </summary>
// --------------------------------------------------------------------------------------------------------------------

#endregion

namespace linear_A.Test
{
    using System;
    using NUnit.Framework;

    /// <summary>
    ///     The matrix initialization tests.
    /// </summary>
    public class MatrixInitializationTests
    {
        /// <summary>
        ///     The Constructs10X10Identity test.
        /// </summary>
        [Test]
        public void Constructs10X10Identity()
        {
            var testMatrix = new IntegerMatrix(10, true);
            for (var i = 0; i < 10; i++)
            {
                for (var j = 0; j < 10; j++)
                {
                    Assert.AreEqual(i == j ? 1 : 0, testMatrix[i, j]);
                }
            }
        }

        /// <summary>
        ///     The Constructs_3X3Identity test.
        /// </summary>
        [Test]
        public void Constructs_3X3Identity()
        {
            var testMatrix = new IntegerMatrix(3, true);
            Assert.AreEqual(1, testMatrix[0, 0]);
            Assert.AreEqual(0, testMatrix[0, 1]);
            Assert.AreEqual(0, testMatrix[0, 2]);
            Assert.AreEqual(0, testMatrix[1, 0]);
            Assert.AreEqual(1, testMatrix[1, 1]);
            Assert.AreEqual(0, testMatrix[1, 2]);
            Assert.AreEqual(0, testMatrix[2, 0]);
            Assert.AreEqual(0, testMatrix[2, 1]);
            Assert.AreEqual(1, testMatrix[2, 2]);
        }

        /// <summary>
        ///     The Constructs_2X2IdentityNoExtraRows.
        /// </summary>
        [Test]
        public void Constructs_2X2IdentityNoExtraRows()
        {
            var testMatrix = new IntegerMatrix(2, true);
            Assert.AreEqual(1, testMatrix[0, 0]);
            Assert.AreEqual(0, testMatrix[0, 1]);
            Assert.AreEqual(0, testMatrix[1, 0]);
            Assert.AreEqual(1, testMatrix[1, 1]);
            // ReSharper disable once UnusedVariable
            Assert.Throws<IndexOutOfRangeException>(() =>
            {
                var test = testMatrix[2, 0];
            });
        }

        /// <summary>
        ///     The Constructs_2X2IdentityNoExtraColumns test.
        /// </summary>
        [Test]
        public void Constructs_2X2IdentityNoExtraColumns()
        {
            var testMatrix = new IntegerMatrix(2, true);
            Assert.AreEqual(1, testMatrix[0, 0]);
            Assert.AreEqual(0, testMatrix[0, 1]);
            Assert.AreEqual(0, testMatrix[1, 0]);
            Assert.AreEqual(1, testMatrix[1, 1]);
            // ReSharper disable once UnusedVariable
            Assert.Throws<IndexOutOfRangeException>(() =>
            {
                var test = testMatrix[0, 2];
            });
        }

        /// <summary>
        ///     The Constructs_1X1IdentityNoExtraRows test.
        /// </summary>
        [Test]
        public void Constructs_1X1IdentityNoExtraRows()
        {
            var testMatrix = new IntegerMatrix(1, true);
            Assert.AreEqual(1, testMatrix[0, 0]);
            // ReSharper disable once UnusedVariable
            Assert.Throws<IndexOutOfRangeException>(() =>
            {
                var test = testMatrix[1, 0];
            });
        }

        /// <summary>
        ///     The Constructs_1X1IdentityNoExtraColumns test.
        /// </summary>
        [Test]
        public void Constructs_1X1IdentityNoExtraColumns()
        {
            var testMatrix = new IntegerMatrix(1, true);
            Assert.AreEqual(1, testMatrix[0, 0]);
            // ReSharper disable once UnusedVariable
            Assert.Throws<IndexOutOfRangeException>(() =>
            {
                var test = testMatrix[0, 1];
            });
        }

        /// <summary>
        ///     The Constructs_1X1MatrixWithZerosNoExtraRows test.
        /// </summary>
        [Test]
        public void Constructs_1X1MatrixWithZerosNoExtraRows()
        {
            var testMatrix = new IntegerMatrix(1, 1);
            Assert.AreEqual(0, testMatrix[0, 0]);
            // ReSharper disable once UnusedVariable
            Assert.Throws<IndexOutOfRangeException>(() =>
            {
                var test = testMatrix[1, 0];
            });
        }

        /// <summary>
        ///     The Constructs_1X1MatrixWithZerosNoExtraColumns test.
        /// </summary>
        [Test]
        public void Constructs_1X1MatrixWithZerosNoExtraColumns()
        {
            var testMatrix = new IntegerMatrix(1, 1);
            Assert.AreEqual(0, testMatrix[0, 0]);
            // ReSharper disable once UnusedVariable
            Assert.Throws<IndexOutOfRangeException>(() =>
            {
                var test = testMatrix[0, 1];
            });
        }

        /// <summary>
        ///     The Constructs_2X3MatrixWithZeros test.
        /// </summary>
        [Test]
        public void Constructs_2X3MatrixWithZeros()
        {
            var testMatrix = new IntegerMatrix(2, 3);
            for (var i = 0; i < 2; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    Assert.AreEqual(0, testMatrix[i, j]);
                }
            }
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeSquareDimensionsCountZero test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeSquareDimensionsCountZero()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new IntegerMatrix(0, false);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeIdentitySquareDimensionsCountZero test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeIdentitySquareDimensionsCountZero()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new IntegerMatrix(0, true);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeRowCountZero test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeRowCountZero()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new IntegerMatrix(0, 1);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeRowCountNegative test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeRowCountNegative()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new IntegerMatrix(-1, 1);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeColumnCountZero test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeColumnCountZero()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new IntegerMatrix(1, 0);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeColumnCountNegative test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeColumnCountNegative()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new IntegerMatrix(1, -1);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeRowAndColumnCountZero test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeRowAndColumnCountZero()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new IntegerMatrix(0, 0);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeRowAndColumnCountNegative test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeRowAndColumnCountNegative()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new IntegerMatrix(-1, -1);
            });
        }
        
        /// <summary>
        ///     The Constructs10X10Identity test.
        /// </summary>
        [Test]
        public void Constructs10X10DoubleIdentity()
        {
            var testMatrix = new IntegerMatrix(10, true);
            for (var i = 0; i < 10; i++)
            {
                for (var j = 0; j < 10; j++)
                {
                    Assert.AreEqual(i == j ? 1 : 0, testMatrix[i, j]);
                }
            }
        }

        /// <summary>
        ///     The Constructs_3X3Identity test.
        /// </summary>
        [Test]
        public void Constructs_3X3DoubleIdentity()
        {
            var testMatrix = new DoubleMatrix(3, true);
            Assert.AreEqual(1, testMatrix[0, 0]);
            Assert.AreEqual(0, testMatrix[0, 1]);
            Assert.AreEqual(0, testMatrix[0, 2]);
            Assert.AreEqual(0, testMatrix[1, 0]);
            Assert.AreEqual(1, testMatrix[1, 1]);
            Assert.AreEqual(0, testMatrix[1, 2]);
            Assert.AreEqual(0, testMatrix[2, 0]);
            Assert.AreEqual(0, testMatrix[2, 1]);
            Assert.AreEqual(1, testMatrix[2, 2]);
        }

        /// <summary>
        ///     The Constructs_2X2IdentityNoExtraRows.
        /// </summary>
        [Test]
        public void Constructs_2X2DoubleIdentityNoExtraRows()
        {
            var testMatrix = new DoubleMatrix(2, true);
            Assert.AreEqual(1, testMatrix[0, 0]);
            Assert.AreEqual(0, testMatrix[0, 1]);
            Assert.AreEqual(0, testMatrix[1, 0]);
            Assert.AreEqual(1, testMatrix[1, 1]);
            // ReSharper disable once UnusedVariable
            Assert.Throws<IndexOutOfRangeException>(() =>
            {
                var test = testMatrix[2, 0];
            });
        }

        /// <summary>
        ///     The Constructs_2X2IdentityNoExtraColumns test.
        /// </summary>
        [Test]
        public void Constructs_2X2DoubleIdentityNoExtraColumns()
        {
            var testMatrix = new DoubleMatrix(2, true);
            Assert.AreEqual(1, testMatrix[0, 0]);
            Assert.AreEqual(0, testMatrix[0, 1]);
            Assert.AreEqual(0, testMatrix[1, 0]);
            Assert.AreEqual(1, testMatrix[1, 1]);
            // ReSharper disable once UnusedVariable
            Assert.Throws<IndexOutOfRangeException>(() =>
            {
                var test = testMatrix[0, 2];
            });
        }

        /// <summary>
        ///     The Constructs_1X1IdentityNoExtraRows test.
        /// </summary>
        [Test]
        public void Constructs_1X1DoubleIdentityNoExtraRows()
        {
            var testMatrix = new DoubleMatrix(1, true);
            Assert.AreEqual(1, testMatrix[0, 0]);
            // ReSharper disable once UnusedVariable
            Assert.Throws<IndexOutOfRangeException>(() =>
            {
                var test = testMatrix[1, 0];
            });
        }

        /// <summary>
        ///     The Constructs_1X1IdentityNoExtraColumns test.
        /// </summary>
        [Test]
        public void Constructs_1X1DoubleIdentityNoExtraColumns()
        {
            var testMatrix = new DoubleMatrix(1, true);
            Assert.AreEqual(1, testMatrix[0, 0]);
            // ReSharper disable once UnusedVariable
            Assert.Throws<IndexOutOfRangeException>(() =>
            {
                var test = testMatrix[0, 1];
            });
        }

        /// <summary>
        ///     The Constructs_1X1MatrixWithZerosNoExtraRows test.
        /// </summary>
        [Test]
        public void Constructs_1X1DoubleMatrixWithZerosNoExtraRows()
        {
            var testMatrix = new DoubleMatrix(1, 1);
            Assert.AreEqual(0, testMatrix[0, 0]);
            // ReSharper disable once UnusedVariable
            Assert.Throws<IndexOutOfRangeException>(() =>
            {
                var test = testMatrix[1, 0];
            });
        }

        /// <summary>
        ///     The Constructs_1X1MatrixWithZerosNoExtraColumns test.
        /// </summary>
        [Test]
        public void Constructs_1X1DoubleMatrixWithZerosNoExtraColumns()
        {
            var testMatrix = new DoubleMatrix(1, 1);
            Assert.AreEqual(0, testMatrix[0, 0]);
            // ReSharper disable once UnusedVariable
            Assert.Throws<IndexOutOfRangeException>(() =>
            {
                var test = testMatrix[0, 1];
            });
        }

        /// <summary>
        ///     The Constructs_2X3MatrixWithZeros test.
        /// </summary>
        [Test]
        public void Constructs_2X3DoubleMatrixWithZeros()
        {
            var testMatrix = new DoubleMatrix(2, 3);
            for (var i = 0; i < 2; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    Assert.AreEqual(0, testMatrix[i, j]);
                }
            }
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeSquareDimensionsCountZero test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeSquareDimensionsCountZeroDouble()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new DoubleMatrix(0, false);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeIdentitySquareDimensionsCountZero test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeIdentitySquareDimensionsCountZeroDouble()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new DoubleMatrix(0, true);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeRowCountZero test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeRowCountZeroDouble()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new DoubleMatrix(0, 1);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeRowCountNegative test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeRowCountNegativeDouble()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new DoubleMatrix(-1, 1);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeColumnCountZero test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeColumnCountZeroDouble()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new DoubleMatrix(1, 0);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeColumnCountNegative test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeColumnCountNegativeDouble()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new DoubleMatrix(1, -1);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeRowAndColumnCountZero test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeRowAndColumnCountZeroDouble()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new DoubleMatrix(0, 0);
            });
        }

        /// <summary>
        ///     The ThrowsArgumentOutOfRangeRowAndColumnCountNegative test.
        /// </summary>
        [Test]
        public void ThrowsArgumentOutOfRangeRowAndColumnCountNegativeDouble()
        {
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                var testMatrix = new DoubleMatrix(-1, -1);
            });
        }
    }
}