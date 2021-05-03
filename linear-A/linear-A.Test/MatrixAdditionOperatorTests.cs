#region Copyright

// --------------------------------------------------------------------------------------------------------------------
// <copyright file="MatrixAdditionTests.cs" company="John-Michael Cummings">
//   John-Michael Cummings 2021
// </copyright>
// <summary>
//   The matrix addition tests.
// </summary>
// --------------------------------------------------------------------------------------------------------------------

#endregion

namespace linear_A.Test
{
    using NUnit.Framework;
    using System;

    /// <summary>
    ///     Test class for the Matrix addition overloaded operator "+" 
    /// </summary>
    public class MatrixAdditionOperatorTests
    {
        /// <summary>
        ///     The AddZeros_1X1MatricesWithOperatorResultIsZero test.
        /// </summary>
        [Test]
        public void AddZeros_1X1MatricesWithOperatorResultIsZero()
        {
            var testMatrix1 = new IntegerMatrix(1, false);
            var testMatrix2 = new IntegerMatrix(1, false);
            var testMatrix3 = testMatrix1 + testMatrix2;
            Assert.AreEqual(testMatrix3[0, 0], 0);
        }

        /// <summary>
        ///     The AddOnes_1X1MatricesWithOperatorResultIsTwo test.
        /// </summary>
        [Test]
        public void AddOnes_1X1MatricesWithOperatorResultIsTwo()
        {
            var testMatrix1 = new IntegerMatrix(1, true);
            var testMatrix2 = new IntegerMatrix(1, true);
            var testMatrix3 = testMatrix1 + testMatrix2;
            Assert.AreEqual(testMatrix3[0, 0], 2);
        }

        /// <summary>
        ///     AddDifferentDimensionMatricesWithOperatorFails test.
        /// </summary>
        [Test]
        public void AddDifferentDimensionMatricesWithOperatorFails()
        {
            var testMatrix1 = new IntegerMatrix(1, false);
            var testMatrix2 = new IntegerMatrix(2, true);
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() => {var m = testMatrix1 + testMatrix2;});
        }

        /// <summary>
        ///     The AddNegative_3X3ToEquivalent_3X3WithOperatorResultAllAreZero test.
        /// </summary>
        [Test]
        public void AddNegative_3X3ToEquivalent_3X3WithOperatorResultAllAreZero()
        {
            var testMatrix1 = new IntegerMatrix(3, false);
            var testMatrix2 = new IntegerMatrix(3, false);

            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    var value = i + j + 1;
                    testMatrix1[i, j] = value;
                    testMatrix2[j, i] = -value;
                }
            }

            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    Assert.AreNotEqual(testMatrix1[i, j], 0);
                    Assert.AreNotEqual(testMatrix2[j, i], 0);
                }
            }

            var testMatrix3 = testMatrix1 + testMatrix2;
            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    Assert.AreEqual(testMatrix3[i, j], 0);
                }
            }
        }

        /// <summary>
        ///     The AddMatrixToNegativeEquivalentWithOperatorResultAllAreZero test.
        /// </summary>
        [Test]
        public void AddMatrixToNegativeEquivalentWithOperatorResultAllAreZero()
        {
            for (var rowCount = 1; rowCount < 50; rowCount++)
            {
                for (var columnCount = 1; columnCount < 50; columnCount++)
                {
                    var testMatrix1 = new IntegerMatrix(rowCount, columnCount);
                    var testMatrix2 = new IntegerMatrix(rowCount, columnCount);

                    for (var i = 0; i < rowCount; i++)
                    {
                        for (var j = 0; j < columnCount; j++)
                        {
                            var value = i + j + 1;
                            testMatrix1[i, j] = value;
                            testMatrix2[i, j] = -value;
                        }
                    }

                    for (var i = 0; i < rowCount; i++)
                    {
                        for (var j = 0; j < columnCount; j++)
                        {
                            Assert.AreNotEqual(testMatrix1[i, j], 0);
                            Assert.AreNotEqual(testMatrix2[i, j], 0);
                        }
                    }

                    var testMatrix3 = testMatrix1 + testMatrix2;
                    for (var i = 0; i < rowCount; i++)
                    {
                        for (var j = 0; j < columnCount; j++)
                        {
                            Assert.AreEqual(testMatrix3[i, j], 0);
                        }
                    }
                }
            }
        }

        /// <summary>
        ///     The hacker rank test case #1 with addition operator.
        /// </summary>
        [Test]
        public void HackerRankTestCase1WithOperator()
        {
            var rowCount = 3;
            var columnCount = 3;
            var testMatrix1 = new IntegerMatrix(rowCount, columnCount);
            var testMatrix2 = new IntegerMatrix(rowCount, columnCount);

            testMatrix1[0, 0] = 1;
            testMatrix1[0, 1] = 2;
            testMatrix1[0, 2] = 3;
            testMatrix1[1, 0] = 2;
            testMatrix1[1, 1] = 3;
            testMatrix1[1, 2] = 4;
            testMatrix1[2, 0] = 1;
            testMatrix1[2, 1] = 1;
            testMatrix1[2, 2] = 1;

            testMatrix2[0, 0] = 4;
            testMatrix2[0, 1] = 5;
            testMatrix2[0, 2] = 6;
            testMatrix2[1, 0] = 7;
            testMatrix2[1, 1] = 8;
            testMatrix2[1, 2] = 9;
            testMatrix2[2, 0] = 4;
            testMatrix2[2, 1] = 5;
            testMatrix2[2, 2] = 7;

            var testMatrix3 = testMatrix1 + testMatrix2;

            Assert.AreEqual(5, testMatrix3[0, 0]);
            Assert.AreEqual(7, testMatrix3[0, 1]);
            Assert.AreEqual(9, testMatrix3[0, 2]);
            Assert.AreEqual(9, testMatrix3[1, 0]);
            Assert.AreEqual(11, testMatrix3[1, 1]);
            Assert.AreEqual(13, testMatrix3[1, 2]);
            Assert.AreEqual(5, testMatrix3[2, 0]);
            Assert.AreEqual(6, testMatrix3[2, 1]);
            Assert.AreEqual(8, testMatrix3[2, 2]);
        }
        
        /// <summary>
        ///     The AddZeros_1X1MatricesWithOperatorResultIsZero test.
        /// </summary>
        [Test]
        public void AddZeros_1X1DoubleMatricesWithOperatorResultIsZero()
        {
            var testMatrix1 = new DoubleMatrix(1, false);
            var testMatrix2 = new DoubleMatrix(1, false);
            var testMatrix3 = testMatrix1 + testMatrix2;
            Assert.AreEqual(testMatrix3[0, 0], 0);
        }

        /// <summary>
        ///     The AddOnes_1X1MatricesWithOperatorResultIsTwo test.
        /// </summary>
        [Test]
        public void AddOnes_1X1DoubleMatricesWithOperatorResultIsTwo()
        {
            var testMatrix1 = new DoubleMatrix(1, true);
            var testMatrix2 = new DoubleMatrix(1, true);
            var testMatrix3 = testMatrix1 + testMatrix2;
            Assert.AreEqual(testMatrix3[0, 0], 2);
        }

        /// <summary>
        ///     AddDifferentDimensionMatricesWithOperatorFails test.
        /// </summary>
        [Test]
        public void AddDifferentDoubleDimensionMatricesWithOperatorFails()
        {
            var testMatrix1 = new DoubleMatrix(1, false);
            var testMatrix2 = new DoubleMatrix(2, true);
            // ReSharper disable once UnusedVariable
            Assert.Throws<ArgumentOutOfRangeException>(() => {var m = testMatrix1 + testMatrix2;});
        }

        /// <summary>
        ///     The AddNegative_3X3ToEquivalent_3X3WithOperatorResultAllAreZero test.
        /// </summary>
        [Test]
        public void AddNegative_3X3DoubleToEquivalent_3X3WithOperatorResultAllAreZero()
        {
            var testMatrix1 = new DoubleMatrix(3, false);
            var testMatrix2 = new DoubleMatrix(3, false);

            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    var value = i + j + 1;
                    testMatrix1[i, j] = value;
                    testMatrix2[j, i] = -value;
                }
            }

            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    Assert.AreNotEqual(testMatrix1[i, j], 0);
                    Assert.AreNotEqual(testMatrix2[j, i], 0);
                }
            }

            var testMatrix3 = testMatrix1 + testMatrix2;
            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    Assert.AreEqual(testMatrix3[i, j], 0);
                }
            }
        }

        /// <summary>
        ///     The AddMatrixToNegativeEquivalentWithOperatorResultAllAreZero test.
        /// </summary>
        [Test]
        public void AddDoubleMatrixToNegativeEquivalentWithOperatorResultAllAreZero()
        {
            for (var rowCount = 1; rowCount < 50; rowCount++)
            {
                for (var columnCount = 1; columnCount < 50; columnCount++)
                {
                    var testMatrix1 = new DoubleMatrix(rowCount, columnCount);
                    var testMatrix2 = new DoubleMatrix(rowCount, columnCount);

                    for (var i = 0; i < rowCount; i++)
                    {
                        for (var j = 0; j < columnCount; j++)
                        {
                            var value = i + j + 1;
                            testMatrix1[i, j] = value;
                            testMatrix2[i, j] = -value;
                        }
                    }

                    for (var i = 0; i < rowCount; i++)
                    {
                        for (var j = 0; j < columnCount; j++)
                        {
                            Assert.AreNotEqual(testMatrix1[i, j], 0);
                            Assert.AreNotEqual(testMatrix2[i, j], 0);
                        }
                    }

                    var testMatrix3 = testMatrix1 + testMatrix2;
                    for (var i = 0; i < rowCount; i++)
                    {
                        for (var j = 0; j < columnCount; j++)
                        {
                            Assert.AreEqual(testMatrix3[i, j], 0);
                        }
                    }
                }
            }
        }

        /// <summary>
        ///     The hacker rank test case #1 with addition operator.
        /// </summary>
        [Test]
        public void HackerRankTestCase1WithOperatorDoubleMatrix()
        {
            const int rowCount = 3;
            const int columnCount = 3;
            var testMatrix1 = new DoubleMatrix(rowCount, columnCount);
            var testMatrix2 = new DoubleMatrix(rowCount, columnCount);

            testMatrix1[0, 0] = 1;
            testMatrix1[0, 1] = 2;
            testMatrix1[0, 2] = 3;
            testMatrix1[1, 0] = 2;
            testMatrix1[1, 1] = 3;
            testMatrix1[1, 2] = 4;
            testMatrix1[2, 0] = 1;
            testMatrix1[2, 1] = 1;
            testMatrix1[2, 2] = 1;

            testMatrix2[0, 0] = 4;
            testMatrix2[0, 1] = 5;
            testMatrix2[0, 2] = 6;
            testMatrix2[1, 0] = 7;
            testMatrix2[1, 1] = 8;
            testMatrix2[1, 2] = 9;
            testMatrix2[2, 0] = 4;
            testMatrix2[2, 1] = 5;
            testMatrix2[2, 2] = 7;

            var testMatrix3 = testMatrix1 + testMatrix2;

            Assert.AreEqual(5, testMatrix3[0, 0]);
            Assert.AreEqual(7, testMatrix3[0, 1]);
            Assert.AreEqual(9, testMatrix3[0, 2]);
            Assert.AreEqual(9, testMatrix3[1, 0]);
            Assert.AreEqual(11, testMatrix3[1, 1]);
            Assert.AreEqual(13, testMatrix3[1, 2]);
            Assert.AreEqual(5, testMatrix3[2, 0]);
            Assert.AreEqual(6, testMatrix3[2, 1]);
            Assert.AreEqual(8, testMatrix3[2, 2]);
        }
    }
}
