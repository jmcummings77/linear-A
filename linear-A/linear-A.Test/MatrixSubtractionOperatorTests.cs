#region Copyright

// --------------------------------------------------------------------------------------------------------------------
// <copyright file="MatrixSubtractionOperatorTests.cs" company="John-Michael Cummings">
//   John-Michael Cummings 2021
// </copyright>
// <summary>
//   The matrix subtraction operator tests.
// </summary>
// --------------------------------------------------------------------------------------------------------------------

#endregion

using Matricks.Core;

namespace linear_A.Test
{
    using System;
    using NUnit.Framework;

    /// <summary>
    ///     The matrix subtraction operator "-" tests.
    /// </summary>
    public class MatrixSubtractionOperatorTests
    {
        /// <summary>
        ///     The SubtractZeros_1X1MatricesWithOperatorAndIsZero test.
        /// </summary>
        [Test]
        public void SubtractZeros_1X1MatricesWithOperatorAndIsZero()
        {
            var testMatrix1 = new IntegerMatrix(1, false);
            var testMatrix2 = new IntegerMatrix(1, false);
            
            var testMatrix3 = testMatrix1 - testMatrix2;
            
            Assert.AreEqual(0, testMatrix3[0, 0]);
        }

        /// <summary>
        ///     The SubtractTwoFromOne_1X1MatricesWithOperatorAndIsZero test.
        /// </summary>
        [Test]
        public void SubtractTwoFromOne_1X1MatricesWithOperatorAndIsZero()
        {
            var testMatrix1 = new IntegerMatrix(1, true);
            var testMatrix2 = new IntegerMatrix(1, true);
            testMatrix1[0, 0]++;
            
            var testMatrix3 = testMatrix1 - testMatrix2;
            
            Assert.AreEqual(1, testMatrix3[0, 0]);
        }

        /// <summary>
        ///     The SubtractDifferentDimensionMatricesFails test.
        /// </summary>
        [Test]
        public void SubtractDifferentDimensionMatricesFails()
        {
            var testMatrix1 = new IntegerMatrix(1, false);
            var testMatrix2 = new IntegerMatrix(2, true);
            
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                // ReSharper disable once UnusedVariable
                var testMatrix3 = testMatrix1 - testMatrix2;
            });
        }

        /// <summary>
        ///     The SubtractNegative_3X3FromEquivalent_3X3WithOperatorAndAllAreDoubled test.
        /// </summary>
        [Test]
        public void SubtractNegative_3X3FromEquivalent_3X3WithOperatorAndAllAreDoubled()
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
                    Assert.AreNotEqual(0, testMatrix1[i, j]);
                    Assert.AreNotEqual(0, testMatrix2[j, i]);
                }
            }

            var testMatrix3 = testMatrix1 - testMatrix2;
            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    var expectedValue = (i + j + 1) * 2;
                    Assert.AreEqual(expectedValue, testMatrix3[i, j], 0);
                }
            }
        }

        /// <summary>
        ///     The SubtractMatrixToNegativeEquivalentWithOperatorAndAllAreZero test.
        /// </summary>
        [Test]
        public void SubtractMatrixToNegativeEquivalentWithOperatorAndAllAreZero()
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
                            Assert.AreNotEqual(0, testMatrix1[i, j]);
                            Assert.AreNotEqual(0, testMatrix2[i, j]);
                        }
                    }

                    var testMatrix3 = testMatrix1 - testMatrix2;
                    for (var i = 0; i < rowCount; i++)
                    {
                        for (var j = 0; j < columnCount; j++)
                        {
                            var expectedValue = (i + j + 1) * 2;
                            Assert.AreEqual(expectedValue, testMatrix3[i, j]);
                        }
                    }
                }
            }
        }

        /// <summary>
        ///     The hacker rank test case 2.
        /// </summary>
        [Test]
        public void HackerRankTestCase2()
        {
            const int rowCount = 3;
            const int columnCount = 3;
            var testMatrix1 = new IntegerMatrix(rowCount, columnCount)
            {
                [0, 0] = 1,
                [0, 1] = 2,
                [0, 2] = 3,
                [1, 0] = 2,
                [1, 1] = 3,
                [1, 2] = 4,
                [2, 0] = 1,
                [2, 1] = 1,
                [2, 2] = 1
            };
            var testMatrix2 = new IntegerMatrix(rowCount, columnCount)
            {
                [0, 0] = 4,
                [0, 1] = 5,
                [0, 2] = 6,
                [1, 0] = 7,
                [1, 1] = 8,
                [1, 2] = 9,
                [2, 0] = 4,
                [2, 1] = 5,
                [2, 2] = 0
            };

            var testMatrix3 = testMatrix1 - testMatrix2;

            Assert.AreEqual(-3, testMatrix3[0, 0]);
            Assert.AreEqual(-3, testMatrix3[0, 1]);
            Assert.AreEqual(-3, testMatrix3[0, 2]);
            Assert.AreEqual(-5, testMatrix3[1, 0]);
            Assert.AreEqual(-5, testMatrix3[1, 1]);
            Assert.AreEqual(-5, testMatrix3[1, 2]);
            Assert.AreEqual(-3, testMatrix3[2, 0]);
            Assert.AreEqual(-4, testMatrix3[2, 1]);
            Assert.AreEqual(1, testMatrix3[2, 2]);
        }
        
        /// <summary>
        ///     The SubtractZeros_1X1MatricesWithOperatorAndIsZero test.
        /// </summary>
        [Test]
        public void SubtractZeros_1X1DoubleMatricesWithOperatorAndIsZero()
        {
            var testMatrix1 = new DoubleMatrix(1, false);
            var testMatrix2 = new DoubleMatrix(1, false);
            
            var testMatrix3 = testMatrix1 - testMatrix2;
            
            Assert.AreEqual(0, testMatrix3[0, 0]);
        }

        /// <summary>
        ///     The SubtractTwoFromOne_1X1MatricesWithOperatorAndIsZero test.
        /// </summary>
        [Test]
        public void SubtractTwoFromOne_1X1DoubleMatricesWithOperatorAndIsZero()
        {
            var testMatrix1 = new DoubleMatrix(1, true);
            var testMatrix2 = new DoubleMatrix(1, true);
            testMatrix1[0, 0]++;
            
            var testMatrix3 = testMatrix1 - testMatrix2;
            
            Assert.AreEqual(1, testMatrix3[0, 0]);
        }

        /// <summary>
        ///     The SubtractDifferentDimensionMatricesFails test.
        /// </summary>
        [Test]
        public void SubtractDifferentDimensionDoubleMatricesFails()
        {
            var testMatrix1 = new DoubleMatrix(1, false);
            var testMatrix2 = new DoubleMatrix(2, true);
            Assert.Throws<ArgumentOutOfRangeException>(() =>
            {
                // ReSharper disable once UnusedVariable
                var testMatrix3 = testMatrix1 - testMatrix2;
            });
        }

        /// <summary>
        ///     The SubtractNegative_3X3FromEquivalent_3X3WithOperatorAndAllAreDoubled test.
        /// </summary>
        [Test]
        public void SubtractNegative_3X3DoubleFromEquivalent_3X3WithOperatorAndAllAreDoubled()
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
                    Assert.AreNotEqual(0, testMatrix1[i, j]);
                    Assert.AreNotEqual(0, testMatrix2[j, i]);
                }
            }

            var testMatrix3 = testMatrix1 - testMatrix2;
            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    var expectedValue = (i + j + 1) * 2;
                    Assert.AreEqual(expectedValue, testMatrix3[i, j], 0);
                }
            }
        }

        /// <summary>
        ///     The SubtractMatrixToNegativeEquivalentWithOperatorAndAllAreZero test.
        /// </summary>
        [Test]
        public void SubtractDoubleMatrixToNegativeEquivalentWithOperatorAndAllAreZero()
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
                            Assert.AreNotEqual(0, testMatrix1[i, j]);
                            Assert.AreNotEqual(0, testMatrix2[i, j]);
                        }
                    }

                    var testMatrix3 = testMatrix1 - testMatrix2;
                    for (var i = 0; i < rowCount; i++)
                    {
                        for (var j = 0; j < columnCount; j++)
                        {
                            var expectedValue = (i + j + 1) * 2;
                            Assert.AreEqual(expectedValue, testMatrix3[i, j]);
                        }
                    }
                }
            }
        }

        /// <summary>
        ///     The hacker rank test case 2.
        /// </summary>
        [Test]
        public void HackerRankTestCase2Double()
        {
            const int rowCount = 3;
            const int columnCount = 3;
            var testMatrix1 = new DoubleMatrix(rowCount, columnCount)
            {
                [0, 0] = 1,
                [0, 1] = 2,
                [0, 2] = 3,
                [1, 0] = 2,
                [1, 1] = 3,
                [1, 2] = 4,
                [2, 0] = 1,
                [2, 1] = 1,
                [2, 2] = 1
            };
            var testMatrix2 = new DoubleMatrix(rowCount, columnCount)
            {
                [0, 0] = 4,
                [0, 1] = 5,
                [0, 2] = 6,
                [1, 0] = 7,
                [1, 1] = 8,
                [1, 2] = 9,
                [2, 0] = 4,
                [2, 1] = 5,
                [2, 2] = 0
            };

            var testMatrix3 = testMatrix1 - testMatrix2;

            Assert.AreEqual(-3, testMatrix3[0, 0]);
            Assert.AreEqual(-3, testMatrix3[0, 1]);
            Assert.AreEqual(-3, testMatrix3[0, 2]);
            Assert.AreEqual(-5, testMatrix3[1, 0]);
            Assert.AreEqual(-5, testMatrix3[1, 1]);
            Assert.AreEqual(-5, testMatrix3[1, 2]);
            Assert.AreEqual(-3, testMatrix3[2, 0]);
            Assert.AreEqual(-4, testMatrix3[2, 1]);
            Assert.AreEqual(1, testMatrix3[2, 2]);
        }
    }
}