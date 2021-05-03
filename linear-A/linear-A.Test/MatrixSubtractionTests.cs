#region Copyright

// --------------------------------------------------------------------------------------------------------------------
// <copyright file="MatrixSubtractionTests.cs" company="John-Michael Cummings">
//   John-Michael Cummings 2021
// </copyright>
// <summary>
//   The matrix subtraction tests.
// </summary>
// --------------------------------------------------------------------------------------------------------------------

#endregion

using Matricks.Core;

namespace linear_A.Test
{
    using NUnit.Framework;

    /// <summary>
    ///     The matrix subtraction tests.
    /// </summary>
    public class MatrixSubtractionTests
    {
        /// <summary>
        ///     The TrySubtractZeros_1X1MatricesResultSuccessAndIsZero test.
        /// </summary>
        [Test]
        public void TrySubtractZeros_1X1MatricesResultSuccessAndIsZero()
        {
            var testMatrix1 = new IntegerMatrix(1, false);
            var testMatrix2 = new IntegerMatrix(1, false);
            
            var success = testMatrix1.TrySubtractMatrix(testMatrix2);
            
            Assert.IsTrue(success);
            Assert.AreEqual(0, testMatrix1[0, 0]);
        }

        /// <summary>
        ///     The TrySubtractTwoFromOne_1X1MatricesResultSuccessAndIsZero test.
        /// </summary>
        [Test]
        public void TrySubtractTwoFromOne_1X1MatricesResultSuccessAndIsZero()
        {
            var testMatrix1 = new IntegerMatrix(1, true);
            var testMatrix2 = new IntegerMatrix(1, true);
            testMatrix1[0, 0]++;
            
            var success = testMatrix1.TrySubtractMatrix(testMatrix2);
            
            Assert.IsTrue(success);
            Assert.AreEqual(1, testMatrix1[0, 0]);
        }

        /// <summary>
        ///     The TrySubtractDifferentDimensionMatricesFails test.
        /// </summary>
        [Test]
        public void TrySubtractDifferentDimensionMatricesFails()
        {
            var testMatrix1 = new IntegerMatrix(1, false);
            var testMatrix2 = new IntegerMatrix(2, true);
            
            var success = testMatrix1.TrySubtractMatrix(testMatrix2);
            
            Assert.IsFalse(success);
            Assert.AreEqual(0, testMatrix1[0, 0]);
        }

        /// <summary>
        ///     The TrySubtractNegative_3X3FromEquivalent_3X3ResultSuccessAndAllAreDoubled test.
        /// </summary>
        [Test]
        public void TrySubtractNegative_3X3FromEquivalent_3X3ResultSuccessAndAllAreDoubled()
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

            var success = testMatrix1.TrySubtractMatrix(testMatrix2);
            Assert.IsTrue(success);
            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    var expectedValue = (i + j + 1) * 2;
                    Assert.AreEqual(expectedValue, testMatrix1[i, j], 0);
                }
            }
        }

        /// <summary>
        ///     The TrySubtractMatrixToNegativeEquivalentResultSuccessAndAllAreZero test.
        /// </summary>
        [Test]
        public void TrySubtractMatrixToNegativeEquivalentResultSuccessAndAllAreZero()
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

                    var success = testMatrix1.TrySubtractMatrix(testMatrix2);
                    Assert.IsTrue(success);
                    for (var i = 0; i < rowCount; i++)
                    {
                        for (var j = 0; j < columnCount; j++)
                        {
                            var expectedValue = (i + j + 1) * 2;
                            Assert.AreEqual(expectedValue, testMatrix1[i, j]);
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

            testMatrix1.TrySubtractMatrix(testMatrix2);

            Assert.AreEqual(-3, testMatrix1[0, 0]);
            Assert.AreEqual(-3, testMatrix1[0, 1]);
            Assert.AreEqual(-3, testMatrix1[0, 2]);
            Assert.AreEqual(-5, testMatrix1[1, 0]);
            Assert.AreEqual(-5, testMatrix1[1, 1]);
            Assert.AreEqual(-5, testMatrix1[1, 2]);
            Assert.AreEqual(-3, testMatrix1[2, 0]);
            Assert.AreEqual(-4, testMatrix1[2, 1]);
            Assert.AreEqual(1,  testMatrix1[2, 2]);
        }
        
        /// <summary>
        ///     The TrySubtractZeros_1X1MatricesResultSuccessAndIsZero test.
        /// </summary>
        [Test]
        public void TrySubtractZeros_1X1DoubleMatricesResultSuccessAndIsZero()
        {
            var testMatrix1 = new DoubleMatrix(1, false);
            var testMatrix2 = new DoubleMatrix(1, false);
            
            var success = testMatrix1.TrySubtractMatrix(testMatrix2);
            
            Assert.IsTrue(success);
            Assert.AreEqual(0, testMatrix1[0, 0]);
        }

        /// <summary>
        ///     The TrySubtractTwoFromOne_1X1MatricesResultSuccessAndIsZero test.
        /// </summary>
        [Test]
        public void TrySubtractTwoFromOne_1X1DoubleMatricesResultSuccessAndIsZero()
        {
            var testMatrix1 = new DoubleMatrix(1, true);
            var testMatrix2 = new DoubleMatrix(1, true);
            testMatrix1[0, 0]++;
            
            var success = testMatrix1.TrySubtractMatrix(testMatrix2);
            
            Assert.IsTrue(success);
            Assert.AreEqual(1, testMatrix1[0, 0]);
        }

        /// <summary>
        ///     The TrySubtractDifferentDimensionMatricesFails test.
        /// </summary>
        [Test]
        public void TrySubtractDifferentDimensionDoubleMatricesFails()
        {
            var testMatrix1 = new DoubleMatrix(1, false);
            var testMatrix2 = new DoubleMatrix(2, true);
            
            var success = testMatrix1.TrySubtractMatrix(testMatrix2);
            
            Assert.IsFalse(success);
            Assert.AreEqual(0, testMatrix1[0, 0]);
        }

        /// <summary>
        ///     The TrySubtractNegative_3X3FromEquivalent_3X3ResultSuccessAndAllAreDoubled test.
        /// </summary>
        [Test]
        public void TrySubtractNegative_3X3DoubleFromEquivalent_3X3ResultSuccessAndAllAreDoubled()
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

            var success = testMatrix1.TrySubtractMatrix(testMatrix2);
            Assert.IsTrue(success);
            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    var expectedValue = (i + j + 1) * 2;
                    Assert.AreEqual(expectedValue, testMatrix1[i, j], 0);
                }
            }
        }

        /// <summary>
        ///     The TrySubtractMatrixToNegativeEquivalentResultSuccessAndAllAreZero test.
        /// </summary>
        [Test]
        public void TrySubtractDoubleMatrixToNegativeEquivalentResultSuccessAndAllAreZero()
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

                    var success = testMatrix1.TrySubtractMatrix(testMatrix2);
                    Assert.IsTrue(success);
                    for (var i = 0; i < rowCount; i++)
                    {
                        for (var j = 0; j < columnCount; j++)
                        {
                            var expectedValue = (i + j + 1) * 2;
                            Assert.AreEqual(expectedValue, testMatrix1[i, j]);
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

            testMatrix1.TrySubtractMatrix(testMatrix2);

            Assert.AreEqual(-3, testMatrix1[0, 0]);
            Assert.AreEqual(-3, testMatrix1[0, 1]);
            Assert.AreEqual(-3, testMatrix1[0, 2]);
            Assert.AreEqual(-5, testMatrix1[1, 0]);
            Assert.AreEqual(-5, testMatrix1[1, 1]);
            Assert.AreEqual(-5, testMatrix1[1, 2]);
            Assert.AreEqual(-3, testMatrix1[2, 0]);
            Assert.AreEqual(-4, testMatrix1[2, 1]);
            Assert.AreEqual(1,  testMatrix1[2, 2]);
        }
    }
}