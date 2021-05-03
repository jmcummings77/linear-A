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

    /// <summary>
    ///     The matrix addition tests.
    /// </summary>
    public class MatrixAdditionTests
    {
        /// <summary>
        ///     The TryAddZeros_1X1MatricesResultSuccessAndIsZero test.
        /// </summary>
        [Test]
        public void TryAddZeros_1X1MatricesResultSuccessAndIsZero()
        {
            var testMatrix1 = new IntegerMatrix(1, false);
            var testMatrix2 = new IntegerMatrix(1, false);
            
            var success = testMatrix1.TryAddMatrix(testMatrix2);
            
            Assert.IsTrue(success);
            Assert.AreEqual(testMatrix1[0, 0], 0);
        }

        /// <summary>
        ///     The TryAddOnes_1X1MatricesResultSuccessAndIsTwo test.
        /// </summary>
        [Test]
        public void TryAddOnes_1X1MatricesResultSuccessAndIsTwo()
        {
            var testMatrix1 = new IntegerMatrix(1, true);
            var testMatrix2 = new IntegerMatrix(1, true);
            
            var success = testMatrix1.TryAddMatrix(testMatrix2);
            
            Assert.IsTrue(success);
            Assert.AreEqual(testMatrix1[0, 0], 2);
        }

        /// <summary>
        ///     The TryAddDifferentDimensionMatricesFails test.
        /// </summary>
        [Test]
        public void TryAddDifferentDimensionMatricesFails()
        {
            var testMatrix1 = new IntegerMatrix(1, false);
            var testMatrix2 = new IntegerMatrix(2, true);
            
            var success = testMatrix1.TryAddMatrix(testMatrix2);
            
            Assert.IsFalse(success);
            Assert.AreEqual(testMatrix1[0, 0], 0);
        }

        /// <summary>
        ///     The TryAddNegative_3X3ToEquivalent_3X3ResultSuccessAndAllAreZero test.
        /// </summary>
        [Test]
        public void TryAddNegative_3X3ToEquivalent_3X3ResultSuccessAndAllAreZero()
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

            var success = testMatrix1.TryAddMatrix(testMatrix2);
            
            Assert.IsTrue(success);
            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    Assert.AreEqual(testMatrix1[i, j], 0);
                }
            }
        }

        /// <summary>
        ///     The TryAddNegative_3X3ToEquivalent_3X3WithOperatorResultAllAreZero test.
        /// </summary>
        [Test]
        public void TryAddNegative_3X3ToEquivalent_3X3WithOperatorResultAllAreZero()
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

            var success = testMatrix1.TryAddMatrix(testMatrix2);
            Assert.IsTrue(success);
            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    Assert.AreEqual(testMatrix1[i, j], 0);
                }
            }
        }

        /// <summary>
        ///     The try_add_matrix_to_negative_equivalent_result_success_and_all_are_zero test.
        /// </summary>
        [Test]
        public void add_matrix_to_negative_equivalent_with_operator_result_all_are_zero()
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
        ///     The try_add_matrix_to_negative_equivalent_result_success_and_all_are_zero test.
        /// </summary>
        [Test]
        public void try_add_matrix_to_negative_equivalent_result_success_and_all_are_zero()
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

                    var success = testMatrix1.TryAddMatrix(testMatrix2);
                    Assert.IsTrue(success);
                    for (var i = 0; i < rowCount; i++)
                    {
                        for (var j = 0; j < columnCount; j++)
                        {
                            Assert.AreEqual(testMatrix1[i, j], 0);
                        }
                    }
                }
            }
        }

        /// <summary>
        ///     The hacker rank test case #1.
        /// </summary>
        [Test]
        public void HackerRankTestCase1()
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
                [2, 2] = 7
            };

            testMatrix1.TryAddMatrix(testMatrix2);

            Assert.AreEqual(5,  testMatrix1[0, 0]);
            Assert.AreEqual(7,  testMatrix1[0, 1]);
            Assert.AreEqual(9,  testMatrix1[0, 2]);
            Assert.AreEqual(9,  testMatrix1[1, 0]);
            Assert.AreEqual(11, testMatrix1[1, 1]);
            Assert.AreEqual(13, testMatrix1[1, 2]);
            Assert.AreEqual(5,  testMatrix1[2, 0]);
            Assert.AreEqual(6,  testMatrix1[2, 1]);
            Assert.AreEqual(8,  testMatrix1[2, 2]);
        }

        /// <summary>
        ///     The hacker rank test case #1 with addition operator.
        /// </summary>
        [Test]
        public void HackerRankTestCase1WithOperator()
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
                [2, 2] = 7
            };

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
        ///     The TryAddZeros_1X1MatricesResultSuccessAndIsZero test.
        /// </summary>
        [Test]
        public void TryAddZeros_1X1DoubleMatricesResultSuccessAndIsZero()
        {
            var testMatrix1 = new DoubleMatrix(1, false);
            var testMatrix2 = new DoubleMatrix(1, false);
            
            var success = testMatrix1.TryAddMatrix(testMatrix2);
            
            Assert.IsTrue(success);
            Assert.AreEqual(testMatrix1[0, 0], 0);
        }

        /// <summary>
        ///     The TryAddOnes_1X1MatricesResultSuccessAndIsTwo test.
        /// </summary>
        [Test]
        public void TryAddOnes_1X1DoubleMatricesResultSuccessAndIsTwo()
        {
            var testMatrix1 = new DoubleMatrix(1, true);
            var testMatrix2 = new DoubleMatrix(1, true);
            
            var success = testMatrix1.TryAddMatrix(testMatrix2);
            
            Assert.IsTrue(success);
            Assert.AreEqual(testMatrix1[0, 0], 2);
        }

        /// <summary>
        ///     The TryAddDifferentDimensionMatricesFails test.
        /// </summary>
        [Test]
        public void TryAddDifferentDimensionDoubleMatricesFails()
        {
            var testMatrix1 = new DoubleMatrix(1, false);
            var testMatrix2 = new DoubleMatrix(2, true);
            
            var success = testMatrix1.TryAddMatrix(testMatrix2);
            
            Assert.IsFalse(success);
            Assert.AreEqual(testMatrix1[0, 0], 0);
        }

        /// <summary>
        ///     The TryAddNegative_3X3ToEquivalent_3X3ResultSuccessAndAllAreZero test.
        /// </summary>
        [Test]
        public void TryAddNegative_3X3DoubleToEquivalent_3X3ResultSuccessAndAllAreZero()
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

            var success = testMatrix1.TryAddMatrix(testMatrix2);
            Assert.IsTrue(success);
            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    Assert.AreEqual(testMatrix1[i, j], 0);
                }
            }
        }

        /// <summary>
        ///     The TryAddNegative_3X3ToEquivalent_3X3WithOperatorResultAllAreZero test.
        /// </summary>
        [Test]
        public void TryAddNegative_3X3DoubleToEquivalent_3X3WithOperatorResultAllAreZero()
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

            var success = testMatrix1.TryAddMatrix(testMatrix2);
            Assert.IsTrue(success);
            for (var i = 0; i < 3; i++)
            {
                for (var j = 0; j < 3; j++)
                {
                    Assert.AreEqual(testMatrix1[i, j], 0);
                }
            }
        }

        /// <summary>
        ///     The try_add_matrix_to_negative_equivalent_result_success_and_all_are_zero test.
        /// </summary>
        [Test]
        public void add_double_matrix_to_negative_equivalent_with_operator_result_all_are_zero()
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
        ///     The try_add_matrix_to_negative_equivalent_result_success_and_all_are_zero test.
        /// </summary>
        [Test]
        public void try_add_double_matrix_to_negative_equivalent_result_success_and_all_are_zero()
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

                    var success = testMatrix1.TryAddMatrix(testMatrix2);
                    Assert.IsTrue(success);
                    for (var i = 0; i < rowCount; i++)
                    {
                        for (var j = 0; j < columnCount; j++)
                        {
                            Assert.AreEqual(testMatrix1[i, j], 0);
                        }
                    }
                }
            }
        }

        /// <summary>
        ///     The hacker rank test case #1.
        /// </summary>
        [Test]
        public void HackerRankTestCase1Double()
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
                [2, 2] = 7
            };

            testMatrix1.TryAddMatrix(testMatrix2);

            Assert.AreEqual(5,  testMatrix1[0, 0]);
            Assert.AreEqual(7,  testMatrix1[0, 1]);
            Assert.AreEqual(9,  testMatrix1[0, 2]);
            Assert.AreEqual(9,  testMatrix1[1, 0]);
            Assert.AreEqual(11, testMatrix1[1, 1]);
            Assert.AreEqual(13, testMatrix1[1, 2]);
            Assert.AreEqual(5,  testMatrix1[2, 0]);
            Assert.AreEqual(6,  testMatrix1[2, 1]);
            Assert.AreEqual(8,  testMatrix1[2, 2]);
        }

        /// <summary>
        ///     The hacker rank test case #1 with addition operator.
        /// </summary>
        [Test]
        public void HackerRankTestCase1WithOperatorDouble()
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
                [2, 2] = 7
            };

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