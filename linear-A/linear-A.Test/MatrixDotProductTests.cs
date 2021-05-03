#region Copyright

// --------------------------------------------------------------------------------------------------------------------
// <copyright file="MatrixDotProductTests.cs" company="John-Michael Cummings">
//   John-Michael Cummings 2021
// </copyright>
// <summary>
//   The matrix dot product tests.
// </summary>
// --------------------------------------------------------------------------------------------------------------------

#endregion

using Matricks.Core;

namespace linear_A.Test
{
    using System.IO;
    using NUnit.Framework;

    /// <summary>
    ///     The matrix dot product tests.
    /// </summary>
    public class MatrixDotProductTests
    {
        /// <summary>
        ///     The hacker rank test case 3.
        /// </summary>
        [Test]
        public void HackerRankTestCase3()
        {
            const int rowCount = 2;
            const int columnCount = 2;
            var testMatrix1 = new IntegerMatrix(rowCount, columnCount)
            {
                [0, 0] = 1, 
                [0, 1] = 2, 
                [1, 0] = 2,
                [1, 1] = 3
            };
            var testMatrix2 = new IntegerMatrix(rowCount, columnCount)
            {
                [0, 0] = 4,
                [0, 1] = 5, 
                [1, 0] = 7, 
                [1, 1] = 8
            };

            var result = testMatrix1.DotProduct(testMatrix2);

            Assert.AreEqual(18, result[0, 0]);
            Assert.AreEqual(21, result[0, 1]);
            Assert.AreEqual(29, result[1, 0]);
            Assert.AreEqual(34, result[1, 1]);
        }

        /// <summary>
        ///     The hacker rank test case 4.
        /// </summary>
        [Test]
        public void HackerRankTestCase4()
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

            var result = testMatrix1.DotProduct(testMatrix2);

            Assert.AreEqual(30, result[0, 0]);
            Assert.AreEqual(36, result[0, 1]);
            Assert.AreEqual(45, result[0, 2]);
            Assert.AreEqual(45, result[1, 0]);
            Assert.AreEqual(54, result[1, 1]);
            Assert.AreEqual(67, result[1, 2]);
            Assert.AreEqual(15, result[2, 0]);
            Assert.AreEqual(18, result[2, 1]);
            Assert.AreEqual(22, result[2, 2]);
        }

        /// <summary>
        ///     The hacker rank test case 5.
        /// </summary>
        [Test]
        public void HackerRankTestCase5()
        {
            const int rowCount = 3;
            const int columnCount = 3;
            var testMatrix1 = new IntegerMatrix(rowCount, columnCount)
            {
                [0, 0] = 1,
                [0, 1] = 1,
                [0, 2] = 0,
                [1, 0] = 0,
                [1, 1] = 1,
                [1, 2] = 0,
                [2, 0] = 0,
                [2, 1] = 0,
                [2, 2] = 1
            };
            var testMatrix2 = new IntegerMatrix(rowCount, columnCount)
            {
                [0, 0] = 1,
                [0, 1] = 1,
                [0, 2] = 0,
                [1, 0] = 0,
                [1, 1] = 1,
                [1, 2] = 0,
                [2, 0] = 0,
                [2, 1] = 0,
                [2, 2] = 1
            };

            for (var i = 1; i < 100; i++)
            {
                testMatrix2 = testMatrix2.DotProduct(testMatrix1);
            }

            Assert.AreEqual(1,   testMatrix2[0, 0]);
            Assert.AreEqual(100, testMatrix2[0, 1]);
            Assert.AreEqual(0,   testMatrix2[0, 2]);
            Assert.AreEqual(0,   testMatrix2[1, 0]);
            Assert.AreEqual(1,   testMatrix2[1, 1]);
            Assert.AreEqual(0,   testMatrix2[1, 2]);
            Assert.AreEqual(0,   testMatrix2[2, 0]);
            Assert.AreEqual(0,   testMatrix2[2, 1]);
            Assert.AreEqual(1,   testMatrix2[2, 2]);
        }

        /// <summary>
        ///     The hacker rank test case 6.
        /// </summary>
        [Test]
        public void HackerRankTestCase6()
        {
            const int rowCount = 3;
            const int columnCount = 3;
            const int x = -2;
            const int y = 1;
            var matrixI = new IntegerMatrix(rowCount, true);
            var matrixA = new IntegerMatrix(rowCount, columnCount)
            {
                [0, 0] = 1,
                [0, 1] = 1,
                [0, 2] = 0,
                [1, 0] = 0,
                [1, 1] = 1,
                [1, 2] = 0,
                [2, 0] = 0,
                [2, 1] = 0,
                [2, 2] = 1
            };

            var matrixASquared = matrixA.DotProduct(matrixA);
            matrixA.Scale(x);
            matrixI.Scale(y);
            var firstSuccess = matrixASquared.TryAddMatrix(matrixA);
            var secondSuccess = matrixASquared.TryAddMatrix(matrixI);

            Assert.IsTrue(firstSuccess);
            Assert.IsTrue(secondSuccess);
            Assert.AreEqual(0, matrixASquared[0, 0]);
            Assert.AreEqual(0, matrixASquared[0, 1]);
            Assert.AreEqual(0, matrixASquared[0, 2]);
            Assert.AreEqual(0, matrixASquared[1, 0]);
            Assert.AreEqual(0, matrixASquared[1, 1]);
            Assert.AreEqual(0, matrixASquared[1, 2]);
            Assert.AreEqual(0, matrixASquared[2, 0]);
            Assert.AreEqual(0, matrixASquared[2, 1]);
            Assert.AreEqual(0, matrixASquared[2, 2]);
        }

        /// <summary>
        ///     The hacker rank test case 7.
        /// </summary>
        [Test]
        public void HackerRankTestCase7()
        {
            const int rowCount = 2;
            const int columnCount = 2;
            var testMatrix1 = new IntegerMatrix(rowCount, columnCount)
            {
                [0, 0] = -2, 
                [0, 1] = -9, 
                [1, 0] = 1, 
                [1, 1] = 4
            };

            var testMatrix2 = testMatrix1;
            for (var i = 1; i < 1000; i++)
            {
                testMatrix2 = testMatrix2.DotProduct(testMatrix1);
            }

            Assert.AreEqual(-2999, testMatrix2[0, 0]);
            Assert.AreEqual(-9000, testMatrix2[0, 1]);
            Assert.AreEqual(1000,  testMatrix2[1, 0]);
            Assert.AreEqual(3001,  testMatrix2[1, 1]);
        }

        /// <summary>
        ///     The hacker rank test case 8.
        /// </summary>
        [Test]
        public void HackerRankTestCase8()
        {
            const int rowCount = 5;
            var matrixA = new IntegerMatrix(rowCount, true)
            {
                [0, 0] = 3,
                [0, 1] = 0,
                [0, 2] = 0,
                [0, 3] = -4,
                [0, 4] = 0,
                [1, 0] = 0,
                [1, 1] = 2,
                [1, 2] = -1,
                [1, 3] = 0,
                [1, 4] = -1,
                [2, 0] = 0,
                [2, 1] = 0,
                [2, 2] = 0,
                [2, 3] = 1,
                [2, 4] = 0,
                [3, 0] = -2,
                [3, 1] = 0,
                [3, 2] = 5,
                [3, 3] = 0,
                [3, 4] = 3,
                [4, 0] = 4,
                [4, 1] = 0,
                [4, 2] = -3,
                [4, 3] = 6,
                [4, 4] = 2
            };

            var success = matrixA.TryGetDeterminant(out var result);

            Assert.IsTrue(success);
            Assert.AreEqual(-114, result);
        }

        /// <summary>
        ///     The hacker rank test case 9.
        /// </summary>
        [Test]
        public void HackerRankTestCase9()
        {
            const int rowCount = 5;
            var matrixA = new IntegerMatrix(rowCount, true)
            {
                [0, 0] = 3,
                [0, 1] = 0,
                [0, 2] = 0,
                [0, 3] = -4,
                [0, 4] = 0,
                [1, 0] = 0,
                [1, 1] = 2,
                [1, 2] = -1,
                [1, 3] = 0,
                [1, 4] = -1,
                [2, 0] = 0,
                [2, 1] = 0,
                [2, 2] = 0,
                [2, 3] = 1,
                [2, 4] = 0,
                [3, 0] = -2,
                [3, 1] = 0,
                [3, 2] = 5,
                [3, 3] = 0,
                [3, 4] = 3,
                [4, 0] = 4,
                [4, 1] = 0,
                [4, 2] = -3,
                [4, 3] = 6,
                [4, 4] = 2
            };

            var success = matrixA.TryGetDeterminant(out var result);
            
            Assert.IsTrue(success);
            Assert.AreEqual(-114, result);
        }
        
         /// <summary>
        ///     The hacker rank test case 3.
        /// </summary>
        [Test]
        public void HackerRankTestCase3Double()
        {
            const int rowCount = 2;
            const int columnCount = 2;
            var testMatrix1 = new DoubleMatrix(rowCount, columnCount)
            {
                [0, 0] = 1, 
                [0, 1] = 2, 
                [1, 0] = 2,
                [1, 1] = 3
            };
            var testMatrix2 = new DoubleMatrix(rowCount, columnCount)
            {
                [0, 0] = 4, 
                [0, 1] = 5, 
                [1, 0] = 7,
                [1, 1] = 8
            };

            var result = testMatrix1.DotProduct(testMatrix2);

            Assert.AreEqual(18, result[0, 0]);
            Assert.AreEqual(21, result[0, 1]);
            Assert.AreEqual(29, result[1, 0]);
            Assert.AreEqual(34, result[1, 1]);
        }

        /// <summary>
        ///     The hacker rank test case 4.
        /// </summary>
        [Test]
        public void HackerRankTestCase4Double()
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

            var result = testMatrix1.DotProduct(testMatrix2);

            Assert.AreEqual(30, result[0, 0]);
            Assert.AreEqual(36, result[0, 1]);
            Assert.AreEqual(45, result[0, 2]);
            Assert.AreEqual(45, result[1, 0]);
            Assert.AreEqual(54, result[1, 1]);
            Assert.AreEqual(67, result[1, 2]);
            Assert.AreEqual(15, result[2, 0]);
            Assert.AreEqual(18, result[2, 1]);
            Assert.AreEqual(22, result[2, 2]);
        }

        /// <summary>
        ///     The hacker rank test case 5.
        /// </summary>
        [Test]
        public void HackerRankTestCase5Double()
        {
            const int rowCount = 3;
            const int columnCount = 3;

            var testMatrix1 = new DoubleMatrix(rowCount, columnCount)
            {
                [0, 0] = 1,
                [0, 1] = 1,
                [0, 2] = 0,
                [1, 0] = 0,
                [1, 1] = 1,
                [1, 2] = 0,
                [2, 0] = 0,
                [2, 1] = 0,
                [2, 2] = 1
            };
            
            var testMatrix2 = new DoubleMatrix(rowCount, columnCount)
            {
                [0, 0] = 1,
                [0, 1] = 1,
                [0, 2] = 0,
                [1, 0] = 0,
                [1, 1] = 1,
                [1, 2] = 0,
                [2, 0] = 0,
                [2, 1] = 0,
                [2, 2] = 1
            };

            for (var i = 1; i < 100; i++)
            {
                testMatrix2 = testMatrix2.DotProduct(testMatrix1);
            }

            Assert.AreEqual(1,   testMatrix2[0, 0]);
            Assert.AreEqual(100, testMatrix2[0, 1]);
            Assert.AreEqual(0,   testMatrix2[0, 2]);
            Assert.AreEqual(0,   testMatrix2[1, 0]);
            Assert.AreEqual(1,   testMatrix2[1, 1]);
            Assert.AreEqual(0,   testMatrix2[1, 2]);
            Assert.AreEqual(0,   testMatrix2[2, 0]);
            Assert.AreEqual(0,   testMatrix2[2, 1]);
            Assert.AreEqual(1,   testMatrix2[2, 2]);
        }

        /// <summary>
        ///     The hacker rank test case 6.
        /// </summary>
        [Test]
        public void HackerRankTestCase6Double()
        {
            const int rowCount = 3;
            const int columnCount = 3;
            const int x = -2;
            const int y = 1;
            
            var matrixI = new DoubleMatrix(rowCount, true);
            var matrixA = new DoubleMatrix(rowCount, columnCount)
            {
                [0, 0] = 1,
                [0, 1] = 1,
                [0, 2] = 0,
                [1, 0] = 0,
                [1, 1] = 1,
                [1, 2] = 0,
                [2, 0] = 0,
                [2, 1] = 0,
                [2, 2] = 1
            };

            var matrixASquared = matrixA.DotProduct(matrixA);
            matrixA.Scale(x);
            matrixI.Scale(y);

            var success = matrixASquared.TryAddMatrix(matrixA);
            
            Assert.IsTrue(success);
         
            success = matrixASquared.TryAddMatrix(matrixI);
            Assert.IsTrue(success);
            Assert.AreEqual(0, matrixASquared[0, 0]);
            Assert.AreEqual(0, matrixASquared[0, 1]);
            Assert.AreEqual(0, matrixASquared[0, 2]);
            Assert.AreEqual(0, matrixASquared[1, 0]);
            Assert.AreEqual(0, matrixASquared[1, 1]);
            Assert.AreEqual(0, matrixASquared[1, 2]);
            Assert.AreEqual(0, matrixASquared[2, 0]);
            Assert.AreEqual(0, matrixASquared[2, 1]);
            Assert.AreEqual(0, matrixASquared[2, 2]);
        }

        /// <summary>
        ///     The hacker rank test case 7.
        /// </summary>
        [Test]
        public void HackerRankTestCase7Double()
        {
            const int rowCount = 2;
            const int columnCount = 2;
            
            var testMatrix1 = new DoubleMatrix(rowCount, columnCount)
            {
                [0, 0] = -2, 
                [0, 1] = -9, 
                [1, 0] = 1, 
                [1, 1] = 4
            };


            var testMatrix2 = testMatrix1;

            for (var i = 1; i < 1000; i++)
            {
                testMatrix2 = testMatrix2.DotProduct(testMatrix1);
            }

            Assert.AreEqual(-2999, testMatrix2[0, 0]);
            Assert.AreEqual(-9000, testMatrix2[0, 1]);
            Assert.AreEqual(1000,  testMatrix2[1, 0]);
            Assert.AreEqual(3001,  testMatrix2[1, 1]);
        }

        /// <summary>
        ///     The hacker rank test case 8.
        /// </summary>
        [Test]
        public void HackerRankTestCase8Double()
        {
            const int rowCount = 5;
            var matrixA = new DoubleMatrix(rowCount, true)
            {
                [0, 0] = 3,
                [0, 1] = 0,
                [0, 2] = 0,
                [0, 3] = -4,
                [0, 4] = 0,
                [1, 0] = 0,
                [1, 1] = 2,
                [1, 2] = -1,
                [1, 3] = 0,
                [1, 4] = -1,
                [2, 0] = 0,
                [2, 1] = 0,
                [2, 2] = 0,
                [2, 3] = 1,
                [2, 4] = 0,
                [3, 0] = -2,
                [3, 1] = 0,
                [3, 2] = 5,
                [3, 3] = 0,
                [3, 4] = 3,
                [4, 0] = 4,
                [4, 1] = 0,
                [4, 2] = -3,
                [4, 3] = 6,
                [4, 4] = 2
            };

            var success = matrixA.TryGetDeterminant(out var result);
           
            Assert.IsTrue(success);
            Assert.AreEqual(-114, result);
        }

        /// <summary>
        ///     The hacker rank test case 9.
        /// </summary>
        [Test]
        public void HackerRankTestCase9Double()
        {
            const int rowCount = 5;
            var matrixA = new DoubleMatrix(rowCount, true)
            {
                [0, 0] = 3,
                [0, 1] = 0,
                [0, 2] = 0,
                [0, 3] = -4,
                [0, 4] = 0,
                [1, 0] = 0,
                [1, 1] = 2,
                [1, 2] = -1,
                [1, 3] = 0,
                [1, 4] = -1,
                [2, 0] = 0,
                [2, 1] = 0,
                [2, 2] = 0,
                [2, 3] = 1,
                [2, 4] = 0,
                [3, 0] = -2,
                [3, 1] = 0,
                [3, 2] = 5,
                [3, 3] = 0,
                [3, 4] = 3,
                [4, 0] = 4,
                [4, 1] = 0,
                [4, 2] = -3,
                [4, 3] = 6,
                [4, 4] = 2
            };

            var success = matrixA.TryGetDeterminant(out var result);
            
            Assert.IsTrue(success);
            Assert.AreEqual(-114, result);
        }

        /// <summary>
        /// The write matrix to file method. Used for debugging.
        /// </summary>
        /// <param name="matrix">
        /// The matrix.
        /// </param>
        /// <param name="header">
        /// The header.
        /// </param>
        /// <param name="fileName">
        /// The file name.
        /// </param>
        private void WriteMatrixToFile(
            IntegerMatrix matrix,
            string header = "Matrix  :   ",
            string fileName = @"hr2.txt")
        {
            using (var file = new StreamWriter(fileName, true))
            {
                file.WriteLine("===============================");
                file.WriteLine(header);

                for (var i = 0; i < matrix.RowCount; i++)
                {
                    var nextLine = string.Empty;
                    for (var j = 0; j < matrix.ColumnCount; j++)
                    {
                        nextLine += matrix[i, j] + ", ";
                    }

                    file.WriteLine(nextLine);
                }

                file.WriteLine("===============================");
                file.WriteLine(string.Empty);
            }
        }
    }
}