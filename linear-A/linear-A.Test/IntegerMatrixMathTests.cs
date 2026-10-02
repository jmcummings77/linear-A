namespace linear_A.Test
{
    using System;
    using System.Collections.Generic;
    using NUnit.Framework;

    public class IntegerMatrixMathTests
    {
        private static IEnumerable<TestCaseData> DeterminantCases()
        {
            yield return new TestCaseData(new[,] { { -7 } }, -7).SetName("Integer determinant: one element");
            yield return new TestCaseData(new[,] { { 1, 2 }, { 3, 4 } }, -2).SetName("Integer determinant: 2x2");
            yield return new TestCaseData(new[,] { { 6, 1, 1 }, { 4, -2, 5 }, { 2, 8, 7 } }, -306)
                .SetName("Integer determinant: dense 3x3");
            yield return new TestCaseData(new[,] { { 3, 2, 0, 1 }, { 4, 0, 1, 2 }, { 3, 0, 2, 1 }, { 9, 2, 3, 1 } }, 24)
                .SetName("Integer determinant: dense 4x4");
            yield return new TestCaseData(new[,] { { 0, 2, 1 }, { 3, 1, 4 }, { 2, 5, 6 } }, -7)
                .SetName("Integer determinant: zero leading entry");
            yield return new TestCaseData(new[,] { { 1, 2, 3 }, { 2, 4, 6 }, { 4, 5, 6 } }, 0)
                .SetName("Integer determinant: linearly dependent rows");
            yield return new TestCaseData(new[,] { { 2, 7, -1 }, { 0, -3, 4 }, { 0, 0, 5 } }, -30)
                .SetName("Integer determinant: triangular diagonal product");
        }

        [TestCaseSource(nameof(DeterminantCases))]
        public void DeterminantMatchesKnownValueAndPreservesEntries(int[,] entries, int expected)
        {
            var matrix = Create(entries);

            Assert.That(matrix.GetDeterminant(), Is.EqualTo(expected));
            Assert.That(matrix.TryGetDeterminant(out var determinant), Is.True);
            Assert.That(determinant, Is.EqualTo(expected));
            AssertMatrix(matrix, entries);
        }

        [TestCase(0, false)]
        [TestCase(1, true)]
        [TestCase(-7, true)]
        public void OneElementInvertibilityDependsOnWhetherEntryIsZero(int value, bool expected)
        {
            Assert.That(Create(new[,] { { value } }).IsInvertible(), Is.EqualTo(expected));
        }

        [Test]
        public void InvertibilityDistinguishesIndependentAndDependentRows()
        {
            Assert.That(Create(new[,] { { 1, 2 }, { 3, 4 } }).IsInvertible(), Is.True);
            Assert.That(Create(new[,] { { 1, 2 }, { 2, 4 } }).IsInvertible(), Is.False);
            Assert.That(new IntegerMatrix(3, true).IsInvertible(), Is.True);
            Assert.That(new IntegerMatrix(2, 3).IsInvertible(), Is.False);
        }

        [Test]
        public void RectangularTransposeSwapsDimensionsAndEntriesAndRoundTrips()
        {
            var entries = new[,] { { 1, -2, 3 }, { 4, 5, -6 } };
            var matrix = Create(entries);

            var transposed = matrix.Transpose();

            AssertMatrix(transposed, new[,] { { 1, 4 }, { -2, 5 }, { 3, -6 } });
            AssertMatrix(transposed.Transpose(), entries);
            transposed[0, 0] = 99;
            AssertMatrix(matrix, entries);
        }

        [Test]
        public void RectangularMultiplicationMatchesHandComputedProduct()
        {
            var leftEntries = new[,] { { 1, 2, 3 }, { -1, 0, 2 } };
            var rightEntries = new[,] { { 4, -2 }, { 0, 5 }, { 3, 1 } };
            var left = Create(leftEntries);
            var right = Create(rightEntries);
            var expected = new[,] { { 13, 11 }, { 2, 4 } };

            AssertMatrix(left.DotProduct(right), expected);
            AssertMatrix(left * right, expected);
            AssertMatrix(left, leftEntries);
            AssertMatrix(right, rightEntries);
        }

        [Test]
        public void MultiplicationPreservesIdentityAndReversesOrderUnderTranspose()
        {
            var entries = new[,] { { 1, 2, 3 }, { -1, 0, 2 } };
            var left = Create(entries);
            var right = Create(new[,] { { 4, -2 }, { 0, 5 }, { 3, 1 } });

            AssertMatrix(new IntegerMatrix(2, true) * left, entries);
            AssertMatrix(left * new IntegerMatrix(3, true), entries);
            AssertMatrix((left * right).Transpose(), right.Transpose() * left.Transpose());
            Assert.That((left * right).GetTrace(), Is.EqualTo((right * left).GetTrace()));
        }

        [Test]
        public void DeterminantObeysTransposeProductAndRowSwapIdentities()
        {
            var left = Create(new[,] { { 1, 2 }, { 3, 4 } });
            var right = Create(new[,] { { 2, 0 }, { 1, 3 } });

            Assert.That(left.Transpose().GetDeterminant(), Is.EqualTo(-2));
            Assert.That((left * right).GetDeterminant(), Is.EqualTo(-12));
            Assert.That(Create(new[,] { { 3, 4 }, { 1, 2 } }).GetDeterminant(), Is.EqualTo(2));
        }

        [Test]
        public void TriangularClassificationChecksTheCorrectSideOfTheDiagonal()
        {
            var upper = Create(new[,] { { 2, 7, -1 }, { 0, -3, 4 }, { 0, 0, 5 } });
            var lower = upper.Transpose();
            var dense = Create(new[,] { { 1, 2 }, { 3, 4 } });

            Assert.That(upper.IsUpperTriangular(), Is.True);
            Assert.That(upper.IsLowerTriangular(), Is.False);
            Assert.That(upper.IsTriangular(), Is.True);
            Assert.That(lower.IsLowerTriangular(), Is.True);
            Assert.That(lower.IsUpperTriangular(), Is.False);
            Assert.That(lower.IsTriangular(), Is.True);
            Assert.That(dense.IsTriangular(), Is.False);
            Assert.That(new IntegerMatrix(3, true).IsUpperTriangular(), Is.True);
            Assert.That(new IntegerMatrix(3, true).IsLowerTriangular(), Is.True);
        }

        [Test]
        public void TraceAndDeterminantReflectSubsequentMutations()
        {
            var matrix = Create(new[,] { { 1, 2 }, { 3, 4 } });
            Assert.That(matrix.GetTrace(), Is.EqualTo(5));
            Assert.That(matrix.GetDeterminant(), Is.EqualTo(-2));

            matrix[0, 0] = 5;
            Assert.That(matrix.GetDeterminant(), Is.EqualTo(14));
            Assert.That(matrix.TryAddMatrix(new IntegerMatrix(2, true)), Is.True);
            Assert.That(matrix.GetTrace(), Is.EqualTo(11));
            Assert.That(matrix.GetDeterminant(), Is.EqualTo(24));
            Assert.That(matrix.TrySubtractMatrix(new IntegerMatrix(2, true)), Is.True);
            matrix.Scale(-2);

            Assert.That(matrix.TryGetTrace(out var trace), Is.True);
            Assert.That(trace, Is.EqualTo(-18));
            Assert.That(matrix.GetDeterminant(), Is.EqualTo(56));
        }

        [Test]
        public void CopyRowsColumnsAndArithmeticResultsDoNotAliasOriginal()
        {
            var entries = new[,] { { 1, 2, 3 }, { 4, 5, 6 } };
            var original = Create(entries);
            var copy = new IntegerMatrix(original);
            var row = original.GetRow(1);
            var column = original.GetColumn(2);
            Assert.That(row, Is.EqualTo(new[] { 4, 5, 6 }));
            Assert.That(column, Is.EqualTo(new[] { 3, 6 }));

            copy[0, 0] = 100;
            row[0] = 100;
            column[0] = 100;
            var sum = original + Create(entries);
            var difference = original - Create(entries);
            sum[0, 0] = 100;
            difference[0, 0] = 100;

            AssertMatrix(original, entries);
        }

        [Test]
        public void NonSquareCharacteristicsRejectUnsupportedShapes()
        {
            var matrix = new IntegerMatrix(2, 3);

            Assert.That(matrix.IsSquare(), Is.False);
            Assert.That(matrix.IsTriangular(), Is.False);
            Assert.That(matrix.IsUpperTriangular(), Is.False);
            Assert.That(matrix.IsLowerTriangular(), Is.False);
            Assert.That(() => matrix.GetTrace(), Throws.TypeOf<NotSquareMatrixException>());
            Assert.That(() => matrix.GetDeterminant(), Throws.TypeOf<NotSquareMatrixException>());
            Assert.That(matrix.TryGetTrace(out var trace), Is.False);
            Assert.That(trace, Is.Zero);
            Assert.That(matrix.TryGetDeterminant(out var determinant), Is.False);
            Assert.That(determinant, Is.Zero);
        }

        [Test]
        public void ShapeMismatchDoesNotPartiallyMutateOperands()
        {
            var entries = new[,] { { 1, 2, 3 }, { 4, 5, 6 } };
            var matrix = Create(entries);
            var mismatch = new IntegerMatrix(3, 2);

            Assert.That(matrix.TryAddMatrix(mismatch), Is.False);
            Assert.That(matrix.TrySubtractMatrix(mismatch), Is.False);
            Assert.That(() => matrix + mismatch, Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => matrix - mismatch, Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => matrix * new IntegerMatrix(2, 2), Throws.TypeOf<ArgumentOutOfRangeException>());
            AssertMatrix(matrix, entries);
        }

        [Test]
        public void DeterminantUsesExactIntermediateArithmetic()
        {
            var maximum = int.MaxValue;
            var matrix = Create(new[,] { { maximum, maximum - 1 }, { maximum - 1, maximum - 2 } });

            // n(n - 2) - (n - 1)^2 = -1, despite both products exceeding Int32.
            Assert.That(matrix.GetDeterminant(), Is.EqualTo(-1));
            Assert.That(matrix.IsInvertible(), Is.True);
        }

        [Test]
        public void UnrepresentableDeterminantIsReportedWithoutLosingInvertibility()
        {
            var matrix = Create(new[,] { { 50000, 0 }, { 0, 50000 } });

            Assert.That(() => matrix.GetDeterminant(), Throws.TypeOf<OverflowException>());
            Assert.That(matrix.TryGetDeterminant(out var determinant), Is.False);
            Assert.That(determinant, Is.Zero);
            Assert.That(matrix.IsInvertible(), Is.True);
        }

        [Test]
        public void DeterminantsMatchIndependentCofactorOracleForSmallRandomMatrices()
        {
            var random = new Random(20261002);
            for (var size = 1; size <= 5; size++)
                for (var sample = 0; sample < 24; sample++)
                {
                    var entries = new int[size, size];
                    for (var row = 0; row < size; row++)
                        for (var column = 0; column < size; column++)
                            entries[row, column] = random.Next(-4, 5);

                    Assert.That(Create(entries).GetDeterminant(), Is.EqualTo(CofactorDeterminant(entries)),
                        $"Size {size}, sample {sample}");
                }
        }

        [Test]
        public void EmptyMatrixObeysEmptyProductAndSumConventions()
        {
            var matrix = new IntegerMatrix();

            Assert.That(matrix.GetDeterminant(), Is.EqualTo(1));
            Assert.That(matrix.GetTrace(), Is.Zero);
            Assert.That(matrix.IsInvertible(), Is.True);
            Assert.That(matrix.Transpose().RowCount, Is.Zero);
            Assert.That(matrix.Transpose().ColumnCount, Is.Zero);
            Assert.That((matrix * matrix).RowCount, Is.Zero);
            Assert.That((matrix * matrix).ColumnCount, Is.Zero);
        }

        [Test]
        public void NullOperandsAreRejectedAndTryOperationsPreserveTheReceiver()
        {
            var entries = new[,] { { 1, 2 }, { 3, 4 } };
            var matrix = Create(entries);

            Assert.That(matrix.TryAddMatrix(null!), Is.False);
            Assert.That(matrix.TrySubtractMatrix(null!), Is.False);
            Assert.That(() => matrix.DotProduct(null!), Throws.TypeOf<ArgumentNullException>());
            Assert.That(() => matrix * null!, Throws.TypeOf<ArgumentNullException>());
            Assert.That(() => null! * matrix, Throws.TypeOf<ArgumentNullException>());
            AssertMatrix(matrix, entries);
        }

        private static long CofactorDeterminant(int[,] entries)
        {
            var size = entries.GetLength(0);
            if (size == 1) return entries[0, 0];

            long determinant = 0;
            for (var excludedColumn = 0; excludedColumn < size; excludedColumn++)
            {
                var minor = new int[size - 1, size - 1];
                for (var row = 1; row < size; row++)
                {
                    var minorColumn = 0;
                    for (var column = 0; column < size; column++)
                        if (column != excludedColumn)
                            minor[row - 1, minorColumn++] = entries[row, column];
                }

                var sign = excludedColumn % 2 == 0 ? 1 : -1;
                determinant += sign * entries[0, excludedColumn] * CofactorDeterminant(minor);
            }

            return determinant;
        }

        private static IntegerMatrix Create(int[,] entries)
        {
            var matrix = new IntegerMatrix(entries.GetLength(0), entries.GetLength(1));
            for (var row = 0; row < matrix.RowCount; row++)
                for (var column = 0; column < matrix.ColumnCount; column++)
                    matrix[row, column] = entries[row, column];
            return matrix;
        }

        private static void AssertMatrix(IntegerMatrix actual, int[,] expected)
        {
            AssertMatrix(actual, Create(expected));
        }

        private static void AssertMatrix(IntegerMatrix actual, IntegerMatrix expected)
        {
            Assert.That(actual.RowCount, Is.EqualTo(expected.RowCount));
            Assert.That(actual.ColumnCount, Is.EqualTo(expected.ColumnCount));
            for (var row = 0; row < expected.RowCount; row++)
                for (var column = 0; column < expected.ColumnCount; column++)
                    Assert.That(actual[row, column], Is.EqualTo(expected[row, column]), $"Entry [{row}, {column}]");
        }
    }
}
