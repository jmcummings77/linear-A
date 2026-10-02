namespace linear_A.Test
{
    using System;
    using System.Collections.Generic;
    using NUnit.Framework;

    public class DoubleMatrixMathTests
    {
        private static IEnumerable<TestCaseData> DeterminantCases()
        {
            yield return new TestCaseData(new[,] { { -0.25 } }, -0.25).SetName("Double determinant: one element");
            yield return new TestCaseData(new[,] { { 0.5, 1.25 }, { -2.0, 0.75 } }, 2.875)
                .SetName("Double determinant: fractional 2x2");
            yield return new TestCaseData(new[,] { { 3.0, 0.5, 0.5 }, { 2.0, -1.0, 2.5 }, { 1.0, 4.0, 3.5 } }, -38.25)
                .SetName("Double determinant: fractional dense 3x3");
            yield return new TestCaseData(new[,] { { 1.5, 1.0, 0.0, 0.5 }, { 2.0, 0.0, 0.5, 1.0 },
                { 1.5, 0.0, 1.0, 0.5 }, { 4.5, 1.0, 1.5, 0.5 } }, 1.5)
                .SetName("Double determinant: fractional dense 4x4");
            yield return new TestCaseData(new[,] { { 0.0, 2.0, 1.0 }, { 3.0, 1.0, 4.0 }, { 2.0, 5.0, 6.0 } }, -7.0)
                .SetName("Double determinant: zero leading entry");
            yield return new TestCaseData(new[,] { { 0.5, 1.0, 1.5 }, { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 } }, 0.0)
                .SetName("Double determinant: linearly dependent rows");
            yield return new TestCaseData(new[,] { { 0.5, 7.0, -1.0 }, { 0.0, -1.5, 4.0 }, { 0.0, 0.0, 2.5 } }, -1.875)
                .SetName("Double determinant: triangular diagonal product");
        }

        [TestCaseSource(nameof(DeterminantCases))]
        public void DeterminantMatchesKnownValueAndPreservesEntries(double[,] entries, double expected)
        {
            var matrix = Create(entries);

            Assert.That(matrix.GetDeterminant(), Is.EqualTo(expected).Within(1e-12));
            Assert.That(matrix.TryGetDeterminant(out var determinant), Is.True);
            Assert.That(determinant, Is.EqualTo(expected).Within(1e-12));
            AssertMatrix(matrix, entries);
        }

        [TestCase(0.0, false)]
        [TestCase(1.0, true)]
        [TestCase(-0.25, true)]
        public void OneElementInvertibilityDependsOnWhetherEntryIsZero(double value, bool expected)
        {
            Assert.That(Create(new[,] { { value } }).IsInvertible(), Is.EqualTo(expected));
        }

        [Test]
        public void InvertibilityDistinguishesIndependentAndDependentRows()
        {
            Assert.That(Create(new[,] { { 0.5, 1.0 }, { 1.5, 2.0 } }).IsInvertible(), Is.True);
            Assert.That(Create(new[,] { { 0.5, 1.0 }, { 1.0, 2.0 } }).IsInvertible(), Is.False);
            Assert.That(new DoubleMatrix(3, true).IsInvertible(), Is.True);
            Assert.That(new DoubleMatrix(2, 3).IsInvertible(), Is.False);
        }

        [Test]
        public void RectangularTransposeSwapsDimensionsAndEntriesAndRoundTrips()
        {
            var entries = new[,] { { 0.5, -1.25, 1.5 }, { 2.0, 2.5, -3.0 } };
            var matrix = Create(entries);

            var transposed = matrix.Transpose();

            AssertMatrix(transposed, new[,] { { 0.5, 2.0 }, { -1.25, 2.5 }, { 1.5, -3.0 } });
            AssertMatrix(transposed.Transpose(), entries);
            transposed[0, 0] = 99;
            AssertMatrix(matrix, entries);
        }

        [Test]
        public void FractionalRectangularMultiplicationMatchesHandComputedProduct()
        {
            var leftEntries = new[,] { { 0.5, 1.0, 1.5 }, { -0.5, 0.0, 1.0 } };
            var rightEntries = new[,] { { 2.0, -1.0 }, { 0.0, 2.5 }, { 1.5, 0.5 } };
            var left = Create(leftEntries);
            var right = Create(rightEntries);
            var expected = new[,] { { 3.25, 2.75 }, { 0.5, 1.0 } };

            AssertMatrix(left.DotProduct(right), expected);
            AssertMatrix(left * right, expected);
            AssertMatrix(left, leftEntries);
            AssertMatrix(right, rightEntries);
        }

        [Test]
        public void MultiplicationPreservesIdentityAndReversesOrderUnderTranspose()
        {
            var entries = new[,] { { 0.5, 1.0, 1.5 }, { -0.5, 0.0, 1.0 } };
            var left = Create(entries);
            var right = Create(new[,] { { 2.0, -1.0 }, { 0.0, 2.5 }, { 1.5, 0.5 } });

            AssertMatrix(new DoubleMatrix(2, true) * left, entries);
            AssertMatrix(left * new DoubleMatrix(3, true), entries);
            AssertMatrix((left * right).Transpose(), right.Transpose() * left.Transpose());
            Assert.That((left * right).GetTrace(), Is.EqualTo((right * left).GetTrace()).Within(1e-12));
        }

        [Test]
        public void DeterminantObeysTransposeProductAndRowSwapIdentities()
        {
            var left = Create(new[,] { { 0.5, 1.0 }, { 1.5, 2.0 } });
            var right = Create(new[,] { { 2.0, 0.0 }, { 1.0, 3.0 } });

            Assert.That(left.Transpose().GetDeterminant(), Is.EqualTo(-0.5).Within(1e-12));
            Assert.That((left * right).GetDeterminant(), Is.EqualTo(-3).Within(1e-12));
            Assert.That(Create(new[,] { { 1.5, 2.0 }, { 0.5, 1.0 } }).GetDeterminant(), Is.EqualTo(0.5).Within(1e-12));
        }

        [Test]
        public void TriangularClassificationChecksTheCorrectSideOfTheDiagonal()
        {
            var upper = Create(new[,] { { 0.5, 7.0, -1.0 }, { 0.0, -1.5, 4.0 }, { 0.0, 0.0, 2.5 } });
            var lower = upper.Transpose();
            var dense = Create(new[,] { { 0.5, 1.0 }, { 1.5, 2.0 } });

            Assert.That(upper.IsUpperTriangular(), Is.True);
            Assert.That(upper.IsLowerTriangular(), Is.False);
            Assert.That(upper.IsTriangular(), Is.True);
            Assert.That(lower.IsLowerTriangular(), Is.True);
            Assert.That(lower.IsUpperTriangular(), Is.False);
            Assert.That(lower.IsTriangular(), Is.True);
            Assert.That(dense.IsTriangular(), Is.False);
            Assert.That(new DoubleMatrix(3, true).IsUpperTriangular(), Is.True);
            Assert.That(new DoubleMatrix(3, true).IsLowerTriangular(), Is.True);
        }

        [Test]
        public void FractionalArithmeticAndCharacteristicsReflectSubsequentMutations()
        {
            var matrix = Create(new[,] { { 0.5, 1.0 }, { 1.5, 2.0 } });
            Assert.That(matrix.GetTrace(), Is.EqualTo(2.5));
            Assert.That(matrix.GetDeterminant(), Is.EqualTo(-0.5).Within(1e-12));

            matrix[0, 0] = 2.5;
            Assert.That(matrix.GetDeterminant(), Is.EqualTo(3.5).Within(1e-12));
            Assert.That(matrix.TryAddMatrix(new DoubleMatrix(2, true)), Is.True);
            Assert.That(matrix.GetTrace(), Is.EqualTo(6.5));
            Assert.That(matrix.GetDeterminant(), Is.EqualTo(9).Within(1e-12));
            Assert.That(matrix.TrySubtractMatrix(new DoubleMatrix(2, true)), Is.True);
            matrix.Scale(-0.5);

            Assert.That(matrix.TryGetTrace(out var trace), Is.True);
            Assert.That(trace, Is.EqualTo(-2.25));
            Assert.That(matrix.GetDeterminant(), Is.EqualTo(0.875).Within(1e-12));
        }

        [Test]
        public void CopyRowsColumnsAndArithmeticResultsDoNotAliasOriginal()
        {
            var entries = new[,] { { 0.5, 1.0, 1.5 }, { 2.0, 2.5, 3.0 } };
            var original = Create(entries);
            var copy = new DoubleMatrix(original);
            var row = original.GetRow(1);
            var column = original.GetColumn(2);
            Assert.That(row, Is.EqualTo(new[] { 2.0, 2.5, 3.0 }));
            Assert.That(column, Is.EqualTo(new[] { 1.5, 3.0 }));

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
            var matrix = new DoubleMatrix(2, 3);

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
            var entries = new[,] { { 0.5, 1.0, 1.5 }, { 2.0, 2.5, 3.0 } };
            var matrix = Create(entries);
            var mismatch = new DoubleMatrix(3, 2);

            Assert.That(matrix.TryAddMatrix(mismatch), Is.False);
            Assert.That(matrix.TrySubtractMatrix(mismatch), Is.False);
            Assert.That(() => matrix + mismatch, Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => matrix - mismatch, Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => matrix * new DoubleMatrix(2, 2), Throws.TypeOf<ArgumentOutOfRangeException>());
            AssertMatrix(matrix, entries);
        }

        [Test]
        public void InvertibilityDoesNotDependOnDeterminantUnderflowOrOverflow()
        {
            var tiny = Create(new[,] { { 1e-200, 0.0 }, { 0.0, 1e-200 } });
            var huge = Create(new[,] { { 1e200, 0.0 }, { 0.0, 1e200 } });

            Assert.That(tiny.GetDeterminant(), Is.Zero);
            Assert.That(tiny.IsInvertible(), Is.True);
            Assert.That(huge.GetDeterminant(), Is.EqualTo(double.PositiveInfinity));
            Assert.That(huge.IsInvertible(), Is.True);
        }

        [Test]
        public void MixedScaleDiagonalDoesNotLoseItsRepresentableDeterminant()
        {
            var matrix = Create(new[,] { { 1e-200, 0.0, 0.0, 0.0 }, { 0.0, 1e-200, 0.0, 0.0 },
                { 0.0, 0.0, 1e200, 0.0 }, { 0.0, 0.0, 0.0, 1e200 } });

            Assert.That(matrix.GetDeterminant(), Is.EqualTo(1.0).Within(1e-12));
            Assert.That(matrix.IsInvertible(), Is.True);
        }

        [Test]
        public void ExtremeRowScalesPreserveTheOffDiagonalContribution()
        {
            var matrix = Create(new[,] { { 1e308, 1e308 }, { 1e-308, 2e-308 } });

            // ad - bc is approximately 2 - 1; eliminating with an underflowed factor loses bc.
            Assert.That(matrix.GetDeterminant(), Is.EqualTo(1.0).Within(1e-12));
            Assert.That(matrix.IsInvertible(), Is.True);
            // Row scaling must also preserve small entries when the matrix is transposed.
            Assert.That(matrix.Transpose().GetDeterminant(), Is.EqualTo(1.0).Within(1e-12));
        }

        [Test]
        public void SmallestNonzeroDoubleIsNotTreatedAsZero()
        {
            var scalar = Create(new[,] { { double.Epsilon } });
            var matrix = Create(new[,] { { 1.0, double.Epsilon }, { double.Epsilon, 1.0 } });

            Assert.That(scalar.GetDeterminant(), Is.EqualTo(double.Epsilon));
            Assert.That(scalar.IsInvertible(), Is.True);
            Assert.That(matrix.IsUpperTriangular(), Is.False);
            Assert.That(matrix.IsLowerTriangular(), Is.False);
            Assert.That(matrix.IsTriangular(), Is.False);
        }

        [TestCase(double.NaN)]
        [TestCase(double.PositiveInfinity)]
        [TestCase(double.NegativeInfinity)]
        public void NonfiniteEntriesAreNotReportedAsAnInvertibleFiniteMatrix(double value)
        {
            var matrix = Create(new[,] { { 1.0, value }, { 0.0, 1.0 } });

            Assert.That(matrix.GetDeterminant(), Is.NaN);
            Assert.That(matrix.IsInvertible(), Is.False);
            Assert.That(matrix.IsLowerTriangular(), Is.False);
        }

        [Test]
        public void EmptyMatrixObeysEmptyProductAndSumConventions()
        {
            var matrix = new DoubleMatrix();

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
            var entries = new[,] { { 0.5, 1.0 }, { 1.5, 2.0 } };
            var matrix = Create(entries);

            Assert.That(matrix.TryAddMatrix(null!), Is.False);
            Assert.That(matrix.TrySubtractMatrix(null!), Is.False);
            Assert.That(() => matrix.DotProduct(null!), Throws.TypeOf<ArgumentNullException>());
            Assert.That(() => matrix * null!, Throws.TypeOf<ArgumentNullException>());
            Assert.That(() => null! * matrix, Throws.TypeOf<ArgumentNullException>());
            AssertMatrix(matrix, entries);
        }

        private static DoubleMatrix Create(double[,] entries)
        {
            var matrix = new DoubleMatrix(entries.GetLength(0), entries.GetLength(1));
            for (var row = 0; row < matrix.RowCount; row++)
                for (var column = 0; column < matrix.ColumnCount; column++)
                    matrix[row, column] = entries[row, column];
            return matrix;
        }

        private static void AssertMatrix(DoubleMatrix actual, double[,] expected)
        {
            AssertMatrix(actual, Create(expected));
        }

        private static void AssertMatrix(DoubleMatrix actual, DoubleMatrix expected)
        {
            Assert.That(actual.RowCount, Is.EqualTo(expected.RowCount));
            Assert.That(actual.ColumnCount, Is.EqualTo(expected.ColumnCount));
            for (var row = 0; row < expected.RowCount; row++)
                for (var column = 0; column < expected.ColumnCount; column++)
                    Assert.That(actual[row, column], Is.EqualTo(expected[row, column]).Within(1e-12), $"Entry [{row}, {column}]");
        }
    }
}
