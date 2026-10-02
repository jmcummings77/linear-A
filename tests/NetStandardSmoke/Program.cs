using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Reflection;
using System.Runtime.Versioning;
using linear_A;

internal static class Program
{
    private static int Main()
    {
        VerifyTarget();
        ArrayOwnershipAndEnumeration();
        BasicArithmeticAndEmptyShapes();
        Determinants();
        CrossProductsAndRotations();
        SymmetricEigenpairs();
        GeneralEigenpairs();
        Console.WriteLine(".NET Standard 2.1 runtime smoke checks passed.");
        return 0;
    }

    private static void Require(bool condition, string message)
    {
        if (!condition) throw new InvalidOperationException(message);
    }

    private static void Near(double actual, double expected, string message, double tolerance = 1e-12)
    {
        Require(double.IsFinite(actual) && Math.Abs(actual - expected) <= tolerance, message);
    }

    private static void Throws<TException>(Action operation) where TException : Exception
    {
        try { operation(); }
        catch (TException) { return; }
        throw new InvalidOperationException($"Expected {typeof(TException).Name}.");
    }

    private static void VerifyTarget()
    {
        var assembly = typeof(DoubleMatrix).Assembly;
        Require(assembly.GetCustomAttribute<TargetFrameworkAttribute>()?.FrameworkName == ".NETStandard,Version=v2.1",
            "The smoke runner must load the .NET Standard 2.1 library, not its .NET 10 build.");
        Require(typeof(IntegerMatrix).Assembly == assembly, "Both legacy types must come from the same assembly.");
        Require(assembly.GetType("linear_A.Matrix`1") is null, "Generic matrices must not appear in the standard target.");
        foreach (var type in new[] { typeof(DoubleMatrix), typeof(IntegerMatrix) })
        {
            Require(type.GetMethod("ToMatrix") is null && type.GetMethod("FromMatrix") is null,
                "Generic conversion methods must not appear in the standard target.");
        }
    }

    private static void ArrayOwnershipAndEnumeration()
    {
        var doubles = new double[,] { { 1, 2, 3 }, { 4, 5, 6 } };
        var matrix = DoubleMatrix.FromArray(doubles);
        doubles[0, 0] = 99;
        var snapshot = matrix.ToArray();
        snapshot[0, 1] = 99;
        Require(matrix[0, 0] == 1 && matrix[0, 1] == 2, "Double array boundaries must copy storage.");
        Require(snapshot.GetLength(0) == 2 && snapshot.GetLength(1) == 3, "ToArray must preserve rectangular shape.");
        Require(((IEnumerable<double>)matrix).SequenceEqual(new double[] { 1, 2, 3, 4, 5, 6 }),
            "Typed double enumeration must use row-major order.");
        IEnumerator originalEnumerator = matrix.GetEnumerator();
        var index = 1;
        while (originalEnumerator.MoveNext()) Require((double)originalEnumerator.Current == index++, "Legacy enumeration order changed.");
        Require(index == 7, "Legacy enumeration omitted values.");

        var integers = new int[,] { { 1, 2, 3 }, { 4, 5, 6 } };
        var integerMatrix = IntegerMatrix.FromArray(integers);
        integers[0, 0] = 99;
        var integerSnapshot = integerMatrix.ToArray();
        integerSnapshot[0, 1] = 99;
        Require(integerMatrix[0, 0] == 1 && integerMatrix[0, 1] == 2, "Integer array boundaries must copy storage.");
        Require(((IEnumerable<int>)integerMatrix).SequenceEqual(new[] { 1, 2, 3, 4, 5, 6 }),
            "Typed integer enumeration must use row-major order.");
        Require(((IEnumerable)integerMatrix).Cast<int>().Sum() == 21, "Nongeneric integer enumeration must remain available.");
        Throws<ArgumentNullException>(() => DoubleMatrix.FromArray(null!));
        Throws<ArgumentNullException>(() => IntegerMatrix.FromArray(null!));
    }

    private static void BasicArithmeticAndEmptyShapes()
    {
        var a = DoubleMatrix.FromArray(new double[,] { { 1, 2, 3 }, { 4, 5, 6 } });
        var gram = a * a.Transpose();
        Require(gram.RowCount == 2 && gram.ColumnCount == 2, "Double product shape changed.");
        Near(gram[0, 0], 14, "Wrong double dot product.");
        Near(gram[0, 1], 32, "Wrong rectangular double product.");
        Near(gram[1, 1], 77, "Wrong double dot product.");
        Near(gram.GetTrace(), 91, "Wrong double trace.");
        var sum = a + a;
        sum.Scale(.5);
        Require(((IEnumerable<double>)(sum - a)).All(value => value == 0), "Double add/subtract/scale disagree.");
        Require(a[0, 0] == 1 && a[1, 2] == 6, "Double operations changed an operand.");

        var integer = IntegerMatrix.FromArray(new[,] { { 1, 2, 3 }, { 4, 5, 6 } });
        var integerGram = integer * integer.Transpose();
        Require(integerGram[0, 0] == 14 && integerGram[0, 1] == 32 && integerGram[1, 1] == 77,
            "Wrong rectangular integer product.");
        Require(integerGram.GetTrace() == 91, "Wrong integer trace.");
        Require(((IEnumerable<int>)((integer + integer) - integer)).SequenceEqual((IEnumerable<int>)integer),
            "Integer add/subtract disagree.");
        var wrapped = IntegerMatrix.FromArray(new[,] { { int.MaxValue } }) + IntegerMatrix.FromArray(new[,] { { 1 } });
        Require(wrapped[0, 0] == int.MinValue, "Legacy integer addition must retain unchecked arithmetic.");

        var emptyLeft = DoubleMatrix.FromArray(new double[2, 0]);
        var emptyRight = DoubleMatrix.FromArray(new double[0, 3]);
        var zero = emptyLeft * emptyRight;
        Require(zero.RowCount == 2 && zero.ColumnCount == 3 && ((IEnumerable<double>)zero).All(value => value == 0),
            "Empty inner dimensions must produce a correctly shaped zero matrix.");
        Require(emptyLeft.Transpose().RowCount == 0 && emptyLeft.Transpose().ColumnCount == 2,
            "Empty rectangular transpose lost its shape.");
        var integerZero = IntegerMatrix.FromArray(new int[2, 0]) * IntegerMatrix.FromArray(new int[0, 3]);
        Require(integerZero.RowCount == 2 && integerZero.ColumnCount == 3 && ((IEnumerable<int>)integerZero).All(value => value == 0),
            "Integer empty shapes must remain usable in arithmetic.");
        Require(new DoubleMatrix().GetDeterminant() == 1 && new DoubleMatrix().GetTrace() == 0,
            "Empty double square conventions changed.");
        Require(new IntegerMatrix().GetDeterminant() == 1 && new IntegerMatrix().GetTrace() == 0,
            "Empty integer square conventions changed.");
        Throws<ArgumentException>(() => new DoubleMatrix(0, 2));
        Throws<ArgumentException>(() => new IntegerMatrix(2, 0));
    }

    private static void Determinants()
    {
        var pivoted = DoubleMatrix.FromArray(new double[,] { { 0, 2 }, { 3, 4 } });
        Near(pivoted.GetDeterminant(), -6, "Pivoted determinant has the wrong sign.");
        var scaled = DoubleMatrix.FromArray(new double[,]
        {
            { 1e200, 0, 0, 0 }, { 0, 1e200, 0, 0 },
            { 0, 0, 1e-200, 0 }, { 0, 0, 0, 1e-200 }
        });
        Near(scaled.GetDeterminant(), 1, "Diagonal product scaling lost a finite determinant.");
        var underflow = DoubleMatrix.FromArray(new double[,] { { 1e-300, 0 }, { 0, 1e-300 } });
        Require(underflow.GetDeterminant() == 0 && underflow.IsInvertible(), "Determinant underflow must not imply singularity.");
        var exact = IntegerMatrix.FromArray(new[,]
        {
            { int.MaxValue, int.MaxValue - 1 }, { int.MaxValue - 1, int.MaxValue - 2 }
        });
        Require(exact.GetDeterminant() == -1, "Integer determinant must retain exact intermediates.");
        var overflow = IntegerMatrix.FromArray(new[,] { { 50000, 0 }, { 0, 50000 } });
        Require(!overflow.TryGetDeterminant(out var value) && value == 0, "Out-of-range integer TryGetDeterminant must fail.");
        Throws<OverflowException>(() => overflow.GetDeterminant());
    }

    private static void CrossProductsAndRotations()
    {
        var x = DoubleMatrix.FromArray(new double[,] { { 1 }, { 0 }, { 0 } });
        var y = DoubleMatrix.FromArray(new double[,] { { 0, 1, 0 } });
        var z = x.GetCrossProduct(y);
        Require(z.RowCount == 3 && z.ColumnCount == 1 && z[2, 0] == 1 && x[0, 0] == 1,
            "Double cross product must follow the right-hand rule and preserve the left shape.");
        var integerX = IntegerMatrix.FromArray(new[,] { { 1, 0, 0 } });
        var integerY = IntegerMatrix.FromArray(new[,] { { 0 }, { 1 }, { 0 } });
        var integerZ = integerX.GetCrossProduct(integerY);
        Require(integerZ.RowCount == 1 && integerZ.ColumnCount == 3 && integerZ[0, 2] == 1,
            "Integer cross product must preserve a row-vector receiver's shape.");
        Require(integerX.CrossProduct(integerY) && integerX[0, 2] == 1, "Mutating integer cross product failed.");
        var tooLarge = IntegerMatrix.FromArray(new[,] { { int.MaxValue, 0, 0 } });
        var twiceY = IntegerMatrix.FromArray(new[,] { { 0, 2, 0 } });
        Require(!tooLarge.CrossProduct(twiceY) && tooLarge[0, 0] == int.MaxValue && tooLarge[0, 2] == 0,
            "A checked cross-product failure must leave the receiver unchanged.");
        Throws<OverflowException>(() => tooLarge.GetCrossProduct(twiceY));

        var axis = DoubleMatrix.FromArray(new double[,] { { 0 }, { 0 }, { 1e300 } });
        var rotation = MatrixRotation.CreateAxisAngle(axis, Math.PI / 2);
        var rotated = rotation * x;
        Near(rotated[0, 0], 0, "Rotation X component is wrong.");
        Near(rotated[1, 0], 1, "Positive Z rotation must map X to Y.");
        Near(rotation.GetDeterminant(), 1, "A rotation must preserve orientation and volume.");
        Near(axis[2, 0] / 1e300, 1, "Axis normalization must preserve its input.");
    }

    private static void SymmetricEigenpairs()
    {
        var a = DoubleMatrix.FromArray(new double[,] { { 2, 1 }, { 1, 2 } });
        var result = a.GetSymmetricEigenDecomposition();
        Near(result.EigenValues[0], 1, "Wrong lower eigenvalue.");
        Near(result.EigenValues[1], 3, "Wrong upper eigenvalue.");
        for (var row = 0; row < 2; row++)
        {
            for (var col = 0; col < 2; col++)
            {
                double av = 0, dot = 0;
                for (var k = 0; k < 2; k++)
                {
                    av += a[row, k] * result.EigenVectors[k, col];
                    dot += result.EigenVectors[k, row] * result.EigenVectors[k, col];
                }
                Near(av, result.EigenVectors[row, col] * result.EigenValues[col], "Eigenpair residual is too large.");
                Near(dot, row == col ? 1 : 0, "Eigenvectors must be orthonormal.");
            }
        }
        var integer = IntegerMatrix.FromArray(new[,] { { 1, 1 }, { 1, 0 } });
        var integerResult = integer.GetSymmetricEigenDecomposition();
        Near(integerResult.EigenValues[1], (1 + Math.Sqrt(5)) / 2, "Integer input must support noninteger eigenvalues.");
        Require(new DoubleMatrix().GetSymmetricEigenDecomposition().EigenValues.Length == 0,
            "Empty symmetric eigenvalue array must be supported.");
        var nonsymmetric = DoubleMatrix.FromArray(new double[,] { { 1, 1 }, { 0, 1 } });
        Require(!nonsymmetric.TryGetEigenValues(out var values) && values.Length == 0, "Nonsymmetric legacy eigen calls must fail explicitly.");
    }

    private static void GeneralEigenpairs()
    {
        var rotation = DoubleMatrix.FromArray(new double[,] { { 0, -1 }, { 1, 0 } });
        var integerRotation = IntegerMatrix.FromArray(new[,] { { 0, -1 }, { 1, 0 } });
        foreach (var result in new[] { rotation.GetEigenDecomposition(), integerRotation.GetEigenDecomposition() })
        {
            Require(result.EigenValues.Length == 2 && result.EigenVectors.GetLength(0) == 2 && result.EigenVectors.GetLength(1) == 2,
                "General eigenpair result has the wrong shape.");
            for (var col = 0; col < 2; col++)
            {
                Near(result.EigenValues[col].Real, 0, "A quarter turn must have purely imaginary eigenvalues.");
                Near(result.EigenValues[col].Imaginary, col == 0 ? -1 : 1, "Complex eigenvalues must sort by real then imaginary part.");
                var normSquared = 0d;
                for (var row = 0; row < 2; row++)
                {
                    var vector = result.EigenVectors[row, col];
                    normSquared += vector.Magnitude * vector.Magnitude;
                    var product = Complex.Zero;
                    for (var k = 0; k < 2; k++) product += rotation[row, k] * result.EigenVectors[k, col];
                    Near((product - result.EigenValues[col] * vector).Magnitude, 0, "Complex right-eigenpair residual is too large.");
                }
                Near(normSquared, 1, "Complex right-eigenvector columns must have unit norm.");
            }
        }
        Require(rotation[0, 1] == -1 && rotation[1, 0] == 1 && integerRotation[0, 1] == -1 && integerRotation[1, 0] == 1,
            "General eigendecomposition must preserve both legacy matrix inputs.");
    }
}
