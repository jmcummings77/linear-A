using System.Diagnostics;
using System.Globalization;
using System.Text.Json;
using linear_A;

internal readonly record struct Result(Matrix<double>? Matrix, double Value = 0, bool? Upper = null, bool? Lower = null, SymmetricEigenDecomposition? Eigen = null, DoubleMatrix? Rotation = null, EigenDecomposition? GeneralEigen = null)
{
    public double Checksum
    {
        get
        {
            if (GeneralEigen is not null)
            {
                var total = 0.0;
                for (var i = 0; i < GeneralEigen.EigenValues.Length; i++) total += (i + 1) * (GeneralEigen.EigenValues[i].Real + Math.Abs(GeneralEigen.EigenValues[i].Imaginary));
                foreach (var value in GeneralEigen.EigenVectors) total += value.Real * value.Real + value.Imaginary * value.Imaginary;
                return total;
            }
            if (Rotation is not null)
            {
                var length = Rotation.RowCount * Rotation.ColumnCount;
                return Rotation[0, 0] + Rotation[length / 2 / Rotation.ColumnCount, length / 2 % Rotation.ColumnCount]
                    + Rotation[Rotation.RowCount - 1, Rotation.ColumnCount - 1];
            }
            if (Eigen is not null)
            {
                var total = 0.0;
                for (var i = 0; i < Eigen.EigenValues.Length; i++) total += (i + 1) * Eigen.EigenValues[i];
                var vectors = Eigen.EigenVectors;
                for (var row = 0; row < vectors.RowCount; row++)
                    for (var col = 0; col < vectors.ColumnCount; col++) total += vectors[row, col] * vectors[row, col];
                return total;
            }
            return Matrix is null ? Value : Matrix.RowCount * Matrix.ColumnCount == 0 ? 0 :
                At(0) + At(Matrix.RowCount * Matrix.ColumnCount / 2) + At(Matrix.RowCount * Matrix.ColumnCount - 1);
        }
    }
    private double At(int index) => Matrix![index / Matrix.ColumnCount, index % Matrix.ColumnCount];
    public void Write()
    {
        if (GeneralEigen is not null)
        {
            var rows = GeneralEigen.EigenValues.Length;
            var real = new double[rows * rows];
            var imaginary = new double[rows * rows];
            for (var row = 0; row < rows; row++)
                for (var col = 0; col < rows; col++)
                {
                    real[row * rows + col] = GeneralEigen.EigenVectors[row, col].Real;
                    imaginary[row * rows + col] = GeneralEigen.EigenVectors[row, col].Imaginary;
                }
            Console.WriteLine(JsonSerializer.Serialize(new
            {
                eigenvalues_real = GeneralEigen.EigenValues.Select(value => value.Real).ToArray(),
                eigenvalues_imag = GeneralEigen.EigenValues.Select(value => value.Imaginary).ToArray(),
                eigenvectors_real = new { rows, cols = rows, values = real },
                eigenvectors_imag = new { rows, cols = rows, values = imaginary }
            }));
        }
        else if (Eigen is not null)
        {
            var vectors = Eigen.EigenVectors;
            var values = new double[vectors.RowCount * vectors.ColumnCount];
            for (var row = 0; row < vectors.RowCount; row++)
                for (var col = 0; col < vectors.ColumnCount; col++) values[row * vectors.ColumnCount + col] = vectors[row, col];
            Console.WriteLine(JsonSerializer.Serialize(new
            {
                eigenvalues = Eigen.EigenValues,
                eigenvectors = new { rows = vectors.RowCount, cols = vectors.ColumnCount, values }
            }));
        }
        else if (Rotation is not null)
        {
            var values = new double[Rotation.RowCount * Rotation.ColumnCount];
            for (var row = 0; row < Rotation.RowCount; row++)
                for (var col = 0; col < Rotation.ColumnCount; col++) values[row * Rotation.ColumnCount + col] = Rotation[row, col];
            Console.WriteLine(JsonSerializer.Serialize(new { rows = Rotation.RowCount, cols = Rotation.ColumnCount, values }));
        }
        else if (Matrix is not null)
        {
            var values = Matrix.ToArray();
            if (values.Any(value => !double.IsFinite(value))) throw new ArithmeticException("Nonfinite result.");
            Console.WriteLine(JsonSerializer.Serialize(new { rows = Matrix.RowCount, cols = Matrix.ColumnCount, values }));
        }
        else if (Upper.HasValue) Console.WriteLine(JsonSerializer.Serialize(new { upper = Upper.Value, lower = Lower!.Value }));
        else
        {
            if (!double.IsFinite(Value)) throw new ArithmeticException("Nonfinite result.");
            Console.WriteLine(JsonSerializer.Serialize(new { value = Value }));
        }
    }
}

internal static class Program
{
    private static int Number(string value, bool positive = false)
    {
        var result = int.Parse(value, CultureInfo.InvariantCulture);
        if (result < (positive ? 1 : 0)) throw new ArgumentException("Invalid dimension or iteration count.");
        return result;
    }
    private static Matrix<double> Read(int rows, int cols, double[] values, ref int offset)
    {
        var matrix = new Matrix<double>(rows, cols);
        for (var row = 0; row < rows; row++)
            for (var col = 0; col < cols; col++) matrix[row, col] = values[offset++];
        return matrix;
    }
    private static Matrix<double> Generate(int size, int seed)
    {
        var matrix = new Matrix<double>(size, size);
        for (var row = 0; row < size; row++)
            for (var col = 0; col < size; col++)
                matrix[row, col] = (((long)row * size * 17 + col * 17L + seed * 13L) % 101 - 50) / 16.0;
        return matrix;
    }
    private static Func<Result> Operation(string op, Matrix<double> a, Matrix<double>? b, double scalar)
    {
        if (op == "rotation2d")
        {
            if (a.RowCount != 0 || a.ColumnCount != 0) throw new ArgumentException("rotation2d expects an empty 0-by-0 input.");
            return () => new(null, Rotation: MatrixRotation.Create2D(scalar));
        }
        if (op == "rotation3d")
        {
            if (!((a.RowCount == 3 && a.ColumnCount == 1) || (a.RowCount == 1 && a.ColumnCount == 3)))
                throw new ArgumentException("Expected a 3-by-1 or 1-by-3 axis.");
            var axis = new DoubleMatrix(a.RowCount, a.ColumnCount);
            for (var row = 0; row < a.RowCount; row++)
                for (var col = 0; col < a.ColumnCount; col++) axis[row, col] = a[row, col];
            return () => new(null, Rotation: MatrixRotation.CreateAxisAngle(axis, scalar));
        }
        if (op.StartsWith("determinant", StringComparison.Ordinal))
        {
            var algorithm = op switch
            {
                "determinant" => DeterminantAlgorithm.Auto,
                "determinant_lu" or "determinant_spd_lu" => DeterminantAlgorithm.Lu,
                "determinant_cholesky" => DeterminantAlgorithm.Cholesky,
                "determinant_cofactor" => DeterminantAlgorithm.Cofactor,
                _ => throw new ArgumentException("Unknown determinant algorithm.")
            };
            return () => new(null, a.GetDeterminant(algorithm));
        }
        return op switch
        {
            "add" => () => new(a + b!),
            "subtract" => () => new(a - b!),
            "multiply" => () => new(a * b!),
            "cross" => () => new(a.CrossProduct(b!)),
            "transpose" => () => new(a.Transpose()),
            "scale" => () => { var result = new Matrix<double>(a); result.Scale(scalar); return new(result); }
            ,
            "trace" => () => new(null, a.GetTrace()),
            "eigen_general" => () => new(null, GeneralEigen: a.GetEigenDecomposition()),
            "eigen_symmetric" => () => new(null, Eigen: a.GetSymmetricEigenDecomposition()),
            "triangular" => () => new(null, 0, a.IsUpperTriangular(), a.IsLowerTriangular()),
            _ => throw new ArgumentException("Unknown operation.")
        };
    }
    public static int Main(string[] args)
    {
        try
        {
            if (args.Length < 4) throw new ArgumentException("Expected check or bench arguments; see benchmarks/PROTOCOL.md.");
            var op = args[1];
            if (args[0] == "check")
            {
                var rows = Number(args[2]);
                var cols = Number(args[3]);
                var binary = op is "add" or "subtract" or "multiply" or "cross";
                var hasScalar = op is "scale" or "rotation2d" or "rotation3d";
                if (args.Length != (binary ? 6 : hasScalar ? 5 : 4)) throw new ArgumentException("Wrong argument count.");
                var br = binary ? Number(args[4]) : 0;
                var bc = binary ? Number(args[5]) : 0;
                var values = Console.In.ReadToEnd().Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries)
                    .Select(value => double.Parse(value, CultureInfo.InvariantCulture)).ToArray();
                if (values.Length != checked(rows * cols + br * bc) || values.Any(value => !double.IsFinite(value)))
                    throw new ArgumentException("Expected exactly the stated number of finite entries.");
                var offset = 0;
                var a = Read(rows, cols, values, ref offset);
                var b = binary ? Read(br, bc, values, ref offset) : null;
                var scalar = hasScalar ? double.Parse(args[4], CultureInfo.InvariantCulture) : 1.25;
                if (!double.IsFinite(scalar)) throw new ArgumentException("Scalar must be finite.");
                Operation(op, a, b, scalar)().Write();
            }
            else if (args[0] == "bench" && args.Length == 5 && op != "triangular")
            {
                var size = Number(args[2], true);
                var iterations = Number(args[3], true);
                var seed = Number(args[4]);
                if (seed == int.MaxValue) throw new ArgumentException("Seed must be <= 2147483646.");
                if ((op is "cross" or "rotation3d") && size != 3 || op == "rotation2d" && size != 2)
                    throw new ArgumentException("Cross and rotation3d require size 3; rotation2d requires size 2.");
                var a = Generate(size, seed);
                var b = Generate(size, seed + 1);
                if (op == "cross")
                {
                    a = new Matrix<double>(new[,] { { a[0, 0] }, { a[0, 1] }, { a[0, 2] } });
                    b = new Matrix<double>(new[,] { { b[0, 0] }, { b[0, 1] }, { -b[0, 2] } });
                }
                else if (op == "rotation2d") a = new Matrix<double>();
                else if (op == "rotation3d") a = new Matrix<double>(new[,] { { 1.0 }, { 2.0 }, { 3.0 } });
                if (op is "determinant_cholesky" or "determinant_spd_lu")
                    for (var row = 0; row < size; row++)
                        for (var col = 0; col <= row; col++)
                            a[row, col] = a[col, row] = (a[row, col] + a[col, row]) / 2;
                if (op.StartsWith("determinant", StringComparison.Ordinal)) for (var i = 0; i < size; i++) a[i, i] += size * 4.0;
                if (op == "eigen_symmetric")
                    for (var row = 0; row < size; row++)
                        for (var col = 0; col < size; col++)
                            a[row, col] = row == col ? 2.0 + (seed % 17) / 16.0 : Math.Abs(row - col) == 1 ? -1.0 : 0.0;
                if (op == "eigen_general")
                    for (var row = 0; row < size; row++)
                        for (var col = 0; col < size; col++)
                        {
                            var block = row / 2;
                            a[row, col] = row == col ? 1 + (seed % 17) / 16.0 + block / 8.0 :
                                row % 2 == 0 && col == row + 1 ? -(0.5 + block / 16.0) :
                                row % 2 == 1 && col == row - 1 ? 0.5 + block / 16.0 :
                                row < col ? ((row * 3L + col * 5L + seed) % 11 - 5) / 32.0 : 0;
                        }
                var operation = Operation(op, a, b, op is "rotation2d" or "rotation3d" ? 0.5 : 1.25);
                for (var i = 0; i < Math.Max(5, Math.Min(iterations, 100)); i++) _ = operation().Checksum;
                var checksum = 0.0;
                var start = Stopwatch.GetTimestamp();
                for (var i = 0; i < iterations; i++) checksum += operation().Checksum;
                var elapsed = Stopwatch.GetElapsedTime(start).TotalNanoseconds;
                if (!double.IsFinite(checksum)) throw new ArithmeticException("Nonfinite checksum.");
                Console.WriteLine(JsonSerializer.Serialize(new { elapsed_ns = elapsed, iterations, checksum }));
            }
            else throw new ArgumentException("Unknown mode or wrong argument count.");
            return 0;
        }
        catch (Exception error)
        {
            Console.Error.WriteLine(error.Message);
            return 2;
        }
    }
}
