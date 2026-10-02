module LinearA.Runner

open System
open System.Diagnostics
open System.Numerics
open System.Globalization
open System.Text.Json
open LinearA

type Result =
    | MatrixResult of Matrix
    | ScalarResult of double
    | GeneralEigenResult of Complex array * Complex[,]
    | EigenResult of double array * Matrix
    | TriangularResult of bool * bool

let number positive (value: string) =
    let parsed = Int32.Parse(value, CultureInfo.InvariantCulture)
    if parsed < (if positive then 1 else 0) then invalidArg "number" "Invalid dimension or iteration count."
    parsed

let checksum result =
    match result with
    | GeneralEigenResult(values, vectors) ->
        let mutable total = 0.0
        for i in 0 .. values.Length - 1 do total <- total + double (i+1) * (values[i].Real + abs values[i].Imaginary)
        for row in 0 .. values.Length - 1 do
            for col in 0 .. values.Length - 1 do
                let value = vectors[row,col]
                total <- total + value.Real * value.Real + value.Imaginary * value.Imaginary
        total
    | EigenResult(values, vectors) ->
        let mutable total = 0.0
        for i in 0 .. values.Length - 1 do total <- total + double (i+1) * values[i]
        for row in 0 .. vectors.Rows - 1 do
            for col in 0 .. vectors.Cols - 1 do total <- total + vectors[row,col] * vectors[row,col]
        total
    | ScalarResult value -> value
    | MatrixResult matrix ->
        let length = matrix.Rows * matrix.Cols
        if length = 0 then 0.0
        else
            let at index = matrix[index / matrix.Cols, index % matrix.Cols]
            at 0 + at (length / 2) + at (length - 1)
    | TriangularResult _ -> invalidOp "Triangular classification is not a benchmark operation."

let operation op (a: Matrix) (b: Matrix) scalar =
    match op with
    | "add" -> MatrixResult(a.Add(b))
    | "subtract" -> MatrixResult(a.Subtract(b))
    | "multiply" -> MatrixResult(a.Multiply(b))
    | "solve" -> MatrixResult(a.Solve(b))
    | "solve_cholesky" -> MatrixResult(a.FactorCholesky().Solve(b))
    | "least_squares" -> MatrixResult(a.LeastSquares(b))
    | "rcond" -> ScalarResult(a.FactorLu().ReciprocalCondition())
    | "cross" -> MatrixResult(a.Cross(b))
    | "rotation2d" ->
        if a.Rows <> 0 || a.Cols <> 0 then invalidArg "matrix" "rotation2d expects an empty 0-by-0 input."
        MatrixResult(Matrix.Rotation2D(scalar))
    | "rotation3d" -> MatrixResult(Matrix.RotationAxisAngle(a, scalar))
    | "scale" -> MatrixResult(a.Scale(scalar))
    | "transpose" -> MatrixResult(a.Transpose())
    | "trace" -> ScalarResult(a.Trace())
    | "eigen_general" -> GeneralEigenResult(a.EigenGeneral())
    | "eigen_symmetric" -> EigenResult(a.EigenSymmetric())
    | "determinant" -> ScalarResult(a.Determinant())
    | "determinant_lu" | "determinant_spd_lu" -> ScalarResult(a.Determinant("lu"))
    | "determinant_cholesky" -> ScalarResult(a.Determinant("cholesky"))
    | "triangular" -> TriangularResult(a.IsUpperTriangular(), a.IsLowerTriangular())
    | _ -> invalidArg "op" "Unknown operation."

let write result =
    let output =
        match result with
        | MatrixResult matrix ->
            let values = matrix.ToArray()
            if Array.exists (Double.IsFinite >> not) values then invalidOp "Nonfinite result."
            box {| rows = matrix.Rows; cols = matrix.Cols; values = values |}
        | GeneralEigenResult(values, vectors) ->
            let n = values.Length
            box {| eigenvalues_real = Array.map (fun (v: Complex) -> v.Real) values
                   eigenvalues_imag = Array.map (fun (v: Complex) -> v.Imaginary) values
                   eigenvectors_real = {| rows = n; cols = n; values = Array.init (n*n) (fun i -> vectors[i/n,i%n].Real) |}
                   eigenvectors_imag = {| rows = n; cols = n; values = Array.init (n*n) (fun i -> vectors[i/n,i%n].Imaginary) |} |}
        | EigenResult(values, vectors) ->
            box {| eigenvalues = values; eigenvectors = {| rows = vectors.Rows; cols = vectors.Cols; values = vectors.ToArray() |} |}
        | ScalarResult value ->
            if not (Double.IsFinite(value)) then invalidOp "Nonfinite result."
            box {| value = value |}
        | TriangularResult(upper, lower) -> box {| upper = upper; lower = lower |}
    Console.WriteLine(JsonSerializer.Serialize(output))

let generate size seed =
    let result = Matrix(size, size)
    for row in 0 .. size - 1 do
        for col in 0 .. size - 1 do
            let index = int64 row * int64 size + int64 col
            result[row, col] <- double ((index * 17L + int64 seed * 13L) % 101L - 50L) / 16.0
    result

[<EntryPoint>]
let main args =
    try
        if args.Length < 4 then invalidArg "args" "Expected check or bench arguments."
        let op = args[1]
        match args[0] with
        | "check" ->
            let rows, cols = number false args[2], number false args[3]
            let binary = op = "add" || op = "subtract" || op = "multiply" || op = "cross" || op = "solve" || op = "solve_cholesky" || op = "least_squares"
            let hasScalar = op = "scale" || op = "rotation2d" || op = "rotation3d"
            if args.Length <> (if binary then 6 elif hasScalar then 5 else 4) then
                invalidArg "args" "Wrong argument count."
            let br, bc = if binary then number false args[4], number false args[5] else 0, 0
            let values = Console.In.ReadToEnd().Split(([||]: char array), StringSplitOptions.RemoveEmptyEntries)
                         |> Array.map (fun value -> Double.Parse(value, CultureInfo.InvariantCulture))
            if int64 values.Length <> int64 rows * int64 cols + int64 br * int64 bc || Array.exists (Double.IsFinite >> not) values then
                invalidArg "values" "Expected exactly the stated number of finite entries."
            let a = Matrix.FromArray(rows, cols, Array.take (rows * cols) values)
            let b = Matrix.FromArray(br, bc, Array.skip (rows * cols) values)
            let scalar = if hasScalar then Double.Parse(args[4], CultureInfo.InvariantCulture) else 1.25
            if not (Double.IsFinite(scalar)) then invalidArg "scalar" "Scalar must be finite."
            operation op a b scalar |> write
        | "bench" when args.Length = 5 && op <> "triangular" ->
            let size, iterations, seed = number true args[2], number true args[3], number false args[4]
            if seed = Int32.MaxValue then invalidArg "seed" "Seed must be <= 2147483646."
            if (op = "cross" || op = "rotation3d") && size <> 3 || op = "rotation2d" && size <> 2 then
                invalidArg "size" "Cross and rotation3d require size 3; rotation2d requires size 2."
            let a, b =
                if op = "cross" then
                    Matrix.FromArray(3,1,(generate 3 seed).ToArray()[0..2]), Matrix.FromArray(3,1,(generate 3 (seed+1)).ToArray()[0..2])
                elif op = "rotation2d" then Matrix(0,0), Matrix(0,0)
                elif op = "rotation3d" then Matrix.FromArray(3,1,[|1.;2.;3.|]), Matrix(0,0)
                else generate size seed, generate size (seed + 1)
            if op = "cross" then b[2,0] <- -b[2,0]
            if op = "determinant_cholesky" || op = "determinant_spd_lu" then
                for row in 0 .. size - 1 do
                    for col in 0 .. row do
                        let value = (a[row, col] + a[col, row]) / 2.0
                        a[row, col] <- value
                        a[col, row] <- value
            if op.StartsWith("determinant", StringComparison.Ordinal) then
                for i in 0 .. size - 1 do a[i, i] <- a[i, i] + double size * 4.0
            if op = "eigen_symmetric" then
                for row in 0 .. size - 1 do
                    for col in 0 .. size - 1 do
                        a[row,col] <- if row = col then 2.0 + double (seed % 17) / 16.0 elif abs (row-col) = 1 then -1.0 else 0.0
            if op = "eigen_general" then
                for row in 0 .. size - 1 do
                    for col in 0 .. size - 1 do
                        let block = row/2
                        a[row,col] <-
                            if row = col then 1.0 + double (seed%17)/16.0 + double block/8.0
                            elif row%2 = 0 && col = row+1 then -(0.5+double block/16.0)
                            elif row%2 = 1 && col = row-1 then 0.5+double block/16.0
                            elif row < col then double ((int64 row*3L+int64 col*5L+int64 seed)%11L-5L)/32.0
                            else 0.0
            let scalar = if op = "rotation2d" || op = "rotation3d" then 0.5 else 1.25
            let run () = operation op a b scalar |> checksum
            for _ in 1 .. max 5 (min iterations 100) do run () |> ignore
            let mutable total = 0.0
            let start = Stopwatch.GetTimestamp()
            for _ in 1 .. iterations do total <- total + run ()
            let elapsed = Stopwatch.GetElapsedTime(start).TotalNanoseconds
            if not (Double.IsFinite(total)) then invalidOp "Nonfinite checksum."
            Console.WriteLine(JsonSerializer.Serialize({| elapsed_ns = elapsed; iterations = iterations; checksum = total |}))
        | _ -> invalidArg "args" "Unknown mode or wrong argument count."
        0
    with error ->
        Console.Error.WriteLine(error.Message)
        2
