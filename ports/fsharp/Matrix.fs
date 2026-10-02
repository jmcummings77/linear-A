namespace LinearA

open System

/// An owned, row-major float64 matrix. Arithmetic operations return independent matrices.
type Matrix private (rows: int, cols: int, values: double array) =
    static member private Length(rows, cols) =
        if rows < 0 || cols < 0 || int64 rows * int64 cols > int64 Array.MaxLength then
            invalidArg "dimensions" "Invalid matrix dimensions."
        rows * cols

    new(rows, cols) = Matrix(rows, cols, Array.zeroCreate (Matrix.Length(rows, cols)))
    member _.Rows = rows
    member _.Cols = cols

    member private _.Index(row, col) =
        if row < 0 || row >= rows || col < 0 || col >= cols then
            invalidArg "index" "Matrix index out of range."
        row * cols + col

    member this.Item
        with get(row, col) = values[this.Index(row, col)]
        and set(row, col) value = values[this.Index(row, col)] <- value

    static member FromArray(rows, cols, data: double array) =
        if isNull (box data) || data.Length <> Matrix.Length(rows, cols) then
            invalidArg "data" "Entries must match the dimensions."
        Matrix(rows, cols, Array.copy data)

    static member Identity(size) =
        let result = Matrix(size, size)
        for i in 0 .. size - 1 do result[i, i] <- 1.0
        result

    member _.Copy() = Matrix(rows, cols, Array.copy values)
    member _.ToArray() = Array.copy values
    member this.GetRow(row) =
        if row < 0 || row >= rows then invalidArg "row" "Row out of range."
        Array.init cols (fun col -> this[row, col])
    member this.GetColumn(col) =
        if col < 0 || col >= cols then invalidArg "col" "Column out of range."
        Array.init rows (fun row -> this[row, col])

    member private this.Combine(other: Matrix, subtract) =
        if isNull (box other) || rows <> other.Rows || cols <> other.Cols then
            invalidArg "other" "Matrix dimensions must match."
        let result = Matrix(rows, cols)
        for row in 0 .. rows - 1 do
            for col in 0 .. cols - 1 do
                result[row, col] <- if subtract then this[row, col] - other[row, col] else this[row, col] + other[row, col]
        result

    member this.Add(other) = this.Combine(other, false)
    member this.Subtract(other) = this.Combine(other, true)
    member _.Scale(scalar) = Matrix(rows, cols, Array.map (fun value -> value * scalar) values)
    member this.Transpose() =
        let result = Matrix(cols, rows)
        for row in 0 .. rows - 1 do
            for col in 0 .. cols - 1 do result[col, row] <- this[row, col]
        result

    member this.Multiply(other: Matrix) =
        if isNull (box other) || cols <> other.Rows then invalidArg "other" "Inner dimensions must match."
        let result = Matrix(rows, other.Cols)
        for row in 0 .. rows - 1 do
            for col in 0 .. other.Cols - 1 do
                let mutable sum = 0.0
                for k in 0 .. cols - 1 do sum <- sum + this[row, k] * other[k, col]
                result[row, col] <- sum
        result

    member this.Trace() =
        if rows <> cols then invalidOp "Trace requires a square matrix."
        let mutable sum = 0.0
        for i in 0 .. rows - 1 do sum <- sum + this[i, i]
        sum

    member private this.Vector3() =
        if not ((rows = 3 && cols = 1) || (rows = 1 && cols = 3)) then
            invalidArg "vector" "Expected a 3-by-1 or 1-by-3 vector."
        if Array.exists (Double.IsFinite >> not) values then invalidArg "vector" "Expected finite vector entries."
        Array.copy values

    /// Right-handed 3D vector cross product, preserving the left vector's shape.
    member this.Cross(other: Matrix) =
        if isNull (box other) then nullArg "other"
        let a, b = this.Vector3(), other.Vector3()
        let result = [| a[1]*b[2] - a[2]*b[1]; a[2]*b[0] - a[0]*b[2]; a[0]*b[1] - a[1]*b[0] |]
        if Array.exists (Double.IsFinite >> not) result then invalidOp "Cross product exceeds the finite float64 range."
        Matrix(rows, cols, result)

    static member private Angle(radians: double) =
        if not (Double.IsFinite(radians)) then invalidArg "radians" "Expected a finite angle in radians."
        cos radians, sin radians

    /// Active counterclockwise rotation of column vectors, in radians.
    static member Rotation2D(radians) =
        let c, s = Matrix.Angle(radians)
        Matrix(2, 2, [|c; -s; s; c|])

    static member RotationX(radians) =
        let c, s = Matrix.Angle(radians)
        Matrix(3, 3, [|1.;0.;0.; 0.;c;-s; 0.;s;c|])

    static member RotationY(radians) =
        let c, s = Matrix.Angle(radians)
        Matrix(3, 3, [|c;0.;s; 0.;1.;0.; -s;0.;c|])

    static member RotationZ(radians) =
        let c, s = Matrix.Angle(radians)
        Matrix(3, 3, [|c;-s;0.; s;c;0.; 0.;0.;1.|])

    /// Right-handed active rotation about a nonzero axis; magnitude is ignored.
    static member RotationAxisAngle(axis: Matrix, radians) =
        if isNull (box axis) then nullArg "axis"
        let v = axis.Vector3()
        let largest = Array.map abs v |> Array.max
        if largest = 0.0 then invalidArg "axis" "Rotation axis must be nonzero."
        let scaled = Array.map (fun value -> value / largest) v
        let norm = sqrt (Array.sumBy (fun value -> value*value) scaled)
        let x, y, z = scaled[0]/norm, scaled[1]/norm, scaled[2]/norm
        let c, s = Matrix.Angle(radians)
        let t = if abs radians < 1.0 then 2.0 * (sin (radians/2.0)) ** 2.0 else 1.0 - c
        Matrix(3, 3, [|c+x*x*t; x*y*t-z*s; x*z*t+y*s;
                      y*x*t+z*s; c+y*y*t; y*z*t-x*s;
                      z*x*t-y*s; z*y*t+x*s; c+z*z*t|])

    /// Auto uses triangular shortcuts; Lu uses partial pivoting; Cholesky requires SPD input.
    member this.Determinant(?algorithm: string) =
        if rows <> cols then invalidOp "Determinant requires a square matrix."
        let method = defaultArg algorithm "auto"
        if method <> "auto" && method <> "lu" && method <> "cholesky" then
            invalidArg "algorithm" "Expected auto, lu, or cholesky."
        if Array.exists (Double.IsFinite >> not) values then invalidOp "Determinant requires finite entries."
        let product (diagonal: double array) sign repeats initialExponent =
            let mutable mantissa = sign
            let mutable exponent = initialExponent
            let mutable zero = false
            for value in diagonal do
                if value = 0.0 then zero <- true
                else
                    let power = Math.ILogB(abs value)
                    let fraction = Math.ScaleB(value, -power)
                    for _ in 1 .. repeats do
                        mantissa <- mantissa * fraction
                        let shift = Math.ILogB(abs mantissa)
                        mantissa <- Math.ScaleB(mantissa, -shift)
                        exponent <- exponent + power + shift
            let result = if zero then 0.0 else Math.ScaleB(mantissa, exponent)
            if not (Double.IsFinite(result)) then invalidOp "Nonfinite determinant."
            result
        if method = "auto" && this.IsTriangular() then
            product (Array.init rows (fun i -> this[i, i])) 1.0 1 0
        elif method = "cholesky" then
            for row in 0 .. rows - 1 do
                for col in 0 .. row - 1 do
                    if this[row, col] <> this[col, row] then invalidArg "matrix" "Cholesky requires exact symmetry."
            let lower = Matrix(rows, cols)
            for row in 0 .. rows - 1 do
                for col in 0 .. row do
                    let mutable value = this[row, col]
                    for k in 0 .. col - 1 do value <- value - lower[row, k] * lower[col, k]
                    if not (Double.IsFinite(value)) then invalidOp "Nonfinite factorization."
                    if row = col then
                        if value <= 0.0 then invalidArg "matrix" "Cholesky requires positive definite input."
                        lower[row, col] <- sqrt value
                    else
                        lower[row, col] <- value / lower[col, col]
                        if not (Double.IsFinite(lower[row, col])) then invalidOp "Nonfinite factorization."
            product (Array.init rows (fun i -> lower[i, i])) 1.0 2 0
        else
            let work = this.Copy()
            let mutable rowExponent = 0
            // Scale only when every entry survives an exact power-of-two round trip.
            for row in 0 .. rows - 1 do
                let largest = Array.init cols (fun col -> abs work[row, col]) |> Array.fold max 0.0
                if largest > 0.0 then
                    let power = Math.ILogB(largest)
                    let safe = seq { for col in 0 .. cols - 1 -> Math.ScaleB(Math.ScaleB(work[row,col], -power), power) = work[row,col] } |> Seq.forall id
                    if safe then
                        for col in 0 .. cols - 1 do work[row,col] <- Math.ScaleB(work[row,col], -power)
                        rowExponent <- rowExponent + power
            let mutable sign = 1.0
            let mutable singular = false
            for k in 0 .. rows - 1 do
                if not singular then
                    let mutable pivotRow = k
                    for row in k + 1 .. rows - 1 do
                        if abs work[row, k] > abs work[pivotRow, k] then pivotRow <- row
                    let pivot = work[pivotRow, k]
                    if pivot = 0.0 then singular <- true
                    else
                        if pivotRow <> k then
                            for col in k .. cols - 1 do
                                let old = work[k, col]
                                work[k, col] <- work[pivotRow, col]
                                work[pivotRow, col] <- old
                            sign <- -sign
                        for row in k + 1 .. rows - 1 do
                            let numerator = work[row, k]
                            let factor = numerator / pivot
                            work[row, k] <- 0.0
                            for col in k + 1 .. cols - 1 do
                                let entry = work[k, col]
                                let update =
                                    if numerator <> 0.0 && entry <> 0.0 && abs factor < Math.ScaleB(1.0, -1022) then
                                        let en, ep, ee = Math.ILogB(abs numerator), Math.ILogB(abs pivot), Math.ILogB(abs entry)
                                        Math.ScaleB((Math.ScaleB(numerator, -en) / Math.ScaleB(pivot, -ep)) * Math.ScaleB(entry, -ee), en - ep + ee)
                                    else factor * entry
                                work[row, col] <- work[row, col] - update
                                if not (Double.IsFinite(work[row, col])) then invalidOp "Nonfinite factorization."
            if singular then 0.0 else product (Array.init rows (fun i -> work[i, i])) sign 1 rowExponent

    /// Complex right eigenpairs of a finite real square matrix, sorted by real then imaginary part.
    /// Columns are normalized; defective matrices need not have an independent eigenbasis.
    /// Hessenberg/QR iteration is bounded per unresolved root; nonconvergence fails.
    member this.EigenGeneral(?maxIterations: int) =
        let input = Array2D.init rows cols (fun row col -> this[row,col])
        GeneralEigen.solve input (defaultArg maxIterations 1000)

    /// All eigenpairs of a finite, exactly symmetric matrix. Eigenvectors are columns.
    /// The tolerance is relative to the matrix Frobenius norm, not each eigenvalue.
    member this.EigenSymmetric(?tolerance: double, ?maxSweeps: int) =
        let tolerance = defaultArg tolerance 1e-12
        let maxSweeps = defaultArg maxSweeps 50
        if not (Double.IsFinite(tolerance)) || tolerance <= 0.0 || tolerance >= 1.0 then
            invalidArg "tolerance" "Expected a finite tolerance strictly between zero and one."
        if maxSweeps < 1 then invalidArg "maxSweeps" "Expected a positive sweep limit."
        if rows <> cols then invalidOp "Eigen decomposition requires a square matrix."
        let mutable diagonal = true
        let mutable largest = 0.0
        for row in 0 .. rows - 1 do
            for col in 0 .. cols - 1 do
                let value = this[row, col]
                if not (Double.IsFinite(value)) then invalidArg "matrix" "Expected finite entries."
                if value <> this[col, row] then invalidArg "matrix" "Expected exact symmetry."
                if row <> col && value <> 0.0 then diagonal <- false
                largest <- max largest (abs value)
        let work = this.Copy()
        let vectors = Matrix.Identity(rows)
        let exponent = if diagonal || largest = 0.0 then 0 else Math.ILogB(largest)
        let hypot x y =
            let large, small = max (abs x) (abs y), min (abs x) (abs y)
            if large = 0.0 then 0.0 else large * sqrt (1.0 + (small / large) * (small / large))
        if not diagonal then
            let mutable norm = 0.0
            for row in 0 .. rows - 1 do
                for col in 0 .. cols - 1 do
                    work[row,col] <- Math.ScaleB(work[row,col], -exponent)
                    norm <- hypot norm work[row,col]
            let limit = tolerance * norm
            let skip = limit / (2.0 * double rows)
            let offNorm () =
                let mutable total = 0.0
                for row in 0 .. rows - 1 do
                    for col in row + 1 .. cols - 1 do
                        total <- hypot total work[row,col]
                sqrt 2.0 * total
            let mutable sweeps = 0
            while offNorm () > limit && sweeps < maxSweeps do
                for p in 0 .. rows - 2 do
                    for q in p + 1 .. rows - 1 do
                        let b = work[p,q]
                        if abs b > skip then
                            let delta = (work[q,q] - work[p,p]) / 2.0
                            let t =
                                if delta = 0.0 then (if b >= 0.0 then 1.0 else -1.0)
                                else b / (delta + Math.CopySign(hypot delta b, delta))
                            let c = 1.0 / sqrt (1.0 + t*t)
                            let s = t*c
                            work[p,p] <- work[p,p] - t*b
                            work[q,q] <- work[q,q] + t*b
                            work[p,q] <- 0.0
                            work[q,p] <- 0.0
                            for k in 0 .. rows - 1 do
                                if k <> p && k <> q then
                                    let x, y = work[k,p], work[k,q]
                                    let left, right = c*x - s*y, s*x + c*y
                                    work[k,p] <- left
                                    work[p,k] <- left
                                    work[k,q] <- right
                                    work[q,k] <- right
                                let x, y = vectors[k,p], vectors[k,q]
                                vectors[k,p] <- c*x - s*y
                                vectors[k,q] <- s*x + c*y
                sweeps <- sweeps + 1
            if offNorm () > limit then invalidOp "Symmetric eigen decomposition did not converge within the sweep limit."
        let eigenvalues = Array.init rows (fun i -> Math.ScaleB(work[i,i], exponent))
        if Array.exists (Double.IsFinite >> not) eigenvalues then invalidOp "Eigenvalues exceed the finite float64 range."
        let order = Array.init rows id |> Array.sortBy (fun i -> eigenvalues[i], i)
        let sorted = Matrix(rows, rows)
        for col in 0 .. rows - 1 do
            let source = order[col]
            let mutable pivot = 0
            let mutable norm = 0.0
            for row in 0 .. rows - 1 do
                norm <- hypot norm vectors[row,source]
                if abs vectors[row,source] > abs vectors[pivot,source] then pivot <- row
            let divisor = if vectors[pivot,source] < 0.0 then -norm else norm
            for row in 0 .. rows - 1 do sorted[row,col] <- vectors[row,source] / divisor
        (Array.map (fun i -> eigenvalues[i]) order, sorted)

    member this.IsUpperTriangular() =
        let mutable result = rows = cols
        if result then
            for row in 0 .. rows - 1 do
                for col in 0 .. row - 1 do
                    if this[row, col] <> 0.0 then result <- false
        result

    member this.IsLowerTriangular() =
        let mutable result = rows = cols
        if result then
            for row in 0 .. rows - 1 do
                for col in row + 1 .. cols - 1 do
                    if this[row, col] <> 0.0 then result <- false
        result

    member this.IsTriangular() = this.IsUpperTriangular() || this.IsLowerTriangular()
