namespace LinearA

open System

/// Reusable double-precision factorization algorithms.
type FactorizationAlgorithm = Lu | Cholesky | Qr

/// Owned factor snapshot; QR requires numerical full column rank.
type Factorization(source: Matrix, algorithm: FactorizationAlgorithm) =
    let finite value =
        if not (Double.IsFinite value) then raise (ArithmeticException("Solver arithmetic exceeds the finite double range."))
        value
    let scaled value scale =
        let result = finite (value / scale)
        if value <> 0.0 && result = 0.0 then raise (ArithmeticException("Solver scaling would discard a nonzero value."))
        result
    let m, n = source.Rows, source.Cols
    let original = source.ToArray()
    let mutable norm = 0.0
    let maximum = original |> Array.fold (fun value entry -> max value (abs (finite entry))) 0.0
    let scale = if maximum = 0.0 then 1.0 else Math.ScaleB(1.0, Math.ILogB maximum)
    let a = original |> Array.map (fun value -> scaled value scale)
    let tau = Array.zeroCreate n
    let permutation = Array.init m id
    do
        if m < n || (algorithm <> Qr && m <> n) then invalidArg "source" "Factorization requires square input, or rows >= columns for QR."
        for i in 0 .. m-1 do
            let mutable sum = 0.0
            for j in 0 .. n-1 do sum <- sum + abs a[i*n+j]
            norm <- max norm sum
        match algorithm with
        | Lu ->
            for k in 0 .. n-1 do
                let mutable pivot = k
                for i in k+1 .. n-1 do
                    if abs a[i*n+k] > abs a[pivot*n+k] then pivot <- i
                if a[pivot*n+k] = 0.0 then invalidArg "source" "Singular matrix: zero computed LU pivot."
                for j in 0 .. n-1 do
                    let saved = a[k*n+j]
                    a[k*n+j] <- a[pivot*n+j]
                    a[pivot*n+j] <- saved
                let saved = permutation[k]
                permutation[k] <- permutation[pivot]
                permutation[pivot] <- saved
                for i in k+1 .. n-1 do
                    a[i*n+k] <- finite (a[i*n+k] / a[k*n+k])
                    for j in k+1 .. n-1 do a[i*n+j] <- finite (a[i*n+j] - a[i*n+k] * a[k*n+j])
        | Cholesky ->
            for i in 0 .. n-1 do
                for j in 0 .. i-1 do
                    if original[i*n+j] <> original[j*n+i] then invalidArg "source" "Cholesky requires exact symmetry."
            for i in 0 .. n-1 do
                for j in 0 .. i do
                    let mutable value = a[i*n+j]
                    for k in 0 .. j-1 do value <- finite (value - a[i*n+k] * a[j*n+k])
                    if i = j then
                        if value <= 0.0 then invalidArg "source" "Cholesky requires positive computed pivots."
                        a[i*n+j] <- sqrt value
                    else a[i*n+j] <- finite (value / a[j*n+j])
        | Qr ->
            let mutable largest = 0.0
            for j in 0 .. n-1 do
                let mutable length = 0.0
                for i in 0 .. m-1 do length <- Double.Hypot(length, a[i*n+j])
                largest <- max largest length
            let threshold = 2.220446049250313e-16 * float m * largest
            for k in 0 .. n-1 do
                let mutable pivot = k
                let mutable length = 0.0
                for j in k .. n-1 do
                    let mutable candidate = 0.0
                    for i in k .. m-1 do candidate <- Double.Hypot(candidate, a[i*n+j])
                    if candidate > length then
                        pivot <- j
                        length <- candidate
                if length <= threshold then invalidArg "source" "QR input is numerically rank deficient."
                for i in 0 .. m-1 do
                    let saved = a[i*n+k]
                    a[i*n+k] <- a[i*n+pivot]
                    a[i*n+pivot] <- saved
                let saved = permutation[k]
                permutation[k] <- permutation[pivot]
                permutation[pivot] <- saved
                let old = a[k*n+k]
                let alpha = -Math.CopySign(length, old)
                let divisor = old-alpha
                tau[k] <- (alpha-old)/alpha
                for i in k+1 .. m-1 do a[i*n+k] <- a[i*n+k]/divisor
                a[k*n+k] <- alpha
                for j in k+1 .. n-1 do
                    let mutable dot = a[k*n+j]
                    for i in k+1 .. m-1 do dot <- dot + a[i*n+k]*a[i*n+j]
                    dot <- dot*tau[k]
                    a[k*n+j] <- finite (a[k*n+j]-dot)
                    for i in k+1 .. m-1 do a[i*n+j] <- finite (a[i*n+j]-a[i*n+k]*dot)

    member private _.SolveInternal(rhs: Matrix, rescale) =
        if rhs.Rows <> m then invalidArg "rhs" "Right-hand side row count must match factorization."
        let p = rhs.Cols
        let work = Array.zeroCreate (m*p)
        for i in 0 .. m-1 do
            for j in 0 .. p-1 do
                let row = if algorithm = Lu then permutation[i] else i
                let value = finite rhs[row,j]
                work[i*p+j] <- if rescale then scaled value scale else value
        if algorithm = Qr then
            for k in 0 .. n-1 do
                for j in 0 .. p-1 do
                    let mutable dot = work[k*p+j]
                    for i in k+1 .. m-1 do dot <- finite (dot + a[i*n+k]*work[i*p+j])
                    dot <- finite (dot*tau[k])
                    work[k*p+j] <- finite (work[k*p+j]-dot)
                    for i in k+1 .. m-1 do work[i*p+j] <- finite (work[i*p+j]-a[i*n+k]*dot)
        else
            for i in 0 .. n-1 do
                for j in 0 .. p-1 do
                    let mutable value = work[i*p+j]
                    for k in 0 .. i-1 do value <- finite (value-a[i*n+k]*work[k*p+j])
                    work[i*p+j] <- if algorithm = Cholesky then finite (value/a[i*n+i]) else value
        for i in n-1 .. -1 .. 0 do
            for j in 0 .. p-1 do
                let mutable value = work[i*p+j]
                for k in i+1 .. n-1 do
                    let coefficient = if algorithm = Cholesky then a[k*n+i] else a[i*n+k]
                    value <- finite (value-coefficient*work[k*p+j])
                work[i*p+j] <- finite (value/a[i*n+i])
        let result = Matrix(n,p)
        for i in 0 .. n-1 do
            for j in 0 .. p-1 do result[(if algorithm = Qr then permutation[i] else i),j] <- work[i*p+j]
        result

    member this.Solve(rhs) = this.SolveInternal(rhs, true)

    /// Computed infinity-norm rcond using an inverse, O(n^3), not a certified bound.
    member this.ReciprocalCondition() =
        if m <> n then invalidOp "Condition diagnostic requires a square matrix."
        if n = 0 then 1.0
        else
            try
                let inverse = this.SolveInternal(Matrix.Identity(n), false)
                let mutable inverseNorm = 0.0
                for i in 0 .. n-1 do
                    let mutable sum = 0.0
                    for j in 0 .. n-1 do sum <- sum + abs inverse[i,j]
                    inverseNorm <- max inverseNorm sum
                min 1.0 ((1.0/norm)/inverseNorm)
            with :? ArithmeticException -> 0.0

[<AutoOpen>]
module MatrixSolvers =
    type Matrix with
        member this.FactorLu() = Factorization(this, Lu)
        member this.FactorCholesky() = Factorization(this, Cholesky)
        member this.FactorQr() = Factorization(this, Qr)
        member this.Solve(rhs) = this.FactorLu().Solve(rhs)
        member this.LeastSquares(rhs) = this.FactorQr().Solve(rhs)
