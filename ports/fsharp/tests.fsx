#load "GeneralEigen.fs"
#load "Matrix.fs"
open System.Numerics
open System
open LinearA

let equal actual expected = if actual <> expected then failwithf "%A <> %A" actual expected
let rejects action =
    let mutable rejected = false
    try action () with :? ArgumentException -> rejected <- true | :? InvalidOperationException -> rejected <- true
    if not rejected then failwith "Invalid operation accepted"
let original = [| 1.; 2.; 3.; 4.; 5.; 6. |]
let a = Matrix.FromArray(2, 3, original)
original[0] <- 99.
equal a[0,0] 1.
let transpose = a.Transpose()
equal (transpose.ToArray()) [|1.;4.;2.;5.;3.;6.|]
equal (a.Multiply(transpose).ToArray()) [|14.;32.;32.;77.|]
equal (a.Multiply(transpose).Determinant()) 54.
equal (Matrix.FromArray(3,3,[|6.;1.;1.;4.;-2.;5.;2.;8.;7.|]).Determinant()) -306.
let copy = a.Copy()
copy[0,0] <- 12.
equal a[0,0] 1.
let row = a.GetRow(0)
row[0] <- 13.
equal a[0,0] 1.
equal (a.GetColumn(1)) [|2.;5.|]
equal ((a.Add(a).Subtract(a)).ToArray()) (a.ToArray())
equal (a.Scale(0.5).ToArray()) [|0.5;1.;1.5;2.;2.5;3.|]
equal (Matrix.Identity(3).Trace()) 3.
equal (Matrix.Identity(3).IsTriangular()) true
equal (Matrix(0,0).Determinant()) 1.
equal (Matrix(2,0).Multiply(Matrix(0,3)).ToArray()) (Array.zeroCreate 6)
rejects (fun () -> a.Trace() |> ignore)
rejects (fun () -> a.Multiply(a) |> ignore)
rejects (fun () -> a[-1,0] |> ignore)
printfn "F# matrix API checks passed"

let near actual expected =
    if abs (actual - expected) > max 1e-12 (abs expected * 1e-12) then failwithf "%g <> %g" actual expected
for algorithm in ["auto"; "lu"; "cholesky"] do
    let spd = Matrix.FromArray(2,2,[|4.;2.;2.;3.|])
    near (spd.Determinant(algorithm)) 8.
    equal (spd.ToArray()) [|4.;2.;2.;3.|]
    equal (Matrix(0,0).Determinant(algorithm)) 1.
    for diagonal in [[|1e200;1e200;1e-200;1e-200|]; [|1e-200;1e-200;1e200;1e200|]] do
        let scaled = Matrix(4,4)
        for i in 0 .. 3 do scaled[i,i] <- diagonal[i]
        near (scaled.Determinant(algorithm)) 1.
for values in [[|1.;2.;0.;1.|]; [|1.;2.;2.;1.|]; [|1.;1.;1.;1.|]] do
    rejects (fun () -> Matrix.FromArray(2,2,values).Determinant("cholesky") |> ignore)
rejects (fun () -> Matrix.Identity(2).Determinant("unknown") |> ignore)
for algorithm in ["auto"; "lu"] do
    equal (Matrix.FromArray(2,2,[|0.;1.;1.;0.|]).Determinant(algorithm)) -1.
    equal (Matrix.FromArray(2,2,[|1.;2.;2.;4.|]).Determinant(algorithm)) 0.
printfn "F# determinant algorithm checks passed"
let mixed = Matrix.FromArray(3,3,[|1e308;1e308;1e-308;1e-308;2e-308;1e308;0.;0.;1e-308|])
for algorithm in ["auto"; "lu"] do near (mixed.Determinant(algorithm) / 1e-308) 1.

let checkEigen (a: Matrix) (expected: double array) =
    let original = a.ToArray()
    let values, vectors = a.EigenSymmetric()
    equal values.Length a.Rows
    equal (a.ToArray()) original
    let scale = Array.fold (fun current value -> max current (abs value)) 0.0 original |> max 1e-300
    for i in 0 .. values.Length - 1 do
        if i > 0 && values[i] < values[i-1] then failwith "Unsorted spectrum"
        near (values[i] / scale) (expected[i] / scale)
        for row in 0 .. a.Rows - 1 do
            let mutable av = 0.0
            for k in 0 .. a.Cols - 1 do av <- av + (a[row,k] / scale) * vectors[k,i]
            near av ((values[i] / scale) * vectors[row,i])
        for j in 0 .. values.Length - 1 do
            let mutable dot = 0.0
            for row in 0 .. a.Rows - 1 do dot <- dot + vectors[row,i] * vectors[row,j]
            near dot (if i = j then 1.0 else 0.0)
checkEigen (Matrix(0,0)) [||]
checkEigen (Matrix(3,3)) [|0.;0.;0.|]
checkEigen (Matrix.Identity(3)) [|1.;1.;1.|]
checkEigen (Matrix.FromArray(2,2,[|2.;1.;1.;2.|])) [|1.;3.|]
checkEigen (Matrix.FromArray(2,2,[|0.;-2.;-2.;0.|])) [|-2.;2.|]
for scale in [1e-300; 1e300] do
    checkEigen (Matrix.FromArray(2,2,[|2.*scale;scale;scale;2.*scale|])) [|scale;3.*scale|]
let values, _ = Matrix.FromArray(2,2,[|1e-300;0.;0.;1e300|]).EigenSymmetric()
equal values [|1e-300;1e300|]
for size in [3;8;16] do
    let a = Matrix(size,size)
    for i in 0 .. size - 1 do
        a[i,i] <- 2.
        if i > 0 then
            a[i,i-1] <- -1.
            a[i-1,i] <- -1.
    checkEigen a (Array.init size (fun i -> 2. - 2. * cos (double (i+1) * Math.PI / double (size+1))))
    if size = 8 then rejects (fun () -> a.EigenSymmetric(maxSweeps=1) |> ignore)
rejects (fun () -> Matrix(2,3).EigenSymmetric() |> ignore)
rejects (fun () -> Matrix.FromArray(2,2,[|1.;1.;0.;1.|]).EigenSymmetric() |> ignore)
rejects (fun () -> Matrix.FromArray(1,1,[|Double.NaN|]).EigenSymmetric() |> ignore)
rejects (fun () -> Matrix.Identity(2).EigenSymmetric(tolerance=0.0) |> ignore)
rejects (fun () -> Matrix.Identity(2).EigenSymmetric(maxSweeps=0) |> ignore)
printfn "F# symmetric eigenpair checks passed"

let vector values = Matrix.FromArray(3,1,values)
let matrixNear (actual: Matrix) (expected: Matrix) =
    equal (actual.Rows,actual.Cols) (expected.Rows,expected.Cols)
    Array.iter2 near (actual.ToArray()) (expected.ToArray())
let vx, vy, vz = vector [|1.;0.;0.|], vector [|0.;1.;0.|], vector [|0.;0.;1.|]
matrixNear (vx.Cross(vy)) vz
matrixNear (vy.Cross(vx)) (vz.Scale(-1.))
equal (vx.Transpose().Cross(vy).Rows) 1
let ca, cb = vector [|1.;2.;3.|], vector [|-4.;5.;-6.|]
let cc = ca.Cross(cb)
equal (cc.ToArray()) [|-27.;-6.;13.|]
near (Array.map2 (*) (ca.ToArray()) (cc.ToArray()) |> Array.sum) 0.
near (Array.map2 (*) (cb.ToArray()) (cc.ToArray()) |> Array.sum) 0.
equal (ca.ToArray()) [|1.;2.;3.|]
matrixNear (ca.Cross(ca)) (Matrix(3,1))
for bad in [Matrix(0,0); Matrix(3,3); Matrix(2,1); vector [|Double.NaN;0.;1.|]] do
    rejects (fun () -> ca.Cross(bad) |> ignore)
    rejects (fun () -> Matrix.RotationAxisAngle(bad, 0.5) |> ignore)
rejects (fun () -> (vector [|1e308;1e308;0.|]).Cross(vector [|0.;1e308;1e308|]) |> ignore)
rejects (fun () -> Matrix.RotationAxisAngle(Matrix(3,1), 0.) |> ignore)
for angle in [Double.NaN;Double.PositiveInfinity;Double.NegativeInfinity] do
    rejects (fun () -> Matrix.Rotation2D(angle) |> ignore)
    rejects (fun () -> Matrix.RotationAxisAngle(vx,angle) |> ignore)
matrixNear (Matrix.RotationX(Math.PI/2.).Multiply(vy)) vz
matrixNear (Matrix.RotationY(Math.PI/2.).Multiply(vz)) vx
matrixNear (Matrix.RotationZ(Math.PI/2.).Multiply(vx)) vy
matrixNear (Matrix.Rotation2D(Math.PI/2.).Multiply(Matrix.FromArray(2,1,[|1.;0.|]))) (Matrix.FromArray(2,1,[|0.;1.|]))
for axis in [vector [|1.;2.;3.|]; vector [|1e300;2e300;3e300|]; vector [|1e-300;2e-300;3e-300|]; vector [|Double.Epsilon;0.;0.|]] do
    let before = axis.ToArray()
    let unit = vector (before |> Array.map (fun v -> v/(Array.map abs before |> Array.max)))
    for angle in [0.;0.5;-1.3;Math.PI;1e300] do
        let rotation = Matrix.RotationAxisAngle(axis,angle)
        matrixNear (rotation.Transpose().Multiply(rotation)) (Matrix.Identity(3))
        near (rotation.Determinant()) 1.
        matrixNear (rotation.Multiply(unit)) unit
        matrixNear (rotation.Multiply(Matrix.RotationAxisAngle(axis,-angle))) (Matrix.Identity(3))
        equal (axis.ToArray()) before
        let ra, rb = rotation.Multiply(ca), rotation.Multiply(cb)
        matrixNear (ra.Cross(rb)) (rotation.Multiply(cc))
matrixNear (Matrix.RotationAxisAngle(vector [|1e300;2e300;3e300|],0.5)) (Matrix.RotationAxisAngle(vector [|1.;2.;3.|],0.5))
equal (Matrix.RotationZ(1e-200)[1,0]) 1e-200
printfn "F# cross product and rotation checks passed"

let verifyGeneral (input: Matrix) =
    let before = input.ToArray()
    let values,vectors = input.EigenGeneral()
    equal (input.ToArray()) before
    equal values.Length input.Rows
    let scale = before |> Array.fold (fun maximum value -> max maximum (abs value)) 1.0
    for col in 0 .. input.Cols-1 do
        let mutable norm,residual = 0.0,0.0
        for row in 0 .. input.Rows-1 do
            norm <- norm + vectors[row,col].Magnitude ** 2.0
            let mutable product = Complex.Zero
            for k in 0 .. input.Cols-1 do product <- product + (input[row,k]/scale)*vectors[k,col]
            residual <- Double.Hypot(residual,(product-(values[col]/scale)*vectors[row,col]).Magnitude)
        near norm 1.0
        if residual > 1e-10 then failwithf "General eigenpair residual %g" residual
    values,vectors
let generalRotation = Matrix.FromArray(2,2,[|0.;-1.;1.;0.|])
let rotationValues,_ = verifyGeneral generalRotation
equal rotationValues [|-Complex.ImaginaryOne;Complex.ImaginaryOne|]
for input in [Matrix(0,0); Matrix(3,3); Matrix.FromArray(2,2,[|2.;1.;0.;2.|]);
              Matrix.FromArray(3,3,[|0.;1.;0.;0.;0.;1.;0.;0.;0.|]);
              Matrix.FromArray(3,3,[|1.;-2.;3.;2.;1.;4.;0.;0.;-3.|]);
              Matrix.FromArray(3,3,[|3.;0.;0.;2.;-1.;0.;1.;4.;2.|])] do
    verifyGeneral input |> ignore
let balancedValues,balancedVectors = verifyGeneral (Matrix.FromArray(2,2,[|0.;1e300;-1e-300;0.|]))
near balancedValues[0].Imaginary -1.0
near balancedValues[1].Imaginary 1.0
if abs (abs balancedVectors[1,0].Imaginary/1e-300-1.0) > 1e-12 then failwith "Lost small balanced eigenvector component"
for scale in [1e-200;1e200] do
    let values,_ = verifyGeneral (Matrix.FromArray(2,2,[|scale;-2.*scale;2.*scale;scale|]))
    near (values[0].Real/scale) 1.0
    near (values[0].Imaginary/scale) -2.0
let diagonalValues,_ = verifyGeneral (Matrix.FromArray(2,2,[|1e300;0.;0.;1e-300|]))
equal diagonalValues [|Complex(1e-300,0.);Complex(1e300,0.)|]
rejects (fun () -> Matrix(2,3).EigenGeneral() |> ignore)
for invalid in [Double.NaN;Double.PositiveInfinity;Double.NegativeInfinity] do
    rejects (fun () -> Matrix.FromArray(1,1,[|invalid|]).EigenGeneral() |> ignore)
for limit in [0;-1;100001] do rejects (fun () -> generalRotation.EigenGeneral(maxIterations=limit) |> ignore)
rejects (fun () -> Matrix.FromArray(3,3,[|1.;2.;3.;4.;5.;6.;7.;8.;10.|]).EigenGeneral(maxIterations=1) |> ignore)
printfn "F# general eigenvalue and complex right-eigenvector checks passed"
