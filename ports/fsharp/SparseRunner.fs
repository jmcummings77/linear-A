module LinearA.SparseRunner
open System
open System.Diagnostics
open System.Globalization
open System.Text.Json
open LinearA
let run (args: string[]) =
    if args.Length<>11 then invalidArg "args" "Invalid sparse protocol."
    let integer i = Int32.Parse(args[i],CultureInfo.InvariantCulture)
    let number i = Double.Parse(args[i],CultureInfo.InvariantCulture)
    let op,rows,cols,nnz,iterations,rtol,atol,limit,jacobi,capture = args[1],integer 2,integer 3,integer 4,integer 5,number 6,number 7,integer 8,integer 9,integer 10
    if not (List.contains op ["spmv";"dense";"cg"]) || min rows (min cols (min nnz iterations))<0 || jacobi<0 || jacobi>1 || capture<0 || capture>1 then invalidArg "args" "Invalid sparse options."
    let raw = Console.In.ReadToEnd().Split(Array.empty<char>,StringSplitOptions.RemoveEmptyEntries)
    let count = if op="cg" then rows else cols
    if raw.Length<>rows+1+2*nnz+count then invalidArg "args" "Incorrect sparse input count."
    let mutable offset = 0
    let take n parse = Array.init n (fun _ -> let value = parse raw[offset] in offset <- offset+1; value)
    let rp = take (rows+1) (fun x -> Int32.Parse(x,CultureInfo.InvariantCulture))
    let ci = take nnz (fun x -> Int32.Parse(x,CultureInfo.InvariantCulture))
    let v = take nnz (fun x -> Double.Parse(x,CultureInfo.InvariantCulture))
    let b = take count (fun x -> Double.Parse(x,CultureInfo.InvariantCulture))
    let a = CSRMatrix(rows,cols,rp,ci,v)
    let dense = if op="dense" then Matrix(rows,cols) else Matrix(0,0)
    let right = if op="dense" then Matrix.FromArray(cols,1,b) else Matrix(0,0)
    if op="dense" then
        for i in 0..rows-1 do
            for p in rp[i]..rp[i+1]-1 do dense[i,ci[p]] <- v[p]
    let compute () =
        if op="spmv" then a.Matvec(b)
        elif op="dense" then dense.Multiply(right).ToArray()
        else
            let r = a.ConjugateGradient(b,relativeTolerance=rtol,absoluteTolerance=atol,maxIterations=limit,jacobi=(jacobi<>0),capture=(capture<>0))
            if iterations>0 then
                if not r.Converged then failwith ("Benchmark CG did not converge: "+r.Reason)
                r.X
            else
                let reason = [|"converged";"iteration_limit";"breakdown";"nonfinite"|] |> Array.findIndex ((=) r.Reason)
                Array.concat [ [|float reason;float r.Iterations;float r.Residuals.Length|];r.X;r.Residuals;Array.concat r.Iterates ]
    if iterations=0 then
        let values = compute()
        printfn "%s" (JsonSerializer.Serialize({|rows=1;cols=values.Length;values=values|}))
    else
        for _ in 1..3 do compute() |> ignore
        let mutable checksum = 0.0
        let timer = Stopwatch.StartNew()
        for _ in 1..iterations do
            for x in compute() do checksum <- checksum+x
        let elapsed = timer.Elapsed.TotalNanoseconds
        if not (Double.IsFinite checksum) then failwith "Nonfinite checksum."
        printfn "%s" (JsonSerializer.Serialize({|elapsed_ns=elapsed;iterations=iterations;checksum=checksum|}))
