module LinearA.SparseRunner
open System
open System.Diagnostics
open System.Globalization
open System.Text.Json
open LinearA
let run (args: string[]) =
    if args.Length<>(if args.Length>1 && args[1]="gmres" then 12 else 11) then invalidArg "args" "Invalid sparse protocol."
    let integer i = Int32.Parse(args[i],CultureInfo.InvariantCulture)
    let number i = Double.Parse(args[i],CultureInfo.InvariantCulture)
    let op,rows,cols,nnz,iterations,rtol,atol,limit,jacobi,capture = args[1],integer 2,integer 3,integer 4,integer 5,number 6,number 7,integer 8,integer 9,integer 10
    if not (List.contains op ["spmv";"dense";"cg";"gmres";"ilu_setup";"ilu_apply";"rcm";"permute";"permutation_check";"rcm_solve";"ilu_solve"]) || min rows (min cols (min nnz iterations))<0 || jacobi<0 || jacobi>(if op="gmres" then 3 else 1) || capture<0 || capture>1 then invalidArg "args" "Invalid sparse options."
    let raw = Console.In.ReadToEnd().Split(Array.empty<char>,StringSplitOptions.RemoveEmptyEntries)
    let count = if op="cg" || op="gmres" then rows else cols
    if raw.Length<>rows+1+2*nnz+count then invalidArg "args" "Incorrect sparse input count."
    let mutable offset = 0
    let take n parse = Array.init n (fun _ -> let value = parse raw[offset] in offset <- offset+1; value)
    let rp = take (rows+1) (fun x -> Int32.Parse(x,CultureInfo.InvariantCulture))
    let ci = take nnz (fun x -> Int32.Parse(x,CultureInfo.InvariantCulture))
    let v = take nnz (fun x -> Double.Parse(x,CultureInfo.InvariantCulture))
    let b = take count (fun x -> Double.Parse(x,CultureInfo.InvariantCulture))
    let a = CSRMatrix(rows,cols,rp,ci,v)
    let factor = if op="ilu_apply" || (op="gmres" && jacobi=2) then Some(ILU0(a)) else None
    let dense = if op="dense" then Matrix(rows,cols) else Matrix(0,0)
    let right = if op="dense" then Matrix.FromArray(cols,1,b) else Matrix(0,0)
    if op="dense" then
        for i in 0..rows-1 do
            for p in rp[i]..rp[i+1]-1 do dense[i,ci[p]] <- v[p]
    let compute () =
        if op="rcm_solve" || op="ilu_solve" then
            let p=if op="rcm_solve" then a.ReverseCuthillMcKee() else [|0..rows-1|]
            let q=if op="rcm_solve" then a.PermuteSymmetric(p) else a
            let rhs=if op="rcm_solve" then CSRMatrix.PermuteVector(p,b) else b
            let result=q.Gmres(rhs,restart=20,relativeTolerance=rtol,absoluteTolerance=atol,maxIterations=limit,preconditioner=ILU0(q))
            if not result.Converged then failwith ("Ordering solve failed: "+result.Reason)
            let x=if op="rcm_solve" then CSRMatrix.PermuteVector(p,result.X,true) else result.X
            Array.append [|float result.Iterations|] x
        elif op="rcm" then a.ReverseCuthillMcKee() |> Array.map float
        elif op="permute" || op="permutation_check" then
            if Array.exists (fun x -> not (Double.IsFinite x) || x<0.0 || x>=float rows || x<>Math.Floor(x)) b then invalidArg "b" "Invalid permutation."
            let p=Array.map int b
            let q=a.PermuteSymmetric(p)
            if op="permute" then Array.concat [Array.map float q.RowOffsets;Array.map float q.ColumnIndices;q.Values]
            else
                let x=Array.init rows (fun i -> float (i+1))
                let y=CSRMatrix.PermuteVector(p,x)
                Array.concat [y;CSRMatrix.PermuteVector(p,y,true);CSRMatrix.PermuteVector(p,q.Matvec(y),true)]
        elif op="ilu_setup" then
            let f=ILU0(a)
            [|float f.Size;float f.NNZ|]
        elif op="ilu_apply" then factor.Value.Apply(b)
        elif op="spmv" then a.Matvec(b)
        elif op="dense" then dense.Multiply(right).ToArray()
        elif op="gmres" then
            let r = a.Gmres(b,restart=integer 11,relativeTolerance=rtol,absoluteTolerance=atol,maxIterations=limit,jacobi=(jacobi=1),capture=(capture<>0),?preconditioner=(if jacobi=3 then Some(ILU0(a)) else factor))
            if iterations>0 then
                if not r.Converged then failwith ("Benchmark GMRES did not converge: "+r.Reason)
                r.X
            else
                let reason = [|"converged";"iteration_limit";"breakdown";"nonfinite";"stagnation"|] |> Array.findIndex ((=) r.Reason)
                Array.concat [[|float reason;float r.Iterations;float r.Residuals.Length|];r.X;r.Residuals;r.EstimatedResiduals;[|float r.Restarts.Length|];Array.map float r.Restarts;Array.concat r.Iterates]
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
