namespace LinearA
open System

type CGResult = { X: double[]; Converged: bool; Iterations: int; Reason: string; Residuals: double[]; Iterates: double[][] }

/// Canonical, owned CSR with zero-based sorted unique column indices.
type CSRMatrix(rows: int, cols: int, offsets: int[], indices: int[], values: double[]) =
    let rp,ci,v = Array.copy offsets,Array.copy indices,Array.copy values
    do
        if rows<0 || rows=Int32.MaxValue || cols<0 || rp.Length<>rows+1 || ci.Length<>v.Length || rp[0]<>0 || rp[rows]<>v.Length then invalidArg "offsets" "Invalid CSR dimensions or arrays."
        if Array.exists (fun x -> x<0 || x>v.Length) rp || Array.exists (Double.IsFinite >> not) v then invalidArg "values" "Invalid CSR offsets or values."
        for i in 0..rows-1 do
            if rp[i]>rp[i+1] then invalidArg "offsets" "CSR offsets must be monotone."
            let mutable previous = -1
            for p in rp[i]..rp[i+1]-1 do
                if ci[p]<=previous || ci[p]>=cols then invalidArg "indices" "CSR columns must be sorted, unique and in range."
                previous <- ci[p]
    let find row col =
        let mutable lo,hi = rp[row],rp[row+1]
        while lo<hi do
            let mid = lo+(hi-lo)/2
            if ci[mid]<col then lo <- mid+1 else hi <- mid
        lo
    member _.Rows = rows
    member _.Cols = cols
    member _.NNZ = v.Length
    member _.RowOffsets = Array.copy rp
    member _.ColumnIndices = Array.copy ci
    member _.Values = Array.copy v
    member _.ReverseCuthillMcKee() =
        if rows<>cols then invalidArg "rows" "RCM requires square matrix."
        let graph = Array.init rows (fun _ -> Collections.Generic.HashSet<int>())
        for i in 0..rows-1 do
            for k in rp[i]..rp[i+1]-1 do
                let j=ci[k]
                if i<>j then
                    graph[i].Add(j) |> ignore
                    graph[j].Add(i) |> ignore
        let key i = graph[i].Count,i
        let seen = Array.create rows false
        let order = ResizeArray<int>()
        for start in Array.sortBy key [|0..rows-1|] do
            if not seen[start] then
                let queue=ResizeArray<int>()
                queue.Add(start)
                seen[start] <- true
                let mutable h=0
                while h<queue.Count do
                    for j in graph[queue[h]] |> Seq.filter (fun j -> not seen[j]) |> Seq.sortBy key |> Seq.toArray do
                        seen[j] <- true
                        queue.Add(j)
                    h <- h+1
                order.AddRange(queue)
        order.ToArray() |> Array.rev
    static member PermuteVector(order:int[],x:double[],?inverse:bool) =
        if order.Length<>x.Length then invalidArg "order" "Invalid permutation length."
        let seen=Array.create x.Length false
        for i in 0..x.Length-1 do
            let j=order[i]
            if j<0 || j>=x.Length || seen[j] || not (Double.IsFinite x[i]) then invalidArg "order" "Invalid permutation or vector."
            seen[j] <- true
        let result=Array.zeroCreate<double> x.Length
        for i in 0..x.Length-1 do
            if defaultArg inverse false then result[order[i]] <- x[i] else result[i] <- x[order[i]]
        result
    member _.PermuteSymmetric(order:int[]) =
        if rows<>cols then invalidArg "rows" "Permutation requires square matrix."
        let inv=CSRMatrix.PermuteVector(order,Array.init rows float,true)
        let offsets,indices,entries=ResizeArray<int>(),ResizeArray<int>(),ResizeArray<double>()
        offsets.Add(0)
        for i in order do
            for k in [|rp[i]..rp[i+1]-1|] |> Array.sortBy (fun k -> inv[ci[k]]) do
                indices.Add(int inv[ci[k]])
                entries.Add(v[k])
            offsets.Add(entries.Count)
        CSRMatrix(rows,cols,offsets.ToArray(),indices.ToArray(),entries.ToArray())
    static member FromDense(a: Matrix) =
        let rp,ci,v = ResizeArray<int>(),ResizeArray<int>(),ResizeArray<double>()
        rp.Add(0)
        for i in 0..a.Rows-1 do
            for j in 0..a.Cols-1 do
                if a[i,j]<>0.0 then
                    ci.Add(j)
                    v.Add(a[i,j])
            rp.Add(v.Count)
        CSRMatrix(a.Rows,a.Cols,rp.ToArray(),ci.ToArray(),v.ToArray())
    member _.Matvec(x: double[]) =
        if x.Length<>cols || Array.exists (Double.IsFinite >> not) x then invalidArg "x" "Invalid vector."
        Array.init rows (fun i ->
            let mutable sum = 0.0
            for p in rp[i]..rp[i+1]-1 do sum <- sum+v[p]*x[ci[p]]
            if not (Double.IsFinite sum) then raise (ArithmeticException("Sparse multiplication outside float64 range."))
            sum)
    member this.ConjugateGradient(b: double[], ?relativeTolerance: double, ?absoluteTolerance: double, ?maxIterations: int, ?jacobi: bool, ?capture: bool) =
        let rtol,atol,limit,jacobi,capture = defaultArg relativeTolerance 1e-10,defaultArg absoluteTolerance 0.0,defaultArg maxIterations 1000,defaultArg jacobi false,defaultArg capture false
        let n = rows
        if cols<>n || b.Length<>n || Array.exists (Double.IsFinite >> not) b then invalidArg "b" "CG requires square matrix and finite matching vector."
        if not (Double.IsFinite rtol) || rtol<0.0 || rtol>=1.0 || not (Double.IsFinite atol) || atol<0.0 || limit<0 || limit>100000 then invalidArg "maxIterations" "Invalid CG options."
        let diagonal = Array.create n 1.0
        for i in 0..n-1 do
            for p in rp[i]..rp[i+1]-1 do
                let j = ci[p]
                let q = find j i
                let other = if q<rp[j+1] && ci[q]=i then v[q] else 0.0
                if v[p]<>other then invalidArg "values" "CG requires exact symmetry."
            if jacobi then
                let q = find i i
                if q=rp[i+1] || ci[q]<>i || v[q]<=0.0 then invalidArg "values" "Jacobi requires positive diagonal."
                diagonal[i] <- v[q]
        let norm a = Array.fold (fun s x -> Double.Hypot(s,x)) 0.0 a
        let dot a b = Array.fold2 (fun s x y -> s+x*y) 0.0 a b
        let mutable x = Array.zeroCreate<double> n
        let mutable r = Array.copy b
        let history,frames = ResizeArray<double>(),ResizeArray<double[]>()
        history.Add(norm r)
        if capture then frames.Add(Array.copy x)
        let threshold = max atol (rtol*history[0])
        let mutable reason = if not (Double.IsFinite history[0]) then "nonfinite" elif history[0]<=threshold then "converged" else "iteration_limit"
        let mutable z = Array.mapi (fun i a -> a/diagonal[i]) r
        let mutable p = Array.copy z
        let mutable rho = dot r z
        let mutable step = 0
        while reason="iteration_limit" && step<limit do
            try
                let q = this.Matvec(p)
                let curvature = dot p q
                if not (Double.IsFinite rho) || not (Double.IsFinite curvature) then reason <- "nonfinite"
                elif rho<=0.0 || curvature<=0.0 then reason <- "breakdown"
                else
                    let alpha = rho/curvature
                    let candidate = Array.mapi (fun i a -> a+alpha*p[i]) x
                    let ax = this.Matvec(candidate)
                    let residual = Array.map2 (-) b ax
                    let length = norm residual
                    if not (Double.IsFinite length) then reason <- "nonfinite"
                    else
                        x <- candidate
                        r <- residual
                        history.Add(length)
                        if capture then frames.Add(Array.copy x)
                        if length<=threshold then reason <- "converged"
                        else
                            z <- Array.mapi (fun i a -> a/diagonal[i]) r
                            let next = dot r z
                            if not (Double.IsFinite next) then reason <- "nonfinite"
                            elif next<=0.0 then reason <- "breakdown"
                            else
                                let beta = next/rho
                                p <- Array.mapi (fun i a -> a+beta*p[i]) z
                                rho <- next
            with
            | :? ArithmeticException -> reason <- "nonfinite"
            | :? ArgumentException -> reason <- "nonfinite"
            step <- step+1
        { X=Array.copy x; Converged=reason="converged"; Iterations=history.Count-1; Reason=reason; Residuals=history.ToArray(); Iterates=frames.ToArray() }
