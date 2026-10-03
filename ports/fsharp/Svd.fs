namespace LinearA

open System

/// Economy A = U diag(Values) Vt. All returned arrays and matrices are owned.
type SingularValueDecomposition = { U: Matrix; Values: double[]; Vt: Matrix }

type SpectralDiagnostics = { Rank: int; ReciprocalCondition: double; RetainedReciprocalCondition: double }

[<AutoOpen>]
module MatrixSvd =
    let rec private decompose (source: Matrix) tolerance maxSweeps =
        if not (Double.IsFinite tolerance) || tolerance <= 0.0 || tolerance >= 1.0 || maxSweeps < 1 || maxSweeps > 10000 then
            invalidArg "tolerance" "Invalid SVD options."
        let m,n = source.Rows,source.Cols
        if m < n then
            let r = decompose (source.Transpose()) tolerance maxSweeps
            { U=r.Vt.Transpose(); Values=r.Values; Vt=r.U.Transpose() }
        else
            let raw = source.ToArray()
            if Array.exists (Double.IsFinite >> not) raw then invalidArg "source" "SVD requires finite input."
            let largest = Array.fold (fun a x -> max a (abs x)) 0.0 raw
            let scaling = if largest=0.0 then 1.0 else largest
            let b = Array.map (fun x -> x/scaling) raw
            if Array.exists2 (fun x y -> x<>0.0 && y=0.0) raw b then raise (ArithmeticException("SVD scaling discards an entry."))
            let v = (Matrix.Identity(n)).ToArray()
            let norm j =
                let mutable sum = 0.0
                for i in 0..m-1 do sum <- Double.Hypot(sum,b[i*n+j])
                sum
            let mutable converged = false
            let mutable sweep = 0
            while sweep <= maxSweeps && not converged do
                let mutable changed = false
                for p in 0..n-1 do
                    for q in p+1..n-1 do
                        let np,nq = norm p,norm q
                        if np<>0.0 && nq<>0.0 then
                            let mutable corr = 0.0
                            for i in 0..m-1 do corr <- corr + (b[i*n+p]/np)*(b[i*n+q]/nq)
                            if abs corr > tolerance then
                                changed <- true
                                if sweep < maxSweeps then
                                    let pair = max np nq
                                    let ap,aq = np/pair,nq/pair
                                    let delta,g = aq*aq-ap*ap,2.0*ap*aq*corr
                                    let t = if delta=0.0 then Math.CopySign(1.0,g) else g/(delta+Math.CopySign(Double.Hypot(delta,g),delta))
                                    if abs t < 2.2250738585072014e-308 then
                                        let small = if np<nq then p else q
                                        for i in 0..m-1 do b[i*n+small] <- 0.0
                                    else
                                        let c = 1.0/Double.Hypot(1.0,t)
                                        let s = c*t
                                        let rotate (data: double[]) count =
                                            for i in 0..count-1 do
                                                let x,y = data[i*n+p],data[i*n+q]
                                                data[i*n+p] <- c*x-s*y
                                                data[i*n+q] <- s*x+c*y
                                        rotate b m
                                        rotate v n
                converged <- not changed
                sweep <- sweep+1
            if not converged then raise (ArithmeticException("SVD did not converge."))
            let norms = Array.init n norm
            let order = Array.init n id |> Array.sortBy (fun j -> -norms[j],j)
            let u,vt,values = Matrix(m,n),Matrix(n,n),Array.zeroCreate<double> n
            for j in 0..n-1 do
                let k = order[j]
                values[j] <- norms[k]*scaling
                if not (Double.IsFinite values[j]) || (norms[k]<>0.0 && values[j]=0.0) then raise (ArithmeticException("Singular value outside double range."))
                for i in 0..n-1 do vt[j,i] <- v[i*n+k]
                if norms[k]<>0.0 then
                    for i in 0..m-1 do u[i,j] <- b[i*n+k]/norms[k]
                else
                    let mutable found = false
                    let mutable axis = 0
                    while axis<m && not found do
                        let candidate = Array.zeroCreate<double> m
                        candidate[axis] <- 1.0
                        for _ in 1..2 do
                            for col in 0..j-1 do
                                let mutable dot = 0.0
                                for i in 0..m-1 do dot <- dot+candidate[i]*u[i,col]
                                for i in 0..m-1 do candidate[i] <- candidate[i]-dot*u[i,col]
                        let length = Array.fold (fun a x -> Double.Hypot(a,x)) 0.0 candidate
                        if length>0.5/sqrt(float m) then
                            for i in 0..m-1 do u[i,j] <- candidate[i]/length
                            found <- true
                        axis <- axis+1
                    if not found then raise (ArithmeticException("Cannot complete SVD null basis."))
            { U=u; Values=values; Vt=vt }


    let private cutoff (a: Matrix) value =
        let x = defaultArg value (float (max a.Rows a.Cols)*2.220446049250313e-16)
        if not (Double.IsFinite x) || x<0.0 || x>1.0 then invalidArg "relativeCutoff" "Invalid relative cutoff."
        x
    let private rank (s: double[]) cutoff =
        s |> Array.takeWhile (fun x -> x>0.0 && (cutoff=0.0 || x/s[0]>cutoff)) |> Array.length
    let private inverseProduct a b c =
        if a=0.0 || c=0.0 then 0.0
        else
            let ae,be,ce = Math.ILogB(a),Math.ILogB(b),Math.ILogB(c)
            let x = Math.ScaleB(Math.ScaleB(a,-ae)/Math.ScaleB(b,-be)*Math.ScaleB(c,-ce),ae-be+ce)
            if not (Double.IsFinite x) then raise (ArithmeticException("SVD inverse outside float64 range."))
            x
    let private applyInverse (a: Matrix) (rhs: Matrix option) relativeCutoff =
        let cutoff = cutoff a relativeCutoff
        if rhs |> Option.exists (fun b -> b.Rows<>a.Rows) then invalidArg "rhs" "Incompatible right-hand side."
        if rhs |> Option.exists (fun b -> b.ToArray() |> Array.exists (Double.IsFinite >> not)) then invalidArg "rhs" "Expected finite right-hand side."
        let r = decompose a 1e-12 100
        let rank = rank r.Values cutoff
        let m,n = a.Rows,a.Cols
        let cols = rhs |> Option.map (fun b -> b.Cols) |> Option.defaultValue m
        let output = Matrix(n,cols)
        for j in 0..cols-1 do
            let scaling =
                match rhs with
                | None -> 1.0
                | Some b -> [0..m-1] |> List.fold (fun s i -> max s (abs b[i,j])) 0.0
            if rank>0 && scaling<>0.0 && (rhs |> Option.exists (fun b -> [0..m-1] |> List.exists (fun i -> b[i,j]<>0.0 && b[i,j]/scaling=0.0))) then
                raise (ArithmeticException("Right-hand side scaling discards an entry."))
            for p in 0..rank-1 do
                let projection =
                    match rhs with
                    | None -> r.U[j,p]
                    | Some b when scaling<>0.0 -> [0..m-1] |> List.sumBy (fun i -> r.U[i,p]*(b[i,j]/scaling))
                    | _ -> 0.0
                let coefficient = if rhs.IsSome then inverseProduct projection r.Values[p] scaling else 0.0
                for i in 0..n-1 do
                    let term = if rhs.IsSome then r.Vt[p,i]*coefficient else inverseProduct projection r.Values[p] r.Vt[p,i]
                    let value = output[i,j]+term
                    if not (Double.IsFinite value) then raise (ArithmeticException("SVD inverse outside float64 range."))
                    output[i,j] <- value
        output

    type Matrix with
        /// One-sided Jacobi SVD; see ports/SVD.md for range and convergence rules.
        member this.Svd(?tolerance: double, ?maxSweeps: int) =
            decompose this (defaultArg tolerance 1e-12) (defaultArg maxSweeps 100)

        member this.Pseudoinverse(?relativeCutoff: double) = applyInverse this None relativeCutoff
        member this.SolveMinimumNorm(rhs: Matrix, ?relativeCutoff: double) = applyInverse this (Some rhs) relativeCutoff
        member this.SpectralDiagnostics(?relativeCutoff: double) =
            let cutoff = cutoff this relativeCutoff
            let s = (this.Svd()).Values
            let rank = rank s cutoff
            { Rank=rank; ReciprocalCondition=(if s.Length>0 && s[0]>0.0 then s[s.Length-1]/s[0] else 0.0);
              RetainedReciprocalCondition=(if rank>0 then s[rank-1]/s[0] else 0.0) }
