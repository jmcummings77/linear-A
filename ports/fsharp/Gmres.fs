namespace LinearA
open System

type GMRESResult = { X:double[]; Converged:bool; Iterations:int; Reason:string; Residuals:double[]; EstimatedResiduals:double[]; Iterates:double[][]; Restarts:int[] }
exception internal GmresStop of string

[<AutoOpen>]
module GmresExtensions =
    type CSRMatrix with
        member a.Gmres(b:double[], ?restart:int, ?relativeTolerance:double, ?absoluteTolerance:double, ?maxIterations:int, ?jacobi:bool, ?capture:bool) =
            let restart,rtol,atol,limit,jacobi,capture = defaultArg restart 30,defaultArg relativeTolerance 1e-10,defaultArg absoluteTolerance 0.,defaultArg maxIterations 1000,defaultArg jacobi false,defaultArg capture false
            let n = a.Rows
            if a.Cols<>n || b.Length<>n || Array.exists (Double.IsFinite >> not) b then invalidArg "b" "GMRES requires square matrix and finite matching vector."
            if restart<1 || restart>1024 || limit<0 || limit>100000 || not(Double.IsFinite rtol) || rtol<0. || rtol>=1. || not(Double.IsFinite atol) || atol<0. then invalidArg "restart" "Invalid GMRES options."
            let diagonal = Array.create n 1.
            if jacobi then
                let rp,ci,v = a.RowOffsets,a.ColumnIndices,a.Values
                for i in 0..n-1 do
                    let mutable found = false
                    for p in rp[i]..rp[i+1]-1 do
                        if ci[p]=i then
                            diagonal[i] <- v[p]
                            found <- true
                    if not found || diagonal[i]=0. then invalidArg "jacobi" "Jacobi requires a nonzero diagonal."
            let norm v = Array.fold (fun s x -> Double.Hypot(s,x)) 0. v
            let finite v = Array.forall Double.IsFinite v
            let mutable x,r = Array.zeroCreate<double> n,Array.copy b
            let history,estimates,frames,restarts = ResizeArray<double>(),ResizeArray<double>(),ResizeArray<double[]>(),ResizeArray<int>()
            history.Add(norm r)
            estimates.Add(history[0])
            if capture then frames.Add(Array.copy x)
            let threshold,m = max atol (rtol*history[0]),min restart (min n limit)
            let stop reason = raise(GmresStop reason)
            let product v =
                try a.Matvec(v) with :? ArithmeticException -> stop "nonfinite" | :? ArgumentException -> stop "nonfinite"
            let mutable reason = "iteration_limit"
            try
                if not(Double.IsFinite history[0]) then stop "nonfinite"
                if history[0]<=threshold then stop "converged"
                while history.Count-1<limit do
                    if history.Count>1 then restarts.Add(history.Count-1)
                    let initial,beta = Array.copy x,norm r
                    let basis = ResizeArray<double[]>()
                    basis.Add(Array.map (fun v -> v/beta) r)
                    let h,cs,sn,g = Array2D.zeroCreate<double> (m+1) m,Array.zeroCreate<double> m,Array.zeroCreate<double> m,Array.zeroCreate<double> (m+1)
                    g[0] <- beta
                    let steps = min m (limit-(history.Count-1))
                    for j in 0..steps-1 do
                        let w = product (Array.mapi (fun i v -> v/diagonal[i]) basis[j])
                        let original = norm w
                        for _ in 1..2 do
                            for k in 0..j do
                                let mutable dot = 0.
                                for i in 0..n-1 do dot <- dot+basis[k][i]*w[i]
                                h[k,j] <- h[k,j]+dot
                                for i in 0..n-1 do w[i] <- w[i]-dot*basis[k][i]
                        let tail = norm w
                        if not(Double.IsFinite original && Double.IsFinite tail) then stop "nonfinite"
                        for k in 0..j do
                            if not(Double.IsFinite h[k,j]) then stop "nonfinite"
                        let happy = tail<=8.*2.220446049250313e-16*original
                        h[j+1,j] <- if happy then 0. else tail
                        if not happy then basis.Add(Array.map (fun v -> v/tail) w)
                        for k in 0..j-1 do
                            let top = cs[k]*h[k,j]+sn[k]*h[k+1,j]
                            h[k+1,j] <- -sn[k]*h[k,j]+cs[k]*h[k+1,j]
                            h[k,j] <- top
                        let pivot = Double.Hypot(h[j,j],h[j+1,j])
                        if not(Double.IsFinite pivot) then stop "nonfinite"
                        if pivot=0. then stop "breakdown"
                        cs[j] <- h[j,j]/pivot
                        sn[j] <- h[j+1,j]/pivot
                        h[j,j] <- pivot
                        h[j+1,j] <- 0.
                        g[j+1] <- -sn[j]*g[j]
                        g[j] <- cs[j]*g[j]
                        let y = g[0..j]
                        for k in j.. -1 ..0 do
                            if h[k,k]=0. then stop "breakdown"
                            let mutable sum = 0.
                            for q in k+1..j do sum <- sum+h[k,q]*y[q]
                            y[k] <- (y[k]-sum)/h[k,k]
                        let candidate = Array.init n (fun i ->
                            let mutable sum = 0.
                            for k in 0..j do sum <- sum+basis[k][i]*y[k]
                            initial[i]+sum/diagonal[i])
                        if not(finite y && finite candidate && Double.IsFinite g[j+1]) then stop "nonfinite"
                        let ax = product candidate
                        let residual = Array.mapi (fun i v -> v-ax[i]) b
                        let length = norm residual
                        if not(Double.IsFinite length) then stop "nonfinite"
                        x <- candidate
                        r <- residual
                        history.Add(length)
                        estimates.Add(abs g[j+1])
                        if capture then frames.Add(Array.copy x)
                        if length<=threshold then stop "converged"
                        if happy then stop "breakdown"
                    if x=initial then stop "stagnation"
            with GmresStop value -> reason <- value
            { X=Array.copy x; Converged=reason="converged"; Iterations=history.Count-1; Reason=reason; Residuals=history.ToArray(); EstimatedResiduals=estimates.ToArray(); Iterates=frames.ToArray(); Restarts=restarts.ToArray() }
