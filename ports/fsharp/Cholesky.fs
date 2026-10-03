namespace LinearA
open System
open System.Collections.Generic

/// Owned lower factor supporting repeated right-hand sides.
type SparseCholesky internal (lower:CSRMatrix) =
    member _.Size=lower.Rows
    member _.NNZ=lower.NNZ
    member _.Lower=CSRMatrix(lower.Rows,lower.Cols,lower.RowOffsets,lower.ColumnIndices,lower.Values)
    member _.Solve(b:double[]) =
        let n,rp,ci,v=lower.Rows,lower.RowOffsets,lower.ColumnIndices,lower.Values
        if b.Length<>n || Array.exists (Double.IsFinite >> not) b then invalidArg "b" "Invalid Cholesky right-hand side."
        let x=Array.copy b
        let finite z = if not(Double.IsFinite z) then raise(ArithmeticException("Nonfinite Cholesky solve."))
        for i in 0..n-1 do
            for p in rp[i]..rp[i+1]-2 do x[i]<-x[i]-v[p]*x[ci[p]]
            x[i]<-x[i]/v[rp[i+1]-1]
            finite x[i]
        for i in n-1.. -1..0 do
            x[i]<-x[i]/v[rp[i+1]-1]
            finite x[i]
            for p in rp[i]..rp[i+1]-2 do
                x[ci[p]]<-x[ci[p]]-v[p]*x[i]
                finite x[ci[p]]
        x

/// Reusable symbolic pattern. FillSteps uses -1 for original entries.
type SparseCholeskySymbolic(a:CSRMatrix) =
    let n,sourceRP,sourceCI=a.Rows,a.RowOffsets,a.ColumnIndices
    let rp=Array.zeroCreate<int> (n+1)
    let columns,births=ResizeArray<int>(),ResizeArray<int>()
    do
        if a.Cols<>n then invalidArg "a" "Cholesky requires square matrix."
        let g=Array.init n (fun _ -> SortedDictionary<int,int>())
        for i in 0..n-1 do
            for p in sourceRP[i]..sourceRP[i+1]-1 do
                let j=sourceCI[p]
                if i<>j then
                    g[i][j]<- -1
                    g[j][i]<- -1
        for k in 0..n-1 do
            let ns=g[k].Keys |> Seq.filter(fun j -> j>k) |> Seq.toArray
            for u in 0..ns.Length-1 do
                for w in 0..u-1 do
                    let i,j=ns[u],ns[w]
                    if not(g[i].ContainsKey j) then
                        g[i][j]<-k
                        g[j][i]<-k
        for i in 0..n-1 do
            for pair in g[i] do
                if pair.Key<i then
                    columns.Add(pair.Key)
                    births.Add(pair.Value)
            columns.Add(i)
            births.Add(-1)
            rp[i+1]<-columns.Count
    let ci,steps=columns.ToArray(),births.ToArray()
    member _.Size=n
    member _.NNZ=ci.Length
    member _.FillCount=steps |> Array.filter(fun k -> k>=0) |> Array.length
    member _.RowOffsets=Array.copy rp
    member _.ColumnIndices=Array.copy ci
    member _.FillSteps=Array.copy steps
    member _.Factorize(a:CSRMatrix) =
        if a.Rows<>n || a.Cols<>n || a.RowOffsets<>sourceRP || a.ColumnIndices<>sourceCI then invalidArg "a" "Cholesky symbolic pattern mismatch."
        let av=a.Values
        let rows=Array.init n (fun i ->
            let row=Dictionary<int,double>()
            for p in sourceRP[i]..sourceRP[i+1]-1 do row[sourceCI[p]]<-av[p]
            row)
        let get i j = match rows[i].TryGetValue(j) with true,v -> v | _ -> 0.
        for i in 0..n-1 do
            for pair in rows[i] do
                if pair.Value<>get pair.Key i then invalidArg "a" "Cholesky requires symmetric values."
        let v=Array.zeroCreate<double> ci.Length
        for i in 0..n-1 do
            for p in rp[i]..rp[i+1]-1 do
                let j=ci[p]
                let mutable s=get i j
                let mutable u,w=rp[i],rp[j]
                while u<p && w<rp[j+1]-1 do
                    if ci[u]=ci[w] then
                        s<-s-v[u]*v[w]
                        u<-u+1
                        w<-w+1
                    elif ci[u]<ci[w] then u<-u+1
                    else w<-w+1
                if not(Double.IsFinite s) then raise(ArithmeticException("Nonfinite Cholesky factor."))
                if i=j then
                    if s<=0. then raise(ArithmeticException("Nonpositive Cholesky pivot."))
                    v[p]<-sqrt s
                else v[p]<-s/v[rp[j+1]-1]
                if not(Double.IsFinite v[p]) then raise(ArithmeticException("Nonfinite Cholesky factor."))
        SparseCholesky(CSRMatrix(n,n,rp,ci,v))
