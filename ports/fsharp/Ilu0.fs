namespace LinearA
open System
/// Owned ILU(0), with no pivoting, fill or diagonal shifts.
type ILU0(a:CSRMatrix) =
    let n,rp,ci,v = a.Rows,a.RowOffsets,a.ColumnIndices,a.Values
    let d = Array.zeroCreate<int> n
    let finite x = if not(Double.IsFinite x) then raise(ArithmeticException("Nonfinite ILU0 arithmetic."))
    let find i j =
        let mutable lo,hi = rp[i],rp[i+1]
        while lo<hi do
            let m = lo+(hi-lo)/2
            if ci[m]<j then lo<-m+1 else hi<-m
        lo
    do
        if a.Cols<>n then invalidArg "a" "ILU0 requires square matrix."
        for i in 0..n-1 do
            d[i]<-find i i
            if d[i]=rp[i+1] || ci[d[i]]<>i then invalidArg "a" "ILU0 requires stored diagonal."
        for i in 0..n-1 do
            for p in rp[i]..d[i]-1 do
                let j=ci[p]
                v[p]<-v[p]/v[d[j]]
                finite v[p]
                for q in d[j]+1..rp[j+1]-1 do
                    let k=find i ci[q]
                    if k<rp[i+1] && ci[k]=ci[q] then
                        v[k]<-v[k]-v[p]*v[q]
                        finite v[k]
            if v[d[i]]=0. then raise(ArithmeticException("Zero ILU0 pivot."))
    member _.Size=n
    member _.NNZ=v.Length
    member _.Apply(b:double[]) =
        if b.Length<>n || Array.exists (Double.IsFinite >> not) b then invalidArg "b" "Invalid ILU0 vector."
        let x=Array.copy b
        for i in 0..n-1 do
            for p in rp[i]..d[i]-1 do x[i]<-x[i]-v[p]*x[ci[p]]
            finite x[i]
        for i in n-1.. -1..0 do
            for p in d[i]+1..rp[i+1]-1 do x[i]<-x[i]-v[p]*x[ci[p]]
            x[i]<-x[i]/v[d[i]]
            finite x[i]
        x
