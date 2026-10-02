namespace LinearA

open System
open System.Numerics

// Hessenberg reduction and real double-shift QR adapted from public-domain NIST
// JAMA (Martin/Wilkinson and EISPACK). https://math.nist.gov/javanumerics/jama/
module internal GeneralEigen =
    let private hessenberg (h: double[][]) (v: double[][]) =
        let n = h.Length
        let ort = Array.zeroCreate<double> n
        for m in 1 .. n-2 do
            let mutable scale = 0.0
            for i in m .. n-1 do scale <- scale + abs h.[i].[m-1]
            if scale <> 0.0 then
                let mutable norm = 0.0
                for i in n-1 .. -1 .. m do
                    ort.[i] <- h.[i].[m-1] / scale
                    norm <- norm + ort.[i]*ort.[i]
                let g = if ort.[m] > 0.0 then -sqrt norm else sqrt norm
                norm <- norm - ort.[m]*g
                ort.[m] <- ort.[m] - g
                for j in m .. n-1 do
                    let mutable f = 0.0
                    for i in n-1 .. -1 .. m do f <- f + ort.[i]*h.[i].[j]
                    f <- f/norm
                    for i in m .. n-1 do h.[i].[j] <- h.[i].[j] - f*ort.[i]
                for i in 0 .. n-1 do
                    let mutable f = 0.0
                    for j in n-1 .. -1 .. m do f <- f + ort.[j]*h.[i].[j]
                    f <- f/norm
                    for j in m .. n-1 do h.[i].[j] <- h.[i].[j] - f*ort.[j]
                ort.[m] <- scale*ort.[m]
                h.[m].[m-1] <- scale*g
        for i in 0 .. n-1 do v.[i].[i] <- 1.0
        for m in n-2 .. -1 .. 1 do
            if h.[m].[m-1] <> 0.0 then
                for i in m+1 .. n-1 do ort.[i] <- h.[i].[m-1]
                for j in m .. n-1 do
                    let mutable g = 0.0
                    for i in m .. n-1 do g <- g + ort.[i]*v.[i].[j]
                    g <- (g/ort.[m])/h.[m].[m-1]
                    for i in m .. n-1 do v.[i].[j] <- v.[i].[j] + g*ort.[i]

    let private divide xr xi yr yi =
        if abs yr > abs yi then
            let r = yi/yr
            let d = yr+r*yi
            (xr+r*xi)/d, (xi-r*xr)/d
        else
            let r = yr/yi
            let d = yi+r*yr
            (r*xr+xi)/d, (r*xi-xr)/d

    let private schur (h: double[][]) (v: double[][]) (d: double[]) (e: double[]) maxIterations =
        let nn = h.Length
        let eps = 2.2204460492503131e-16
        let mutable n = nn-1
        let mutable exshift = 0.0
        let mutable p = 0.0
        let mutable q = 0.0
        let mutable r = 0.0
        let mutable s = 0.0
        let mutable z = 0.0
        let mutable t = 0.0
        let mutable w = 0.0
        let mutable x = 0.0
        let mutable y = 0.0
        let mutable norm = 0.0
        for i in 0 .. nn-1 do
            for j in max (i-1) 0 .. nn-1 do norm <- norm + abs h.[i].[j]
        let mutable iter = 0
        while n >= 0 do
            let mutable l = n
            let mutable searching = true
            while l > 0 && searching do
                s <- abs h.[l-1].[l-1] + abs h.[l].[l]
                if s = 0.0 then s <- norm
                if abs h.[l].[l-1] <= eps*s then searching <- false else l <- l-1
            if l = n then
                h.[n].[n] <- h.[n].[n] + exshift
                d.[n] <- h.[n].[n]
                e.[n] <- 0.0
                n <- n-1
                iter <- 0
            elif l = n-1 then
                w <- h.[n].[n-1]*h.[n-1].[n]
                p <- (h.[n-1].[n-1]-h.[n].[n])/2.0
                q <- p*p+w
                z <- sqrt (abs q)
                h.[n].[n] <- h.[n].[n]+exshift
                h.[n-1].[n-1] <- h.[n-1].[n-1]+exshift
                x <- h.[n].[n]
                if q >= 0.0 then
                    z <- if p >= 0.0 then p+z else p-z
                    d.[n-1] <- x+z
                    d.[n] <- if z <> 0.0 then x-w/z else d.[n-1]
                    e.[n-1] <- 0.0
                    e.[n] <- 0.0
                    x <- h.[n].[n-1]
                    s <- abs x+abs z
                    p <- x/s
                    q <- z/s
                    r <- sqrt (p*p+q*q)
                    p <- p/r
                    q <- q/r
                    for j in n-1 .. nn-1 do
                        z <- h.[n-1].[j]
                        h.[n-1].[j] <- q*z+p*h.[n].[j]
                        h.[n].[j] <- q*h.[n].[j]-p*z
                    for i in 0 .. n do
                        z <- h.[i].[n-1]
                        h.[i].[n-1] <- q*z+p*h.[i].[n]
                        h.[i].[n] <- q*h.[i].[n]-p*z
                    for i in 0 .. nn-1 do
                        z <- v.[i].[n-1]
                        v.[i].[n-1] <- q*z+p*v.[i].[n]
                        v.[i].[n] <- q*v.[i].[n]-p*z
                else
                    d.[n-1] <- x+p
                    d.[n] <- x+p
                    e.[n-1] <- z
                    e.[n] <- -z
                n <- n-2
                iter <- 0
            else
                x <- h.[n].[n]
                y <- 0.0
                w <- 0.0
                if l < n then
                    y <- h.[n-1].[n-1]
                    w <- h.[n].[n-1]*h.[n-1].[n]
                if iter = 10 then
                    exshift <- exshift+x
                    for i in 0 .. n do h.[i].[i] <- h.[i].[i]-x
                    s <- abs h.[n].[n-1]+abs h.[n-1].[n-2]
                    x <- 0.75*s
                    y <- x
                    w <- -0.4375*s*s
                if iter = 30 then
                    s <- (y-x)/2.0
                    s <- s*s+w
                    if s > 0.0 then
                        s <- sqrt s
                        if y < x then s <- -s
                        s <- x-w/((y-x)/2.0+s)
                        for i in 0 .. n do h.[i].[i] <- h.[i].[i]-s
                        exshift <- exshift+s
                        x <- 0.964
                        y <- x
                        w <- x
                if iter >= maxIterations then invalidOp "General eigendecomposition did not converge within the iteration limit."
                iter <- iter+1
                let mutable m = n-2
                let mutable searching = true
                while m >= l && searching do
                    z <- h.[m].[m]
                    r <- x-z
                    s <- y-z
                    p <- (r*s-w)/h.[m+1].[m]+h.[m].[m+1]
                    q <- h.[m+1].[m+1]-z-r-s
                    r <- h.[m+2].[m+1]
                    s <- abs p+abs q+abs r
                    p <- p/s
                    q <- q/s
                    r <- r/s
                    if m = l || abs h.[m].[m-1]*(abs q+abs r) < eps*(abs p*(abs h.[m-1].[m-1]+abs z+abs h.[m+1].[m+1])) then
                        searching <- false
                    else m <- m-1
                for i in m+2 .. n do
                    h.[i].[i-2] <- 0.0
                    if i > m+2 then h.[i].[i-3] <- 0.0
                for k in m .. n-1 do
                    let notlast = k <> n-1
                    let mutable skip = false
                    if k <> m then
                        p <- h.[k].[k-1]
                        q <- h.[k+1].[k-1]
                        r <- if notlast then h.[k+2].[k-1] else 0.0
                        x <- abs p+abs q+abs r
                        if x = 0.0 then skip <- true
                        else
                            p <- p/x
                            q <- q/x
                            r <- r/x
                    if not skip then
                        s <- sqrt (p*p+q*q+r*r)
                        if p < 0.0 then s <- -s
                        if s <> 0.0 then
                            if k <> m then h.[k].[k-1] <- -s*x
                            elif l <> m then h.[k].[k-1] <- -h.[k].[k-1]
                            p <- p+s
                            x <- p/s
                            y <- q/s
                            z <- r/s
                            q <- q/p
                            r <- r/p
                            for j in k .. nn-1 do
                                p <- h.[k].[j]+q*h.[k+1].[j]
                                if notlast then
                                    p <- p+r*h.[k+2].[j]
                                    h.[k+2].[j] <- h.[k+2].[j]-p*z
                                h.[k].[j] <- h.[k].[j]-p*x
                                h.[k+1].[j] <- h.[k+1].[j]-p*y
                            for i in 0 .. min n (k+3) do
                                p <- x*h.[i].[k]+y*h.[i].[k+1]
                                if notlast then
                                    p <- p+z*h.[i].[k+2]
                                    h.[i].[k+2] <- h.[i].[k+2]-p*r
                                h.[i].[k] <- h.[i].[k]-p
                                h.[i].[k+1] <- h.[i].[k+1]-p*q
                            for i in 0 .. nn-1 do
                                p <- x*v.[i].[k]+y*v.[i].[k+1]
                                if notlast then
                                    p <- p+z*v.[i].[k+2]
                                    v.[i].[k+2] <- v.[i].[k+2]-p*r
                                v.[i].[k] <- v.[i].[k]-p
                                v.[i].[k+1] <- v.[i].[k+1]-p*q
        if norm <> 0.0 then
            for n in nn-1 .. -1 .. 0 do
                p <- d.[n]
                q <- e.[n]
                if q = 0.0 then
                    let mutable l = n
                    h.[n].[n] <- 1.0
                    for i in n-1 .. -1 .. 0 do
                        w <- h.[i].[i]-p
                        r <- 0.0
                        for j in l .. n do r <- r+h.[i].[j]*h.[j].[n]
                        if e.[i] < 0.0 then
                            z <- w
                            s <- r
                        else
                            l <- i
                            if e.[i] = 0.0 then h.[i].[n] <- -r/(if w <> 0.0 then w else eps*norm)
                            else
                                x <- h.[i].[i+1]
                                y <- h.[i+1].[i]
                                q <- (d.[i]-p)*(d.[i]-p)+e.[i]*e.[i]
                                t <- (x*s-z*r)/q
                                h.[i].[n] <- t
                                h.[i+1].[n] <- if abs x > abs z then (-r-w*t)/x else (-s-y*t)/z
                            t <- abs h.[i].[n]
                            if (eps*t)*t > 1.0 then
                                for j in i .. n do h.[j].[n] <- h.[j].[n]/t
                elif q < 0.0 then
                    let mutable l = n-1
                    if abs h.[n].[n-1] > abs h.[n-1].[n] then
                        h.[n-1].[n-1] <- q/h.[n].[n-1]
                        h.[n-1].[n] <- -(h.[n].[n]-p)/h.[n].[n-1]
                    else
                        let re,im = divide 0.0 (-h.[n-1].[n]) (h.[n-1].[n-1]-p) q
                        h.[n-1].[n-1] <- re
                        h.[n-1].[n] <- im
                    h.[n].[n-1] <- 0.0
                    h.[n].[n] <- 1.0
                    for i in n-2 .. -1 .. 0 do
                        let mutable ra = 0.0
                        let mutable sa = 0.0
                        for j in l .. n do
                            ra <- ra+h.[i].[j]*h.[j].[n-1]
                            sa <- sa+h.[i].[j]*h.[j].[n]
                        w <- h.[i].[i]-p
                        if e.[i] < 0.0 then
                            z <- w
                            r <- ra
                            s <- sa
                        else
                            l <- i
                            if e.[i] = 0.0 then
                                let re,im = divide (-ra) (-sa) w q
                                h.[i].[n-1] <- re
                                h.[i].[n] <- im
                            else
                                x <- h.[i].[i+1]
                                y <- h.[i+1].[i]
                                let mutable vr = (d.[i]-p)*(d.[i]-p)+e.[i]*e.[i]-q*q
                                let vi = (d.[i]-p)*2.0*q
                                if vr = 0.0 && vi = 0.0 then vr <- eps*norm*(abs w+abs q+abs x+abs y+abs z)
                                let re,im = divide (x*r-z*ra+q*sa) (x*s-z*sa-q*ra) vr vi
                                h.[i].[n-1] <- re
                                h.[i].[n] <- im
                                if abs x > abs z+abs q then
                                    h.[i+1].[n-1] <- (-ra-w*h.[i].[n-1]+q*h.[i].[n])/x
                                    h.[i+1].[n] <- (-sa-w*h.[i].[n]-q*h.[i].[n-1])/x
                                else
                                    let re,im = divide (-r-y*h.[i].[n-1]) (-s-y*h.[i].[n]) z q
                                    h.[i+1].[n-1] <- re
                                    h.[i+1].[n] <- im
                            t <- max (abs h.[i].[n-1]) (abs h.[i].[n])
                            if (eps*t)*t > 1.0 then
                                for j in i .. n do
                                    h.[j].[n-1] <- h.[j].[n-1]/t
                                    h.[j].[n] <- h.[j].[n]/t
            for j in nn-1 .. -1 .. 0 do
                for i in 0 .. nn-1 do
                    z <- 0.0
                    for k in 0 .. j do z <- z+v.[i].[k]*h.[k].[j]
                    v.[i].[j] <- z

    let private logAdd a b =
        let maximum = max a b
        maximum + Math.Log2(2.0 ** (a-maximum) + 2.0 ** (b-maximum))

    let private balance (h: double[][]) =
        let n = h.Length
        let exponents = Array.zeroCreate<int> n
        let mutable changed = true
        let mutable sweep = 0
        while changed && sweep < 64 do
            changed <- false
            sweep <- sweep + 1
            for i in 0 .. n-1 do
                let mutable rowMax,colMax = 0.0,0.0
                for j in 0 .. n-1 do
                    if j <> i then
                        rowMax <- max rowMax (abs h.[i].[j])
                        colMax <- max colMax (abs h.[j].[i])
                if rowMax <> 0.0 && colMax <> 0.0 then
                    let mutable rowSum,colSum = 0.0,0.0
                    for j in 0 .. n-1 do
                        if j <> i then
                            rowSum <- rowSum + abs h.[i].[j]/rowMax
                            colSum <- colSum + abs h.[j].[i]/colMax
                    let rowLog,colLog = Math.Log2(rowMax)+Math.Log2(rowSum),Math.Log2(colMax)+Math.Log2(colSum)
                    let shift = int (max -512.0 (min 512.0 (floor ((rowLog-colLog)/2.0+0.5))))
                    if shift <> 0 && logAdd (rowLog-double shift) (colLog+double shift) < logAdd rowLog colLog + Math.Log2(0.95) then
                        let factor = 2.0 ** double shift
                        let mutable safe = true
                        for j in 0 .. n-1 do
                            if j <> i then
                                let row,col = h.[i].[j]/factor,h.[j].[i]*factor
                                if not(Double.IsFinite(row) && Double.IsFinite(col)) ||
                                   (h.[i].[j] <> 0.0 && (row = 0.0 || row*factor <> h.[i].[j])) ||
                                   (h.[j].[i] <> 0.0 && (col = 0.0 || col/factor <> h.[j].[i])) then safe <- false
                        if safe then
                            for j in 0 .. n-1 do
                                if j <> i then
                                    h.[i].[j] <- h.[i].[j]/factor
                                    h.[j].[i] <- h.[j].[i]*factor
                            exponents.[i] <- exponents.[i]+shift
                            changed <- true
        exponents

    let solve (source: double[,]) maxIterations =
        if maxIterations < 1 || maxIterations > 100000 then invalidArg "maxIterations" "Expected an iteration limit in 1..100000."
        let n = source.GetLength(0)
        if n <> source.GetLength(1) then invalidArg "source" "General eigendecomposition requires a square matrix."
        let mutable scale = 0.0
        let mutable upper = true
        let mutable lower = true
        for i in 0 .. n-1 do
            for j in 0 .. n-1 do
                let value = source.[i,j]
                if not (Double.IsFinite(value)) then invalidArg "source" "Expected finite entries."
                scale <- max scale (abs value)
                if i > j && value <> 0.0 then upper <- false
                if i < j && value <> 0.0 then lower <- false
        let reverse = lower && not upper
        let h = Array.init n (fun i -> Array.init n (fun j -> source.[(if reverse then n-1-i else i),(if reverse then n-1-j else j)]))
        let balanceExponents = if upper || lower then Array.zeroCreate<int> n else balance h
        scale <- 0.0
        for row in h do
            for value in row do scale <- max scale (abs value)
        let v = Array.init n (fun _ -> Array.zeroCreate<double> n)
        let d,e = Array.zeroCreate<double> n,Array.zeroCreate<double> n
        if upper && lower then
            for i in 0 .. n-1 do
                d.[i] <- source.[i,i]
                v.[i].[i] <- 1.0
            scale <- 1.0
        else
            let exponentBits = (BitConverter.DoubleToInt64Bits(scale) >>> 52) &&& 0x7ffL
            if exponentBits <> 0L then scale <- BitConverter.Int64BitsToDouble(exponentBits <<< 52)
            for i in 0 .. n-1 do
                for j in 0 .. n-1 do h.[i].[j] <- h.[i].[j]/scale
            hessenberg h v
            schur h v d e maxIterations
            if upper || lower then
                for i in 0 .. n-1 do
                    d.[i] <- source.[(if reverse then n-1-i else i),(if reverse then n-1-i else i)]
                    e.[i] <- 0.0
                scale <- 1.0
        let eigenvalues = Array.init n (fun i -> Complex(d.[i]*scale,e.[i]*scale))
        if eigenvalues |> Array.exists (fun value -> not(Double.IsFinite(value.Real) && Double.IsFinite(value.Imaginary))) then
            invalidOp "An eigenvalue exceeds the finite double range."
        let vectors = Array2D.zeroCreate<Complex> n n
        for col in 0 .. n-1 do
            let mutable exponent = Int32.MinValue
            let mutable maximum = 0.0
            let mutable largest = 0
            for row in 0 .. n-1 do
                let i = if reverse then n-1-row else row
                let re = if e.[col] < 0.0 then v.[i].[col-1] else v.[i].[col]
                let im = if e.[col] > 0.0 then v.[i].[col+1] elif e.[col] < 0.0 then -v.[i].[col] else 0.0
                if not(Double.IsFinite(re) && Double.IsFinite(im)) then invalidOp "Nonfinite eigenvector."
                vectors.[row,col] <- Complex(re,im)
                let entryMagnitude = max (abs re) (abs im)
                if entryMagnitude <> 0.0 then exponent <- max exponent (Math.ILogB(entryMagnitude)+balanceExponents.[row])
            if exponent = Int32.MinValue then invalidOp "Zero eigenvector."
            for row in 0 .. n-1 do
                let re,im = vectors.[row,col].Real,vectors.[row,col].Imaginary
                let shift = balanceExponents.[row]-exponent
                vectors.[row,col] <- Complex(Math.ScaleB(re,shift),Math.ScaleB(im,shift))
                let magnitude = vectors.[row,col].Magnitude
                if magnitude > maximum then
                    maximum <- magnitude
                    largest <- row
            if maximum = 0.0 then invalidOp "Zero eigenvector."
            let mutable norm = 0.0
            for row in 0 .. n-1 do
                vectors.[row,col] <- vectors.[row,col]/maximum
                norm <- Double.Hypot(norm,vectors.[row,col].Magnitude)
            let pivot = vectors.[largest,col]
            let phase = Complex.Conjugate(pivot)/pivot.Magnitude
            for row in 0 .. n-1 do vectors.[row,col] <- (vectors.[row,col]/norm)*phase
        let order = Array.init n id |> Array.sortBy (fun i -> eigenvalues.[i].Real,eigenvalues.[i].Imaginary,i)
        Array.map (fun i -> eigenvalues.[i]) order, Array2D.init n n (fun row col -> vectors.[row,order.[col]])
