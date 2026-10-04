namespace LinearA
open System
/// Symmetric V-cycle for the unit five-point Dirichlet Laplacian.
type GeometricMultigrid(width:int) =
    do if width<1 || width>255 || (width &&& (width+1))<>0 then invalidArg "width" "Grid width must be 2^k-1 in 1..255."
    let multiply w (x:double[]) = Array.init x.Length (fun i -> 4.*x[i]-(if i%w>0 then x[i-1] else 0.)-(if i%w+1<w then x[i+1] else 0.)-(if i>=w then x[i-w] else 0.)-(if i+w<x.Length then x[i+w] else 0.))
    let rec cycle w (b:double[]) =
        if w=1 then [|b[0]/4.|] else
            let x=Array.zeroCreate<double> b.Length
            let smooth () =
                for _ in 1..2 do
                    let ax=multiply w x
                    for i in 0..x.Length-1 do x[i]<-x[i]+(b[i]-ax[i])/6.
            smooth()
            let r=Array.map2 (-) b (multiply w x)
            let c=w/2
            let bc=Array.zeroCreate<double> (c*c)
            for y in 0..c-1 do
                for j in 0..c-1 do
                    for dy in -1..1 do
                        for dx in -1..1 do bc[y*c+j]<-bc[y*c+j]+(if dy=0 then 1. else 0.5)*(if dx=0 then 1. else 0.5)*r[(2*y+1+dy)*w+2*j+1+dx]
            let ec=cycle c bc
            for y in 0..c-1 do
                for j in 0..c-1 do
                    for dy in -1..1 do
                        for dx in -1..1 do x[(2*y+1+dy)*w+2*j+1+dx]<-x[(2*y+1+dy)*w+2*j+1+dx]+(if dy=0 then 1. else 0.5)*(if dx=0 then 1. else 0.5)*ec[y*c+j]
            smooth()
            if Array.exists (Double.IsFinite >> not) x then raise(ArithmeticException("Nonfinite multigrid cycle."))
            x
    member _.Width=width
    member _.Size=width*width
    member _.Levels=int(Math.Log2(float(width+1)))
    member _.Matrix=
        let rp,ci,v=ResizeArray<int>(),ResizeArray<int>(),ResizeArray<double>()
        rp.Add(0)
        for i in 0..width*width-1 do
            let add j=ci.Add(j);v.Add(if j=i then 4. else -1.)
            if i>=width then add(i-width)
            if i%width>0 then add(i-1)
            add i
            if i%width+1<width then add(i+1)
            if i+width<width*width then add(i+width)
            rp.Add(v.Count)
        CSRMatrix(width*width,width*width,rp.ToArray(),ci.ToArray(),v.ToArray())
    member this.Apply(b:double[])=
        if b.Length<>this.Size || Array.exists (Double.IsFinite >> not) b then invalidArg "b" "Invalid multigrid RHS."
        cycle width b
    interface ISymmetricPreconditioner with
        member this.Size=this.Size
        member this.Apply(b)=this.Apply(b)
