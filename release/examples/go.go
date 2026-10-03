package main
import (
 "fmt"
 "math"
 matrix "github.com/jmcummings77/linear-A/ports/go"
)
func main() {
 a, err := matrix.New(2,2,[]float64{4,1,2,3}); if err != nil { panic(err) }
 b, err := matrix.New(2,1,[]float64{6,8}); if err != nil { panic(err) }
 x, err := a.Solve(b); if err != nil { panic(err) }
 r, err := a.Multiply(x); if err != nil { panic(err) }
 for i,v := range x.Values() { if math.IsNaN(v) || math.Abs(v-float64(i+1))>1e-12 || math.IsNaN(r.Values()[i]) || math.Abs(r.Values()[i]-b.Values()[i])>1e-12 { panic("incorrect solution") } }
 decomposition, err := a.SVD(); if err != nil { panic(err) }; if len(decomposition.Values)!=2 || decomposition.Values[1]<=0 { panic("incorrect SVD") }
 inverse,err:=a.Pseudoinverse();if err!=nil{panic(err)}
 minimum,err:=a.SolveMinimumNorm(b);if err!=nil{panic(err)}
 diagnostics,err:=a.SpectralDiagnostics();if err!=nil{panic(err)}
 if math.Abs(inverse.Values()[0]-.3)>1e-12 || math.Abs(minimum.Values()[1]-2)>1e-12 || diagnostics.Rank!=2{panic("incorrect SVD inverse")}
 ridge,err:=a.SolveRidge(b,1);if err!=nil{panic(err)}
 if math.Abs(ridge.Values()[0]-140.0/131)>1e-12 || math.Abs(ridge.Values()[1]-230.0/131)>1e-12{panic("incorrect ridge solution")}
 sparse,err:=matrix.NewCSR(2,2,[]int{0,2,4},[]int{0,1,0,1},[]float64{4,1,1,3});if err!=nil{panic(err)}
 product,err:=sparse.Matvec([]float64{1,2});if err!=nil{panic(err)}
 options:=matrix.DefaultCGOptions();options.Jacobi=true;options.Capture=true
 goOpts:=matrix.DefaultGMRESOptions();goOpts.Restart=2;goOpts.Jacobi=true;goOpts.Capture=true
 gm,err:=sparse.GMRES([]float64{6,7},goOpts);if err!=nil{panic(err)}
 if !gm.Converged||math.Abs(gm.X[0]-1)>1e-12||math.Abs(gm.X[1]-2)>1e-12{panic("incorrect GMRES")}
 cg,err:=sparse.ConjugateGradient([]float64{6,7},options);if err!=nil{panic(err)}
 if product[0]!=6||product[1]!=7||!cg.Converged||math.Abs(cg.X[0]-1)>1e-12||math.Abs(cg.X[1]-2)>1e-12{panic("incorrect sparse solver")}
 ilu,err:=matrix.NewILU0(sparse);if err!=nil{panic(err)}
 first,err:=ilu.Apply([]float64{6,7});if err!=nil{panic(err)}
 second,err:=ilu.Apply([]float64{11,13});if err!=nil{panic(err)}
 if math.Abs(first[0]-1)>1e-12||math.Abs(second[0]-20.0/11)>1e-12{panic("incorrect ILU reuse")}
 iluOptions:=matrix.DefaultGMRESOptions();iluOptions.Preconditioner=ilu
 solved,err:=sparse.GMRES([]float64{11,13},iluOptions);if err!=nil||!solved.Converged{panic("incorrect ILU GMRES")}
 fmt.Println("solution: 1, 2")
}
