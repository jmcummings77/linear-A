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
 fmt.Println("solution: 1, 2")
}
