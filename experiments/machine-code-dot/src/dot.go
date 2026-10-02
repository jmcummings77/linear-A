package main
import ("fmt"; "os"; "strconv"; "time"; "math")
var sink float64
//go:noinline
func dot(a, b []float64) float64 {
    sum := 0.0
    for i := range a { sum += float64(a[i] * b[i]) }
    return sum
}
func main() {
    if len(os.Args) != 4 { panic("expected n iterations seed") }
    args := make([]int, 3)
    for i := range args { v,e := strconv.Atoi(os.Args[i+1]); if e != nil { panic(e) }; args[i] = v }
    n, iterations, seed := args[0], args[1], args[2]
    if n < 0 || n > 1048576 || iterations < 1 || iterations > 100000000 || seed < 0 || seed > 1000000 { panic("out of range") }
    count := n*32; if count == 0 { count=32 }
    a, b := make([]float64,count), make([]float64,count)
    for i := range a { a[i] = float64((i*17+seed*13)%101-50)/16; b[i] = float64((i*29+(seed+1)*7)%103-51)/16 }
    warm := iterations; if warm > 128 { warm=128 }; if warm < 8 { warm=8 }
    for i:=0; i<warm; i++ { off := (i%32)*n; sink = dot(a[off:off+n],b[off:off+n]) }
    checksum := 0.0; start := time.Now()
    for i:=0; i<iterations; i++ { off := (i%32)*n; checksum += dot(a[off:off+n],b[off:off+n]) }
    elapsed := time.Since(start).Nanoseconds()
    if math.IsNaN(checksum) || math.IsInf(checksum,0) { panic("nonfinite checksum") }
    fmt.Printf("{\"elapsed_ns\":%d,\"iterations\":%d,\"checksum\":%.17g}\n",elapsed,iterations,checksum)
}
