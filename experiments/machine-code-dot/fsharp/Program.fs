open System
open System.Diagnostics
open System.Globalization
open System.Runtime.CompilerServices

[<MethodImpl(MethodImplOptions.NoInlining ||| MethodImplOptions.AggressiveOptimization)>]
let dot (a: float[]) (b: float[]) offset n =
    let mutable sum = 0.0
    for i = offset to offset + n - 1 do
        sum <- sum + a.[i] * b.[i]
    sum

[<EntryPoint>]
let main args =
    let n, iterations, seed = int args.[0], int args.[1], int args.[2]
    if n < 0 || n > 1048576 || iterations < 1 || iterations > 100000000 || seed < 0 || seed > 1000000 then
        invalidArg "args" "out of range"
    let count = max n 1 * 32
    let a = Array.init count (fun i -> float ((i * 17 + seed * 13) % 101 - 50) / 16.0)
    let b = Array.init count (fun i -> float ((i * 29 + (seed + 1) * 7) % 103 - 51) / 16.0)
    let mutable sink = 0.0
    for i = 0 to max 8 (min iterations 128) - 1 do
        sink <- dot a b ((i % 32) * n) n
    let mutable checksum = 0.0
    let start = Stopwatch.GetTimestamp()
    for i = 0 to iterations - 1 do
        checksum <- checksum + dot a b ((i % 32) * n) n
    let elapsed = Stopwatch.GetElapsedTime(start).TotalNanoseconds
    if not (Double.IsFinite checksum && Double.IsFinite sink) then failwith "nonfinite result"
    printfn "{\"elapsed_ns\":%s,\"iterations\":%d,\"checksum\":%s}" (elapsed.ToString("R", CultureInfo.InvariantCulture)) iterations (checksum.ToString("R", CultureInfo.InvariantCulture))
    0
