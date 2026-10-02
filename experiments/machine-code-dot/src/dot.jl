@noinline function dot_kernel(a::Vector{Float64}, b::Vector{Float64}, offset::Int, n::Int)
    total = 0.0
    @inbounds for i in (offset+1):(offset+n)
        total += a[i] * b[i]
    end
    total
end
function main()
    n, iterations, seed = parse.(Int, ARGS)
    @assert 0 <= n <= 1048576 && 1 <= iterations <= 100000000 && 0 <= seed <= 1000000
    count = max(n, 1) * 32
    a = [((i * 17 + seed * 13) % 101 - 50) / 16 for i in 0:count-1]
    b = [((i * 29 + (seed + 1) * 7) % 103 - 51) / 16 for i in 0:count-1]
    sink = 0.0
    for i in 0:max(8,min(iterations,128))-1
        sink += dot_kernel(a,b,(i%32)*n,n)
    end
    checksum = 0.0
    start = time_ns()
    for i in 0:iterations-1
        checksum += dot_kernel(a,b,(i%32)*n,n)
    end
    elapsed = time_ns()-start
    @assert isfinite(checksum) && isfinite(sink)
    println("{\"elapsed_ns\":",elapsed,",\"iterations\":",iterations,",\"checksum\":",checksum,"}")
end
main()
