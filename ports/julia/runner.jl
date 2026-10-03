include(joinpath(@__DIR__, "src", "LinearAMatrices.jl"))
using .LinearAMatrices

const OPERATIONS = Set(["pseudoinverse", "spectral_diagnostics", "solve_ridge", "solve_minimum_norm", "svd", "svd_one_sweep","add", "subtract", "scale", "transpose", "multiply", "solve", "solve_cholesky", "least_squares", "rcond", "cross", "rotation2d", "rotation3d", "trace", "eigen_symmetric", "eigen_general", "determinant", "determinant_lu", "determinant_spd_lu", "determinant_cholesky", "triangular"])

function integer(text, minimum=0, maximum=typemax(Int))
    value = parse(Int, text)
    minimum <= value <= maximum || throw(ArgumentError("integer argument is outside the allowed range"))
    value
end

function operation(name, a, b=nothing, scalar=1.25)
    name == "add" && return a + b
    name == "subtract" && return a - b
    name == "multiply" && return a * b
    name == "pseudoinverse" && return pseudoinverse(a;relative_cutoff=scalar)
    name == "solve_ridge" && return solve_ridge(a,b,scalar)
    name == "solve_minimum_norm" && return solve_minimum_norm(a,b)
    if name == "spectral_diagnostics"
        d=spectral_diagnostics(a;relative_cutoff=scalar)
        return Matrix64(1,3,[d.rank,d.reciprocal_condition,d.retained_reciprocal_condition])
    end
    if name in ("svd","svd_one_sweep")
        r=svd(a;max_sweeps=name=="svd" ? 100 : 1)
        return Matrix64(a.rows+1+a.cols,min(a.rows,a.cols),vcat(rowmajor(r.u),r.values,rowmajor(transpose(r.vt))))
    end
    name == "solve" && return solve(a,b)
    name == "solve_cholesky" && return solve(factor_cholesky(a),b)
    name == "least_squares" && return least_squares(a,b)
    name == "rcond" && return reciprocal_condition(factor_lu(a))
    name == "cross" && return cross(a, b)
    if name == "rotation2d"
        size(a) == (0,0) || throw(ArgumentError("rotation2d expects an empty 0-by-0 input"))
        return rotation2d(scalar)
    end
    name == "rotation3d" && return rotation_axis_angle(a, scalar)
    name == "scale" && return scale(a, scalar)
    name == "transpose" && return transpose(a)
    name == "trace" && return trace(a)
    name == "eigen_symmetric" && return eigen_symmetric(a)
    name == "eigen_general" && return eigen_general(a)
    name == "determinant" && return determinant(a)
    name in ("determinant_lu", "determinant_spd_lu") && return determinant(a, :lu)
    name == "determinant_cholesky" && return determinant(a, :cholesky)
    name == "triangular" && return triangular(a)
    throw(ArgumentError("unsupported operation"))
end

function check(args)
    (length(args) >= 3 && args[1] in OPERATIONS) || throw(ArgumentError("usage: check OP ROWS COLS [BROWS BCOLS | SCALAR]"))
    name, rows, cols = args[1], integer(args[2]), integer(args[3])
    binary = name in ("add", "subtract", "multiply", "cross", "solve", "solve_cholesky", "least_squares", "solve_minimum_norm", "solve_ridge")
    has_scalar = name in ("scale", "rotation2d", "rotation3d", "pseudoinverse", "spectral_diagnostics")
    length(args) == (name == "solve_ridge" ? 6 : binary ? 5 : has_scalar ? 4 : 3) || throw(ArgumentError("incorrect number of operation arguments"))
    brows, bcols = binary ? (integer(args[4]), integer(args[5])) : (0, 0)
    scalar = name == "solve_ridge" ? finite(parse(Float64,args[6])) : has_scalar ? finite(parse(Float64, args[4])) : 1.25
    values = Float64[finite(parse(Float64, token)) for token in split(read(stdin, String))]
    count = Base.Checked.checked_mul(rows, cols)
    expected = Base.Checked.checked_add(count, Base.Checked.checked_mul(brows, bcols))
    length(values) == expected || throw(ArgumentError("input value count does not match matrix dimensions"))
    a = Matrix64(rows, cols, values[1:count])
    b = binary ? Matrix64(brows, bcols, values[count + 1:end]) : nothing
    operation(name, a, b, scalar)
end

function preparebenchmark(args)
    (length(args) == 4 && args[1] in OPERATIONS && args[1] != "triangular") || throw(ArgumentError("usage: bench OP SIZE ITERATIONS SEED"))
    name, size, iterations, seed = args[1], integer(args[2], 1), integer(args[3], 1), integer(args[4], 0, 2147483646)
    ((name in ("cross", "rotation3d") && size != 3) || (name == "rotation2d" && size != 2)) &&
        throw(ArgumentError("cross and rotation3d require size 3; rotation2d requires size 2"))
    a, b = Matrix64(size, size), Matrix64(size, size)
    for index in 0:length(a.data) - 1
        a.data[index + 1] = (mod(index * 17 + seed * 13, 101) - 50) / 16.0
        b.data[index + 1] = (mod(index * 17 + (seed + 1) * 13, 101) - 50) / 16.0
    end
    if name == "cross"
        a, b = Matrix64(3,1,a.data[1:3]), Matrix64(3,1,b.data[1:3])
        b[3,1] = -b[3,1]
    elseif name == "rotation2d"
        a, b = Matrix64(0,0), Matrix64(0,0)
    elseif name == "rotation3d"
        a, b = Matrix64(3,1,[1.,2.,3.]), Matrix64(0,0)
    end
    if name in ("determinant_cholesky", "determinant_spd_lu")
        for row in 1:size, col in 1:row
            value = (a[row, col] + a[col, row]) / 2.0
            a[row, col] = a[col, row] = value
        end
    end
    if startswith(name, "determinant")
        for i in 1:size
            a[i, i] += size * 4
        end
    end
    if name == "eigen_general"
        for i in 0:size-1,j in 0:size-1
            k=i÷2;a0=1+(seed%17)/16+k/8;b0=.5+k/16
            a[i+1,j+1]=i==j ? a0 : iseven(i)&&j==i+1 ? -b0 : isodd(i)&&j==i-1 ? b0 : i<j ? (mod(i*3+j*5+seed,11)-5)/32 : 0.0
        end
    end
    if name == "eigen_symmetric"
        for row in 1:size, col in 1:size
            a[row,col] = row == col ? 2.0 + (seed % 17)/16.0 : abs(row-col) == 1 ? -1.0 : 0.0
        end
    end
    name, a, b, iterations
end

function runiterations(name, a, b, iterations)
    total = 0.0
    for _ in 1:iterations
        result = operation(name, a, b, name in ("rotation2d", "rotation3d") ? 0.5 : 1.25)
        if name == "eigen_general"
            total += sum(i*(real(v)+abs(imag(v))) for (i,v) in enumerate(result.values)) + sum(abs2,result.vectors)
        elseif name == "eigen_symmetric"
            value = 0.0
            for i in eachindex(result.values)
                value += i * result.values[i]
            end
            for entry in result.vectors.data
                value += entry * entry
            end
            total += value
        else
            total += result isa Matrix64 ? checksum(result) : result
        end
    end
    finite(total)
end

function benchmark(args)
    name, a, b, iterations = preparebenchmark(args)
    runiterations(name, a, b, max(5, min(iterations, 100)))
    start = time_ns()
    total = runiterations(name, a, b, iterations)
    elapsed = time_ns() - start
    (elapsed_ns=elapsed, iterations=iterations, checksum=total)
end

function writejson(io, result::Matrix64)
    print(io, "{\"rows\":", result.rows, ",\"cols\":", result.cols, ",\"values\":[")
    for (i, value) in enumerate(result.data)
        i > 1 && print(io, ',')
        print(io, finite(value))
    end
    print(io, "]}")
end

writejson(io, result::Real) = print(io, "{\"value\":", finite(result), "}")
function writejson(io, result::NamedTuple)
    if haskey(result,:vectors) && result.vectors isa Matrix{ComplexF64}
        print(io,"{\"eigenvalues_real\":[")
        for (i,v) in enumerate(result.values);i>1&&print(io,',');print(io,finite(real(v)));end
        print(io,"],\"eigenvalues_imag\":[")
        for (i,v) in enumerate(result.values);i>1&&print(io,',');print(io,finite(imag(v)));end
        n=length(result.values)
        print(io,"],\"eigenvectors_real\":")
        writejson(io,Matrix64(n,n,[real(result.vectors[i,j]) for i in 1:n for j in 1:n]))
        print(io,",\"eigenvectors_imag\":")
        writejson(io,Matrix64(n,n,[imag(result.vectors[i,j]) for i in 1:n for j in 1:n]))
        print(io,'}')
        return
    end
    if haskey(result, :vectors)
        print(io, "{\"eigenvalues\":[")
        for (i, value) in enumerate(result.values)
            i > 1 && print(io, ',')
            print(io, finite(value))
        end
        print(io, "],\"eigenvectors\":")
        writejson(io, result.vectors)
        print(io, '}')
        return
    end
    print(io, '{')
    for (i, key) in enumerate(keys(result))
        i > 1 && print(io, ',')
        print(io, '"', key, "\":", result[key])
    end
    print(io, '}')
end

function main(args=ARGS)
    (!isempty(args) && args[1] in ("check", "bench")) || throw(ArgumentError("usage: runner.jl check|bench ..."))
    result = args[1] == "check" ? check(args[2:end]) : benchmark(args[2:end])
    writejson(stdout, result)
    println()
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main()
    catch error
        showerror(stderr, error)
        println(stderr)
        exit(1)
    end
end
