module LinearAMatrices

export Factorization, factor_lu, factor_cholesky, factor_qr, solve, least_squares, reciprocal_condition,
       eigen_general, eigen_symmetric, cross, rotation2d, rotationx, rotationy, rotationz, rotation_axis_angle,
       Matrix64, identitymatrix, row, column, scale, trace, determinant, triangular, rowmajor, checksum, finite

function finite(value::Real)
    result = Float64(value)
    isfinite(result) || throw(ArgumentError("matrix arithmetic requires finite float64 values"))
    result
end

"""A mutable matrix's values in owned row-major storage; public indices are one-based."""
struct Matrix64
    rows::Int
    cols::Int
    data::Vector{Float64}

    function Matrix64(rows::Integer=0, cols::Integer=0, values=nothing)
        (rows isa Bool || cols isa Bool || rows < 0 || cols < 0) && throw(ArgumentError("dimensions must be nonnegative integers"))
        count = Base.Checked.checked_mul(Int(rows), Int(cols))
        data = values === nothing ? zeros(Float64, count) : Float64[finite(value) for value in values]
        length(data) == count || throw(ArgumentError("value count does not match dimensions"))
        new(Int(rows), Int(cols), data)
    end
end

Base.size(a::Matrix64) = (a.rows, a.cols)
Base.copy(a::Matrix64) = Matrix64(a.rows, a.cols, a.data)
rowmajor(a::Matrix64) = copy(a.data)

function identitymatrix(size::Integer)
    result = Matrix64(size, size)
    for i in 1:size
        result.data[(i - 1) * size + i] = 1.0
    end
    result
end

function index(a::Matrix64, row::Integer, col::Integer)
    (1 <= row <= a.rows && 1 <= col <= a.cols) || throw(BoundsError(a, (row, col)))
    (row - 1) * a.cols + col
end

Base.getindex(a::Matrix64, row::Integer, col::Integer) = a.data[index(a, row, col)]
function Base.setindex!(a::Matrix64, value::Real, row::Integer, col::Integer)
    position = index(a, row, col)
    a.data[position] = finite(value)
end

function row(a::Matrix64, rowindex::Integer)
    1 <= rowindex <= a.rows || throw(BoundsError(a, rowindex))
    a.data[(rowindex - 1) * a.cols + 1:rowindex * a.cols]
end

function column(a::Matrix64, colindex::Integer)
    1 <= colindex <= a.cols || throw(BoundsError(a, colindex))
    Float64[a.data[(i - 1) * a.cols + colindex] for i in 1:a.rows]
end

function vector3(a::Matrix64)
    ((a.rows == 3 && a.cols == 1) || (a.rows == 1 && a.cols == 3)) ||
        throw(DimensionMismatch("expected a 3-by-1 or 1-by-3 vector"))
    finite.(a.data)
end

"""Right-handed vector cross product; the result retains the left vector's shape."""
function cross(a::Matrix64, b::Matrix64)
    x, y = vector3(a), vector3(b)
    Matrix64(a.rows, a.cols, [x[2]*y[3]-x[3]*y[2], x[3]*y[1]-x[1]*y[3], x[1]*y[2]-x[2]*y[1]])
end

"""Active counterclockwise 2D column-vector rotation, with the angle in radians."""
function rotation2d(radians::Real)
    angle = finite(radians)
    c, s = cos(angle), sin(angle)
    Matrix64(2,2,[c,-s,s,c])
end

function rotationx(radians::Real)
    angle = finite(radians)
    c, s = cos(angle), sin(angle)
    Matrix64(3,3,[1,0,0, 0,c,-s, 0,s,c])
end

function rotationy(radians::Real)
    angle = finite(radians)
    c, s = cos(angle), sin(angle)
    Matrix64(3,3,[c,0,s, 0,1,0, -s,0,c])
end

function rotationz(radians::Real)
    angle = finite(radians)
    c, s = cos(angle), sin(angle)
    Matrix64(3,3,[c,-s,0, s,c,0, 0,0,1])
end

"""Right-handed active rotation around a finite, nonzero axis, in radians."""
function rotation_axis_angle(axis::Matrix64, radians::Real)
    v, angle = vector3(axis), finite(radians)
    largest = maximum(abs, v)
    largest > 0 || throw(ArgumentError("rotation axis must be nonzero"))
    v ./= largest
    v ./= sqrt(sum(x*x for x in v))
    x,y,z = v
    c,s = cos(angle), sin(angle)
    t = abs(angle) < 1 ? 2*sin(angle/2)^2 : 1-c
    Matrix64(3,3,[c+x*x*t, x*y*t-z*s, x*z*t+y*s,
                  y*x*t+z*s, c+y*y*t, y*z*t-x*s,
                  z*x*t-y*s, z*y*t+x*s, c+z*z*t])
end

function combine(a::Matrix64, b::Matrix64, subtract::Bool)
    size(a) == size(b) || throw(DimensionMismatch("matrix dimensions must match"))
    result = Matrix64(a.rows, a.cols)
    for i in eachindex(a.data)
        result.data[i] = finite(subtract ? a.data[i] - b.data[i] : a.data[i] + b.data[i])
    end
    result
end

Base.:+(a::Matrix64, b::Matrix64) = combine(a, b, false)
Base.:-(a::Matrix64, b::Matrix64) = combine(a, b, true)

function scale(a::Matrix64, scalar::Real)
    scalar = finite(scalar)
    result = Matrix64(a.rows, a.cols)
    for i in eachindex(a.data)
        result.data[i] = finite(a.data[i] * scalar)
    end
    result
end

Base.:*(a::Matrix64, scalar::Real) = scale(a, scalar)
Base.:*(scalar::Real, a::Matrix64) = scale(a, scalar)

function Base.transpose(a::Matrix64)
    result = Matrix64(a.cols, a.rows)
    for row in 1:a.rows, col in 1:a.cols
        result.data[(col - 1) * a.rows + row] = a.data[(row - 1) * a.cols + col]
    end
    result
end

function Base.:*(a::Matrix64, b::Matrix64)
    a.cols == b.rows || throw(DimensionMismatch("left columns must equal right rows"))
    result = Matrix64(a.rows, b.cols)
    for row in 1:a.rows, col in 1:b.cols
        total = 0.0
        for k in 1:a.cols
            total += a.data[(row - 1) * a.cols + k] * b.data[(k - 1) * b.cols + col]
        end
        result.data[(row - 1) * b.cols + col] = finite(total)
    end
    result
end

function square(a::Matrix64)
    a.rows == a.cols || throw(DimensionMismatch("operation requires a square matrix"))
end

function trace(a::Matrix64)
    square(a)
    result = 0.0
    for i in 1:a.rows
        result += a.data[(i - 1) * a.cols + i]
    end
    finite(result)
end

# Keep the diagonal product in a normalized mantissa and binary exponent.
function diagonalproduct(values, sign=1.0, repeats=1, exponent=0)
    mantissa = sign
    zero = false
    for value in values
        if value == 0.0
            zero = true
            continue
        end
        fraction, power = frexp(value)
        for _ in 1:repeats
            mantissa, shift = frexp(mantissa * fraction)
            exponent += power + shift
        end
    end
    zero ? 0.0 : finite(ldexp(mantissa, exponent))
end

"""Auto shortcuts triangular inputs; :lu uses partial pivoting; :cholesky requires SPD."""
function determinant(a::Matrix64, algorithm::Symbol=:auto)
    square(a)
    algorithm in (:auto, :lu, :cholesky) || throw(ArgumentError("expected auto, lu, or cholesky"))
    all(isfinite, a.data) || throw(ArgumentError("determinant requires finite entries"))
    n = a.rows
    if algorithm == :auto
        shape = triangular(a)
        (shape.upper || shape.lower) && return diagonalproduct((a[i, i] for i in 1:n))
    elseif algorithm == :cholesky
        for row in 1:n, col in 1:row-1
            a[row, col] == a[col, row] || throw(ArgumentError("Cholesky requires exact symmetry"))
        end
        lower = zeros(Float64, n*n)
        for row in 1:n, col in 1:row
            value = a[row, col]
            for k in 1:col-1
                value -= lower[(row-1)*n+k] * lower[(col-1)*n+k]
            end
            finite(value)
            if row == col
                value > 0.0 || throw(ArgumentError("Cholesky requires positive definite input"))
                lower[(row-1)*n+col] = sqrt(value)
            else
                lower[(row-1)*n+col] = finite(value / lower[(col-1)*n+col])
            end
        end
        return diagonalproduct((lower[(i-1)*n+i] for i in 1:n), 1.0, 2)
    end
    work = copy(a.data)
    rowexponent = 0
    # Preserve every entry exactly; never scale a tiny entry out of the matrix.
    for row in 1:n
        largest = maximum(abs, @view(work[(row-1)*n+1:row*n]); init=0.0)
        if largest > 0.0
            _, power = frexp(largest)
            if all(ldexp(ldexp(work[(row-1)*n+col], -power), power) == work[(row-1)*n+col] for col in 1:n)
                for col in 1:n
                    work[(row-1)*n+col] = ldexp(work[(row-1)*n+col], -power)
                end
                rowexponent += power
            end
        end
    end
    sign = 1.0
    for k in 1:n
        pivotrow = k
        for row in k + 1:n
            if abs(work[(row - 1) * n + k]) > abs(work[(pivotrow - 1) * n + k])
                pivotrow = row
            end
        end
        pivot = work[(pivotrow - 1) * n + k]
        pivot == 0.0 && return 0.0
        if pivotrow != k
            for col in k:n
                left, right = (k - 1) * n + col, (pivotrow - 1) * n + col
                work[left], work[right] = work[right], work[left]
            end
            sign = -sign
        end
        for row in k + 1:n
            numerator = work[(row - 1) * n + k]
            factor = numerator / pivot
            work[(row - 1) * n + k] = 0.0
            for col in k + 1:n
                index = (row - 1) * n + col
                entry = work[(k - 1) * n + col]
                update = factor * entry
                if numerator != 0.0 && entry != 0.0 && abs(factor) < floatmin(Float64)
                    fn, en = frexp(numerator)
                    fp, ep = frexp(pivot)
                    fe, ee = frexp(entry)
                    update = ldexp((fn / fp) * fe, en - ep + ee)
                end
                work[index] = finite(work[index] - update)
            end
        end
    end
    diagonalproduct((work[(i-1)*n+i] for i in 1:n), sign, 1, rowexponent)
end

"""Real eigenpairs of exactly symmetric input; eigenvectors are columns.
Tolerance is relative to the Frobenius norm of the whole matrix.
"""
function eigen_symmetric(a::Matrix64; tolerance=1e-12, max_sweeps=50)
    (tolerance isa Real && !(tolerance isa Bool) && isfinite(tolerance) && 0 < tolerance < 1) ||
        throw(ArgumentError("expected a finite tolerance strictly between zero and one"))
    (max_sweeps isa Integer && !(max_sweeps isa Bool) && max_sweeps > 0) ||
        throw(ArgumentError("expected a positive sweep limit"))
    square(a)
    n = a.rows
    largest, diagonal = 0.0, true
    for row in 1:n, col in 1:n
        value = finite(a[row,col])
        value == a[col,row] || throw(ArgumentError("expected exact symmetry"))
        diagonal &= row == col || value == 0.0
        largest = max(largest, abs(value))
    end
    work, vectors = copy(a), identitymatrix(n)
    exponent = diagonal || largest == 0.0 ? 0 : last(frexp(largest))
    if !diagonal
        norm = 0.0
        for i in eachindex(work.data)
            work.data[i] = ldexp(work.data[i], -exponent)
            norm = hypot(norm, work.data[i])
        end
        limit = tolerance * norm
        skip = limit / (2*n)
        function offnorm()
            total = 0.0
            for row in 1:n, col in row+1:n
                total = hypot(total, work[row,col])
            end
            sqrt(2.0) * total
        end
        sweeps = 0
        while offnorm() > limit && sweeps < max_sweeps
            for p in 1:n-1, q in p+1:n
                b = work[p,q]
                if abs(b) > skip
                    delta = (work[q,q] - work[p,p]) / 2
                    t = delta == 0.0 ? copysign(1.0, b) : b / (delta + copysign(hypot(delta, b), delta))
                    c = 1 / sqrt(1 + t*t)
                    s = t*c
                    work[p,p] -= t*b
                    work[q,q] += t*b
                    work[p,q] = work[q,p] = 0.0
                    for k in 1:n
                        if k != p && k != q
                            x, y = work[k,p], work[k,q]
                            work[k,p] = work[p,k] = c*x - s*y
                            work[k,q] = work[q,k] = s*x + c*y
                        end
                        x, y = vectors[k,p], vectors[k,q]
                        vectors[k,p], vectors[k,q] = c*x - s*y, s*x + c*y
                    end
                end
            end
            sweeps += 1
        end
        offnorm() <= limit || throw(ErrorException("symmetric eigen decomposition did not converge within the sweep limit"))
    end
    values = [finite(ldexp(work[i,i], exponent)) for i in 1:n]
    order = sortperm(values)
    sorted = Matrix64(n,n)
    for col in 1:n
        source, pivot, norm = order[col], 1, 0.0
        for row in 1:n
            norm = hypot(norm, vectors[row,source])
            abs(vectors[row,source]) > abs(vectors[pivot,source]) && (pivot = row)
        end
        divisor = vectors[pivot,source] < 0.0 ? -norm : norm
        for row in 1:n
            sorted[row,col] = vectors[row,source] / divisor
        end
    end
    (values=values[order], vectors=sorted)
end

function triangular(a::Matrix64)
    a.rows == a.cols || return (upper=false, lower=false)
    upper, lower = true, true
    for row in 1:a.rows, col in 1:a.cols
        if a.data[(row - 1) * a.cols + col] != 0.0
            if row > col
                upper = false
            elseif row < col
                lower = false
            end
        end
    end
    (upper=upper, lower=lower)
end

function checksum(a::Matrix64)
    n = length(a.data)
    n == 0 ? 0.0 : finite(a.data[1] + a.data[n ÷ 2 + 1] + a.data[n])
end

include("eigen_general.jl")
include("solve.jl")

include("svd.jl")
export solve_ridge, svd, pseudoinverse, solve_minimum_norm, spectral_diagnostics

include("sparse.jl")
include("gmres.jl")
export gmres

end
