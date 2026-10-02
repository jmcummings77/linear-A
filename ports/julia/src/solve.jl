"""Reusable factor snapshot. QR requires numerical full column rank and rows ≥ columns."""
struct Factorization
    rows::Int
    cols::Int
    algorithm::Symbol
    scale::Float64
    norm::Float64
    data::Vector{Float64}
    tau::Vector{Float64}
    permutation::Vector{Int}
end

function solver_scaled(value, scale)
    result = finite(value / scale)
    value != 0 && result == 0 && throw(ArgumentError("solver scaling would discard a nonzero value"))
    result
end

function Factorization(source::Matrix64, algorithm::Symbol)
    algorithm in (:lu, :cholesky, :qr) || throw(ArgumentError("unknown factorization"))
    m, n = source.rows, source.cols
    (m < n || (algorithm != :qr && m != n)) && throw(DimensionMismatch("factorization requires square input, or rows ≥ columns for QR"))
    original = rowmajor(source)
    maximum_value = isempty(original) ? 0.0 : maximum(abs, original)
    scaling = maximum_value == 0 ? 1.0 : ldexp(1.0, last(frexp(maximum_value)) - 1)
    a = [solver_scaled(value, scaling) for value in original]
    norm = 0.0
    for i in 1:m
        row_sum = 0.0
        for j in 1:n; row_sum += abs(a[(i-1)*n+j]); end
        norm = max(norm, row_sum)
    end
    p, tau = collect(1:m), zeros(n)
    at(i,j) = (i-1)*n+j
    if algorithm == :lu
        for k in 1:n
            pivot = k
            for i in k+1:n
                abs(a[at(i,k)]) > abs(a[at(pivot,k)]) && (pivot = i)
            end
            a[at(pivot,k)] == 0 && throw(ArgumentError("singular matrix: zero computed LU pivot"))
            for j in 1:n; a[at(k,j)], a[at(pivot,j)] = a[at(pivot,j)], a[at(k,j)]; end
            p[k], p[pivot] = p[pivot], p[k]
            for i in k+1:n
                a[at(i,k)] = finite(a[at(i,k)] / a[at(k,k)])
                for j in k+1:n; a[at(i,j)] = finite(a[at(i,j)] - a[at(i,k)] * a[at(k,j)]); end
            end
        end
    elseif algorithm == :cholesky
        for i in 1:n, j in 1:i-1
            original[at(i,j)] == original[at(j,i)] || throw(ArgumentError("Cholesky requires exact symmetry"))
        end
        for i in 1:n, j in 1:i
            value = a[at(i,j)]
            for k in 1:j-1; value = finite(value - a[at(i,k)] * a[at(j,k)]); end
            if i == j
                value > 0 || throw(ArgumentError("Cholesky requires positive computed pivots"))
                a[at(i,j)] = sqrt(value)
            else
                a[at(i,j)] = finite(value / a[at(j,j)])
            end
        end
    else
        largest = 0.0
        for j in 1:n
            length = 0.0
            for i in 1:m; length = hypot(length, a[at(i,j)]); end
            largest = max(largest, length)
        end
        threshold = eps(Float64) * m * largest
        for k in 1:n
            pivot, length = k, 0.0
            for j in k:n
                candidate = 0.0
                for i in k:m; candidate = hypot(candidate, a[at(i,j)]); end
                if candidate > length; pivot, length = j, candidate; end
            end
            length > threshold || throw(ArgumentError("QR input is numerically rank deficient"))
            for i in 1:m; a[at(i,k)], a[at(i,pivot)] = a[at(i,pivot)], a[at(i,k)]; end
            p[k], p[pivot] = p[pivot], p[k]
            old = a[at(k,k)]; alpha = -copysign(length, old); divisor = old - alpha
            tau[k] = (alpha - old) / alpha
            for i in k+1:m; a[at(i,k)] /= divisor; end
            a[at(k,k)] = alpha
            for j in k+1:n
                dot = a[at(k,j)]
                for i in k+1:m; dot += a[at(i,k)] * a[at(i,j)]; end
                dot *= tau[k]; a[at(k,j)] = finite(a[at(k,j)] - dot)
                for i in k+1:m; a[at(i,j)] = finite(a[at(i,j)] - a[at(i,k)] * dot); end
            end
        end
    end
    Factorization(m, n, algorithm, scaling, norm, a, tau, p)
end

factor_lu(a::Matrix64) = Factorization(a, :lu)
factor_cholesky(a::Matrix64) = Factorization(a, :cholesky)
factor_qr(a::Matrix64) = Factorization(a, :qr)
solve(a::Matrix64, b::Matrix64) = solve(factor_lu(a), b)
least_squares(a::Matrix64, b::Matrix64) = solve(factor_qr(a), b)
solve(f::Factorization, b::Matrix64) = solve_internal(f, b, true)

function solve_internal(f::Factorization, rhs::Matrix64, rescale::Bool)
    m, n, p = f.rows, f.cols, rhs.cols
    rhs.rows == m || throw(DimensionMismatch("right-hand side row count must match factorization"))
    a, work = f.data, zeros(length(rhs.data))
    at(i,j) = (i-1)*n+j
    bt(i,j) = (i-1)*p+j
    for i in 1:m, j in 1:p
        value = finite(rhs.data[bt(f.algorithm == :lu ? f.permutation[i] : i, j)])
        work[bt(i,j)] = rescale ? solver_scaled(value, f.scale) : value
    end
    if f.algorithm == :qr
        for k in 1:n, j in 1:p
            dot = work[bt(k,j)]
            for i in k+1:m; dot = finite(dot + a[at(i,k)] * work[bt(i,j)]); end
            dot = finite(dot * f.tau[k]); work[bt(k,j)] = finite(work[bt(k,j)] - dot)
            for i in k+1:m; work[bt(i,j)] = finite(work[bt(i,j)] - a[at(i,k)] * dot); end
        end
    else
        for i in 1:n, j in 1:p
            value = work[bt(i,j)]
            for k in 1:i-1; value = finite(value - a[at(i,k)] * work[bt(k,j)]); end
            work[bt(i,j)] = f.algorithm == :cholesky ? finite(value / a[at(i,i)]) : value
        end
    end
    for i in n:-1:1, j in 1:p
        value = work[bt(i,j)]
        for k in i+1:n
            coefficient = f.algorithm == :cholesky ? a[at(k,i)] : a[at(i,k)]
            value = finite(value - coefficient * work[bt(k,j)])
        end
        work[bt(i,j)] = finite(value / a[at(i,i)])
    end
    result = Matrix64(n,p)
    for i in 1:n, j in 1:p
        result.data[bt(f.algorithm == :qr ? f.permutation[i] : i, j)] = work[bt(i,j)]
    end
    result
end

"""Computed infinity-norm rcond using an inverse (O(n³)), not a certified error bound."""
function reciprocal_condition(f::Factorization)
    n = f.cols
    f.rows == n || throw(DimensionMismatch("condition diagnostic requires a square matrix"))
    n == 0 && return 1.0
    inverse = try
        solve_internal(f, identitymatrix(n), false)
    catch error
        error isa ArgumentError || rethrow()
        return 0.0
    end
    norm = maximum(sum(abs(inverse[i,j]) for j in 1:n) for i in 1:n)
    min(1.0, (1.0 / f.norm) / norm)
end
