# Julia float64 matrix

Julia 1.10 or newer; standard library only. `Matrix64` uses a flat **row-major**
`Vector{Float64}`, matching the runner protocol. Its public indexing is idiomatic
Julia **one-based** indexing. It does not use native column-major matrix arithmetic,
LinearAlgebra, BLAS, or third-party numeric libraries. Its `data` field is an
implementation detail; use `rowmajor` to obtain an independent copy.

From the repository root:

```sh
julia --startup-file=no --project=ports/julia ports/julia/test/runtests.jl
printf '1 2 3 4' | julia --startup-file=no --project=ports/julia ports/julia/runner.jl check determinant 2 2
julia --startup-file=no --project=ports/julia ports/julia/runner.jl bench multiply 32 10 42
```

The package exposes `Matrix64`, `identitymatrix`, `copy`, `row`, `column`, `scale`,
`transpose`, `cross`, `trace`, `determinant`, `eigen_symmetric`, `eigen_general`,
`triangular`, and the `+`, `-`, `*` operators.
Indexing uses `matrix[row, column]`; dimensions are available through `size` and
`rows`/`cols`. Copies and arithmetic results own their storage. Invalid shapes,
indices, nonfinite input, and nonfinite results throw explicit exceptions.
`determinant(a)` defaults to `:auto`, which shortcuts triangular inputs. Use
`determinant(a, :lu)` for partial-pivot elimination or `determinant(a, :cholesky)`
for exact symmetric positive-definite input. Both use cubic arithmetic work and
scaled determinant products; Cholesky rejects nonpositive computed pivots.
Rounding and extreme dynamic-range limitations apply.

`eigen_symmetric(a; tolerance=1e-12, max_sweeps=50)` returns a named tuple with
`values::Vector{Float64}` and `vectors::Matrix64`, both independently owned.
Values are ascending and the matching unit eigenvectors are **columns** of
`vectors`. For `Matrix64(2, 2, [2., 1., 1., 2.])`, the values are `[1., 3.]`.
Inputs must be finite, exactly symmetric, and square. Empty input returns empty
values and a 0×0 matrix. Use `eigen_general` for nonsymmetric matrices and
complex eigenpairs.

The cyclic Jacobi iteration uses a Frobenius norm tolerance, finite and strictly
between zero and one, and a positive integer sweep limit. Nonconvergence and
nonfinite eigenvalues throw errors. Tiny eigenvalues alongside huge entries may
have poor relative accuracy. Eigenvector signs and bases within repeated
eigenspaces are not unique; compare residuals and orthogonality, not exact vectors.

`cross(a, b)` computes the right-handed 3D vector cross product. Each input must
be 3×1 or 1×3; different orientations are accepted and the result keeps `a`'s
shape. Inputs remain unchanged and results own their storage. Nonfinite input
or a nonfinite result raises an error.

`rotation2d(radians)`, `rotationx`, `rotationy`, `rotationz`, and
`rotation_axis_angle(axis, radians)` construct active, right-handed rotation
matrices. Angles are finite radians. Use `R*v` for column vectors and
`v*transpose(R)` for row vectors. Positive Z rotation takes X toward Y; positive
2D rotation is counterclockwise. An arbitrary axis must be a finite nonzero
three-component row or column vector; it is normalized without changing the
input, including when its entries are huge or subnormal. Compose `B*A` to apply
rotation A and then B. These functions do not call the standard `LinearAlgebra`
module; qualify `LinearAMatrices.cross` if both modules are imported.

The runner follows [the shared protocol](../../benchmarks/PROTOCOL.md). `time_ns()`
measures the warmed loop, including allocation, consumption, and any garbage
collection triggered within it. Startup and warmup compilation are excluded.
No fast-math annotations are used. `runner.jl` can be included without running
its `main` function, allowing external profiling wrappers.

For actual Julia stack samples from the warmed computation loop:

```sh
julia --startup-file=no --threads=1 --project=ports/julia ports/julia/profile.jl multiply 64 2000 42 /tmp/julia-profile.json
```

This uses the standard `Profile` sampler at a nominal 1 ms interval and emits
root-to-leaf frame stacks plus millisecond weights. Samples are marked
`chronological: false` because profiler buffers can contain multiple runtime
threads; weights represent sample counts times the requested interval, not exact
per-frame durations. The output is separate from benchmark timing. If a workload
finishes before a sample is captured, the command reports an error; increase the
iteration count.

`eigen_general(a; max_iterations=1000)` accepts any finite real square `Matrix64`.
It returns a named tuple `(values, vectors)`: `values` is a sorted
`Vector{ComplexF64}`, ordered by real part then imaginary part; `vectors` is a
complex matrix whose corresponding columns are unit **right** eigenvectors,
`A v = λ v`. The input remains unchanged. Empty, singular, and defective inputs
are supported; defective inputs may have dependent columns, and general
vectors are not promised to be orthogonal or a complete basis.

The native implementation uses power-of-two similarity balancing, scaled
Householder Hessenberg reduction, and real double-shift QR with machine-epsilon
deflation. `max_iterations` is an integer from 1 through 100000 and applies per
unconverged root. Nonconvergence and nonfinite results throw. Small components
in mixed-scale matrices and ill-conditioned eigenpairs can lose relative
accuracy. `eigen_symmetric` retains its separate orthonormal-vector contract.
The QR reduction follows the [public-domain JAMA implementation](https://math.nist.gov/javanumerics/jama/)
of the EISPACK algorithms; no numerical packages are required.

The runner accepts `check eigen_general ROWS COLS` and
`bench eigen_general SIZE ITERATIONS SEED`. Its portable JSON separates real
and imaginary parts and exports both vector matrices in row-major order.
Benchmark inputs use upper block-triangular rotation blocks with known spectra;
the checksum consumes weighted eigenvalues and every vector component.

Exact triangular eigenvalues retain the original diagonal entries, including
mixed huge/tiny values that cannot coexist in a uniformly scaled workspace.
