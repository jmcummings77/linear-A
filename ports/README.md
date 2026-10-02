# Matrix implementations

These ports are small, reusable implementations of the same float64 matrix API.
The comparison includes the existing C# library and ten ports. They provide
matrix construction, identity matrices, dimensions, indexed access, independent
copying, row/column extraction, addition, subtraction, scaling, transpose,
multiplication, trace, determinant, real symmetric eigendecomposition, 3D vector
cross products, rotation matrix construction, and triangular classification.

The shared contract covers IEEE 754 double precision. It does not reproduce the
C# generic library's integer, decimal, or `BigInteger` behavior.

| Implementation | Reusable API and storage | Toolchain |
| --- | --- | --- |
| [C# baseline](../README.md) | Existing `Matrix<double>`; rectangular managed arrays | .NET SDK from `global.json` |
| [F#](fsharp/README.md) | `LinearA.Matrix`; flat managed row-major array | .NET 10 SDK |
| [Rust](rust/README.md) | `linear_a::Matrix`; owned row-major `Vec<f64>` | Rust/Cargo 1.69+ |
| [Go](go/README.md) | `matrix.Matrix`; row-major `[]float64` | Go 1.22+ |
| [TypeScript](typescript/README.md) | `Matrix`; row-major `Float64Array` | Node.js 22+, npm; pinned TypeScript compiler |
| [Python](python/README.md) | `Matrix`; row-major list of Python floats | Python 3.9+ |
| [C++](cpp/README.md) | `linear_a::Matrix`; owned row-major `std::vector<double>` | C++17 compiler |
| [C](c/README.md) | `matrix` and `m_*` functions; owned row-major allocation | C11 compiler |
| [ARM64 assembly/C](assembly/README.md) | C matrix API with ARM64 arithmetic kernels | ARM64 host and C11 compiler/assembler |
| [Julia](julia/README.md) | `Matrix64`; row-major `Vector{Float64}` | Julia 1.10+ |
| [WebAssembly](wasm/README.md) | C matrix kernels in WebAssembly memory; JavaScript `Matrix` API | Emscripten SDK and Node.js 22+ for the runner |

The assembly implementation uses assembly for the arithmetic kernels. C provides
allocation, validation, ownership, indexing, row/column extraction, triangular
checks, and determinant pivot control. Symmetric eigendecomposition uses the
shared C Jacobi solver. Cross products use an ARM64 component kernel; rotation
factories and axis normalization use C. Its measurements include these C support
routines.
No implementation calls BLAS, NumPy, or another numeric library.

The WebAssembly target compiles the existing C implementation with Emscripten.
It supplies a reusable asynchronous JavaScript module factory for Node.js and
browsers, with explicit `dispose()` for matrix allocations. It compares the same
C kernels in a WebAssembly runtime, including the JavaScript bridge, rather than
introducing another independent arithmetic algorithm.

## Contracts and differences

The [runner protocol](../benchmarks/PROTOCOL.md) uses row-major values and finite
float64 input/output. Addition/subtraction require equal shapes; multiplication
requires matching inner dimensions. Trace and determinant require square
matrices. Empty rectangular shapes are supported; the 0×0 determinant is one and
trace is zero. Triangular checks use exact zeros.

Arithmetic results returned by the ports do not mutate their inputs. Copies and
row/column extraction own independent storage. Use each language's explicit
copy operation: copying a C struct or assigning a Go struct can still share its
underlying storage. Julia's public indices are one-based; the other APIs use
zero-based indices. Each port's README describes its errors and ownership rules.
The original C# library additionally retains its documented in-place methods;
the comparison runner copies before invoking those methods.

Each port offers Auto, partial-pivot LU, and Cholesky determinant selection.
Auto may shortcut tiny or triangular matrices; it never assumes a matrix is
positive definite. Cholesky explicitly requires exact symmetry and positive
computed pivots. Both factorizations cost cubic arithmetic work; Cholesky uses
roughly half as many factorization operations on its narrower input domain.
Separate pivot mantissas and exponents prevent premature product overflow or
underflow. Floating-point rounding, cancellation, and extreme dynamic-range
limitations remain; no conditioning estimate is provided.

The C# runner now calls `Matrix<double>.GetDeterminant` directly. Generic C# also
provides exact Bareiss elimination for integers and decimal, and an explicit
cofactor baseline. See [determinant comparisons](../benchmarks/README.md#determinant-algorithm-comparisons)
for identical-input comparisons between algorithms.

## Eigenvalues and eigenvectors

All implementations solve finite, exactly symmetric square matrices with cyclic
Jacobi rotations. They return ascending real eigenvalues and corresponding
orthonormal **columns** of an eigenvector matrix: `A Q = Q diag(values)` and
`Qᵀ Q = I`. Results own their storage and leave the input unchanged. Indefinite,
singular, repeated-eigenvalue, and empty matrices are supported. Repeated
eigenvalues can have any orthonormal basis of their eigenspace.

The default relative Frobenius-norm tolerance is `1e-12`, with at most 50 sweeps.
APIs offer a finite tolerance strictly between zero and one and a positive sweep
limit; see each port for its option syntax and limits. Nonconvergence, invalid
input, and eigenvalues outside the finite float64 range fail explicitly.
Accuracy is normwise: tiny eigenvalues beside much larger entries may have large
relative errors. General nonsymmetric matrices and complex spectra are unsupported.

| Implementation | Symmetric eigenpair API |
| --- | --- |
| [C#](../README.md#eigenvalues-and-eigenvectors) | `GetSymmetricEigenDecomposition()` → `EigenValues`, `EigenVectors` |
| [F#](fsharp/README.md) | `EigenSymmetric()` → `(values, vectors)` |
| [Rust](rust/README.md) | `eigen_symmetric()` → `eigenvalues`, `eigenvectors` |
| [Go](go/README.md) | `EigenSymmetric()` → `Eigenvalues`, `Eigenvectors` |
| [TypeScript](typescript/README.md) | `eigenSymmetric()` → `{ values, vectors }` |
| [Python](python/README.md) | `eigen_symmetric()` → named tuple `values`, `vectors` |
| [C++](cpp/README.md) | `eigen_symmetric()` → `values`, `vectors` |
| [C](c/README.md), [ARM64/C](assembly/README.md) | `m_eigen_symmetric(&a, &values, &vectors)`; free both outputs with `m_free` |
| [Julia](julia/README.md) | `eigen_symmetric(a)` → named tuple `values`, `vectors` |
| [WebAssembly](wasm/README.md) | `eigenSymmetric()` → `{ values, vectors }`; dispose the returned vector matrix |

C# converts generic entries to double for this operation. Its old integer-array
eigenvalue/eigenvector signatures are obsolete because general integer matrices
do not have integer spectra or normalized integer eigenvectors. Use the new
decomposition API on `IntegerMatrix` or `Matrix<T>` instead.

## Cross products and rotations

Cross products take two finite three-component vectors, each independently 3×1
or 1×3, and return an independent vector with the **left operand's shape**.
The convention is right-handed: X cross Y is Z. Invalid shapes and nonfinite
arithmetic fail explicitly. This is a vector operation; no general matrix cross
product is defined.

Rotation factories construct 2×2 or 3×3 matrices for right-handed **active
column-vector** rotations. Angles are finite radians. Positive 2D angles rotate
counterclockwise. Arbitrary axes must be nonzero finite three-component row or
column vectors; normalization is stable for extreme magnitudes and does not
change the axis. Apply a rotation with `R * v` to a column, or `v * transpose(R)`
to a row. Every port also offers rotations about the principal X, Y, and Z axes;
these APIs do not expose Euler angles or quaternions.

| Implementation | Cross product | 2D / arbitrary-axis rotation factories |
| --- | --- | --- |
| [C#](../README.md#cross-products-and-rotations) | `a.CrossProduct(b)` on `Matrix<T>` | `MatrixRotation.Create2D(angle)` / `CreateAxisAngle(axis, angle)` |
| [F#](fsharp/README.md) | `a.Cross(b)` | `Matrix.Rotation2D(angle)` / `RotationAxisAngle(axis, angle)` |
| [Rust](rust/README.md) | `a.cross(&b)` | `Matrix::rotation_2d(angle)` / `rotation_axis_angle(&axis, angle)` |
| [Go](go/README.md) | `a.Cross(b)` | `matrix.Rotation2D(angle)` / `RotationAxisAngle(axis, angle)` |
| [TypeScript](typescript/README.md), [WASM](wasm/README.md) | `a.cross(b)` | `Matrix.rotation2D(angle)` / `rotationAxisAngle(axis, angle)` |
| [Python](python/README.md) | `a.cross(b)` | `Matrix.rotation2d(angle)` / `rotation_axis_angle(axis, angle)` |
| [C++](cpp/README.md) | `a.cross(b)` | `Matrix::rotation_2d(angle)` / `rotation_axis_angle(axis, angle)` |
| [C](c/README.md), [ARM64/C](assembly/README.md) | `m_cross(&a, &b, &out)` | `m_rotation_2d(angle, &out)` / `m_rotation_axis_angle(&axis, angle, &out)` |
| [Julia](julia/README.md) | `cross(a, b)` | `rotation2d(angle)` / `rotation_axis_angle(axis, angle)` |

Returned matrices follow each language's existing ownership rules; free C outputs
and dispose WASM outputs. C# rotation factories return `DoubleMatrix`, while its
generic cross product retains `T` and uses checked numeric operators. The legacy
C# classes have an independent `GetCrossProduct` result and a boolean in-place
`CrossProduct` that returns false without mutation on invalid input or overflow.

## Build and check

Each linked README contains standalone build, test, and runner commands. For an
integrated check from the repository root:

```sh
python3 benchmarks/run.py --verify-only --require-all
```

All eleven implementations are selected by default. On a host without every
toolchain, choose an explicit subset:

```sh
python3 benchmarks/run.py --verify-only --implementations csharp python typescript
```

Verification covers 99 shared fixtures: 39 arithmetic, 23 eigenvalue, and 37
cross-product/rotation cases.
The eigenvalue tests compare independent spectra, eigenpair residuals, and
orthogonality without requiring a particular sign or repeated-eigenspace basis.
Create a separate eigendecomposition performance report with:

```sh
python3 benchmarks/run.py --operations eigen_symmetric --output benchmarks/reports/eigen --require-all
python3 benchmarks/run.py --operations cross rotation2d rotation3d --output benchmarks/reports/vectors --require-all
```

The [benchmark guide](../benchmarks/README.md) explains executable overrides,
reports, timing, and profiling. Toolchain minimums are compatibility floors;
actual compiler/runtime versions are recorded in each report. Comparing an older
local runtime with a newer one measures that version difference as well as the
implementation.
