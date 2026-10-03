# Matrix implementations and shared contract

Eleven reusable implementations provide the same float64 matrix operations:
construction, identity matrices, dimensions, indexed access, independent copying,
row/column extraction, addition, subtraction, scaling, transpose, multiplication,
trace, determinant, symmetric and general real eigendecomposition, 3D vector
cross products, rotation factories, linear-system solves, reusable factors,
QR least squares, economy SVD, pseudoinverses, minimum-norm least squares, and triangular classification. Each has its
own public API, tests, and runner for the [shared protocol](../benchmarks/PROTOCOL.md).

The comparison contract uses IEEE 754 double precision. Language-specific APIs
may offer additional numeric types or compatibility methods outside that contract.

| Implementation | Reusable API and storage | Toolchain |
| --- | --- | --- |
| [ARM64 assembly/C](assembly/README.md) | C matrix API with ARM64 arithmetic kernels | ARM64 host and C11 compiler/assembler |
| [C](c/README.md) | `matrix` and `m_*` functions; owned row-major allocation | C11 compiler |
| [C#](../linear-A/README.md) | `Matrix<double>`; rectangular managed arrays | .NET SDK from `global.json` |
| [C++](cpp/README.md) | `linear_a::Matrix`; owned row-major `std::vector<double>` | C++17 compiler |
| [F#](fsharp/README.md) | `LinearA.Matrix`; flat managed row-major array | .NET 10 SDK |
| [Go](go/README.md) | `matrix.Matrix`; row-major `[]float64` | Go 1.22+ |
| [Julia](julia/README.md) | `Matrix64`; row-major `Vector{Float64}` | Julia 1.10+ |
| [Python](python/README.md) | `Matrix`; row-major list of Python floats | Python 3.9+ |
| [Rust](rust/README.md) | `linear_a::Matrix`; owned row-major `Vec<f64>` | Rust/Cargo 1.69+ |
| [TypeScript](typescript/README.md) | `Matrix`; row-major `Float64Array` | Node.js 22+, npm; pinned TypeScript compiler |
| [WebAssembly](wasm/README.md) | C kernels in WebAssembly memory; JavaScript `Matrix` API | Emscripten SDK and Node.js 22+ for the runner |

No implementation calls BLAS, NumPy, or another numeric library. Some numerical
code is deliberately shared:

- ARM64 uses assembly arithmetic kernels with C allocation, validation, ownership,
  indexing, row/column extraction, triangular checks, and determinant pivot
  control. Both eigensolvers, SVD, linear-system solvers, and rotation factories execute in C; cross products
  use an ARM64 component kernel. Timings include this support work.
- C and C++ share the general Hessenberg/QR eigensolver, SVD, and system-factorization kernels.
  Their public APIs and storage management remain language-specific.
- WebAssembly compiles the C implementation with Emscripten. Its asynchronous
  JavaScript factory supports Node.js and browsers, with explicit `dispose()`
  for allocations. Timings include JavaScript dispatch and the WASM bridge.

## Shapes, ownership, and arithmetic

The runner protocol uses row-major values and finite float64 input/output.
Addition/subtraction require equal shapes; multiplication requires matching inner
dimensions. Trace and determinant require square matrices. Empty rectangular
shapes are supported; the 0×0 determinant is one and trace is zero. Triangular
checks use exact zeros.

Results in the shared contract do not mutate their inputs. Copies and row/column
extraction own independent storage. Use each API's explicit copy operation:
copying a C struct or assigning a Go struct can still share underlying storage.
Julia's public indices are one-based; the other APIs use zero-based indices.
Each implementation guide describes errors and ownership. The C# API also has
in-place methods; its comparison runner copies before invoking them.

Each implementation offers Auto, partial-pivot LU, and Cholesky determinant
selection. Auto may shortcut tiny or triangular matrices and never assumes
positive definiteness. Cholesky requires exact symmetry and positive computed
pivots. Both factorizations cost cubic arithmetic work; Cholesky uses roughly
half as many factorization operations on its narrower input domain. Separate
pivot mantissas and exponents prevent premature product overflow or underflow.
Rounding, cancellation, and extreme dynamic-range limitations remain. Square
solver factors expose a computed reciprocal condition diagnostic; see
[linear systems and numerical accuracy](SOLVING.md) for all eleven APIs, shapes,
rank decisions, ownership, and limitations.

C# additionally offers `Matrix<T>` with integer, decimal, and `BigInteger`
behavior, exact Bareiss elimination, and explicit cofactor expansion. Those
extensions are documented in the [C# guide](../linear-A/README.md). The
[determinant suite](../benchmarks/README.md#determinant-algorithm-comparisons)
compares algorithms on identical inputs, including a small cofactor comparison.

## Eigenvalues and eigenvectors

Both eigensolver families return independent results and preserve the input.
Empty matrices, singular matrices, and repeated eigenvalues are supported.
Eigenvectors occupy **columns**, satisfying `A Q = Q diag(values)` up to
floating-point error. Numeric accuracy is normwise: tiny eigenvalues beside
much larger entries may have large relative errors.

### Symmetric real matrices

The symmetric solver requires finite, square, exactly symmetric input and uses
cyclic Jacobi rotations. Real eigenvalues ascend and their corresponding columns
are orthonormal: `Qᵀ Q = I`. Indefinite matrices are valid. Repeated eigenvalues
can have any orthonormal basis of their eigenspace; individual signs are arbitrary.

The default relative Frobenius-norm tolerance is `1e-12`, with at most 50 sweeps.
APIs accept a finite tolerance strictly between zero and one and a positive sweep
limit; see each implementation for option syntax and limits. Invalid input,
nonconvergence, and eigenvalues outside the finite float64 range fail explicitly.

| Implementation | Symmetric eigenpair API and result |
| --- | --- |
| [ARM64/C](assembly/README.md) | `m_eigen_symmetric(&a, &values, &vectors)` |
| [C](c/README.md) | `m_eigen_symmetric(&a, &values, &vectors)` |
| [C#](../linear-A/README.md) | `GetSymmetricEigenDecomposition()` → `EigenValues`, `EigenVectors` |
| [C++](cpp/README.md) | `eigen_symmetric()` → `values`, `vectors` |
| [F#](fsharp/README.md) | `EigenSymmetric()` → `(values, vectors)` |
| [Go](go/README.md) | `EigenSymmetric()` → `Eigenvalues`, `Eigenvectors` |
| [Julia](julia/README.md) | `eigen_symmetric(a)` → named tuple `values`, `vectors` |
| [Python](python/README.md) | `eigen_symmetric()` → named tuple `values`, `vectors` |
| [Rust](rust/README.md) | `eigen_symmetric()` → `eigenvalues`, `eigenvectors` |
| [TypeScript](typescript/README.md) | `eigenSymmetric()` → `{ values, vectors }` |
| [WebAssembly](wasm/README.md) | `eigenSymmetric()` → `{ values, vectors }` |

Free both C/ARM64 outputs with `m_free`; dispose the WASM vector matrix.
C# converts generic numeric entries to double for eigendecomposition.

### General real matrices and complex eigenpairs

The general solver accepts any finite real square input, including nonsymmetric
and defective matrices. It returns eigenvalues sorted by real part, then imaginary
part, and unit complex right eigenvector columns. These columns need not be
orthogonal or independent; a defective matrix does not have a complete eigenbasis.
Complex matrix inputs are outside the shared contract.

Power-of-two similarity balancing precedes scaled Hessenberg reduction and real
double-shift QR, followed by complex back-substitution and normalization. The
iteration limit defaults to 1000 steps between root deflations and accepts
1 through 100000. Deflation uses machine epsilon. Invalid options, nonconvergence,
and nonfinite results fail explicitly. Exact triangular eigenvalues are retained
before workspace scaling, while ill-conditioned eigenpairs can still be sensitive
to rounding. The algorithms follow the public-domain JAMA/EISPACK routines cited
in each implementation.

| Implementation | General eigenpair API and result |
| --- | --- |
| [ARM64/C](assembly/README.md) | `m_eigen_general(&a, &real, &imag, &vectors_real, &vectors_imag)` |
| [C](c/README.md) | `m_eigen_general(&a, &real, &imag, &vectors_real, &vectors_imag)` |
| [C#](../linear-A/README.md) | `GetEigenDecomposition()` → complex `EigenValues`, `EigenVectors` |
| [C++](cpp/README.md) | `eigen_general()` → `values_real`, `values_imag`, `vectors_real`, `vectors_imag` |
| [F#](fsharp/README.md) | `EigenGeneral()` → complex `(values, vectors)` arrays |
| [Go](go/README.md) | `EigenGeneral()` → `EigenvaluesReal`, `EigenvaluesImag`, `EigenvectorsReal`, `EigenvectorsImag` |
| [Julia](julia/README.md) | `eigen_general(a)` → complex named tuple `values`, `vectors` |
| [Python](python/README.md) | `eigen_general()` → `values_real`, `values_imag`, `vectors_real`, `vectors_imag` |
| [Rust](rust/README.md) | `eigen_general()` → `eigenvalues_real`, `eigenvalues_imag`, `eigenvectors_real`, `eigenvectors_imag` |
| [TypeScript](typescript/README.md) | `eigenGeneral()` → `valuesReal`, `valuesImag`, `vectorsReal`, `vectorsImag` |
| [WebAssembly](wasm/README.md) | `eigenGeneral()` → `valuesReal`, `valuesImag`, `vectorsReal`, `vectorsImag` |

The protocol always splits real and imaginary parts into separate arrays/matrices,
regardless of the reusable API's complex representation. Free all four C/ARM64
outputs; dispose both WASM vector matrices. See each guide for option syntax.

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
to a row. Every implementation also offers rotations about the principal X, Y,
and Z axes; these APIs do not expose Euler angles or quaternions.

| Implementation | Cross product | 2D / arbitrary-axis rotation factories |
| --- | --- | --- |
| [ARM64/C](assembly/README.md) | `m_cross(&a, &b, &out)` | `m_rotation_2d(angle, &out)` / `m_rotation_axis_angle(&axis, angle, &out)` |
| [C](c/README.md) | `m_cross(&a, &b, &out)` | `m_rotation_2d(angle, &out)` / `m_rotation_axis_angle(&axis, angle, &out)` |
| [C#](../linear-A/README.md) | `a.CrossProduct(b)` on `Matrix<T>` | `MatrixRotation.Create2D(angle)` / `CreateAxisAngle(axis, angle)` |
| [C++](cpp/README.md) | `a.cross(b)` | `Matrix::rotation_2d(angle)` / `rotation_axis_angle(axis, angle)` |
| [F#](fsharp/README.md) | `a.Cross(b)` | `Matrix.Rotation2D(angle)` / `RotationAxisAngle(axis, angle)` |
| [Go](go/README.md) | `a.Cross(b)` | `matrix.Rotation2D(angle)` / `RotationAxisAngle(axis, angle)` |
| [Julia](julia/README.md) | `cross(a, b)` | `rotation2d(angle)` / `rotation_axis_angle(axis, angle)` |
| [Python](python/README.md) | `a.cross(b)` | `Matrix.rotation2d(angle)` / `rotation_axis_angle(axis, angle)` |
| [Rust](rust/README.md) | `a.cross(&b)` | `Matrix::rotation_2d(angle)` / `rotation_axis_angle(&axis, angle)` |
| [TypeScript](typescript/README.md) | `a.cross(b)` | `Matrix.rotation2D(angle)` / `rotationAxisAngle(axis, angle)` |
| [WebAssembly](wasm/README.md) | `a.cross(b)` | `Matrix.rotation2D(angle)` / `rotationAxisAngle(axis, angle)` |

Returned matrices follow each language's ownership rules; free C outputs and
dispose WASM outputs. C# rotation factories return `DoubleMatrix`, while its
generic cross product retains `T` and uses checked numeric operators. Its
compatibility classes also have an independent `GetCrossProduct` result and a
boolean in-place `CrossProduct` that returns false without mutation on failure.

## Build and check

Each implementation guide contains standalone build, test, and runner commands.
For an integrated check from the repository root:

```sh
python3 benchmarks/run.py --verify-only --require-all --output benchmarks/reports/verification
```

All eleven implementations are selected by default. On a host without every
toolchain, choose an explicit subset:

```sh
python3 benchmarks/run.py --verify-only --implementations c go python --require-all --output benchmarks/reports/verification
```

Verification covers 166 shared fixtures: 39 arithmetic, 23 symmetric eigenvalue,
20 general eigenvalue, 37 cross-product/rotation, and 47 solver cases. Eigenvalue tests
compare independent spectra and residuals, plus orthogonality for symmetric
results or unit complex column norms for general results. They do not require
a particular sign, complex phase, or basis within a repeated eigenspace.
Create separate performance reports with:

```sh
python3 benchmarks/run.py --operations eigen_symmetric eigen_general --output benchmarks/reports/eigen --require-all
python3 benchmarks/run.py --operations cross rotation2d rotation3d --output benchmarks/reports/vectors --require-all
```

The [benchmark guide](../benchmarks/README.md) explains executable overrides,
reports, timing, and profiling. Toolchain minimums are compatibility floors;
actual compiler/runtime versions are recorded in each report. Comparing an older
runtime with a newer one measures that version difference as well as the
implementation.

## Singular value decomposition

All eleven ports provide an economy SVD for finite rectangular real matrices,
plus pseudoinverses, minimum-norm solving and cutoff-based spectral diagnostics.
See [SVD APIs and numerical limits](SVD.md) for shapes, ownership, convergence,
and shared-kernel boundaries. The image playground uses the WebAssembly API.
