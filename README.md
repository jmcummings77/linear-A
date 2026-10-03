# linear-A

A hobby linear algebra library with eleven implementations in ARM64 assembly,
C, C#, C++, F#, Go, Julia, Python, Rust, TypeScript, and WebAssembly. Each provides
a reusable matrix API and participates in the same correctness checks and
performance comparisons. This is a learning project with no intended production
use case.

[View live report](https://jmcummings77.github.io/linear-A/latest/) ·
[Browse all reports](https://jmcummings77.github.io/linear-A/) ·
[Shared API and implementation differences](ports/README.md) ·
[Release kits and installation](release/README.md) ·
[Versioned API contract](docs/api/0.3.0/README.md)

## Choose an implementation

Start with the guide for the language you want to use. Each guide covers its
API, ownership and error rules, build steps, and tests; you only need that
implementation's toolchain to work with it.

| Implementation | Library location | Getting started |
| --- | --- | --- |
| ARM64 assembly | [ports/assembly](ports/assembly) | [Assembly guide](ports/assembly/README.md) · ARM64 host and C11 compiler/assembler |
| C | [ports/c](ports/c) | [C guide](ports/c/README.md) · C11 compiler |
| C# | [linear-A/linear-A](linear-A/linear-A) | [C# guide](linear-A/README.md) · .NET SDK selected by `global.json` |
| C++ | [ports/cpp](ports/cpp) | [C++ guide](ports/cpp/README.md) · C++17 compiler |
| F# | [ports/fsharp](ports/fsharp) | [F# guide](ports/fsharp/README.md) · .NET 10 SDK |
| Go | [ports/go](ports/go) | [Go guide](ports/go/README.md) · Go 1.22+ |
| Julia | [ports/julia](ports/julia) | [Julia guide](ports/julia/README.md) · Julia 1.10+ |
| Python | [ports/python](ports/python) | [Python guide](ports/python/README.md) · Python 3.9+ |
| Rust | [ports/rust](ports/rust) | [Rust guide](ports/rust/README.md) · Rust/Cargo 1.69+ |
| TypeScript | [ports/typescript](ports/typescript) | [TypeScript guide](ports/typescript/README.md) · Node.js 22+ and npm |
| WebAssembly | [ports/wasm](ports/wasm) | [WebAssembly guide](ports/wasm/README.md) · Emscripten; Node.js or a browser |

All eleven implementations cover the shared IEEE 754 float64 contract. Their
interfaces follow their language's conventions: Julia uses one-based indexing,
the others use zero-based indexing, and allocation, copying, mutation, and
disposal rules are documented per implementation. Additional numeric types and
compatibility APIs are described in the relevant language guide.

Some numerical kernels are shared. ARM64 assembly uses C for allocation,
validation, and selected algorithms; WebAssembly compiles the C implementation;
C++ shares the general eigensolver, linear-system factorization, and SVD kernels
with C. These boundaries are included
in the documentation and benchmark notes. No implementation calls BLAS, NumPy,
or another numerical library.

## Shared mathematics

The implementations provide construction and identity matrices, indexed access,
independent copies, row/column extraction, addition, subtraction, scaling,
transpose, matrix multiplication, trace, determinants, triangular classification,
eigenvalues and eigenvectors, singular value decomposition, linear-system
solving, least squares, 3D vector
cross products, and rotation matrices.
Rectangular matrices are supported wherever the operation permits them.

For example, multiplying this rectangular matrix by its transpose produces the
same result in every implementation:

```text
A = [1 2 3]       A Aᵀ = [14 32]       det(A Aᵀ) = 54
    [4 5 6]              [32 77]
```

The shared shape and mathematical conventions are:

- Addition and subtraction require equal shapes; multiplication requires matching
  inner dimensions. Trace and determinant require square matrices. Empty
  rectangular shapes are supported; the 0×0 trace is zero and determinant is one.
- Determinants offer automatic selection, partial-pivot LU, and Cholesky.
  Cholesky requires finite, exactly symmetric positive-definite input. LU and
  Cholesky use cubic arithmetic work; automatic selection can take safe shortcuts.
- Linear systems offer reusable pivoted LU and Cholesky factors, plus
  column-pivoted Householder QR least squares. Multiple right-hand-side columns
  share one factorization. See the [solver guide](ports/SOLVING.md).
- Economy SVD returns U, descending singular values, and Vᵀ for rectangular
  real matrices using one-sided Jacobi rotations. Pseudoinverses and minimum-norm least squares support wide and rank-deficient systems, with configurable singular-value cutoffs and rank/condition diagnostics. Explore the fit/stability tradeoff in the report’s accuracy playground. See the [SVD guide](ports/SVD.md).
- Symmetric eigendecomposition uses cyclic Jacobi rotations and returns ascending
  real eigenvalues with orthonormal eigenvector columns. General real square
  matrices use balanced Hessenberg reduction and double-shift QR and can return
  complex right eigenpairs. Both satisfy `A v = λ v`; general eigenvectors need
  not be orthogonal, and defective matrices can have dependent columns.
- Cross products accept three-component row or column vectors, use the right-hand
  rule, and retain the left operand's shape. Rotations use finite radian angles
  and active column-vector conventions: `R v` rotates a column, while `v Rᵀ`
  rotates a row. Positive 2D angles rotate counterclockwise.

Floating-point rounding, cancellation, overflow, and underflow still matter.
Square factors provide a computed reciprocal condition diagnostic, not a
certified error bound. Triangular and symmetry checks are exact;
eigensolver tolerances and iteration limits are explicit. See the
[shared API guide](ports/README.md) for detailed contracts and language-specific
differences, and the [runner protocol](benchmarks/PROTOCOL.md) for portable inputs,
outputs, and timing boundaries.

## Build, test, and compare

Each language guide contains standalone build and test commands. The common
harness additionally verifies every selected implementation against independent
references before measuring it. From the repository root:

```sh
# Verify all eleven implementations with their toolchains installed.
python3 benchmarks/run.py --verify-only --require-all --output benchmarks/reports/verification

# Or choose a subset by runner ID.
python3 benchmarks/run.py --verify-only --require-all --implementations go python rust --output benchmarks/reports/verification

# Save a new quick comparison and available profiles separately.
python3 benchmarks/run.py --suite quick --profiles --require-all --output benchmarks/reports/quick
```

The harness requires Python 3.9+ and the selected implementations' toolchains.
ARM64 assembly requires a native ARM64 host. Verification uses 166 shared
fixtures: 39 arithmetic, 23 symmetric eigenvalue, 20 general eigenvalue, and 37
cross-product/rotation cases, plus 47 solver cases. The references use worked examples, exact small
determinants and linear solves, analytic spectra, eigenpair residuals, and geometric identities;
no language implementation serves as the correctness oracle for the others.

The [benchmark guide](benchmarks/README.md) covers executable overrides, timing,
profiling, and separate determinant, eigenvalue, and vector comparisons. Results
compare particular implementations, algorithms, and toolchain versions on a
particular host; they are not a universal language ranking.

## Reports and browser exploration

The [published report](https://jmcummings77.github.io/linear-A/latest/) includes
operation timings, language filters, sample variation, and flamecharts. Its
embedded WebAssembly library runs correctness checks and bounded benchmarks on
your device. The geometry view animates matrix vector fields, shows real and
complex eigendirections, and steps through captured multiplication updates.
Matrices are editable or randomized, and geometric overlays can be hidden.
An accuracy playground lets you perturb a linear system and compare solution
sensitivity, residuals, backward error, and reciprocal condition.

Browser runs and illustrative animation frames stay separate from the saved
measurements and sampled profiles. Every report records its hardware, toolchain
versions, timestamp, and Git provenance. The
[saved HTML](benchmarks/reports/latest/index.html) also works offline, and
[results.json](benchmarks/reports/latest/results.json) contains its recorded data.

Committed report changes on `main` publish automatically to
[GitHub Pages](https://jmcummings77.github.io/linear-A/), without rerunning
benchmarks. The [machine-code dot-product experiment](experiments/machine-code-dot/README.md)
is a separate comparison of a small scalar kernel. The
[multiplication locality study](experiments/matmul-locality/README.md) isolates
loop order, cache blocking, and explicit ARM64 SIMD with independent correctness
checks, allocation accounting, and native stack profiles.

[Generated numerical checks and native fuzzing](tests/bugfinding/README.md) exercise
all eleven ports with reproducible inputs, independent references, and bounded
counterexample reduction.

Explore the [numerical playground](https://jmcummings77.github.io/linear-A/applications/):
draggable PCA, QR curve fitting, and a low-rank image explorer.
[Source and numerical limits](applications/README.md).

Compare saved runs in the [snapshot comparison](https://jmcummings77.github.io/linear-A/compare/),
with matched workloads, sample variation, and provenance caveats.
[File and Git-ref workflows](benchmarks/comparison/README.md).

Contributions to any implementation, the shared math checks, reports, and docs
are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for the development workflow.
