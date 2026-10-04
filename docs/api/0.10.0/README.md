# linear-A 0.10.0 API contract

This version describes release kits built from the Git revision recorded in their
`release.json`. It does not imply publication to any package registry. A tagged
release fixes that revision; development artifacts also record whether the
checkout was dirty. See [compatibility](compatibility.md) for every port.

## Shared numerical interface

All eleven ports provide owned, rectangular, row-major float64 matrices. Empty
shapes are supported. Construction copies supplied values; arithmetic returns
new matrices. Julia indices start at one; other ports start at zero. Dimensions
must agree: an m×k matrix times a k×n matrix produces m×n. Elementwise addition
and subtraction require equal shapes. A scalar scales each element.

| Capability | Operations and result |
| --- | --- |
| Arithmetic | Identity, copy, indexing, row/column extraction, add, subtract, scale, transpose, multiply, trace, triangular checks |
| Determinants | Square input; auto, partial-pivot LU, or Cholesky. Empty determinant is 1. Cholesky requires exact symmetry and positive computed pivots. |
| Eigenpairs | Symmetric real eigenvalues with corresponding unit eigenvector columns; general real matrices return complex eigenpairs as separate real/imaginary components. Iterative convergence can fail. |
| SVD | Economy U, descending nonnegative singular values, and Vᵀ for rectangular real inputs; one-sided Jacobi with explicit convergence/range failures. |
| Ridge | Rectangular regularized least squares, minimizing ‖AX − B‖F² + λ‖X‖F²; SVD filters, finite nonnegative λ, zero uses the default minimum-norm cutoff. |
| SVD inverse | Pseudoinverse, minimum-norm least squares for any rectangular rank, and spectral rank / reciprocal 2-norm conditions with an explicit relative cutoff. |
| Geometry | Three-dimensional vector cross product; 2D, axis, and normalized axis-angle rotations. Right-handed active rotations of column vectors, angles in radians. |
| Linear solves | `A X = B`, including multiple right-hand sides. Reusable LU and Cholesky factors own a snapshot. LU uses partial pivoting. |
| Least squares | Tall or square, numerically full-column-rank systems via column-pivoted QR. Reusable QR factors and reciprocal-condition diagnostics. |

No solver promises exact answers or certified error bounds. Rank tests are
numerical; nonfinite, incompatible, or out-of-range inputs can fail. LU and Cholesky reject singular systems; QR rejects numerical rank deficiency. SVD minimum-norm solving accepts deficient rank. Inspect residuals and conditioning when interpreting a result.
The source kit's `ports/SOLVING.md` specifies solver scaling, rank thresholds,
conditioning costs, and failure rules. `ports/README.md` describes the complete
shared interface, including differences in exception/status conventions.

## A small consumer in every language

Each kit contains `example/` with an executable example that solves
`[[4,1],[2,3]] x = [6,8]`, checks `x = [1,2]` to 1e-12, then checks multiplication
reconstructs the right-hand side. These are installation checks, not a substitute
for the analytical and generated conformance tests.

| Port | Construction | Solve | Storage / lifetime |
| --- | --- | --- | --- |
| ARM64 assembly | `m_from_array` | `m_solve` | C ownership; selected ARM64 arithmetic kernels |
| C | `m_from_array` | `m_solve` | Initialize outputs with `{0}`, check status, `m_free` |
| C# | `new Matrix<double>(double[,])` | `a.Solve(b)` | Managed; `linear_A` namespace |
| C++ | `linear_a::Matrix(rows, cols, values)` | `a.solve(b)` | RAII values; shared C headers must remain adjacent |
| F# | `Matrix.FromArray(rows, cols, values)` | `a.Solve(b)` | Managed; `LinearA` namespace |
| Go | `matrix.New(rows, cols, values)` | `a.Solve(b)` | Check returned errors; assignment is a shallow struct copy |
| Julia | `Matrix64(rows, cols, values)` | `solve(a,b)` | Managed; one-based indexing |
| Python | `Matrix(rows, cols, values)` | `a.solve(b)` | Managed; imports `matrix`, `general_eigen`, `solve` |
| Rust | `Matrix::new(rows, cols, &values)` | `a.solve(&b)` | Owned values; check `Result` |
| TypeScript | `new Matrix(rows, cols, values)` | `a.solve(b)` | Owned `Float64Array`; returned values are copies |
| WebAssembly | `await createMatrixAPI()`, then `new Matrix(...)` | `a.solve(b)` | Explicit `dispose()` on matrices and factorizations |

C# generic numeric types beyond double and netstandard2.1 legacy classes are
additional APIs, outside this shared float64 contract. Python's existing
unqualified module names can conflict with other installed packages; use a
virtual environment. No stable native ABI or backwards-compatibility promise is
made for this hobby library. Keep matching headers and sources together.

Version 0.4.0 added SVD-based pseudoinverses, minimum-norm solving, and rank/condition diagnostics to every port. See the source kit’s `ports/SVD.md` for
entry points, tolerance semantics, ownership, and precision limits.

Ridge API names, numerical limits, and examples are specified in [the shared SVD contract](../../../ports/SVD.md#ridge-regularization).

The shared API includes canonical CSR matrices, sparse matrix-vector multiplication, and conjugate gradient with optional Jacobi preconditioning to all eleven ports. See [the sparse contract](../../../ports/SPARSE.md) for construction, ownership, stopping rules, and numerical limits. Installed-consumer examples exercise these APIs too.

## GMRES and ILU(0)

All ports provide restarted, right-preconditioned GMRES with true residual stopping and optional captured iterates. Reusable ILU(0) factors own a snapshot and apply to multiple finite right-hand sides. No fill, pivoting or diagonal shifts are performed; zero computed pivots and nonfinite factorization/solve arithmetic fail explicitly. Stored zeros belong to the retained pattern. See `ports/SPARSE.md` in each kit for ownership, method names, stopping rules and limitations.

## Sparse ordering

All ports support deterministic Reverse Cuthill–McKee, symmetric CSR permutation, and forward/inverse vector permutation. See [the sparse contract](../../../ports/SPARSE.md#reverse-cuthillmckee-and-permutations) for indexing, graph semantics, failure behavior, and API names.

## Sparse Cholesky

Since version 0.8.0, the library provides reusable symbolic analysis, positive-definite numerical factorization and repeated triangular solves in every port. See [the sparse Cholesky contract](../../../ports/CHOLESKY.md) for APIs, fill histories, ownership and numerical limits.

## Approximate minimum degree

Version 0.10.0 adds approximate minimum degree ordering in every port; see [the AMD contract](../../../ports/AMD.md).

Version 0.10.0 adds reusable unshifted IC(0) and preconditioned conjugate gradient in every port. See [the shared IC(0) contract](../../../ports/IC0.md) for API names, failure semantics, ownership and ordering.
