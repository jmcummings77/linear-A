# Design and numerical scope

[Back to the project](../README.md) · [Shared API](../ports/README.md)

`linear-A` is a hobby and learning project: implement the same numerical
operations in eleven implementations, check their results, and study how algorithms,
storage, and runtimes affect performance. It has no intended production use
case. The reusable APIs let an implementation stand on its own; the common
runner makes their behavior and measurements comparable.

## Two engineering decisions

**One binary64 contract, checked against independent references.** The shared
comparison uses IEEE 754 double precision, portable shapes, and explicit
mathematical conventions. Verification uses worked examples, exact small
determinants and solves, analytic spectra, eigenpair residuals, and geometric
identities. No language implementation is the correctness oracle for the
others. Language-specific numeric types and compatibility APIs sit outside this
contract. The [runner protocol](../benchmarks/PROTOCOL.md) defines serialization,
errors, and timing boundaries; the [benchmark guide](../benchmarks/README.md)
explains verification and measurement.

**Implement the numerical work without external numerical libraries, and make
sharing explicit.** No implementation calls BLAS, NumPy, or another numerical
library. ARM64 assembly uses C for allocation, validation, ownership, and
selected algorithms. WebAssembly compiles the C implementation and exposes it
through a JavaScript API. C++ shares C kernels for general eigendecomposition,
system factorization, and SVD; newer sparse algorithms document their sharing
in their contracts. These boundaries matter when interpreting performance:
assembly timings include C support work, and WebAssembly timings include
JavaScript dispatch and the bridge.

## Find the mathematical contract

The guides below specify algorithms, shapes, options, failures, and per-language
APIs. Rectangular matrices are supported wherever the operation permits them.

| Area | Operations and detailed guide |
| --- | --- |
| Dense matrices | Construction, identity, copies, row/column extraction, arithmetic, transpose, multiplication, trace, determinant, and triangular classification: [shared API](../ports/README.md#shapes-ownership-and-arithmetic). |
| Dense systems and decompositions | Reusable pivoted LU and Cholesky factors, multiple right-hand sides, and column-pivoted Householder QR: [solving](../ports/SOLVING.md). Economy SVD, pseudoinverses, minimum-norm least squares, and ridge regularization: [SVD](../ports/SVD.md). |
| Eigenpairs | Cyclic Jacobi for symmetric matrices; balanced Hessenberg reduction and double-shift QR for general real matrices, including complex right eigenpairs: [eigenvalues and eigenvectors](../ports/README.md#eigenvalues-and-eigenvectors). |
| Sparse systems | Canonical CSR, matrix-vector products, CG, restarted GMRES, Jacobi, ILU(0), and RCM: [sparse guide](../ports/SPARSE.md). Specialized contracts cover [AMD ordering](../ports/AMD.md), [sparse Cholesky](../ports/CHOLESKY.md), and [IC(0)](../ports/IC0.md). |
| Geometry | Right-handed 3D cross products and active column-vector rotations with finite radian angles: [cross products and rotations](../ports/README.md#cross-products-and-rotations). |

## Ownership and indexing

Shared operations preserve inputs; copies and extracted rows/columns have
independent storage. Assigning a C or Go struct can still share storage, so use
the API's explicit copy operation. C# also provides in-place methods; its
comparison runner copies before calling them. Julia's dense API uses one-based
indices; its CSR and permutation indices are zero-based, as in the other ports.
Consult the
[implementation guides](../ports/README.md) for allocation, disposal, mutation,
and error rules.

## Numerical limits

Floating-point rounding, cancellation, overflow, and underflow remain possible.
Triangular and symmetry checks are exact. Cholesky requires positive-definite
input; choosing it does not make an arbitrary matrix suitable. Computed
reciprocal condition diagnostics are not certified error bounds. Eigensolver
tolerances and iteration limits are explicit; general eigenvectors need not
be orthogonal, and defective matrices can have dependent columns. Solver and
decomposition guides explain convergence, rank decisions, and residual checks.

For release-specific behavior and API evolution, use the
[versioned API contracts](api/); the root README links to the current version.
Installation and package contents belong in
the [release guide](../release/README.md).
