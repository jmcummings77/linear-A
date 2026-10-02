# C float64 matrices

`matrix.h` provides an owned, contiguous row-major matrix API. Initialize each
output with `matrix result = {0}`, check the returned `matrix_status`, and call
`m_free` when finished. Outputs must be empty and distinct from their inputs.
Copying, arithmetic, transpose, and row/column extraction allocate independent
storage. Errors leave outputs unchanged. Use `m_copy` instead of copying the
struct: a plain struct copy would share the allocation and cannot be freed twice.
Rows and columns are zero-based. Empty dimensions are supported.

`m_get`/`m_set` provide indexed access; dimensions are `rows` and `cols`. Row and
column extraction return a 1-by-cols or rows-by-1 matrix. `m_identity` creates a
square identity matrix. All numeric operations use float64, reject nonfinite
inputs/results, and preserve normal floating-point rounding and underflow.
`m_determinant` selects Auto: guarded formulas for sizes 0 through 3 and an exact
triangular diagonal shortcut, with partial-pivot LU as the general fallback.
`m_determinant_with_algorithm(&a, M_DETERMINANT_LU, &value)` selects LU explicitly;
`M_DETERMINANT_CHOLESKY` selects Cholesky for an exactly symmetric positive
definite matrix and rejects asymmetry or nonpositive computed pivots. Auto never
assumes positive definiteness. LU and Cholesky take cubic time; the triangular
shortcut scans quadratic storage. Products use scaled factors to avoid
intermediate overflow or underflow, but factorization can still overflow, and
the final double can underflow or exceed its finite range. A tiny computed
determinant is not a numerical conditioning test. Triangular checks use exact
zeros. All determinant methods leave the input unchanged.

`m_eigen_symmetric(&a, &values, &vectors)` decomposes an exactly symmetric real
matrix as `A V = V diag(values)`. Initialize both outputs with `{0}` and free
both with `m_free`. Eigenvalues are sorted ascending in an n-by-1 matrix;
corresponding unit eigenvectors are columns of the n-by-n vector matrix.
`m_eigen_symmetric_with_options(&a, tolerance, max_sweeps, &values, &vectors)`
accepts `0 < tolerance < 1` and an integer sweep limit from 1 through 10000.
Errors leave both outputs unchanged; outputs must be distinct and empty.

The solver uses cyclic Jacobi rotations with power-of-two scaling. The default
relative tolerance is `1e-12`, measured by the off-diagonal Frobenius norm against
the original matrix norm, with at most 50 sweeps. It preserves the input. The
largest-magnitude component of each eigenvector is made nonnegative; repeated
eigenvalues can have any orthonormal basis of their eigenspace. Accuracy is
relative to the matrix norm, so very small eigenvalues beside much larger ones
need not have small relative error. Exact diagonal inputs retain their original
values even across extreme scales. Empty input returns empty eigenvalues and a
0-by-0 eigenvector matrix. Nonsymmetric matrices, invalid options, failure to
converge, and eigenvalues outside the finite double range are rejected.

`m_cross(&a, &b, &result)` creates the vector cross product. Rotation factories
are `m_rotation_2d(radians, &result)`, `m_rotation_x`, `m_rotation_y`,
`m_rotation_z`, and `m_rotation_axis_angle(&axis, radians, &result)`. Outputs
follow the usual empty/distinct allocation rule; errors leave them unchanged.

The cross product is defined for three-dimensional vectors, represented as
3-by-1 columns or 1-by-3 rows. The operands may use different orientations; the
result keeps the left operand's shape and owns independent storage. Components
use the right-hand rule: `x × y = z`. This is a vector operation, not a matrix
product. Arithmetic rejects nonfinite results, including product overflow.

Rotations are right-handed **active transformations of column vectors**, with
angles in radians. In 2D, `R = [[cos θ, -sin θ], [sin θ, cos θ]]`; a positive
quarter-turn maps x to y. In 3D, positive X rotates y toward z, positive Y
rotates z toward x, and positive Z rotates x toward y. An arbitrary axis may be
a row or column three-vector and must be finite and nonzero. Axis normalization
uses its largest component before computing its norm, so huge, tiny, and
subnormal axes work without magnitude overflow. Rodrigues' formula uses a
stable `1 - cos θ` expression for small angles. Rotation matrices can be applied
with ordinary matrix multiplication. Both operands and axes stay unchanged.

Build and test from the repository root with a C11 compiler:

```sh
clang -std=c11 -O3 -ffp-contract=off ports/c/matrix.c ports/c/runner.c -lm -o /tmp/linear-a-c
clang -std=c11 -O3 -ffp-contract=off ports/c/matrix.c ports/c/tests.c -lm -o /tmp/linear-a-c-tests
/tmp/linear-a-c-tests
/tmp/linear-a-c bench multiply 32 100 42
```

The runner implements [the shared protocol](../../benchmarks/PROTOCOL.md). Its
parser and timer are also used by the C++ and ARM64 runners; the C++ arithmetic
implementation is separate. No third-party libraries or fast-math are used.


`m_eigen_general(&a, &real, &imag, &vectors_real, &vectors_imag)` also accepts
nonsymmetric real square matrices. It returns real and imaginary eigenvalues in
separate n-by-1 matrices, sorted lexicographically by `(real, imag)`, and the
corresponding complex unit eigenvectors as columns of two n-by-n matrices:
`A (Vr + i Vi) = (Vr + i Vi) diag(real + i imag)`. All four outputs must be empty,
distinct, and freed with `m_free`. Errors leave every output unchanged.
`m_eigen_general_with_options` takes `max_iterations` immediately after the
input; the default is 1000 and the allowed range is 1 through 100000.

The general solver balances by power-of-two diagonal similarity, reduces to
Hessenberg form, and applies implicit real double-shift QR with machine-epsilon
deflation. The iteration bound resets after each one- or two-root deflation.
Nonconvergence returns `M_NO_CONVERGENCE`; nonfinite input or output returns
`M_NONFINITE`. Power-of-two work scaling and exact diagonal/triangular spectra
improve range handling, but rounding and conditioning still limit accuracy.
Defective matrices may have dependent eigenvectors, and general eigenvectors
are not guaranteed to be orthogonal. A largest component fixes the arbitrary
complex phase; repeated eigenvalues need not select a unique basis.

The shared numerical kernel in `general_eigen.h` is adapted from the
[public-domain JAMA orthes/hqr2 implementation](https://math.nist.gov/javanumerics/jama/),
based on Martin/Wilkinson and EISPACK. C++ uses the same general QR arithmetic
with its own storage and public API; WebAssembly and ARM64 use this C API.

## Linear systems

See [the shared solver guide](../SOLVING.md) for this implementation’s reusable
LU, Cholesky, and column-pivoted QR APIs, multiple right-hand sides, least squares,
condition diagnostics, ownership, and numerical limits.
