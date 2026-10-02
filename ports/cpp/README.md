# C++ float64 matrices

`matrix.hpp` provides a C++17 implementation; its general eigensolver also
includes the shared `../c/general_eigen.h` numerical kernel. `linear_a::Matrix` owns a
`std::vector<double>` and copies by value. Moves transfer ownership and leave the
source as a 0-by-0 matrix. Constructors accept dimensions and an optional vector
of row-major values; `Matrix::identity` creates identity matrices.

The public API provides `rows`, `cols`, `get`, `set`, `row`, `column`, `values`,
`scale`, `transpose`, `cross`, `trace`, `determinant`, `eigen_symmetric`,
`eigen_general`, `triangular`, and `+`, `-`, `*`.
Arithmetic returns independent matrices. Invalid dimensions throw
`std::invalid_argument`, invalid indices throw `std::out_of_range`, oversized
allocations throw `std::length_error` or `std::bad_alloc`, and nonfinite values
throw `std::overflow_error`. Empty dimensions are supported. `data()` permits
bulk interoperation; normal callers can use `set` for finite-value validation.

`a.determinant()` defaults to `linear_a::DeterminantAlgorithm::Auto`: guarded
formulas for sizes 0 through 3 and an exact triangular shortcut, with
partial-pivot LU as the general fallback. Pass `DeterminantAlgorithm::Lu` to
select LU, or `DeterminantAlgorithm::Cholesky` for an exactly symmetric positive
definite matrix. Cholesky rejects asymmetry and nonpositive computed pivots with
`std::domain_error`; invalid algorithm values throw `std::invalid_argument`.
Auto never assumes positive definiteness. LU and Cholesky take cubic time;
triangular detection scans quadratic storage. Scaled determinant products avoid
intermediate product overflow or underflow, but factorization and final double
rounding retain floating-point limitations. This is not an exact determinant or
conditioning calculation. Triangular checks use exact zeros.

`a.eigen_symmetric()` returns `SymmetricEigenResult` with a `std::vector<double>`
named `values`, sorted ascending, and an independently owned `Matrix` named
`vectors`, whose columns are the corresponding unit eigenvectors. The result
satisfies `A V = V diag(values)`. Optional positional arguments are the tolerance
and maximum sweeps: `a.eigen_symmetric(1e-12, 50)`. Tolerance must be finite and
between zero and one; sweeps must be from 1 through 10000. Asymmetry throws
`NotSymmetricError` (an `std::invalid_argument`); nonconvergence throws
`EigenConvergenceError` (an `std::runtime_error`).

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

`a.cross(b)` returns the vector cross product. Static rotation factories are
`Matrix::rotation_2d(radians)`, `rotation_x`, `rotation_y`, `rotation_z`, and
`rotation_axis_angle(axis, radians)`. Invalid vector shapes and zero axes throw
`std::invalid_argument`; nonfinite inputs or results throw `std::overflow_error`.

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

Build and test from the repository root:

```sh
clang++ -std=c++17 -O3 -ffp-contract=off ports/cpp/runner.cpp -o /tmp/linear-a-cpp
clang++ -std=c++17 -O3 -ffp-contract=off ports/cpp/tests.cpp -o /tmp/linear-a-cpp-tests
/tmp/linear-a-cpp-tests
/tmp/linear-a-cpp bench multiply 32 100 42
```

The runner uses the shared C protocol parser/timer through `bridge.hpp`, which
adapts calls to this independent C++ arithmetic implementation. It does not link
the C matrix library. See [the shared protocol](../../benchmarks/PROTOCOL.md).


`a.eigen_general(1000)` accepts any finite real square matrix and returns
`GeneralEigenResult`: vectors `values_real` and `values_imag`, plus independently
owned matrices `vectors_real` and `vectors_imag`. Complex eigenvalues are sorted
by `(real, imag)`; the corresponding complex eigenvectors occupy columns and
have unit norm. The optional argument bounds QR steps between root deflations
and must be from 1 through 100000. Invalid options or shape throw
`std::invalid_argument`, failure to converge throws `EigenConvergenceError`, and
nonfinite input or output throws `std::overflow_error`. Empty matrices work.

The Hessenberg/double-shift QR arithmetic is shared with the C implementation in
`general_eigen.h`, adapted from [public-domain JAMA](https://math.nist.gov/javanumerics/jama/).
The C++ API owns its `std::vector` workspace and results. Power-of-two similarity
balancing and work scaling improve range handling; diagonal/triangular spectra
retain their exact input entries. General eigenvectors need not be orthogonal
or independent for defective matrices, and ill-conditioned problems retain
floating-point accuracy limitations. The symmetric Jacobi API is unchanged.

## Linear systems

See [the shared solver guide](../SOLVING.md) for this implementation’s reusable
LU, Cholesky, and column-pivoted QR APIs, multiple right-hand sides, least squares,
condition diagnostics, ownership, and numerical limits.
