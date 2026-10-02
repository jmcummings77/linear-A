# Rust port

A dependency-free `f64` matrix library and comparison runner. Requires Rust and
Cargo 1.69 or newer; uses edition 2021. See the shared
[runner protocol](../../benchmarks/PROTOCOL.md) for input and timing details.

From this directory:

```sh
cargo test
cargo build --release
printf '1 2 3 4 5 6\n' | target/release/linear-a-rust check transpose 2 3
target/release/linear-a-rust bench multiply 32 100 7
```

To use the library from another Cargo package, add a path dependency:

```toml
[dependencies]
linear-a-rust = { path = "../linear-A/ports/rust" }
```

```rust
use linear_a::Matrix;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let a = Matrix::new(2, 2, &[1.0, 2.0, 3.0, 4.0])?;
    let identity = Matrix::identity(2)?;
    let product = a.multiply(&identity)?;
    assert_eq!(product.values(), a.values());
    assert_eq!(a.determinant()?, -2.0);
    Ok(())
}
```

`Matrix` exposes `rows`, `cols`, `get`, `set`, `values`, `row`, `column`, `add`,
`subtract`, `scale`, `transpose`, `multiply`, `cross`, `trace`, `determinant`, `determinant_with`, and
`triangular`. Use `clone()` for an independent copy. Construction copies its
input, and `values()` returns a read-only slice. Operations return owned results
without changing their inputs. Shapes, indices and nonfinite arithmetic are
checked through `Result<_, MatrixError>`.

Matrices use contiguous row-major storage. Multiplication uses the row/inner/column
loop order. `determinant()` uses `DeterminantAlgorithm::Auto`: triangular inputs
multiply their diagonal; other inputs use partial-pivot Gaussian elimination.
Use `determinant_with(DeterminantAlgorithm::Lu)` to force elimination or
`DeterminantAlgorithm::Cholesky` for exactly symmetric positive-definite matrices.
Cholesky rejects nonpositive computed pivots, including rounding failures.
Elimination takes O(n³) arithmetic operations and O(n²) scratch storage and rejects
nonfinite intermediates. Binary scaling preserves mixed large/small diagonal
products and avoids underflowing elimination factors. This is ordinary
floating-point elimination: small/large determinants can underflow/overflow, and
ill-conditioned matrices can lose accuracy. The 0x0 determinant is 1; its trace
is 0. Triangular classification uses exact zeros and returns `(upper, lower)`;
rectangular matrices are neither.

The runner also accepts `determinant_lu`, `determinant_spd_lu`, and
`determinant_cholesky`. The two SPD benchmarks use the same symmetric input,
`(A + Aᵀ)/2 + 4nI`; ordinary determinant input generation is unchanged.

`a.cross(&b)` computes the right-handed cross product of three-component vectors.
Either operand can be 3-by-1 or 1-by-3; the result preserves the left shape and
owns independent storage. Invalid shapes and nonfinite arithmetic return errors.

`Matrix::rotation_2d(radians)`, `rotation_x`, `rotation_y`, `rotation_z`, and
`rotation_axis_angle(&axis, radians)` create right-handed **active rotations of
column vectors**. Angles are finite radians; positive 2D angles rotate
counterclockwise. The arbitrary axis must be a finite, nonzero three-component
row or column vector. Stable normalization accepts huge, tiny, and subnormal
axis scales without changing the input. The result is a 2-by-2 or 3-by-3 matrix.
For a row-vector application, use the transpose of the rotation matrix.

The runner accepts `check cross AR AC BR BC`, `check rotation2d 0 0 ANGLE`, and
`check rotation3d AR AC ANGLE` (the last reads the axis). Their benchmark sizes are
fixed at 3, 2, and 3 respectively. Rotation benchmarks construct matrices at
0.5 radians; the 3D axis is `[1, 2, 3]`. Cross benchmarks use the first three
usual generated values for each vector, negating the right vector's third entry
so the ordinary three-entry checksum is informative.

`eigen_symmetric()` returns `SymmetricEigenDecomposition { eigenvalues,
eigenvectors }`: ascending real eigenvalues and corresponding unit eigenvectors
as columns of a matrix. Input must be exactly symmetric; nonsymmetric matrices
return an error, including matrices whose eigenvalues would be complex. The
cyclic Jacobi solver supports indefinite, singular, repeated-eigenvalue and empty
matrices. Repeated eigenspaces may have any orthonormal basis.
`eigen_symmetric_with(tolerance, max_sweeps)` controls the relative Frobenius-norm
tolerance and positive sweep limit (defaults `1e-12`, `50`). Nonconvergence and
nonfinite eigenvalues return errors. Accuracy is normwise: very small eigenvalues
among much larger entries can have large relative error. The runner's
`eigen_symmetric` benchmark uses a symmetric tridiagonal matrix with a known
analytic spectrum, consuming both eigenvalues and eigenvectors in its checksum.

The runner warms each workload before timing, consumes each result, and includes
result allocation and drop in the timed loop. Release builds retain line-table
debug information for external profilers; ordinary benchmarks enable no profiler.

`eigen_general()` accepts any finite real square matrix and returns
`GeneralEigenDecomposition` with `eigenvalues_real`, `eigenvalues_imag`,
`eigenvectors_real`, and `eigenvectors_imag`. Eigenvalues are sorted by real part,
then imaginary part. Matching columns form unit complex **right** eigenvectors:
`A v = λ v`. The input remains unchanged. Empty, singular, and defective matrices
are supported; defective matrices may yield dependent eigenvector columns, and
no orthogonality or complete-basis guarantee is made.

The solver uses power-of-two similarity balancing, a scaled Householder
Hessenberg reduction, and real double-shift QR with machine-epsilon deflation.
`eigen_general_with(max_iterations)` bounds iterations per unconverged root;
the default is 1000 and the accepted range is 1 through 100000. Nonconvergence
or nonfinite results return `MatrixError`. Results are approximate: tiny
components in mixed-scale matrices and ill-conditioned eigenpairs can lose
relative accuracy. Use `eigen_symmetric()` for the stronger orthogonality
contract on exactly symmetric input. The QR reduction follows the
[public-domain JAMA implementation](https://math.nist.gov/javanumerics/jama/)
of the EISPACK algorithms; it adds no dependencies.

The runner accepts `check eigen_general ROWS COLS` and
`bench eigen_general SIZE ITERATIONS SEED`. Its JSON separates eigenvalue and
eigenvector real/imaginary parts; both vector matrices are row-major. Benchmark
inputs have known real/complex spectra from upper block-triangular rotation
blocks. The checksum consumes weighted eigenvalues and every vector component.

Exact triangular eigenvalues retain the original diagonal entries, including
mixed huge/tiny values that cannot coexist in a uniformly scaled workspace.
