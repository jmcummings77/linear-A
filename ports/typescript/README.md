# TypeScript float64 matrix

Node.js 22 or newer; TypeScript is a build-time dependency. `Matrix` owns a flat
row-major `Float64Array` and uses zero-based indices. Construction, copies, rows,
columns, and exported values have independent storage. Arithmetic returns new
matrices. Invalid dimensions, shapes, nonfinite input, and nonfinite results throw
`RangeError`.

`matrix.determinant()` defaults to `"auto"`: bounded, exactly representable small
matrices and triangular matrices have shortcuts; other inputs use partial-pivot
LU. Pass `"lu"` to select LU explicitly, or `"cholesky"` for an exactly symmetric
positive-definite matrix. Cholesky rejects nonsymmetric inputs and nonpositive
computed pivots. All methods preserve the input and return one for an empty
matrix. LU uses power-of-two row scaling where entries round-trip unchanged;
determinant products accumulate exponents separately. Rounding and conditioning
still matter: final underflow may return zero, and final overflow throws a
`RangeError`. A computed zero is not an exact singularity test.

`matrix.eigenSymmetric(tolerance = 1e-12, maxSweeps = 50)` returns
`{ values: Float64Array, vectors: Matrix }`: eigenvalues are ascending, and
column `j` of `vectors` is the corresponding unit eigenvector. It accepts exactly
symmetric square matrices, including empty matrices, and preserves the input.
Cyclic Jacobi rotations stop when the off-diagonal Frobenius norm is within
`tolerance` times the original norm; nonconvergence throws `RangeError`. Repeated
eigenvalues can have different valid orthonormal bases. Accuracy is relative to
the matrix norm, so a tiny eigenvalue in a mixed-scale matrix may have large
relative error. Use `eigenGeneral()` below for nonsymmetric inputs and complex eigenpairs.

`a.cross(b)` computes the right-handed 3D vector cross product. Each operand may
be a `3 × 1` column or `1 × 3` row; the result retains the left operand's shape.
It preserves both inputs and rejects other shapes or a nonfinite result.

`Matrix.rotation2D(radians)`, `rotationX`, `rotationY`, and `rotationZ` construct
active rotations acting on column vectors. Positive 3D angles follow the
right-hand rule; positive 2D angles are counterclockwise. All angles are radians.
`Matrix.rotationAxisAngle(axis, radians)` accepts a nonzero `3 × 1` or `1 × 3`
axis and normalizes it safely even at very large or subnormal magnitudes.
All rotation constructors return new matrices; a zero axis or nonfinite angle
throws `RangeError`.

From the repository root:

```sh
npm --prefix ports/typescript ci
npm --prefix ports/typescript test
printf '1 2 3 4' | node ports/typescript/dist/runner.js check determinant 2 2
node ports/typescript/dist/runner.js bench multiply 32 10 42
node ports/typescript/dist/runner.js bench determinant_cholesky 16 10 42
```

Import `Matrix` from `dist/matrix.js`. It exposes `rows`, `cols`, `get`, `set`,
`identity`, `copy`, `row`, `column`, `add`, `subtract`, `scale`, `transpose`,
`multiply`, `cross`, `trace`, `determinant`, `eigenSymmetric`, `eigenGeneral`, and
`triangular`. Declaration files are emitted alongside JavaScript. The CLI follows
[the shared protocol](../../benchmarks/PROTOCOL.md)
and uses `process.hrtime.bigint()`. Allocation and garbage-collection work triggered
during the timed loop are included; compilation and process startup are excluded.
No numeric libraries or BLAS are used.

The runner also accepts `determinant_lu`, `determinant_cholesky`, and
`determinant_spd_lu`. The two SPD benchmark operations use the same symmetrized,
strictly diagonally dominant input so Cholesky and LU can be compared fairly.
The ordinary determinant and explicit LU use the original general input.
`determinant_cofactor` is unsupported by this implementation.
`eigen_symmetric` checks return both eigenvalues and eigenvectors; its benchmark
uses a shared symmetric tridiagonal matrix with an analytic spectrum and
consumes the weighted spectrum and squared norm of every eigenvector entry.

The runner also accepts `check cross AR AC BR BC`, `check rotation2d 0 0 ANGLE`
(empty input), and `check rotation3d AR AC ANGLE` (axis input). Their benchmark
sizes are fixed at 3, 2, and 3 respectively. Rotation benchmarks use angle `0.5`;
3D rotation uses axis `[1, 2, 3]`. Cross inputs use the first three seeded entries,
with the third entry of the right operand negated to avoid a trivial checksum.

## General real eigenpairs

`matrix.eigenGeneral(maxIterations=1000)` accepts any finite real square matrix, including
nonsymmetric, singular and defective inputs. It returns separate real/imaginary
value arrays and real/imaginary matrix objects containing unit right eigenvector
columns. Values are sorted by real part then imaginary part. The result fields are
`valuesReal`, `valuesImag`, `vectorsReal`, `vectorsImag`.

The solver uses power-of-two similarity balancing and real Hessenberg/double-shift
QR adapted from public-domain [NIST/MathWorks JAMA](https://math.nist.gov/javanumerics/jama/).
The iteration limit bounds steps between deflations (1–100000). Input storage is
unchanged and outputs are independent. Invalid inputs, nonconvergence and
nonfinite outputs fail explicitly. General columns need not be orthogonal or
independent: defective matrices do not have a complete eigenbasis. Ill-conditioned
roots can be sensitive despite small residuals; complex matrix inputs are not
part of this API. The existing symmetric solver retains its stricter contract.

The runner's `eigen_general` operation supports both correctness checks and the
shared block-triangular benchmark described in the
[runner protocol](../../benchmarks/PROTOCOL.md).

## Linear systems

See [the shared solver guide](../SOLVING.md) for this implementation’s reusable
LU, Cholesky, and column-pivoted QR APIs, multiple right-hand sides, least squares,
condition diagnostics, ownership, and numerical limits.
