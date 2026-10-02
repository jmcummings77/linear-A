# Go port

A standard-library-only `float64` matrix package and comparison runner. Requires
Go 1.22 or newer. See the shared [runner protocol](../../benchmarks/PROTOCOL.md)
for input and timing details.

From this directory:

```sh
go test ./...
go build -o runner ./cmd/runner
printf '1 2 3 4 5 6\n' | ./runner check transpose 2 3
./runner bench multiply 32 100 7
```

Use the package as follows (handle all returned errors in application code):

```go
import (
    "fmt"
    matrix "github.com/jmcummings77/linear-A/ports/go"
)

a, err := matrix.New(2, 2, []float64{1, 2, 3, 4})
if err != nil { panic(err) }
identity, err := matrix.Identity(2)
if err != nil { panic(err) }
product, err := a.Multiply(identity)
if err != nil { panic(err) }
fmt.Println(product.Values()) // [1 2 3 4]
```

`Matrix` exposes `Rows`, `Cols`, `At`, `Set`, `Copy`, `Values`, `Row`, `Column`,
`Add`, `Subtract`, `Scale`, `Transpose`, `Multiply`, `Cross`, `Trace`, `Determinant`, `DeterminantWith`, and
`Triangular`. Constructors and `Values` copy their input/output slices. Use
`Copy()` for an independent matrix; ordinary Go struct assignment copies a slice
header and shares storage. Operations return independent results and leave inputs
unchanged. Invalid dimensions, indices and nonfinite arithmetic return errors.
The zero value of `Matrix` is a valid empty 0x0 matrix.

Storage is contiguous row-major order. Multiplication uses row/inner/column loops.
`Determinant()` uses `DeterminantAuto`: triangular inputs multiply their diagonal;
other inputs use partial-pivot Gaussian elimination. Select `DeterminantLU` or
`DeterminantCholesky` with `DeterminantWith(algorithm)`. Cholesky requires exact
symmetry and finite positive computed pivots, returning an error if these checks
fail. Elimination costs O(n³) arithmetic operations and O(n²) scratch storage,
with exact-zero singularity checks. Binary scaling preserves mixed large/small
diagonal products and avoids underflowing elimination factors. Floating-point limitations
apply: tiny/large determinants can underflow/overflow, and ill-conditioned inputs
can lose accuracy. Nonfinite intermediates return errors. Empty determinants are
1; empty traces are 0. `Triangular` returns `(upper, lower)` using exact zeros;
rectangular matrices are neither.

The runner also accepts `determinant_lu`, `determinant_spd_lu`, and
`determinant_cholesky`. The two SPD benchmarks use the same symmetric input,
`(A + Aᵀ)/2 + 4nI`; ordinary determinant input generation is unchanged.

`a.Cross(b)` computes the right-handed cross product of three-component vectors.
Either operand can be 3-by-1 or 1-by-3; the result preserves the left shape and
owns independent storage. Invalid shapes and nonfinite arithmetic return errors.

Package functions `Rotation2D(radians)`, `RotationX`, `RotationY`, `RotationZ`,
and `RotationAxisAngle(axis, radians)` create right-handed **active rotations of
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

`EigenSymmetric()` returns a `SymmetricEigenDecomposition` with ascending
`Eigenvalues` and corresponding unit eigenvectors as columns of `Eigenvectors`.
Input must be exactly symmetric; nonsymmetric matrices return an error, including
matrices whose eigenvalues would be complex. The cyclic Jacobi solver supports
indefinite, singular, repeated-eigenvalue and empty matrices. Repeated eigenspaces
may have any orthonormal basis. `EigenSymmetricWith(tolerance, maxSweeps)` controls
the relative Frobenius-norm tolerance and positive sweep limit (defaults `1e-12`,
`50`). Nonconvergence and nonfinite eigenvalues return errors. Accuracy is
normwise: very small eigenvalues among much larger entries can have large relative
error. The runner's `eigen_symmetric` benchmark uses a symmetric tridiagonal matrix
with a known analytic spectrum and consumes eigenvalues and eigenvectors in its
checksum.

The benchmark loop includes allocation and any garbage collection that occurs
naturally during measurement; it does not force a collection after every result.
For a separate profiling run, opt into Go's standard CPU profiler:

```sh
LINEAR_A_CPU_PROFILE=cpu.pprof ./runner bench multiply 64 5000 7
go tool pprof -top ./runner cpu.pprof
```

Profiling covers warmup and timed operations, excludes input preparation, and
flushes before JSON serialization. The profiler is disabled unless the environment
variable is set. Profiling runs incur sampling overhead and should be kept separate
from the baseline timing report. Use a long enough run to collect CPU samples.

`a.EigenGeneral()` accepts any finite real square matrix and returns
`GeneralEigenDecomposition`: `EigenvaluesReal`, `EigenvaluesImag`,
`EigenvectorsReal`, and `EigenvectorsImag`. Eigenvalues are sorted by real part,
then imaginary part; corresponding matrix columns form unit complex **right**
eigenvectors satisfying `A v = λ v`. Results own their storage and input is
unchanged. Empty, singular, and defective inputs are supported. Eigenvectors
need not be orthogonal; defective inputs may produce dependent columns.

`a.EigenGeneralWith(maxIterations)` sets a per-root iteration limit from 1 to
100000 (default 1000). Power-of-two similarity balancing precedes scaled
Householder Hessenberg reduction and real double-shift QR with machine-epsilon
deflation. Nonconvergence or nonfinite output returns an error. Small components
in mixed-scale inputs and ill-conditioned eigenpairs may lose relative accuracy.
Use `EigenSymmetric` when its symmetric-input and orthogonality contract applies.
The QR reduction follows the [public-domain JAMA implementation](https://math.nist.gov/javanumerics/jama/)
of the EISPACK algorithms, without adding dependencies.

The runner supports `check eigen_general ROWS COLS` and
`bench eigen_general SIZE ITERATIONS SEED`. JSON has separate real/imaginary
arrays and row-major vector matrices. Benchmarks use upper block-triangular
rotation blocks with known spectra, consuming weighted eigenvalues and all
vector components in the checksum.

Exact triangular eigenvalues retain the original diagonal entries, including
mixed huge/tiny values that cannot coexist in a uniformly scaled workspace.
