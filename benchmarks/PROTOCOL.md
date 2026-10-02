# Matrix runner protocol

The eleven [implementations](../ports/README.md) use this shared float64 contract
for correctness and timing requests. Matrix values use contiguous
row-major order in the protocol. Implementations own their storage. Matrix results
must not mutate either input. Empty shapes are allowed; a 0x0 determinant is 1 and
trace is 0. Invalid shapes/indices/input fail explicitly. Determinant defaults to
Auto (safe shortcuts and partial-pivot elimination); explicit LU and Cholesky are
also available. Cholesky requires exact symmetry and positive computed pivots.
Symmetric eigendecomposition returns ascending real eigenvalues and orthonormal
eigenvector columns. General real eigendecomposition accepts nonsymmetric input
and returns complex right eigenpairs as separate real and imaginary arrays.

## Correctness requests

`runner check OP ROWS COLS [BROWS BCOLS | SCALAR]`

Read exactly ROWS*COLS finite numbers from standard input, followed by BROWS*BCOLS
numbers for binary operations. Separate values with whitespace. Operations:

- `add`, `subtract`, `multiply`, `cross`: additional BROWS BCOLS arguments.
- `scale`: additional SCALAR argument.
- `rotation2d`: `ROWS COLS` must be `0 0`; additional ANGLE argument in radians,
  with no input values. Returns a 2-by-2 rotation matrix.
- `rotation3d`: `ROWS COLS` describes a 3-by-1 or 1-by-3 axis; additional ANGLE
  argument in radians. Read its three components; return a 3-by-3 rotation matrix.
- `transpose`, `trace`, `determinant`, `eigen_symmetric`, `eigen_general`, `triangular`: no additional arguments.
- `determinant_lu`, `determinant_spd_lu`, `determinant_cholesky`: no additional
  arguments. The two LU names are identical for correctness requests.
- The [C# implementation](../linear-A/README.md) also exposes `determinant_cofactor`
  for the small algorithm comparison; this is outside the shared operation set.

Write one JSON object to stdout:

- Matrix result: `{"rows":2,"cols":2,"values":[1,2,3,4]}`.
- Scalar result: `{"value":42}`.
- Triangular classification: `{"upper":true,"lower":false}`. Non-square -> both false.
- Symmetric eigendecomposition:
  `{"eigenvalues":[1,3],"eigenvectors":{"rows":2,"cols":2,"values":[1,0,0,1]}}`.
  This example represents `diag(1,3)`. Vector entries are row-major, with one
  eigenvector per **column**; column j corresponds to `eigenvalues[j]`.

Use enough digits to round-trip float64. Invalid requests exit nonzero with a useful
message on stderr. No diagnostics on stdout. Trace/determinant require square
matrices. Addition/subtraction require identical shapes. Matrix multiplication
requires left columns == right rows. NaN/infinity results fail rather than emit
invalid JSON. Triangular classification uses exact zeros.

`cross` takes two 3-by-1 or 1-by-3 vectors, independently of orientation, and
returns the left shape with components
`[a1*b2-a2*b1, a2*b0-a0*b2, a0*b1-a1*b0]` (zero-based subscripts).
Inputs and results must be finite; invalid shapes or nonfinite arithmetic fail.
Zero and parallel vectors are valid. There is no general matrix cross product.

Rotations are right-handed **active rotations of column vectors** with finite
radian angles. `rotation2d` returns `[[cos θ, -sin θ], [sin θ, cos θ]]`.
`rotation3d` uses a finite, nonzero axis normalized stably by its maximum absolute
component before taking its norm. It returns Rodrigues' rotation matrix.
Positive rotation about Z maps X toward Y. A column vector is transformed by
`R v`; a row vector uses `v Rᵀ`. Factories construct matrices without applying
them to other inputs. Reusable APIs also expose principal X/Y/Z-axis factories;
Euler-angle and quaternion APIs are outside this contract.

`eigen_symmetric` requires finite, square, exactly symmetric input. Use cyclic
Jacobi rotations with relative Frobenius tolerance `1e-12` and a maximum of 50
sweeps; nonconvergence must fail. Stop when the working off-diagonal Frobenius norm
is at most the tolerance times the original matrix Frobenius norm, using scaling
and stable norm evaluation to avoid unnecessary overflow or underflow. Diagonal
input preserves its eigenvalues, including mixed extreme scales. Empty input
returns `[]` and a 0-by-0 eigenvector matrix. Indefinite, singular, and repeated
spectra are valid; repeated eigenspaces may have any orthonormal basis. The
accuracy target is normwise, not a relative-error guarantee for each eigenvalue.
Reusable APIs may expose tolerance and sweep-limit options; invalid options and
eigenvalues outside the finite float64 range fail explicitly.

Reusable library APIs must expose dimensions/indexed access, copying, add/subtract,
scale, transpose, multiply, cross, rotation factories, row/column extraction,
trace, determinant, symmetric/general real eigendecomposition, and triangular
classification. Use idiomatic error handling and no third-party numeric libraries.
Document shared numerical code and include its work in timings: ARM64 uses C
allocation/validation, both C eigensolvers, and C rotation factories alongside
assembly arithmetic and cross-product kernels. WebAssembly compiles the C kernels;
C++ shares the general eigensolver kernel with C and owns its storage separately.

## Benchmark requests

`runner bench OP SIZE ITERATIONS SEED`

Supported operations: add, subtract, scale, transpose, multiply, cross, rotation2d,
rotation3d, trace, determinant, determinant_lu, determinant_spd_lu,
determinant_cholesky, eigen_symmetric, and eigen_general. The C# runner also exposes
determinant_cofactor for the small algorithm comparison.
SIZE and ITERATIONS must be positive. SEED is an integer in 0..2147483646.
Except for the vector/rotation workloads below, create square A and B before timing:

```
A[index] = (((index * 17 + SEED * 13) % 101) - 50) / 16.0
B[index] = (((index * 17 + (SEED + 1) * 13) % 101) - 50) / 16.0
```

For all determinant variants, add SIZE*4 to each diagonal entry of A to avoid
trivial singular workloads. For `determinant_spd_lu` and `determinant_cholesky`,
first replace A with `(A + transpose(A))/2`, before that diagonal shift. Both
variants must receive identical symmetric positive-definite inputs. `determinant_lu`
uses the same general input as `determinant`; C# `determinant_cofactor` does too.

For `eigen_symmetric`, replace A with the symmetric tridiagonal Toeplitz matrix:

```
A[i,j] = 2 + (SEED % 17)/16.0  if i == j
         -1.0                 if abs(i-j) == 1
          0.0                 otherwise
```

Do not add the determinant diagonal shift. For size n its ascending eigenvalues
are `2 + (SEED % 17)/16.0 - 2*cos(k*pi/(n+1))`, for k from 1 through n.
The harness verifies the full eigenbasis for each exact timed input before
calibration: check this independent spectrum, `A Q = Q diag(values)`, and
`Qᵀ Q = I`. The quick suite uses eigenvalue sizes 8 and 16; the full suite uses
16, 32, and 48. Shared correctness verification includes 119 fixtures: 39
arithmetic, 23 symmetric eigenvalue, 20 general eigenvalue, and 37
cross-product/rotation cases, independently of selected timed operations.

Vector and rotation workloads use fixed sizes in both suites:

| Operation | Required SIZE | Prepared input and operation |
| --- | --- | --- |
| `cross` | 3 | A and B are 3-by-1 vectors from the first three usual generated entries above; **negate B[2]** before computing A cross B. |
| `rotation2d` | 2 | A is empty 0-by-0; construct a 2-by-2 rotation at 0.5 radians. |
| `rotation3d` | 3 | A is the 3-by-1 axis `[1,2,3]`; construct a 3-by-3 rotation at 0.5 radians. |

Reject any other SIZE for these operations. Before calibration or timing, the
harness probes the complete result for the exact timed input against an
independent reference. Cross references use rational arithmetic and rotation
references act on basis vectors with quaternion multiplication. These references
are verification code, not alternate public rotation APIs. Use the ordinary
matrix checksum for all three operations; cross's negated B[2] avoids a
trivially zero sum for typical consecutive generated vectors.

Scale uses 1.25. Each iteration operates on the same original
inputs, creating independent matrix results. Time includes output allocation and
per-iteration disposal where deterministic; garbage collection occurs when the
runtime schedules it. Input creation, parsing, compilation, process startup, and
JSON serialization are excluded. Use the platform's monotonic high-resolution timer.

Warm up with `max(5, min(ITERATIONS, 100))` operations before starting the timer.
During the timed loop accumulate a checksum from each result: for a matrix of L
values use `values[0] + values[L/2] + values[L-1]` (integer division; empty -> 0);
for a scalar use its value. For symmetric eigenpairs use
`sum((i+1)*eigenvalues[i]) + sum(Q[i,j]*Q[i,j])`, with zero-based i and every
vector entry included. The independent reference is the weighted analytic
spectrum plus n, since an orthonormal n-by-n basis has squared Frobenius norm n.
Runners must consume the computed vectors rather than substitute n. The
[general eigenpair checksum](#general-real-eigendecomposition) consumes both real
and imaginary parts. Include this consumption in elapsed time. Write:

`{"elapsed_ns":123456,"iterations":10,"checksum":42.0}`

The harness checks checksums independently; do not hardcode checksums or optimize
away operations. Avoid fast-math. Report nonfinite outputs and overflow as errors.
Tests and benchmark runners must compile with the standard language toolchains.

## General real eigendecomposition

`check eigen_general ROWS COLS` accepts any finite real square input. It returns
exactly `eigenvalues_real`, `eigenvalues_imag` (length n arrays),
`eigenvectors_real`, `eigenvectors_imag` (n-by-n row-major matrix objects).
Column j represents the unit complex right eigenvector for eigenvalue j.
Sort lexicographically by real then imaginary part. Eigenvectors may be
nonorthogonal or dependent; defective matrices do not imply a full basis.
Empty inputs yield empty arrays and 0-by-0 matrices. Input storage is preserved.

Use power-of-two similarity balancing, Hessenberg reduction and double-shift QR
with machine-epsilon deflation. Bound QR steps between deflations by
`maxIterations=1000` by default, with allowed range 1–100000. Nonconvergence and
nonfinite results are errors. The public-domain NIST/MathWorks JAMA orthes/hqr2
routines are the numerical reference. Exact triangular eigenvalues are retained
before global workspace scaling. Accuracy is normwise; ill-conditioned roots
can be sensitive even when the residual is small.

For `bench eigen_general SIZE ITERATIONS SEED`, fill each entry A[i,j] with
zero-based k=floor(i/2), a=1+(seed%17)/16+k/8, b=0.5+k/16:

- i=j: a.
- i even and j=i+1: -b.
- i odd and j=i-1: b.
- Other i<j: ((3*i+5*j+seed)%11-5)/32.
- Otherwise: 0.

The block-triangular spectrum is a±ib per complete 2-by-2 block, with a real
last root for odd n. No determinant diagonal shift is added. Before timing,
verify the complete spectrum with multiplicities and every complex residual
`A v − λ v`, plus unit column norms. Consume all results with
`sum((i+1)*(real[i]+abs(imag[i]))) + sum(Vreal²+Vimag²)` per iteration.
The matrix sizes match the symmetric eigen workload.
