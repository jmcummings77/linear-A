# Singular value decomposition

Every port provides an economy SVD of a finite real m×n matrix:

A = U diag(s) Vᵀ, with k = min(m,n), U of shape m×k, k descending nonnegative
singular values, and Vᵀ of shape k×n. Columns of U and rows of Vᵀ are orthonormal.
The input is unchanged and each result owns independent storage. For an empty
input, k=0 and the factors retain the corresponding empty rectangular shapes.

| Port | Entry point | Result |
| --- | --- | --- |
| ARM64 assembly / C | `m_svd(a, tolerance, max_sweeps, &u, &s, &vt)` | Three owned matrices; `s` is k×1 |
| C++ | `a.svd(tolerance, max_sweeps)` | `Matrix::Svd { u, values, vt }` |
| C# | `a.Svd(tolerance, maxSweeps)` | `SingularValueDecomposition(U, Values, Vt)` for `Matrix<double>` |
| F# | `a.Svd(tolerance=..., maxSweeps=...)` | Record `{ U; Values; Vt }` |
| Go | `a.SVD()` / `a.SVDWith(tolerance, maxSweeps)` | `*SingularValueDecomposition`, error |
| Julia | `svd(a; tolerance=..., max_sweeps=...)` | Named tuple `(u, values, vt)` |
| Python | `a.svd(tolerance, max_sweeps)` | Named tuple `(u, values, vt)` |
| Rust | `a.svd()` / `a.svd_with(tolerance, max_sweeps)` | `Result<SingularValueDecomposition, MatrixError>` |
| TypeScript | `a.svd(tolerance, maxSweeps)` | `{ u, values: Float64Array, vt }` |
| WebAssembly | `a.svd(tolerance, maxSweeps)` | `{ u, values: Float64Array, vt }`; dispose both matrices |

Defaults are tolerance=1e-12 and maxSweeps=100. C requires both options explicitly.
Tolerance must be finite and strictly between 0 and 1; the integral sweep limit
must be in 1..10000. A convergence check after the final permitted rotation
sweep avoids reporting a failure merely because convergence occurred on that
last sweep. Nonconvergence is an error, never an apparently successful partial
result. The tolerance bounds normalized column correlations, not an absolute
error in singular values or a certified accuracy estimate.

## Algorithm and numerical limits

Cyclic one-sided Jacobi rotations orthogonalize the columns of a scaled copy of
A while accumulating V. Column norms give singular values and normalized
columns give U. Wide inputs use the transposed problem, then exchange the
factors. Norms use `hypot`; pair scaling avoids squaring large or tiny column
norms. Exact zero columns receive an orthonormal completion by twice applying
modified Gram–Schmidt to coordinate vectors. No numerical rank cutoff is
applied to ordinary small singular values. Repeated singular values and zero
singular directions can have different valid bases, and signs may differ.

Unlike the previous image demo, this does not construct AᵀA or AAᵀ, which would
square the condition number. It is a small handwritten implementation of the
one-sided Jacobi approach, not a port of LAPACK or a claim to LAPACK's robustness
or speed. [LAPACK DGESVJ](https://www.netlib.org/lapack/explore-html/d9/deb/group__gesvj_ga7aec05d2a1523bbeee77ece21b12187c.html)
provides a reference description of that algorithm family.

Scaling rejects an input if dividing by its largest entry would turn a nonzero
entry into zero. Singular values that overflow, or nonzero computed norms that
underflow when rescaled, also fail explicitly. If a required rotation angle is
below the smallest normal float64 value, the smaller working column is deflated
to zero. Its norm is negligible relative to the larger column, but this limits
relative accuracy at extreme dynamic ranges. There is no guarantee of relative
accuracy for arbitrarily small singular values. Reconstruction and orthogonality
should be considered together; an apparently tiny reconstruction residual alone
does not establish accurate singular vectors.

C outputs must be distinct, empty `{0}` matrices and must not alias the input;
on failure they remain empty. Free each successful result using `m_free`.
WebAssembly copies singular values into JavaScript memory and returns two owned
matrix handles; call `u.dispose()` and `vt.dispose()`. A failed call frees its
intermediate handles and leaves the input usable.

C++ shares the SVD arithmetic header with C, as it already does for solving and
general eigenpairs. ARM64 uses the C SVD alongside its assembly arithmetic
kernels; WebAssembly compiles that C implementation. The other seven ports have
native implementations. Agreement among shared kernels is not independent
numerical evidence.

## Validation

The shared protocol check `svd` returns the usual matrix JSON envelope containing
`[U; sᵀ; V]`, an (m+1+n)×k row-major matrix, so one process invocation preserves
matching factors. This packing belongs to the test/benchmark adapter only; the
public APIs return separate factors. `svd_one_sweep` exercises bounded-iteration
failure. It is not included in timed benchmark suites.

The conformance harness checks dimensions, finite results, descending values,
both orthogonality identities and scaled reconstruction residuals. Independent
analytic cases include rotated rational orthogonal factors with known spectra,
signed diagonals, repeated values, vectors, zeros and rank-one matrices. Uniform
scales span 1e-300 through 1e300, and isolated tiny singular values are checked
relatively. Seeded rectangular inputs exercise reconstruction without assuming
particular vectors. Negative tests cover scaling loss, overflow and exhausted
iterations. C ownership/error paths also run with address and undefined-behavior
sanitizers. The image playground uses this same public WebAssembly SVD API.
