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

## Pseudoinverse, minimum-norm least squares, and spectral diagnostics

Every port also supplies the following operations. They use the default SVD
convergence settings above; the rank cutoff is a separate parameter.

| Port | Pseudoinverse | Minimum-norm solve | Diagnostics |
| --- | --- | --- | --- |
| ARM64 / C | `m_pseudoinverse(a, cutoff, &out)` | `m_solve_minimum_norm(a, b, cutoff, &out)` | `m_spectral_diagnostics(a, cutoff, &d)` |
| C++ | `a.pseudoinverse(cutoff)` | `a.solve_minimum_norm(b, cutoff)` | `a.spectral_diagnostics(cutoff)` |
| C# / F# | `a.Pseudoinverse(relativeCutoff)` | `a.SolveMinimumNorm(b, relativeCutoff)` | `a.SpectralDiagnostics(relativeCutoff)` |
| Go | `a.Pseudoinverse()` / `PseudoinverseWith(cutoff)` | `a.SolveMinimumNorm(b)` / `SolveMinimumNormWith(b, cutoff)` | `a.SpectralDiagnostics()` / `SpectralDiagnosticsWith(cutoff)` |
| Julia | `pseudoinverse(a; relative_cutoff=...)` | `solve_minimum_norm(a,b; relative_cutoff=...)` | `spectral_diagnostics(a; relative_cutoff=...)` |
| Python | `a.pseudoinverse(relative_cutoff)` | `a.solve_minimum_norm(b, relative_cutoff)` | `a.spectral_diagnostics(relative_cutoff)` |
| Rust | `a.pseudoinverse()` / `pseudoinverse_with(cutoff)` | `a.solve_minimum_norm(&b)` / `solve_minimum_norm_with(&b, cutoff)` | `a.spectral_diagnostics()` / `spectral_diagnostics_with(cutoff)` |
| TypeScript / WASM | `a.pseudoinverse(cutoff)` | `a.solveMinimumNorm(b, cutoff)` | `a.spectralDiagnostics(cutoff)` |

The default relative cutoff is `max(m,n) * 2^-52`. C uses `-1` to request this
default, and C++ accepts the same sentinel. Other ports use omitted/optional
arguments or the methods without `With`. Explicit cutoffs must be finite and
in [0,1]. A singular value is retained exactly when it is positive and
`s / s_max > cutoff`. Equality is discarded; zero retains every computed
positive value, and one retains none. This is a computed numerical rank, not an
exact algebraic rank. Comparisons use ratios to avoid underflow in a product
of the cutoff and a small largest value; zero cutoff tests positivity directly.

The n×m pseudoinverse is `V diag(s⁺) Uᵀ`, where retained values are reciprocated
and other entries of `s⁺` are zero. If a nonzero singular value is discarded,
this is the Moore–Penrose inverse of the *truncated matrix* Aτ. In general it
no longer satisfies `A A⁺ A = A` for the original matrix. Setting cutoff zero
does not guarantee that rounding preserves an exact null space.

For an m×q right-hand side B, the n×q solve returns the minimum Euclidean-norm
least-squares solution for Aτ, column by column. This handles underdetermined,
overdetermined and rank-deficient systems, including multiple right-hand sides.
It applies the factors directly instead of materializing the pseudoinverse.
Each RHS column is scaled before projection; scaling that would discard a
nonzero entry is rejected when the retained rank is nonzero. Binary exponent arithmetic
avoids an unnecessary overflow from forming `1/s` first. An unrepresentable
coefficient, term, or partial sum still fails even if later cancellation would
produce a finite answer. These algorithms do not use arbitrary precision.
The original SVD's convergence and extreme-dynamic-range limits still apply.

Diagnostics return three fields (using each language's naming convention):

- `rank`: number of retained singular values.
- `reciprocal_condition`: `s_min/s_max` of the complete economy spectrum,
  independent of the cutoff. It is the reciprocal 2-norm condition for a
  full-rank rectangular matrix, not the infinity-norm estimate from LU/QR.
- `retained_reciprocal_condition`: smallest retained value divided by the
  largest. This describes the retained subspace and can improve as rank falls.

Both condition fields are zero for an empty or all-zero spectrum; the retained
field is also zero if rank is zero. Ratios can underflow to zero. A zero original
ratio corresponds to an infinite/unrepresentable condition number, not a
well-conditioned system. A wide full-row-rank matrix can have a positive
reciprocal condition despite having a null space. These computed diagnostics
are not certified bounds.

Zero and empty inputs produce zero inverse/solution matrices with the natural
rectangular shapes. Inputs remain unchanged. C output matrices must be distinct
and empty; failures leave them unchanged. C diagnostic outputs are written only
on success. WASM inverse/solution matrices need `dispose()`; diagnostics are
ordinary copied JavaScript numbers and do not need disposal.

Shared analytic tests check the four Penrose identities for untruncated cases,
known minimum-norm answers, residual orthogonality, multiple RHS columns,
cutoff equality/adjacent floating-point boundaries, zeros, empty shapes and
extreme scales. Truncated cases are checked against the inverse of Aτ. A tiny
1×1 case verifies that direct solving succeeds when an explicit inverse cannot
be represented. The live accuracy playground uses these public WASM APIs.

## Ridge regularization

Every port supports minimizing **‖AX − B‖F² + λ‖X‖F²**, independently for each
right-hand side. A is m×n, B is m×r, and the owned result is n×r. This is ridge
regression / Tikhonov regularization with the identity penalty. λ is the penalty
coefficient, **not its square root**. Every variable is penalized, including an
intercept column if supplied. General penalty matrices and weighted fits are
not part of this API.

| Port | API |
| --- | --- |
| C / ARM64 assembly | `m_solve_ridge(&a, &b, lambda, &out)` |
| C++ | `a.solve_ridge(b, lambda)` |
| C# / F# | `a.SolveRidge(b, lambda)` |
| Go | `a.SolveRidge(b, lambda)` → matrix, error |
| Rust | `a.solve_ridge(&b, lambda)` → Result |
| Python | `a.solve_ridge(b, regularization)` |
| Julia | `solve_ridge(a, b, lambda)` |
| TypeScript / WebAssembly | `a.solveRidge(b, lambda)` |

For example, A = [2], B = [3], λ = 1 produces X = [1.2]. With λ = 0 it
produces [1.5]. C outputs must be empty and distinct from inputs. WASM results
must be disposed. Inputs are not modified. Empty shapes and zero matrices are
supported; a zero A returns zero X.

λ must be finite and nonnegative. At exactly zero, behavior matches the default
minimum-norm solver, including its numerical rank cutoff. For positive λ,
every computed positive singular value participates with filter
σ / (σ² + λ); zero singular directions contribute zero. Thus the limit as λ
approaches zero may differ from the λ = 0 result if the default cutoff removes
a small positive value. There is no implicit truncation for positive λ.

The implementation applies economy SVD factors directly, without forming a
Gram matrix or explicit inverse. Let d = max(σ, √λ). The denominator is
represented as d² [(σ/d)² + (√λ/d)²]; binary mantissas and exponents combine the
projection, σ, right-hand-side scale, and denominator before final rounding.
This avoids overflow from σ² and avoids rounding a tiny filter to zero before
multiplying a large right-hand side. The usual float64 underflow rounding still
applies to final coefficients.

The existing SVD convergence and input-scaling limitations remain. Right-hand
sides are scaled per column before projection; nonzero entries lost during
that scaling cause an error when positive singular directions participate.
An unrepresentable coefficient or accumulated result causes an error, even if
later cancellation could produce a finite result. No certified error bound or
automatic choice of λ is provided. Scaling both A and B by c requires scaling
λ by c² to preserve the objective's minimizer.

The live accuracy report compares ordinary least squares, truncated SVD, and
ridge on identical noisy inputs. Its tradeoff curve plots residual versus
solution size over positive λ. Sensitivity uses a separate perturbation Δ;
these metrics do not claim error against an unknown true solution. Saved
benchmark timings and profiles remain separate from these browser solves.
