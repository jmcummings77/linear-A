# Sparse Cholesky

All eleven ports support symbolic analysis, numerical factorization, and repeated
solves of finite real symmetric positive-definite CSR systems. C++, ARM64 and
WebAssembly share the C scalar sparse kernel; ARM64 does not use assembly for
this factorization. The other ports implement the same algorithm natively.

## Contract

Analyze the stored pattern before supplying values to the numerical phase.
Analysis uses the undirected union of the pattern and its transpose, including
stored zeros. Every diagonal belongs to the lower factor. At pivot k, connecting
all surviving neighbors creates fill. The lower CSR pattern is sorted; each
row's diagonal is last. `fillSteps` aligns with column indices: -1 denotes an
original edge or diagonal, otherwise the zero-based pivot that created it.
The C representation uses `SIZE_MAX` instead of -1. Fill count excludes inserted
diagonals and includes entries that later cancel numerically to zero.

The symbolic plan owns the original pattern and lower pattern. Reuse requires
exactly the same dimensions, row offsets and column indices, including explicit
zeros. Values may change. Numerical factorization requires exact numeric
symmetry (a missing entry equals zero), positive computed pivots, and finite
arithmetic. Rectangular inputs, pattern mismatches, asymmetry, nonpositive
pivots, and nonfinite arithmetic fail explicitly. Empty square matrices work.
There are no shifts, pivoting, tolerance cutoffs, hidden reordering, or fallbacks.
Tiny positive finite pivots are accepted. Floating-point positive pivots do not
certify exact-real positive definiteness under arbitrary rounding error.

A factor owns its lower matrix independently of the source and plan. Its solve
method accepts a finite vector of matching length; forward substitution followed
by backward scatter solves LLᵀx=b. Returned vectors are owned. Reuse the factor
for any number of right-hand sides; a failed solve does not invalidate it.
Accessors return copies or immutable views. Julia and C expose fields that are
read-only by contract. C outputs start zero-initialized; failed construction
leaves them unchanged. C matrix solves also leave the output unchanged on
failure; the internal array kernel may partially write its destination.

## APIs

| Port | Analyze | Factorize | Solve / lower factor |
| --- | --- | --- | --- |
| C#, F# | `new SparseCholeskySymbolic(a)` | `s.Factorize(a)` | `f.Solve(b)`, `f.Lower` |
| Rust | `SparseCholeskySymbolic::new(&a)?` | `s.factorize(&a)?` | `f.solve(&b)?`, `f.lower()` |
| Go | `NewSparseCholeskySymbolic(a)` | `s.Factorize(a)` | `f.Solve(b)`, `f.Lower()` |
| Python | `SparseCholeskySymbolic(a)` | `s.factorize(a)` | `f.solve(b)`, `f.lower` |
| TypeScript, WASM | `new SparseCholeskySymbolic(a)` | `s.factorize(a)` | `f.solve(b)`, `f.lower` |
| Julia | `SparseCholeskySymbolic(a)` | `cholesky_factorize(s,a)` | `cholesky_solve(f,b)`, `f.lower` |
| C++ | `SparseCholeskySymbolic s(a)` | `s.factorize(a)` | `f.solve(b)`, `f.lower()` |
| C, ARM64 | `m_cholesky_analyze(&a,&s)` | `m_cholesky_factorize(&s,&a,&f)` | `m_cholesky_solve(&f,&b,&x)`, `f.lower` |

Free C plans with `m_cholesky_symbolic_free`, factors with `m_cholesky_free`, and
returned matrices normally. C error statuses distinguish argument/pattern,
asymmetry, nonpositive definiteness, allocation, and arithmetic range errors.
C++ uses RAII and noncopyable plans/factors. WASM plans, factors, and exported
lower matrices require `dispose()`; disposal is idempotent and operations on
disposed handles fail. Exporting a lower factor creates a new owned CSR matrix
(except the directly exposed C/Julia field).

Plans expose size, lower NNZ, row offsets, column indices and fill steps, using
each port's existing naming conventions. Fill count is `FillCount`,
`fill_count()` (Rust/C++), `fill_count` (Python), `fillCount` (TS/WASM), or the
number of nonnegative `fill_steps` (Julia); C stores `fill_count` directly.
Julia's CSR indices, offsets and fill steps are zero-based despite one-based
array indexing.

Apply `ReverseCuthillMcKee` and `PermuteSymmetric` before analysis to compare
orders. Permute the right-hand side and restore the solution with the same
new-to-old permutation. RCM reduces bandwidth heuristically; it need not reduce
fill, factorization time, or total time. AMD is not implemented.

## Measurement and visualization

Run `python3 benchmarks/cholesky_bench.py --require-all` for the cross-port
suite. The [interactive report](https://jmcummings77.github.io/linear-A/cholesky/)
shows actual symbolic fill births and independently checks LLᵀ reconstruction
and original-system solution residuals in WebAssembly.

Symbolic timing includes pattern export/checksum; numerical timing reuses a plan
and includes factor export/checksum; solve timing reuses a factor. Total timing
measures the complete path independently, including RCM and permutations when
selected. Logical factor storage assumes 8-byte values and indices; it excludes
symbolic plans, temporaries and runtime overhead. Native index sizes vary.
The algorithm is a simple row-oriented implementation with sparse intersections,
not a supernodal solver. The animation is an elimination trace, not a profiler.

Reference: [Cornell CS 6210 sparse factorization notes](https://www.cs.cornell.edu/courses/cs6210/2025fa/lec/2025-09-22.html).
