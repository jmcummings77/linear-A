# CSR matrices and iterative solvers

All eleven ports expose the same float64 sparse contract. The C++, ARM64 and
WebAssembly ports share the C sparse solver; the other implementations use native
language loops. This is a small educational solver, not a replacement for a
specialized sparse numerical package.

## Storage and ownership

A CSR matrix owns three arrays: row offsets (rows + 1), column indices (nnz),
and finite values (nnz). Offsets start at zero, end at nnz, and never decrease.
Columns within each row must be strictly increasing, unique, and in range.
Explicit zero values are allowed. Empty rows and zero-sized dimensions work.
All CSR offsets and column indices are **zero-based, including Julia**. This
is distinct from Julia's ordinary array indexing.

Construction copies inputs; it never silently sorts or merges duplicates.
Dense conversion omits zero entries. Accessors return copies or immutable
borrows, except for Julia's publicly accessible array fields and C's struct
fields: do not mutate their invariants. Julia and C public operations revalidate
the structure. Matvec accepts a finite matching vector and returns a new vector;
out-of-range arithmetic fails. Constructors and output allocations may fail.

| Port | Construction / dense conversion | Product / solver |
| --- | --- | --- |
| C / ARM64 | `m_csr_create`, `m_csr_from_dense` | `m_csr_matvec`, `m_csr_cg` |
| C++ | `linear_a::CSRMatrix`, `CSRMatrix::from_dense` | `matvec`, `conjugate_gradient` |
| C# (.NET 10) | `linear_A.CSRMatrix`, `CSRMatrix.FromDense` | `Matvec`, `ConjugateGradient` |
| F# | `LinearA.CSRMatrix`, `CSRMatrix.FromDense` | `Matvec`, `ConjugateGradient` |
| Go | `matrix.NewCSR`, `matrix.CSRFromDense` | `Matvec`, `ConjugateGradient` |
| Rust | `CSRMatrix::new`, `CSRMatrix::from_dense` | `matvec`, `conjugate_gradient` |
| TypeScript / WASM | `new CSRMatrix`, `CSRMatrix.fromDense` | `matvec`, `conjugateGradient` |
| Python | `CSRMatrix`, `CSRMatrix.from_dense` (from `matrix`) | `matvec`, `conjugate_gradient` |
| Julia | `CSRMatrix`, `csr_from_dense` | `matvec(a,x)`, `conjugate_gradient(a,b)` |

C outputs must be zero-initialized and empty. Caller-provided pointer lengths
must match the declared dimensions/counts. `m_csr_matvec` takes a dense n×1
matrix as its vector and writes a new dense column. Free storage with
`m_csr_free`, `m_free`, and `m_cg_free`. Never shallow-copy owning C structs.
WASM `CSRMatrix` instances require `dispose()`; returned vectors and CG result
arrays are independent JavaScript copies that survive disposal.

## Solver contract

CG accepts a square **symmetric positive-definite** matrix and a single finite
right-hand-side vector. Exact symmetry is checked sparsely; a missing entry
compares as zero. Positive definiteness is a caller precondition, not certified
by the implementation. A nonpositive computed search curvature produces
`breakdown`; passing this check does not prove positive definiteness. Optional
Jacobi preconditioning requires positive stored diagonal entries.

The starting guess is zero. Defaults are relative tolerance 1e-10, absolute
tolerance 0, maximum 1000 iterations, no preconditioner, no frame capture.
Relative tolerance must be finite in [0,1), absolute tolerance finite and
nonnegative, and the limit an integer in [0,100000]. Zero iterations returns
the initial guess unless it already satisfies the tolerance. In Go start from
`DefaultCGOptions()`; a zero-valued options struct explicitly requests a zero
limit and zero tolerances. Rust uses `CGOptions::default()`.

Stop when `norm(b - A*x, 2) <= max(atol, rtol * norm(b, 2))`.
The true residual is recomputed at every accepted update and replaces the
recursive residual. This costs **two SpMVs per iteration** and makes recorded
histories directly interpretable. It is a residual-replacement variant of
preconditioned CG; floating-point convergence is not guaranteed within n steps.
Residual 2-norms need not decrease monotonically. Tiny tolerances may be
unattainable, and poorly conditioned systems can stagnate.

The result contains `x`, `converged`, accepted `iterations`, `reason`,
`residuals` (initial norm plus one per accepted update), and optional `iterates`
(initial zero vector plus each accepted x). C# / F# / Go use capitalized field
names. C uses `matrix_cg_result`, with row-major flattened frames and reason
codes 0=converged, 1=iteration_limit, 2=breakdown, 3=nonfinite; convergence is
`reason == 0`, and history length is `iterations + 1`.

Invalid shapes, arrays, symmetry, options and Jacobi diagonals are input errors.
Runtime `iteration_limit`, `breakdown` and `nonfinite` are explicit result
states, preserving the last accepted estimate. Dot products use ordinary
float64 accumulation: extreme scaling can overflow or underflow even when an
exact solution exists. Norms use repeated hypot; an unrepresentable initial
norm returns `nonfinite` with an infinite initial history entry. Do not serialize
that as a successful finite result. There are no certified error bounds.

Without capture, storage is O(nnz + n + iterations). Capturing frames costs
O(n × iterations). The C kernel reserves history and optional frames up to the
requested limit before solving; the other implementations may grow them as
needed. WASM then copies accepted frames into JavaScript. Choose modest limits
for interactive work.

```python
from matrix import CSRMatrix
A = CSRMatrix(2, 2, [0, 2, 4], [0, 1, 0, 1], [4, 1, 1, 3])
result = A.conjugate_gradient([6, 7], jacobi=True, capture=True)
assert result.converged
# result.x ≈ [1, 2]; inspect result.residuals and result.reason.
```

## Comparison and live report

`python3 benchmarks/sparse.py --require-all` builds and validates every port,
then measures dense/CSR products and CG with/without Jacobi on the same
five-point SPD diffusion systems. `--no-build --verify-only` runs just the 39
shared checks. The checks include analytic answers, invalid storage, zero and
empty cases, breakdown/limit/range cases, and independent residual checks of
every captured frame. Timed solutions are checked outside the timers too.

The report uses logical storage (float64 coefficients, 64-bit indices), **not
measured RSS**. Object overhead, solver vectors, frame capture and allocators
are excluded. Timers include per-call validation, allocation, copies and
checksums; construction, dense conversion and process startup are excluded.
Fresh-process samples have three warmup calls. Interpret JIT/GC/compiler effects
and small-workload measurements accordingly; this is not a language ranking.

`python3 benchmarks/sparse.py --render-only benchmarks/reports/sparse/results.json`
regenerates HTML without changing measurements. The portable HTML embeds a
sanitized WASM bundle, reruns all shared sparse fixtures in a worker, and solves
a bounded heat-equilibrium problem. Animation shows solver iterates, not
physical time. A low iteration limit is displayed as an unconverged estimate.

## Restarted GMRES

GMRES accepts a square real matrix and finite matching right-hand side, including
nonsymmetric and indefinite systems. The initial guess is zero. Restart length
(default 30) must be an integer from 1 to 1024; it is capped by matrix dimension
and the total iteration limit. Tolerances, iteration limits and capture defaults
match CG. Optional **right Jacobi** preconditioning requires a stored nonzero
diagonal; negative diagonal entries are allowed. Missing or zero diagonal
entries are input errors.

The implementation applies two-pass modified Gram–Schmidt to `A * (v / diag)`,
uses Givens rotations to update the least-squares problem, and computes each
candidate as `x_base + (V * y) / diag`. It recomputes `norm(b - A*x)` after every
accepted update. Only that true residual can establish convergence. The
projected residual estimate is also recorded; it is in the original system's
norm because preconditioning is on the right. Finite precision can separate
these two histories.

The result contains `x`, `converged`, `iterations`, `reason`, true `residuals`,
`estimated_residuals`, optional `iterates`, and `restarts` (with each port's
normal casing). Histories include the initial residual, and captured iterates
include the initial zero vector. Restart indices identify the accepted-iteration
count at the start of each subsequent cycle, excluding the initial cycle. A
failed cycle can therefore start at the final accepted iteration count.

Stop reasons are `converged`, `iteration_limit`, `breakdown`, `nonfinite`, and
`stagnation`. A negligible Arnoldi remainder is a happy breakdown only when the
true residual satisfies the requested tolerance. An unchanged solution after a
completed restart cycle reports stagnation, including at the iteration limit.
A nonsingular matrix does not guarantee convergence with a short restart.
Runtime arithmetic failure returns the last accepted result; an unrepresentable
initial norm reports `nonfinite` and can leave infinity in the initial history.
Invalid inputs throw or return an input-error status instead.

| Port | Method | Options / ownership |
| --- | --- | --- |
| C / ARM64 | `m_csr_gmres` | `m_gmres_free` releases a zero-initialized result |
| C++ | `CSRMatrix::gmres` | returned vectors own their storage |
| C# / F# | `CSRMatrix.Gmres` | named arguments |
| Go | `CSRMatrix.GMRES` | start with `DefaultGMRESOptions()` |
| Rust | `CSRMatrix::gmres` | `GMRESOptions::default()` |
| Python | `CSRMatrix.gmres` | keyword arguments |
| TypeScript / WASM | `CSRMatrix.gmres` | options object; copied result arrays |
| Julia | `gmres(a, b)` | keyword arguments |

Restarted GMRES stores a basis of at most `restart + 1` vectors and a small
Hessenberg matrix. The report's logical workspace model counts float64 basis,
Hessenberg, rotations, and work-vector entries. It excludes CSR storage, result
histories, frame capture, objects, allocator overhead, and transient copies; it
is not measured resident memory. Capture adds one full solution per accepted
iteration and is disabled in timed benchmarks.

The transport example discretizes diffusion plus advection using a five-point
stencil and first-order upwinding. Direction is clockwise from right on the
screen; the top boundary is 100 and the other boundaries are zero. Animation
shows solver iterates, not physical time. The recorded benchmark instead uses
zero boundaries and a manufactured `x = ones` solution, with `b = A*x`.

The sparse runner protocol appends `RESTART` after the existing arguments for
operation `gmres`. Verification output is one packed row:
`reason, iterations, history_length, x..., true_residuals..., estimated_residuals..., restart_count, restart_indices..., captured_frames...`.
Reason codes follow the five stop reasons above in order, starting at zero.

## Reusable ILU(0) preconditioning

ILU(0) performs incomplete LU elimination on the **stored CSR pattern**. It drops
fill outside that pattern and retains explicit zero entries as potential factor
entries. It never pivots, reorders, shifts the diagonal, or silently substitutes
another preconditioner. A stored diagonal is required in every row. An exactly
zero computed pivot or nonfinite factor arithmetic fails construction, even if
the original matrix is nonsingular. Negative pivots are allowed. Very small
nonzero pivots are accepted, but can make application overflow or convergence
poor. This is not a general-purpose robust sparse factorization.

The factor owns a snapshot and can be applied repeatedly to different finite
vectors. `apply(b)` performs a unit-lower forward solve followed by an upper
backward solve, returning a new vector. It does not solve the original system
exactly unless the retained pattern admits the full LU factorization. Use it as
a **right preconditioner** in GMRES; true residuals still refer to the original
matrix. Factors need matching dimension, but can deliberately approximate a
different matrix of that dimension. Jacobi and ILU cannot both be selected.

| Port | Factor construction / application | GMRES integration |
| --- | --- | --- |
| C / ARM64 | `m_ilu0_create`, `m_ilu0_apply`, `m_ilu0_free` | `m_csr_gmres_preconditioned` |
| C++ | `linear_a::ILU0(a)`, `apply(b)` | final `const ILU0*` argument to `gmres` |
| C# / F# | `ILU0(a)`, `Apply(b)` | named `preconditioner` argument |
| Go | `NewILU0(a)`, `Apply(b)` | `GMRESOptions.Preconditioner` |
| Rust | `ILU0::new(&a)`, `apply(&b)` | `gmres_preconditioned(&b, options, &factor)` |
| Python | `ILU0(a)` from `matrix`, `apply(b)` | `preconditioner=factor` |
| TypeScript / WASM | `new ILU0(a)`, `apply(b)` | `{preconditioner: factor}` |
| Julia | `ILU0(a)`, `ilu_apply(factor,b)` | `preconditioner=factor` |

C output factors must be zero-initialized and freed before reuse; failures leave
output ownership unchanged. C reports `M_SINGULAR` for a zero pivot and
`M_SOLVER_RANGE` for nonfinite arithmetic. Do not shallow-copy owning C structs.
C++ factors are noncopyable and free their storage on destruction. WASM factors
require `dispose()`; disposing the original matrix does not invalidate a factor.
Returned vectors remain independent after disposal. Julia's array fields and
C's struct fields must not be mutated; public application revalidates their
invariants. Go's zero-valued factor is not initialized; construct it with
`NewILU0`. During GMRES, application overflow returns `nonfinite` with the last
accepted result. Standalone application fails explicitly.

Factorization uses binary searches in sorted CSR rows; no dense n×n workspace
is allocated. Application is linear in the stored entry count, including
validation. The report uses a portable logical factor-storage model with
float64 values and 64-bit indices: `16*nnz + 16*n + 8` bytes (values, column
indices, row offsets, and diagonal positions). Actual index widths and runtime
object overhead differ. Setup time, reused solve time and measured one-shot
time are reported separately. Their independently sampled medians need not add
up. The repeated-RHS chart estimates `setup + count * reused_solve`; it does not
claim to measure a multi-RHS batch or predict every new RHS's convergence.

The live demo caches one ILU factor. Boundary-temperature changes modify only
the RHS; restart changes modify only the solver. Both reuse the factor. Matrix
coefficient changes rebuild it. Live setup/solve timings describe this browser,
include captured frames, and are separate from the saved benchmark timings.

The sparse benchmark protocol accepts `ilu_setup` and `ilu_apply` with the
ordinary ten arguments. Setup returns `[dimension, nnz]`; application returns
the preconditioned vector. For `gmres`, the preconditioner field formerly named
`JACOBI` is 0 (none), 1 (Jacobi), 2 (saved ILU), or 3 (ILU built for each solve).
CG retains only 0 and 1. Reused setup occurs outside the timed region. The
versioned package examples exercise factor reuse across distinct RHS vectors.

## Reverse Cuthill–McKee and permutations

All eleven ports expose RCM ordering, symmetric CSR permutation, and vector
permutation. These operations are independent of ILU and can be composed with
any solver. They require square matrices; numeric symmetry is not required.

The returned permutation is **new-to-old**: `p[new] = old`. Form
`B[i,j] = A[p[i],p[j]]`, `c[i] = b[p[i]]`, solve `By=c`, and restore
`x[p[i]] = y[i]`. The inverse vector operation performs that last step. Both
matrix permutations and vectors own their outputs; the input is unchanged.
Stored zeros are preserved, rows are sorted canonically, and `nnz` is unchanged.
Permutation inputs must contain each index exactly once and match the dimension.
Vector values must be finite. Empty square matrices and empty permutations work.

RCM constructs the undirected **union of the stored pattern and its transpose**,
including explicit off-diagonal zeros, ignoring the diagonal and duplicate edges.
Start each component at the unvisited vertex minimizing `(degree, original index)`;
visit unvisited neighbors in that same order with a FIFO queue. Reverse the whole
traversal, including the component order. This deterministic minimum-degree seed
variant does not search for pseudo-peripheral vertices. It can differ from other
RCM implementations. It uses sparse adjacency storage, never an n-by-n graph.
See the [Cuthill–McKee description](https://www.boost.org/doc/libs/1_70_0/libs/graph/doc/cuthill_mckee_ordering.html)
for the degree-ordered breadth-first traversal.

| Port | RCM | Symmetric permutation | Vector permutation / inverse |
| --- | --- | --- | --- |
| C#, F# | `a.ReverseCuthillMcKee()` | `a.PermuteSymmetric(p)` | `CSRMatrix.PermuteVector(p, x, inverse)` |
| Rust | `a.reverse_cuthill_mckee()` | `a.permute_symmetric(&p)` | `CSRMatrix::permute_vector(&p, &x, inverse)` |
| Go | `a.ReverseCuthillMcKee()` | `a.PermuteSymmetric(p)` | `matrix.PermuteVector(p, x, inverse)` |
| Python | `a.reverse_cuthill_mckee()` | `a.permute_symmetric(p)` | `CSRMatrix.permute_vector(p, x, inverse)` |
| TypeScript, WASM | `a.reverseCuthillMcKee()` | `a.permuteSymmetric(p)` | `CSRMatrix.permuteVector(p, x, inverse)` |
| Julia | `reverse_cuthill_mckee(a)` | `permute_symmetric(a, p)` | `permute_vector(p, x; inverse=false)` |
| C++ | `a.reverse_cuthill_mckee()` | `a.permute_symmetric(p)` | `CSRMatrix::permute_vector(p, x, inverse)` |
| C, ARM64 | `m_csr_rcm(a, p, count)` | `m_csr_permute(a, p, count, out)` | `m_permute_vector(p, count, x, inverse, out)` |

Julia permutation values are zero-based, matching its CSR column indices.
C output arrays must have `count` elements. Vector input and output may alias;
invalid permutations leave the output unchanged. CSR output must be zero-initialized
and freed with `m_csr_free`. The C++ and ARM64 ports use the shared C ordering
kernel; ARM64 does not claim a handwritten assembly RCM implementation. WASM
ordering and permutations execute in the compiled C kernel; its adapter copies
CSR arrays back to JavaScript to construct the owned result. Dispose WASM CSR
results when finished. Julia validates its exposed arrays again at entry.

The report defines structural bandwidth as `max(abs(row-column))` over stored
entries, or zero for an empty pattern. Diagonal-only matrices have bandwidth zero.

RCM targets bandwidth, not optimal fill, stability, or convergence. ILU(0) can
become better or worse after reordering, or encounter a zero pivot. There is no
fallback to a different ordering. The [ordering report](https://jmcummings77.github.io/linear-A/ordering/)
measures these effects and includes preprocessing in its total-time workloads.

Sparse runner operations `rcm` and `permute` return respectively the permutation
and packed CSR offsets/indices/values. `permute` takes the permutation in the vector
slot. `permutation_check` additionally exercises forward/inverse vector mapping
and the matvec identity. `ilu_solve` and `rcm_solve` return `[iterations, ...x]` in
original coordinates, using ILU(0), restart 20 and the protocol tolerances/limit.
Their timed paths include setup and solve; `rcm_solve` also includes ordering,
CSR/RHS permutation and restoring the solution. All use the existing ten-argument
sparse protocol, with Jacobi and capture set to zero for ordering benchmarks.

## Sparse Cholesky

For direct SPD solves, all ports provide reusable symbolic analysis, numerical Cholesky factors and triangular solves. See [Sparse Cholesky](CHOLESKY.md) for fill semantics, APIs and numerical limits.

## Reusable IC(0) for conjugate gradient

Every port also supports a fixed incomplete Cholesky preconditioner for CG.
See [IC(0)](IC0.md) for its zero-fill pattern, positive-pivot failures, ownership,
API names and ordering semantics. Jacobi and IC(0) are mutually exclusive.

For `cg`, the benchmark protocol's `JACOBI` field accepts 0 (none), 1 (Jacobi),
2 (prepared IC(0)) and 3 (IC(0) built for each solve). `ic0_factor` returns packed
lower CSR offsets, indices and values; `ic0_apply` reuses a prepared factor.
The original CG flags and entry points retain their behavior.
