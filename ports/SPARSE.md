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
