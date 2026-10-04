# Geometric multigrid

All eleven ports provide a reusable symmetric V-cycle for the unit,
constant-coefficient five-point Dirichlet Laplacian. This is a structured-grid
preconditioner, not an algebraic multigrid implementation.

The interior grid width must be `2^k - 1`, from 1 through 255. There are `width²`
unknowns, stored row by row. The boundary outside the interior grid is zero.
The generated CSR operator has diagonal 4 and neighbor entries −1. Every level
uses this same **unscaled** stencil; no hidden mesh-spacing factor is included.
The hierarchy, transfer weights and operators are implicit. Setup retains the
width; each application allocates its workspaces and starts from zero.

| Port | Create | Generated CSR | Apply | CG integration |
|---|---|---|---|---|
| C# | `new GeometricMultigrid(width)` | `.Matrix` | `.Apply(b)` | `preconditioner: mg` |
| F# | `GeometricMultigrid(width)` | `.Matrix` | `.Apply(b)` | `preconditioner=mg` |
| Rust | `GeometricMultigrid::new(width)?` | `.matrix()?` | `.apply(&b)?` | `conjugate_gradient_preconditioned(..., Some(&mg))` |
| Go | `NewGeometricMultigrid(width)` | `.Matrix()` | `.Apply(b)` | `CGOptions{Preconditioner: mg}` |
| TypeScript | `new GeometricMultigrid(width)` | `.matrix` | `.apply(b)` | `{preconditioner: mg}` |
| Python | `GeometricMultigrid(width)` | `.matrix` | `.apply(b)` | `preconditioner=mg` |
| Julia | `GeometricMultigrid(width)` | `multigrid_matrix(mg)` | `multigrid_apply(mg,b)` | `preconditioner=mg` |
| C / ARM64 | `m_multigrid_create(width,&mg)` | `m_multigrid_matrix(&mg,&a)` | `m_multigrid_apply(&mg,&b,&x)` | `m_csr_cg_multigrid(...)` |
| C++ | `GeometricMultigrid(width)` | `.matrix()` | `.apply(b)` | `conjugate_gradient_multigrid(...,&mg)` |
| WebAssembly | `new api.GeometricMultigrid(width)` | `.matrix` | `.apply(b)` | `{preconditioner: mg}` |

The ARM64 port uses the shared C multigrid kernel; it does not have a separate
handwritten assembly V-cycle. C++ and WebAssembly also use that kernel.
Other ports implement the cycle in their own languages.

C objects must be zero-initialized before creation. The multigrid object itself
contains no heap allocation and needs no free call. Its generated CSR matrices
and apply outputs are owned and use the existing free functions. WebAssembly
objects have an idempotent `dispose()`; subsequent use fails. Generated matrices
and returned vectors are independent results. Applications do not modify the
RHS. Preserve the width for the lifetime of a preconditioner.

## One V-cycle

1. Start with `x = 0` at the current level.
2. On the one-point grid, return `x = b / 4`.
3. Run two weighted-Jacobi sweeps: `x ← x + (b − A x) / 6`.
4. Restrict the residual: `b_coarse = Pᵀ(b − A x)`.
5. Recursively solve the coarse correction equation.
6. Add its bilinear interpolation: `x ← x + P e_coarse`.
7. Run the same two Jacobi sweeps and return `x`.

Coarse nodes sit at odd zero-based fine-grid indices. Interpolation uses tensor
products of weights `[1/2, 1, 1/2]` around each coarse node. Restriction is exactly
the transpose of that interpolation. It is **not** normalized full weighting
`Pᵀ/4`: with the unscaled operator, that normalization would require a matching
coarse-equation rescaling. Coarse operators are rediscretized five-point
operators, **not Galerkin products** `Pᵀ A P`.

The fixed V-cycle is symmetric positive definite. To see why, let
`S = I − A/6` be the one-sweep error operator on a level. Two pre-sweeps and two
post-sweeps give the inverse-like map

`B = A⁻¹(I − S⁴) + S² P B_coarse Pᵀ S²`.

For the Dirichlet stencil, all eigenvalues of A lie strictly between 0 and 8.
Thus the eigenvalues of S lie strictly between −1/3 and 1, making the first
term positive definite. The second term is symmetric positive semidefinite
when the coarse map is SPD. The one-point inverse is SPD, completing the
induction. This argument does not require Galerkin coarse operators.

CG still requires an SPD system matrix. A matching-sized structured
preconditioner does not certify a caller's arbitrary matrix. Supplying both
Jacobi and an explicit preconditioner is rejected. C#, F#, Rust and Go expose
an SPD preconditioner interface/trait; custom implementations must be fixed,
linear, positive definite, preserve their input and return the matching size.

## Trace and report

Python, TypeScript and WebAssembly expose `trace(b)`, returning `{x, frames}`.
Each frame stores `level`, `width`, `phase`, `x`, `b`, and the actual residual
`b − A_level x`. Phases are `enter`, `pre_smooth`, `coarse_solve`, `correct`,
and `post_smooth`. C's header-level `la_mg_apply_trace` callback drives the
WebAssembly frames. A callback error cancels the operation without publishing
a partial output. The ordinary apply path does not store frames.

On the fine grid, the report plots error against a known exact solution.
On coarse grids it plots the correction equation's current x. Coarse correction
is not the physical solution or simply a downsampled fine-grid error. A captured
V-cycle is one application to the initial RHS, separate from the full CG run.
Animation time does not represent execution time or a sampled flamechart.

Run local conformance and scaling checks with:

```sh
python3 benchmarks/multigrid_bench.py --sizes 7 15 31 --samples 3
```

For a public snapshot, commit the source first and use the clean-capture workflow:

```sh
python3 benchmarks/reproduce.py --runner multigrid_bench.py --output .build/publishable/multigrid -- --sizes 7 15 31 --samples 3
```

The standalone report is saved under `benchmarks/reports/multigrid/` and linked
from the report directory. The suite compares plain CG, Jacobi, IC(0) and MG,
with prepared and directly measured total solves. Jacobi is scalar scaling on
this operator and does not reduce CG iterations. Saved runs use a deterministic
broad-spectrum solution; the live sine examples can converge in very few plain
CG steps because a single sine mode is an eigenvector.

Reported memory is a **logical array model**, not measured peak resident memory.
MG retained state is one width; apply workspaces are allocated per call.
IC(0) retained storage counts factor values, indices and row offsets.
Language objects, temporary arrays, allocator overhead and source CSR storage
are excluded. Read the embedded methodology for the exact accounting.

Tests use independent dense transfer matrices with exact rational arithmetic,
construct the entire preconditioner map to check symmetry and positive pivots,
validate generated stencils, exercise reuse and failures, and check the true
residual of every captured level and CG iterate.
