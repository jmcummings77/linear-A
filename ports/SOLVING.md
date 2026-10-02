# Linear systems and numerical accuracy

All eleven implementations provide double-precision system solving, reusable
factorizations, and least squares. Factors copy their input; editing or disposing
the source afterward does not change subsequent solves. A right-hand side is a
matrix, so one factorization can solve several columns at once.

For `A = [[4,1],[1,3]]` and `b = [[6],[7]]`, the solution is `x = [[1],[2]]`.
Factor A once, then solve again when only b changes. This saves the cubic
factorization work; each additional right-hand-side column costs quadratic work.

| Implementation | LU factor / solve | Cholesky / QR factors | Least squares / condition diagnostic |
| --- | --- | --- | --- |
| ARM64 assembly | C API below; solver control and arithmetic execute in C | C API | C API |
| C | `m_factorize(&a, M_LU, &f)`; `m_factor_solve(f, &b, &x)`; convenience `m_solve` | `M_CHOLESKY`, `M_QR` | `m_least_squares`; `m_factor_rcond` |
| C# | `a.FactorLu()`; `factor.Solve(b)`; `a.Solve(b)` | `a.FactorCholesky()`, `a.FactorQr()` | `a.LeastSquares(b)`; `factor.ReciprocalCondition()` |
| C++ | `a.factor_lu()`; `factor.solve(b)`; `a.solve(b)` | `a.factor_cholesky()`, `a.factor_qr()` | `a.least_squares(b)`; `factor.reciprocal_condition()` |
| F# | `a.FactorLu()`; `factor.Solve(b)`; `a.Solve(b)` | `a.FactorCholesky()`, `a.FactorQr()` | `a.LeastSquares(b)`; `factor.ReciprocalCondition()` |
| Go | `a.FactorLU()`; `factor.Solve(b)`; `a.Solve(b)` | `a.FactorCholesky()`, `a.FactorQR()` | `a.LeastSquares(b)`; `factor.ReciprocalCondition()` |
| Julia | `factor_lu(a)`; `solve(factor,b)`; `solve(a,b)` | `factor_cholesky(a)`, `factor_qr(a)` | `least_squares(a,b)`; `reciprocal_condition(factor)` |
| Python | `a.factor_lu()`; `factor.solve(b)`; `a.solve(b)` | `a.factor_cholesky()`, `a.factor_qr()` | `a.least_squares(b)`; `factor.reciprocal_condition()` |
| Rust | `a.factor_lu()`; `factor.solve(&b)`; `a.solve(&b)` | `a.factor_cholesky()`, `a.factor_qr()` | `a.least_squares(&b)`; `factor.reciprocal_condition()` |
| TypeScript | `a.factorLU()`; `factor.solve(b)`; `a.solve(b)` | `a.factorCholesky()`, `a.factorQR()` | `a.leastSquares(b)`; `factor.reciprocalCondition()` |
| WebAssembly | `a.factorLU()`; `factor.solve(b)`; `a.solve(b)` | `a.factorCholesky()`, `a.factorQR()` | `a.leastSquares(b)`; `factor.reciprocalCondition()` |

C outputs start empty (`matrix x = {0}` and `matrix_factor *f = NULL`). Check each
status, free x with `m_free(&x)`, and release factors with `m_factor_free(&f)`;
the latter clears the pointer and permits repeated cleanup. Failed operations
leave outputs unchanged. C++ factors are move-only; calling a moved factor throws
`std::logic_error`. WebAssembly factors own allocations: use `try/finally` and
`factor.dispose()` as with matrices. Factors and RHS matrices must belong to the
same WASM instance. Go returns `(value, error)`; Rust returns `Result`.

C# solvers are extensions on `Matrix<double>` in .NET 10. They do not implicitly
convert integers or decimals, and are not added to the legacy .NET Standard API.
F# users loading sources in scripts load `Solve.fs` after `Matrix.fs`; compiled
consumers use the `LinearA` namespace. Treat factor internals as read-only,
including Julia's array fields and Python's underscore-prefixed storage.

C and C++ share `ports/c/solve_core.h` arithmetic. ARM64 uses this C solver and
WebAssembly compiles it. Other ports implement the algorithms in their own
languages. Cross-port agreement alone is therefore not an independent proof.

## Shapes and algorithms

For A with m rows and n columns, b must have m rows and any number of columns p.
The result has n rows and p columns. Inputs and results must be finite float64.
Zero RHS columns and empty 0×0 factors are supported. QR also accepts m×0 input
and returns a 0×p result.

- **LU with partial row pivoting:** square A; rejects a zero computed pivot.
  `solve` uses LU. It does not reject every ill-conditioned matrix or estimate
  rank using a tolerance. A nonzero pivot is not a guarantee of accuracy.
- **Cholesky:** square A, exact symmetry, strictly positive computed pivots.
  Indefinite, semidefinite, and nonsymmetric inputs fail. Independently perturbing
  one off-diagonal entry may break symmetry and correctly cause rejection.
- **Column-pivoted Householder QR:** m ≥ n and numerical full column rank.
  `least_squares` minimizes the Euclidean residual independently for each RHS
  column; it does not form normal equations. Remaining column norms are
  recomputed to choose each pivot. A pivot norm at or below
  `epsilon * max(m,n) * largest_initial_column_norm` is rejected, where
  `epsilon = 2^-52` and norms are Euclidean. There is no rank-deficient or
  underdetermined minimum-norm solution API yet.

Solvers normalize A by a power of two and apply the same scale to b. This supports
uniformly tiny or huge systems, but does not remove all overflow, underflow,
cancellation, or growth limitations. Scaling that rounds a nonzero entry to zero
is rejected; nonfinite intermediate arithmetic fails. Other rounding or
underflow may still occur. LU and QR pivot decisions near thresholds can differ
slightly between runtimes. No iterative refinement is performed.

## Interpreting the diagnostics

`reciprocal_condition` (using the spelling in each API) computes
`1 / (||A||∞ ||A^-1||∞)` for a **square** factor. The matrix infinity norm is the
maximum absolute row sum. The normalized inverse is formed by solving against
identity with the reusable factors; the implementation avoids multiplying the
two potentially large norms. This diagnostic costs O(n³) time and O(n²) extra
storage. It returns one for empty input, clamps rounding above one, and returns
zero when inverse arithmetic overflows. Rectangular factors reject it.

This is a computed diagnostic, not a cheap condition estimator or a certified
error bound. An rcond close to zero indicates sensitivity; zero may also reflect
range limitations. The inverse and norm computations have rounding error.
Use factors to solve systems rather than constructing an inverse yourself.

The report's accuracy playground solves editable 2×2 systems in WASM and shows:

- The baseline and perturbed solution intersections, with automatic view fitting.
- Infinity-norm residual `||b - Ax||∞`, independently recomputed in JavaScript.
- Normwise backward error `||b - Ax||∞ / (||A||∞ ||x||∞ + ||b||∞)`.
  When both numerator and denominator are zero, it reports zero.
- Reciprocal condition and relative infinity-norm **solution change** between the
  two computed answers. This is not forward error against an exact reference.

A tiny residual can coexist with a large response to a tiny input perturbation.
The nearly parallel preset demonstrates this. Inputs recompute after a 300 ms
pause; a requested perturbation that rounds away is explicitly reported.
The playground is separate from saved benchmark timings and profiles.

For background on conditioning, forward error, and backward error, see the
[LAPACK Users' Guide](https://www.netlib.org/lapack/lug/node75.html).

## Verification

The shared suite adds 47 solver fixtures: exact rational Gauss-Jordan references,
multiple right-hand sides, row/column pivoting, independent exact-rational
least-squares references, empty shapes, invalid inputs, numerical rank rejection,
subnormal input, uniform extreme scales, and known condition numbers. Normal
equations are used only by the exact rational reference, never by production QR.
Native tests additionally exercise factor reuse, snapshots, and disposal.
The C and C++ tests can run with address and undefined-behavior sanitizers.

The shared runner supports `check solve`, `check solve_cholesky`,
`check least_squares`, and `check rcond`. Solver performance is not included in
the existing saved measurements or default timing suite. See the
[protocol](../benchmarks/PROTOCOL.md) for argument and result formats.
