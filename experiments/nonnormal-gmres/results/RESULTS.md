# Nonnormal GMRES experiment

Local correctness and convergence diagnostics; no performance measurements.

All effective matrices B = [[1, 2s], [0, 1]] have eigenvalues 1, 1. Their numerical ranges are disks centered at 1 with radius |s|. The SVD independently checks ||I-B||₂ = 2|s|, attaining the Crouzeix constant for nonzero s.

| s | SVD norm of I−B | Supremum of ∣1−z∣ on W(B) |
| --- | ---: | ---: |
| 0 | 0 | 0 |
| 0.25 | 0.5 | 0.25 |
| 0.75 | 1.5 | 0.75 |
| 2 | 4 | 2 |
| 4 | 8 | 4 |

| s | Right Jacobi | Restart | Steps | Stop | Final true residual / initial |
| --- | --- | ---: | ---: | --- | ---: |
| 0 | no | 1 | 1 | converged | 0 |
| 0 | yes | 1 | 1 | converged | 0 |
| 0 | no | 2 | 1 | converged | 0 |
| 0 | yes | 2 | 1 | converged | 0 |
| 0.25 | no | 1 | 13 | converged | 4.76416e-11 |
| 0.25 | yes | 1 | 13 | converged | 4.76416e-11 |
| 0.25 | no | 2 | 2 | converged | 0 |
| 0.25 | yes | 2 | 2 | converged | 0 |
| 0.75 | no | 1 | 24 | iteration_limit | 1.08411e-07 |
| 0.75 | yes | 1 | 24 | iteration_limit | 1.08411e-07 |
| 0.75 | no | 2 | 2 | converged | 4.96507e-16 |
| 0.75 | yes | 2 | 2 | converged | 4.96507e-16 |
| 2 | no | 1 | 24 | iteration_limit | 0.970023 |
| 2 | yes | 1 | 24 | iteration_limit | 0.970023 |
| 2 | no | 2 | 2 | converged | 0 |
| 2 | yes | 2 | 2 | converged | 0 |
| 4 | no | 1 | 11 | stagnation | 0.992276 |
| 4 | yes | 1 | 11 | stagnation | 0.992276 |
| 4 | no | 2 | 2 | converged | 1.79018e-15 |
| 4 | yes | 2 | 2 | converged | 1.79018e-15 |

## Interpretation

- Restart 2 retains the degree-two annihilating polynomial (1−z)²; the bounded test systems solve in at most two steps.
- Restart 1 discards that Krylov space. A common spectrum does not imply common convergence. Iteration-limit and stagnation results are retained.
- Right Jacobi uses A = B·diag(1,8), so the effective operator is A·diag(1,8)⁻¹ = B. Residuals are checked in the original system b−Ax.
- Each captured step is checked against the per-cycle candidate bound ||r_start||·min(1, 2|s|ʲ), where j is the step within that cycle. This is an exact-arithmetic comparison with an explicit scale-dependent rounding allowance, not a certified error bound or stopping criterion. For |s| ≥ 1 the disk contains zero and this envelope is uninformative.
- The first projections are independently checked by exact rational minimal-residual line searches. Full trajectories, bounds, allowances, inputs and source hashes are in results.json.

Source state: `3c581d3c188e6039f6798fafc91579eb60fc8cfc`, dirty: `False`. Convergence diagnostics with no timings. Individual hashes do not archive uncommitted source; published captures also require provenance.json.
