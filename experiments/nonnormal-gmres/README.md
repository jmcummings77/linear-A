# Nonnormal GMRES and graph workloads

This experiment turns two connections from OpenAI's October 2026 mathematics
collection into small, independently checkable numerical examples. It uses the
existing Python numerical kernels and adds cross-port conformance fixtures.
It does not implement a new solver or the paper's large graph construction.

[Open the HTML report](https://jmcummings77.github.io/linear-A/nonnormal-gmres/) ·
[Saved offline report](results/index.html) · [Captured data](results/results.json)

The report presents the results, lets you inspect residual histories by restart
length and right preconditioner, and explains the numerical-range and graph
checks. It is self-contained and works offline.

## Run the convergence experiment

Python 3.9+ and the standard library are sufficient, from the repository root:

```sh
python3 experiments/nonnormal-gmres/run.py
python3 experiments/nonnormal-gmres/run.py --verify-only
python3 -m unittest discover -s tests/benchmark_harness -p 'test_nonnormal_gmres.py' -v
```

The default local output is `.build/nonnormal-gmres/`: open `index.html` for the
interactive report or `RESULTS.md` for its text companion. All captured iterates,
true and projected residuals, restart boundaries, analytic envelopes, inputs
and source hashes are in `results.json`. Use `--output PATH` for a separate
run and `--limit 2` to compare the first two steps. The limit is bounded to 2–100
and defaults to 24. Verification-only mode writes no report.

These are convergence diagnostics, with no timings or language ranking.
The generated report records the actual source revision and dirty state. A
source hash does not archive uncommitted files. To replace the published report,
commit the source and start from a clean checkout, then capture the study in a
fresh clone:

```sh
python3 benchmarks/reproduce.py --runner nonnormal-gmres \
  --output .build/publishable/nonnormal-gmres -- --limit 24
```

Review the capture, replace the contents of `experiments/nonnormal-gmres/results/`
with the complete capture including `provenance.json`, and run
`python3 benchmarks/check_provenance.py`. The saved provenance records the source
commit, Git tree, invocation, and artifact hashes. Local development runs with
working changes remain separate from this published capture.

## Same eigenvalues, different GMRES behavior

Set

```text
N = [0 2]       B = I + sN = [1 2s]       b = [0]
    [0 0]                    [0  1]           [1]
```

The experiment uses `s = 0, 1/4, 3/4, 2, 4`. Every B has eigenvalues 1 and 1.
Its complex numerical range is the disk centered at 1 with radius |s|: for a
unit vector v, `v* N v = 2 conjugate(v1) v2`, whose possible values fill the unit
disk. For `p(z) = 1-z`, the exact operator norm `||p(B)||₂ = 2|s|`, while the
maximum of `|p(z)|` on that disk is |s|. The existing SVD computes the first
quantity and is checked against the formula. The general eigensolver is checked
against the repeated eigenvalue. Multiplication verifies `(I-B)² = 0`.

This is the sharp example for the scalar Crouzeix inequality
`||p(B)||₂ ≤ 2 max(z in W(B)) |p(z)|`. The release's complete inequality extends
to matrix-valued polynomials; this single-right-hand-side experiment only uses
its scalar specialization. It demonstrates the connection, not a newly enabled
GMRES optimization or an independent verification of the general theorem.

GMRES(2) retains the degree-two polynomial that annihilates B and solves these
bounded examples in at most two steps. GMRES(1) starts over after every step;
it can take many more iterations or stagnate even though the spectrum is
unchanged. The report preserves iteration-limit and stagnation outcomes.

For the worked example `s=1/4`, the first iterate is `x1=(0,4/5)` and the
residual is `r1=(-2/5,1/5)`. Restarting gives `x2=(-32/65,68/65)` and
`r2=(-2/65,-3/65)`. Without that restart, the second iterate solves
`x=(-1/2,1)`. The experiment's independent oracle obtains the first two
restarted iterates with exact rational line minimizations:
`alpha = (rᵀBr)/||Br||₂²`, `x_next=x+alpha*r`, `r_next=r-alpha*Br`.
It does not use Arnoldi, Givens rotations, or another language's solver.

Each captured step is compared with the candidate-polynomial bound
`||r_cycle_start||₂ min(1, 2|s|^j)` for the step j **within that cycle**.
The no-progress bound comes from the constant polynomial 1. The disk bound is
deliberately loose at j=2, where the actual operator polynomial vanishes.
When |s|≥1 the disk contains zero and this envelope is uninformative.
The comparison uses the measured cycle-start residual; it is not a global
unrestarted degree-k optimality claim.

Checks allow `256 * epsilon * (1 + |x0| + |2s D x1| + |D x1|)` absolute rounding
error in residual comparisons, where D is 1 or 8. This conservative allowance
is scoped to these tiny, bounded systems and the independent expression for
`b-Ax`. It is not a certified general finite-precision bound. The solver's
existing true-residual stopping rule remains authoritative.

## Right preconditioning

Every run is repeated with `A = B diag(1,8)` and right Jacobi. The actual Krylov
operator is `A diag(1,8)⁻¹ = B`, so its numerical range supplies the same
envelope. Residuals are always checked in the original coordinates `b-Ax`;
the physical solution is `(-2s,1/8)`. No numerical range of the unpreconditioned
A is substituted for the one of B.

The shared GMRES fixtures also check the identity and `s=1/4,1,4`, at restart
lengths 1 and 2, both without preconditioning and with right Jacobi. They test
the first projections, true/projected residuals and termination counts using
closed formulas. Negative harness tests ensure that a wrong iterate with a
self-consistent residual history is still rejected.

## Sparse graph controls

The companion graph fixtures use the classical Petersen graph and triangular
prism. Both are cubic, connected, nonbipartite Ramanujan graphs; their complete
adjacency spectra are known exactly. For either graph adjacency H, use the
shifted Laplacian `A=(3+mu)I-H` with `mu=1/4`. Its positive energy is
`mu sum(x_i²) + sum_edges (x_i-x_j)²`, and its constant-mode eigenvalue is mu.
The unshifted Laplacian would be singular. Spectral-gap claims on the zero-mean
subspace must not be substituted for the full-space condition number.

The scalable control workload is the prism `C_m × K_2`, with 2m vertices and
degree 3. Its adjacency eigenvalues are `2 cos(2 pi k/m) ± 1`. Larger prisms
are **not** claimed to be Ramanujan graphs; they provide a known-spectrum
comparison as the nonconstant gap closes. Graph edge formulas independently
manufacture the RHS, and unit tests verify graph energy and analytic spectra.

```sh
# Shared sparse correctness checks, including the new graph and GMRES cases.
python3 benchmarks/sparse.py --verify-only --languages python c --require-all

# Local graph timings; sizes are vertices per prism ring, not grid widths.
python3 benchmarks/sparse.py --workload prism --sizes 8 16 32 \
  --languages python c --require-all --output .build/prism-benchmark
```

The existing diffusion workload remains the default. Published timing captures
still follow the repository's [clean-source capture workflow](../../benchmarks/README.md#publishing-reproducible-source-records).
No saved public measurements are replaced by these commands. These graphs are
general CSR inputs, outside the structured-grid multigrid contract.

## Sources and scope

- [A direct proof of the complete Crouzeix inequality, family 325](https://github.com/openai/math/blob/main/preprints/A-direct-proof-of-the-complete-Crouzeix-inequality-September-26-2026/paper.pdf)
  and its [formalization scope](https://github.com/openai/math/blob/main/lean/docs/325.md).
- [Deterministic nonbipartite Ramanujan graphs in every fixed degree, family 178](https://github.com/openai/math/blob/main/preprints/Deterministic-nonbipartite-Ramanujan-graphs-in-every-fixed-degree-September-23-2026/paper.pdf).
  The new paper motivates controlled spectral workloads. Our classical small
  graphs and prism controls do not depend on its new construction or proof.
- [Local GMRES contract](../../ports/SPARSE.md#restarted-gmres) and
  [independent graph references](../../benchmarks/graph_reference.py).
