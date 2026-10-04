# Reports and browser exploration

[Project overview](../README.md) · [Getting started](getting-started.md)

Use the [report directory](https://jmcummings77.github.io/linear-A/) to explore
recorded results and numerical demonstrations. Reports compare particular
implementations, algorithms, toolchains, and workloads on a particular host;
they do not establish a universal language ranking.

## Dense calculations

| Report | What to explore |
| --- | --- |
| [Main comparison](https://jmcummings77.github.io/linear-A/latest/) | Operation timings, language filters, sample variation, and sampled profiles. |
| [Determinants](https://jmcummings77.github.io/linear-A/determinants/) | Automatic selection, LU, and Cholesky algorithm costs. |
| [Eigenvalues](https://jmcummings77.github.io/linear-A/eigen/) | Symmetric eigendecomposition timings and geometric interpretation. |
| [Vectors](https://jmcummings77.github.io/linear-A/vectors/) | Cross products and rotation matrix construction. |
| [Snapshot comparison](https://jmcummings77.github.io/linear-A/compare/) | Match saved workloads, inspect individual samples, and compare baseline and candidate runs. |

The [snapshot comparison guide](../benchmarks/comparison/README.md) explains
local file loading, Git-revision comparisons, and provenance warnings. Its timing
envelope describes observed samples, not a confidence interval or a statistical
regression test. Check workload, toolchain, and host differences before attributing
a ratio to a code change.

## Sparse solvers and ordering

These reports separate solver behavior from preparation costs and show how
matrix structure affects work and storage. Their linked contracts define
supported inputs, stopping rules, and numerical limits.

| Report | What to explore | Contract |
| --- | --- | --- |
| [Sparse systems](https://jmcummings77.github.io/linear-A/sparse/) | Heat equilibrium, conjugate gradient, runtime and storage scaling. | [CSR and iterative solvers](../ports/SPARSE.md) |
| [GMRES](https://jmcummings77.github.io/linear-A/gmres/) | Restart lengths and an interactive advection–diffusion solver. | [GMRES](../ports/SPARSE.md#restarted-gmres) |
| [ILU setup and reuse](https://jmcummings77.github.io/linear-A/ilu/) | Preparation costs versus repeated solves with ILU(0). | [ILU(0)](../ports/SPARSE.md#reusable-ilu0-preconditioning) |
| [Sparse ordering](https://jmcummings77.github.io/linear-A/ordering/) | Natural versus Reverse Cuthill–McKee ordering, bandwidth, and ILU/GMRES costs. | [Ordering](../ports/SPARSE.md#reverse-cuthillmckee-and-permutations) |
| [Sparse Cholesky](https://jmcummings77.github.io/linear-A/cholesky/) | Fill creation, analysis, factorization, repeated solves, and factor storage. | [Cholesky](../ports/CHOLESKY.md) |
| [AMD ordering](https://jmcummings77.github.io/linear-A/amd/) | Natural, RCM, and AMD orderings, reference fill, and an elimination graph. | [AMD](../ports/AMD.md) |
| [IC(0)](https://jmcummings77.github.io/linear-A/ic0/) | Plain CG, Jacobi-CG, IC(0), and full Cholesky with live residual fields. | [IC(0)](../ports/IC0.md) |
| [Geometric multigrid](https://jmcummings77.github.io/linear-A/multigrid/) | Structured V-cycle, IC(0), Jacobi, and plain CG scaling with live coarse-grid corrections. | [Multigrid](../ports/MULTIGRID.md) |

## Controlled experiments

- The [machine-code dot-product experiment](../experiments/machine-code-dot/README.md)
  compares a manually encoded 36-byte ARM64 function with scalar implementations.
  It is independent of the matrix API benchmarks.
- The [multiplication locality study](../experiments/matmul-locality/README.md)
  isolates loop order, cache blocking, and explicit ARM64 SIMD, with allocation
  accounting and separate native profiles. Its controlled compiler settings and
  optional vendor-library baseline matter when interpreting the results.

## Browser tools and offline use

The main report embeds WebAssembly for correctness checks and bounded benchmarks
on your device. Its editable geometry view animates matrix vector fields,
eigendirections, and captured multiplication updates. The accuracy playground
lets you perturb a system and inspect solution sensitivity, residuals, backward
error, and reciprocal condition. Browser runs and illustrative frames remain
separate from saved measurements and sampled profiles.

The [numerical playground](https://jmcummings77.github.io/linear-A/applications/)
adds draggable PCA, QR curve fitting, and low-rank image exploration. Read its
[source and numerical limits](../applications/README.md).

Open the [saved main report](../benchmarks/reports/latest/index.html) offline, or
inspect its [recorded data](../benchmarks/reports/latest/results.json).

## Reading and publishing evidence

Reports retain hardware, toolchain, and timestamp metadata. Historical captures
made from dirty working trees did not archive their exact source; a commit and
source hash cannot reconstruct it. The original dot-product capture omitted Git
provenance. Those reports retain their disclosed limitations.

New published captures require committed source and a clean checkout, build in
a fresh clone, and preserve a recoverable commit, invocation, and artifact hashes.
The [saved multigrid report](../benchmarks/reports/multigrid/provenance.json) records
its clean source commit, source tree, command, and artifact hashes.
See the [capture workflow](../benchmarks/README.md#publishing-reproducible-source-records).
The [benchmark guide](../benchmarks/README.md) contains commands and profiling
details; the [runner protocol](../benchmarks/PROTOCOL.md) defines portable inputs,
outputs, and timing boundaries. Committed reports publish to GitHub Pages without
rerunning measurements.
