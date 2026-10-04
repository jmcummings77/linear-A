# linear-A

## One-minute introduction

Different inputs, correctness checks, and timing rules make numerical code hard
to compare. linear-A puts eleven matrix implementations through the same
mathematical contract and benchmark protocol so you can study how algorithms
and implementation choices affect the results. It is a learning project.

For example, every implementation can multiply a rectangular matrix by its
transpose and compute the determinant:

```text
A = [1 2 3]       A Aᵀ = [14 32]       det(A Aᵀ) = 54
    [4 5 6]              [32 77]
```

Two engineering decisions shape the project:

1. **One shared float64 contract, independent correctness references.** The
   implementations use the same inputs and operation rules; no language serves
   as the correctness oracle for another.
2. **Explicit numerical kernels and dependencies.** The library implementations
   do not call BLAS or NumPy. Shared C kernels and language-specific ownership
   rules are documented so comparisons account for what actually runs.

Run the example with **Python 3.9+**, from the repository root; no packages are
needed:

```sh
PYTHONPATH=ports/python python3 - <<'PY'
from matrix import Matrix

a = Matrix(2, 3, [1, 2, 3, 4, 5, 6])
gram = a.multiply(a.transpose())
print(gram.row(0), gram.row(1))  # [14.0, 32.0] [32.0, 77.0]
print(gram.determinant())       # 54.0
PY
```

For another language or shell, follow [Getting started](docs/getting-started.md).
To explore without installing a toolchain, [open the reports](https://jmcummings77.github.io/linear-A/).

## Where to go next

| I want to… | Guide |
| --- | --- |
| Use an implementation or run its tests | [Getting started](docs/getting-started.md): all eleven languages and their toolchains |
| Understand the algorithms and tradeoffs | [Design and numerical scope](docs/design.md) |
| Look up shapes, ownership, errors, and APIs | [Shared API contract](ports/README.md) |
| Explore results and interactive examples | [Reports and exploration](docs/reports.md) |
| Run comparisons, profile, or capture benchmarks from committed source | [Benchmark guide](benchmarks/README.md) |
| Install a release kit | [Release kits and installation](release/README.md) |
| Contribute a change | [Contributing](CONTRIBUTING.md) |

[Versioned API contract](docs/api/0.10.0/README.md) ·
[Earlier API contracts](docs/api/) ·
[License](LICENSE)
