# Approximate minimum degree ordering

All eleven ports expose a deterministic ordering of a square CSR matrix's stored
pattern. Values are irrelevant, including stored zeros. The graph is the union
of the pattern and its transpose, without diagonal edges. The result uses the
same convention as RCM: `p[new] = old`, with zero-based indices in every port.
Inputs are not modified. Empty square matrices return an empty permutation;
rectangular matrices are rejected. Ordering alone does not require symmetry or
positive definiteness; Cholesky still requires both its numerical conditions.

| Port | Entry point |
|---|---|
| C#, F# | `a.ApproximateMinimumDegree()` |
| Rust | `a.approximate_minimum_degree()` |
| Go | `a.ApproximateMinimumDegree()` |
| TypeScript, WebAssembly | `a.approximateMinimumDegree()` |
| Python | `a.approximate_minimum_degree()` |
| Julia | `approximate_minimum_degree(a)` |
| C++, C, ARM64 | `a.approximate_minimum_degree()` (C++), `m_csr_amd(&a, p, n)` (C/ARM64) |

The C, C++, ARM64 and WebAssembly sparse ordering implementations share the C
kernel. ARM64 assembly is used for other arithmetic kernels, not this graph
routine. Other ports implement the algorithm in their own language. Returned
permutations are owned copies. In C the caller supplies `n` output entries and
checks the returned status; allocation failure can leave a partial permutation.
Do not alias the output with the input CSR arrays.

For example, the Python API composes ordering with a reusable factor:

```python
p = a.approximate_minimum_degree()
q = a.permute_symmetric(p)
plan = SparseCholeskySymbolic(q)
factor = plan.factorize(q)
y = factor.solve(CSRMatrix.permute_vector(p, b))
x = CSRMatrix.permute_vector(p, y, inverse=True)
```

Reuse `p` and `factor` for new right-hand sides; rebuild the numerical factor
when values change, and rebuild the symbolic plan when the stored pattern changes.

## Quotient graph and degree bound

Direct adjacency sets store original edges. Elements store implicit cliques
created by previous eliminations. Initially every degree is exact. Select the
active vertex with the smallest bound, capped by the number of other active
vertices; break ties by the smallest original vertex index.

At pivot `k`, form its exact remaining neighbor set `L` by uniting direct
neighbors with the members of every incident element. Absorb those incident
elements into one new element on `L`, remove `k`, and remove direct edges within
`L` because the new element represents them. For each `i` in `L`, update:

```
bound(i) = |L| - 1 + |direct(i)|
           + sum(|element(e) \ L| for other elements e incident on i)
```

Cap the result at the number of other active vertices. External sets can overlap,
so this is an upper bound rather than an exact degree. Bounds of unaffected
vertices remain valid. This removes repeated set unions from degree updates;
it does not explicitly insert every fill edge. Selection and element lookup use
scans, so this educational implementation is not optimized for huge graphs.

This is a transparent AMD variant, not a reproduction of SuiteSparse AMD.
It omits supervariable merging, aggressive absorption, dense-vertex handling and
assembly-tree postordering. It does not guarantee minimal fill, or a faster
end-to-end solve. See the [SuiteSparse AMD contract](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/v7.10.1/AMD/Include/amd.h)
for that implementation's heuristics and permutation convention.

## Verification and comparison

`benchmarks/amd_reference.py` materializes the filled graph independently,
reconstructs uncovered original edges, and checks every estimate against the
exact active degree at each step. Seeded random tests exercise strict upper
bounds, deterministic ties, disconnected graphs and asymmetric patterns. Shared
fixtures verify every port and solve the original system after AMD permutation.
The browser independently replays every fill edge and compares it with the
WebAssembly symbolic factor's history; displayed steps are algorithm frames,
not measured execution-time frames or a flamechart.

Run `python3 benchmarks/amd_bench.py --require-all --sizes 8 --samples 3` for grids, scrambled grids,
binary trees and irregular graphs. Add `--amd-library /path/to/libamd` to compare
an installed SuiteSparse AMD library. The library version is recorded, but the
local path is not exported. Reference ordering timings include Python ctypes
conversion, allocation, permutation validation and checksum; they are not pure
C-kernel timings. Reference samples run separately in one Python process; port samples use fresh runner processes. Reference fill is measured by explicit elimination of its
returned permutation, not by its approximate internal statistics. No SuiteSparse
source or binary is bundled into the library.

The report separately measures ordering, symbolic analysis, numerical
factorization, a prepared solve, and total solve cost (including permutations).
Natural order needs no ordering call. Logical factor storage uses 64-bit values
and indices, excluding temporary graphs and runtime overhead. Existing Cholesky
reports remain historical artifacts; the expanded comparison lives at
[`amd/`](https://jmcummings77.github.io/linear-A/amd/).
