# Zero-fill incomplete Cholesky and preconditioned CG

All eleven ports expose a reusable `IC0` factor for symmetric positive-definite
systems. Construction snapshots the matrix. `Apply(b)` solves `L Lᵀ z = b`
through forward and backward substitution; it does not solve `A z = b` exactly
unless the incomplete factor happens to be exact. CG accepts this factor as a
fixed preconditioner and always checks the true residual `b − A x`.

## Pattern and arithmetic

The lower pattern contains the union of the stored input pattern and its
transpose, plus every diagonal. Stored zeros count as structural entries,
including a zero stored on only one side. The numeric input must be exactly
symmetric, treating missing entries as zero. Factorization performs Cholesky
updates only on this pattern; it never inserts fill, shifts the diagonal, pivots,
or reorders implicitly. Exporting the lower factor returns an owned copy.

Missing/nonpositive computed pivots and nonfinite arithmetic fail explicitly.
An SPD matrix can still produce a nonpositive incomplete pivot: for

```
A = [1.5  1   0   1]
    [1   1.5  1   0]
    [0    1  1.5 -1]
    [1    0  -1  1.5]
```

`(A − 1.5 I)² = 2 I`, so the eigenvalues `1.5 ± sqrt(2)` are positive.
Natural-order incomplete LDLᵀ has pivots `3/2, 5/6, 3/10, −5/2` and fails.
Full Cholesky succeeds. The independent test oracle uses dense masked LDLᵀ
with exact rational arithmetic, avoiding the sparse square-root implementation. Larger timing inputs use the same independent dense LDLᵀ oracle in floating point.

## APIs

| Port | Construct / apply | CG preconditioner argument |
| --- | --- | --- |
| C# | `new IC0(a)` / `f.Apply(b)` | `a.ConjugateGradient(b, preconditioner: f)` |
| F# | `IC0(a)` / `f.Apply(b)` | `a.ConjugateGradient(b, preconditioner=f)` |
| Rust | `IC0::new(&a)?` / `f.apply(&b)?` | `a.conjugate_gradient_preconditioned(&b, options, Some(&f))?` |
| Go | `NewIC0(a)` / `f.Apply(b)` | `CGOptions.Preconditioner = f` |
| TypeScript | `new IC0(a)` / `f.apply(b)` | `a.conjugateGradient(b, {preconditioner:f})` |
| Python | `IC0(a)` / `f.apply(b)` | `a.conjugate_gradient(b, preconditioner=f)` |
| Julia | `IC0(a)` / `ic0_apply(f,b)` | `conjugate_gradient(a,b;preconditioner=f)` |
| C++ | `linear_a::IC0 f(a)` / `f.apply(b)` | Last optional argument `&f` to `a.conjugate_gradient(...)` |
| C / ARM64 | `m_ic0_create(&a,&f)` / `m_ic0_apply(&f,&b,&out)` | `m_csr_cg_preconditioned(...,&f,&result)` |
| WebAssembly | `new api.IC0(a)` / `f.apply(b)` | `a.conjugateGradient(b,{preconditioner:f})` |

Size and lower nonzero count are exposed using each port's normal naming.
Julia uses `ic0_size`, `ic0_nnz`, `ic0_lower`; C exposes `f.lower` under the
same read-only-by-contract ownership rules as its existing Cholesky factor.
Use `m_ic0_free` for C/ARM64 and `dispose()` for WebAssembly. C outputs must
start empty; failed construction does not replace an existing output. C++
factors use RAII and are noncopyable. Julia factor fields are read-only by
contract. ARM64 and WebAssembly use the shared C sparse factorization kernel.

Passing both Jacobi and IC(0), or a factor with incompatible dimensions, is an
argument error, even for a zero RHS. Overflow during factor application inside
CG returns `nonfinite`, retaining the last accepted iterate and its history.
No successful result is fabricated. Reuse with different RHS vectors is supported.
A factor of a different same-sized matrix is allowed, but the caller must ensure
it is an appropriate fixed SPD preconditioner. Positive IC(0) pivots do not prove
that the original matrix is SPD. F# exposes `ISymmetricPreconditioner` to avoid a
compilation dependency cycle; custom implementations must obey that fixed SPD,
matching-size, finite-result contract.

Apply RCM or AMD to `A`, and the same new-to-old permutation to `b`, before
constructing the factor. Restore solution coordinates afterward. Ordering may
reduce full-Cholesky fill while weakening IC(0); neither heuristic guarantees
fewer CG iterations or faster total time.

## Reproducing the comparison

```
python3 benchmarks/ic0_bench.py --sizes 8 --samples 3
```

The [public report](https://jmcummings77.github.io/linear-A/ic0/) compares all
ports with natural, RCM and AMD orderings. Setup includes factor construction
and lower export. Prepared solves reuse a factor across repeated calls with
the same RHS. Total rows measure fresh setup plus solve directly. Every recorded
row starts with an already-permuted system; ordering and permutation costs are
excluded, rather than folded into a projected total. Three serial samples and
one small workload are descriptive, not a general performance ranking.

The live report runs four methods in WebAssembly and animates accepted solution
and true-residual fields. Those are algorithm iterates, not stack frames,
physical time, or sampled flamecharts. The breakdown preset demonstrates a
failure alongside successful full Cholesky. Historical reports retain their
original measurements and embedded builds.
