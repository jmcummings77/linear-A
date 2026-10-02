# ARM64 float64 matrix kernels

This is an ARM64 assembly/C implementation with an explicit division of work:

- `kernels.S` performs addition, subtraction, scalar multiplication, transpose,
  matrix multiplication, vector cross products, trace, and the row-update arithmetic in determinant
  elimination. It uses scalar double-precision instructions and separate
  multiplication/addition to match the other implementations' rounding order.
- The [C matrix API](../c/matrix.h) supplies allocation, ownership, copying,
  indexed access, row/column extraction, shape/finite-value validation, and
  triangular classification. C also controls Gaussian pivot selection, row
  swaps, elimination factors, and the determinant's signed diagonal product.
- The [shared C runner](../c/runner.c) supplies input parsing, JSON output,
  benchmark generation, monotonic timing, and checksums.

This is not an all-assembly program. Benchmark timings include the C support
work as well as the assembly kernels. `MATRIX_USE_ASM` dispatches the listed
matrix arithmetic and ordinary LU row updates to assembly; unsupported
architectures fail to compile. Determinant Auto's tiny/triangular shortcuts and
explicit Cholesky execute in C. Symmetric eigendecomposition also executes in C:
its cyclic Jacobi rotations and eigenvector updates use the shared C solver;
there is no separate assembly eigen kernel. General eigendecomposition also
uses the C balancing/Hessenberg/double-shift QR solver, including complex
eigenpairs, normalization, and sorting; its timings include that C work. Rotation factories, trigonometric
functions, axis normalization, and Rodrigues' formula execute in C. Cross-product
components use the dedicated `m_asm_cross` kernel; C validates vector shapes,
checks finite results, and owns allocation. A rare LU update with a subnormal elimination
factor also uses C's exponent-scaled product to retain a representable result.
These boundaries are included in the reported timings.

Apple Silicon macOS and ARM64 Linux symbols are supported through the `.S`
preprocessor. Build on an ARM64 host from the repository root:

```sh
clang -std=c11 -O3 -ffp-contract=off -DMATRIX_USE_ASM ports/c/matrix.c ports/c/runner.c ports/assembly/kernels.S -lm -o /tmp/linear-a-assembly
clang -std=c11 -O3 -ffp-contract=off -DMATRIX_USE_ASM ports/c/matrix.c ports/assembly/tests.c ports/assembly/kernels.S -lm -o /tmp/linear-a-assembly-tests
/tmp/linear-a-assembly-tests
/tmp/linear-a-assembly bench multiply 32 100 42
```

The reusable API and its memory ownership are identical to the C implementation.
Empty dimensions are supported. Float64 rounding, cancellation, underflow, and
intermediate overflow limitations also match the C implementation. Tests exercise
the same API with the assembly kernels linked. See the
[shared runner protocol](../../benchmarks/PROTOCOL.md) for check/bench requests.
