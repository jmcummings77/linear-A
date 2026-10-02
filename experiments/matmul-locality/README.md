# Multiplication locality study

A controlled native experiment: take the same row-major float64 product through
scalar i-j-k loops, reordered i-k-j loops, 32×32 blocking, and explicit two-lane
ARM64 NEON. An optional Apple Accelerate BLAS implementation supplies context.
This isolates algorithm changes in one language and runtime; it is not another
language ranking and does not change any of the eleven library implementations.

[Published dashboard](https://jmcummings77.github.io/linear-A/matmul-locality/) ·
[Saved HTML](results/index.html) · [Measurements](results/results.json) ·
[Recorded findings](RESULTS.md)

## Hypotheses and controls

| Variant | Change | Hypothesis and possible cost |
| --- | --- | --- |
| Scalar i-j-k | One dot product per output entry | Minimal output writes, but B is read with a row-length stride |
| Scalar i-k-j | Stream contiguous rows of B and C | Better spatial locality, but repeatedly update C |
| 32×32 tiles | Limit active A/B/C regions | More reuse may help larger problems; loop overhead and shorter spans may lose |
| Tiles + NEON | Two adjacent output columns together | Fewer arithmetic/load instructions per output; tails and memory traffic remain |
| Accelerate BLAS | Vendor library | Context only; packing, threading, instructions, reductions and allocations are opaque |

All handwritten variants preserve ascending reduction-index order per output.
They use separate multiplication and addition, with no fast-math or contraction.
Scalar auto-vectorization is explicitly **disabled** in this controlled study:
SIMD is a distinct step, not an accidental compiler difference. This also means
the scalar results do **not** represent the best code a default optimizing
compiler could produce. Functions are compiled in a separate translation unit
from the runner, without LTO. Distinct buffers satisfy the kernels' `restrict`
contract. Output initialization is included in each kernel. NEON is available
only on ARM64; the three scalar variants also build on other 64-bit Clang targets.

The tile size is fixed, not selected from the measurements. Negative results are
retained. Source is in [src/kernels.c](src/kernels.c), with the bounded test and
timing interface in [src/runner.c](src/runner.c).

## Run and reproduce

Requires Python 3.9+ and Clang. On this ARM64 Mac:

```sh
# Correctness and native memory/undefined-behavior checks; no timing output.
python3 experiments/matmul-locality/run.py --verify-only --sanitize --blas

# Fresh measurements and separate native profiles, including optional BLAS.
python3 experiments/matmul-locality/run.py --blas --profiles

# Save a repeat separately, preserving the committed measurement snapshot.
python3 experiments/matmul-locality/run.py --blas --profiles \
  --output /tmp/linear-a-locality-repeat

# Render the saved measurements again without compiling, profiling or timing.
python3 experiments/matmul-locality/run.py \
  --render experiments/matmul-locality/results/results.json
```

Omit `--blas` outside macOS. Native sampling currently uses macOS `sample` and
attaches only to the process launched for that profile. If it is unavailable,
the report records that fact; it never substitutes invented stacks. `--samples`
selects 3–30 independent process samples per workload (default 5). `--shapes`
accepts M×K×N, such as `--shapes 32x32x32 96x257x65`; dimensions are capped at 1024.
ASan/UBSan builds are correctness-only and cannot generate performance results.
CI exercises the handwritten kernels with sanitizers, without a performance gate.

## Correctness

The Python oracle uses exact `Fraction` arithmetic independently of the kernels.
Twelve cases per variant cover rectangular and empty shapes, one element,
odd SIMD tails, both sides of tile edges, decimal cancellation, and mixed scales.
For non-exact arithmetic, every entry must satisfy a forward-error tolerance of
`32 * max(1,K) * epsilon * sum(abs(a*b)) + 1e-300`. This deliberately conservative
bound scales with cancellation rather than dividing by a nearly zero answer.
Output storage is poisoned with NaNs before the check, detecting missing writes.

Timed inputs are small integers divided by 16. Their products and sums are exact
in binary64 within the bounded shapes, even with a different reduction order.
The native runner independently checks every output entry using an integer oracle
before and after every batch. The Python harness separately checks a three-entry
checksum consumed on every timed call. Both allocation modes are correctness-tested.
Sanitizer runs exercise the same bounds and tails; their timings are discarded.

## Timing and allocation boundaries

Each sample launches a fresh process, checks its inputs/results, warms up three
calls, then times a calibrated batch near 50 ms with a monotonic clock. Two
calibration batches are excluded. Zero-duration probes increase iterations with
a bounded retry. The recorded iteration count is the count actually measured.
Variants run serially in seeded shuffled order each round. The dashboard retains
all samples and order positions, and reports median, MAD, minimum and maximum.
Five fresh processes do not provide five independent machines or a confidence
interval; the same CPU state can affect successive samples.

Two output lifecycles expose different costs:

- **Reused:** all three buffers are allocated before timing. Each call still
  initializes C, performs the product and consumes the checksum.
- **Allocated:** input buffers are reused; a fresh output buffer is allocated and
  freed inside each call. The harness counts these calls and requested bytes.

Allocation accounting describes harness output allocations, not peak resident
memory, allocator bookkeeping, or BLAS's internal allocations. All handwritten
kernels make zero internal allocations. The existing library's validation and
ownership wrappers are intentionally excluded. Accelerate is requested to use
one thread through a child-process setting; actual internal behavior is not
measured. BLAS may use fused arithmetic and other instructions unlike the
handwritten kernels. These limitations are visible in the report.

Input buffers remain resident throughout a batch. There is no cache flush or
streaming data set. No hardware cache-miss or bandwidth counters are collected;
locality explanations are hypotheses supported by access patterns, not direct
measurements of cache behavior. The host is not isolated from other activity.

## Profiles and publication

Profiles come from separate roughly six-second workloads at 256×256×256, using
reused outputs. macOS sampling collects a two-second window of main-thread stacks.
The flame graph uses aggregate counts times the nominal interval, not timestamps
or exact elapsed function times. Startup, warmup, checks and allocation may appear;
BLAS worker threads are excluded. Click a frame to zoom. These are flame graphs,
not chronological flamecharts. Emitted disassembly is also saved when `xcrun`
is available so compiler claims can be inspected.

Exports reuse the repository's path and raw-profile sanitizers. They retain
function names, compiler version, source context, basic hardware, timestamps,
source hash and Git provenance, while omitting local paths, environment dumps,
credentials and machine identifiers. The source hash covers the C sources and
Python harness. A report may truthfully record an uncommitted source state; the
hash still identifies the measured source contents. Rendering alone preserves
every saved measurement. Existing comparison reports and raw profiles remain
separate and untouched.

References: [Arm NEON intrinsics](https://developer.arm.com/architectures/instruction-sets/intrinsics/)
and [Apple's double-precision GEMM interface](https://developer.apple.com/documentation/accelerate/cblas_dgemm(_:_:_:_:_:_:_:_:_:_:_:_:_:_:)).
