# Recorded findings

Snapshot: October 2, 2026, macOS ARM64, eight logical CPUs. These are observations
from one unisolated machine, not performance guarantees. The [dashboard](results/index.html)
contains all 250 process samples, ranges, MAD, allocation modes, compiler settings,
source hash, Git provenance, and five separately collected native profiles.

## Reused-output medians

All entries are **microseconds per multiplication**; lower is faster.
Scalar auto-vectorization and floating-point contraction were disabled.

| M × K × N | i-j-k | i-k-j | Tiles | Tiles + NEON | Accelerate |
| --- | ---: | ---: | ---: | ---: | ---: |
| 32 × 32 × 32 | 13.626 | 12.494 | 12.839 | 7.809 | 0.736 |
| 128 × 128 × 128 | 1,641.409 | 907.829 | 773.600 | 609.186 | 16.226 |
| 256 × 256 × 256 | 14,891.333 | 7,416.500 | 6,236.143 | 5,506.889 | 115.103 |
| 512 × 512 × 512 | 190,052.000 | 56,921.000 | 54,714.000 | 65,764.000 | 1,228.368 |
| 96 × 257 × 65 | 1,214.026 | 704.569 | 678.366 | 546.157 | 16.785 |

## What changed

- **Loop reordering:** at 256³, i-k-j took 7.417 ms versus 14.891 ms for
  i-j-k, a 2.01× median speedup. The contiguous inner-loop access pattern is a
  plausible explanation; this experiment did not measure cache misses.
- **Blocking:** 32×32 tiles reduced the 256³ median to 6.236 ms, 1.19× as fast
  as i-k-j. At 32³, the tiled median was slightly worse (12.839 versus
  12.494 µs), with overlapping ranges. There is no convincing small-size win.
- **Explicit SIMD:** NEON reduced the 256³ median to 5.507 ms, 1.13× as fast as
  scalar tiles and 2.70× as fast as the initial loop. The saved disassembly
  contains `fmul.2d` and `fadd.2d`, with no fused multiply-add in the handwritten
  kernels. The scalar loops contain scalar floating-point arithmetic.
- **A negative result:** at 512³, the NEON median was 65.764 ms versus 54.714 ms
  for scalar tiles—about 20% slower. It also lost in allocate/free mode.
  Simply widening this inner loop is insufficient to establish a better large
  matrix kernel. Packing, register blocking and tile selection were not tested;
  neither cache behavior nor scheduling effects can be inferred from this run.
- **Vendor context:** Accelerate was fastest at every measured shape. At 256³
  it took 115.103 µs, about 47.8× faster than this simple NEON kernel. This gap
  shows how much these deliberately constrained examples leave unexplored;
  it does not isolate any one vendor optimization. Internal allocations and
  actual worker-thread activity were not measured.

## Allocation and uncertainty

Allocate/free mode issues one output allocation per multiplication, versus zero
inside timing for reused output. At 512³ that requests 2,097,152 bytes per call.
Those are requested bytes, not memory footprint or an estimate of allocator cost.
Some allocate/free medians are lower than reused medians: the modes were measured
in separate batches, so subtracting them would not yield a reliable allocation
latency. Frequency, scheduling, thermal state and allocator reuse can differ.

Noise matters here: the reused 128³ i-j-k MAD was 18.8% of its median, and the
512³ Accelerate MAD was 15.8%. Five processes are not enough to turn small
ratios into strong conclusions. No performance threshold or default-library
change is based on these observations.

## Evidence and correctness

All 60 shape/numerical cases passed the independent rational reference checks
in both optimized and AddressSanitizer/UndefinedBehaviorSanitizer builds. Every
timed batch passed full output verification against the independent integer
oracle before and after timing, plus the consumed per-call checksum. Profiles
are real main-thread aggregate samples, not chronological traces or exact
function durations. The handwritten profiles are dominated by their respective
multiplication functions; they do not explain cache misses or the 512³ reversal.

The saved snapshot records an uncommitted source state honestly; its source hash
identifies the measured C files and Python harness. Historical cross-language
measurements and profiles were not modified. See [methodology](README.md) for
commands, timing boundaries and publication sanitization.
