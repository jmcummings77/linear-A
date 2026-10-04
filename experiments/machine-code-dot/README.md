# 36-byte dot-product experiment

One deliberately small function written directly as ARM64 instruction words,
compared with scalar implementations in C#, F#, Rust, Go, TypeScript, Python,
Julia, C++, C, and assembly. This experiment is independent of the matrix ports
and their shared harness.

[View the recorded report](https://jmcummings77.github.io/linear-A/machine-code-dot/).

The function computes `sum(a[i] * b[i])` for two arrays of float64 values. It
accepts a count, handles zero length, and returns a double. It does not allocate,
validate pointers, or check the arrays' lengths: the caller must provide at least
`n` readable doubles in each array. The native runner validates its inputs and
allocates those arrays before calling it.

## What is actually machine code here?

[src/raw_dot.h](src/raw_dot.h) contains nine manually selected 32-bit instruction
encodings, totaling **36 bytes**. The kernel was not produced by compiling C,
assembly, or another source language. The C runner copies these literal bytes
into a macOS JIT memory page, flushes the instruction cache, enables execution,
and calls the resulting function pointer. C is used for the test harness and
memory management, not for generating the kernel instructions.

An independently assembled version of the same nine instructions is included as
a control. The harness disassembles it and requires every instruction word to
match the manual encoding exactly. It also exports the raw 36-byte binary.

| Offset | Instruction word | Meaning |
| --- | --- | --- |
| 0x00 | `9e6703e0` | Set the double accumulator to zero |
| 0x04 | `b40000e2` | Return immediately if length is zero |
| 0x08 | `fc408401` | Load a double from the first vector; advance pointer |
| 0x0c | `fc408422` | Load a double from the second vector; advance pointer |
| 0x10 | `1e620821` | Multiply the loaded values |
| 0x14 | `1e612800` | Add the product to the accumulator |
| 0x18 | `f1000442` | Decrement the remaining count |
| 0x1c | `54ffff61` | Repeat while the count is nonzero |
| 0x20 | `d65f03c0` | Return the accumulator |

Words are stored in little-endian byte order on the target Mac. This version is
intentionally scalar: it does not use SIMD, multiple accumulators, or fused
multiply-add. Such optimizations would introduce another experimental variable.

## Run

Requires **ARM64 macOS**, Clang/Xcode command-line tools, Rust, Go, Node.js with
TypeScript type stripping (the measured run uses Node 25), Python 3.9+, .NET 10,
and Julia. No third-party numerical libraries or Python packages are used.

```sh
python3 experiments/machine-code-dot/run.py \
  --dotnet /path/to/dotnet \
  --julia /path/to/julia
```

Omit the executable overrides when the matching tools are on `PATH`.

Options include `--sizes 16 256 4096 65536`, `--samples 5`, `--no-build`, and
`--output /tmp/linear-a-dot-repeat` for a separate local run.
The harness uses a private build directory and private runtime caches under this
experiment. The F# dependency is restored from the .NET SDK's bundled package
when available. It does not modify the main matrix projects or their harness.

To replace the published measurements, first commit the source and start from a
clean checkout, then use the shared capture wrapper. It builds in a fresh clone
and records the recoverable commit, source tree, command and artifact hashes.
Place the documented toolchains on `PATH` before running:

```sh
python3 benchmarks/reproduce.py --runner machine-code-dot \
  --output .build/publishable/machine-code-dot -- --samples 5
```

Review the capture, copy its complete contents into `experiments/machine-code-dot/results/`,
and run `python3 benchmarks/check_provenance.py` before committing the report.
See the [publication workflow](../../benchmarks/README.md#publishing-reproducible-source-records).
The existing snapshot predates Git provenance recording. Its exact measured
source cannot be reconstructed, so its frozen measurements and disassembly are
historical evidence; they are not a recoverable source snapshot.

Output in `results/`:

- `index.html`: interactive comparison dashboard, with a vector-length selector.
- `results.json`: metadata, iteration counts, and every timing sample.
- `raw_dot.bin`: the 36-byte function.
- `disassembly.txt`: verification of the assembly control.

## Method and limits

All implementations use the same deterministic binary-exact inputs and a
left-to-right scalar sum. An independent integer oracle verifies checksums
exactly, including zero-length, one-element, odd-sized, and longer inputs.
Every timed result is also verified. All four native harness variants use the
same indirect function-call mechanism.

The benchmark rotates among 32 vector pairs so repeatedly calling the function
is not simply an invariant expression. Allocation, input generation, startup,
ahead-of-time compilation, and JSON output occur outside the timed region.
JIT optimization may still occur during a short trial, especially for Node. Each process warms
up for 8–128 calls, then times a calibrated batch. Five independent samples run
in shuffled serial order by default. .NET tiered compilation is disabled, and
its kernels request optimized code without inlining. Python uses `array('d')`;
TypeScript uses `Float64Array`; other languages use their normal double arrays
or slices. No runtime consumes another language's numeric implementation.

These are scalar-loop experiments, not comparisons with NumPy, BLAS, Julia's
standard `dot`, or highly tuned vectorized libraries. Normal language/runtime
bounds-checking differences remain; Julia explicitly uses `@inbounds`, while
native pointer kernels rely on the caller. Compilers may auto-vectorize products and unroll loops within strict
floating-point semantics. In the recorded run, Clang generated eight-element
blocks with SIMD multiplication and ordered scalar additions. The manual kernel
uses one scalar multiply per iteration; this explains part of its disadvantage
on short vectors. See `results/c-disassembly.txt`. Native C/C++ builds disable
floating-point contraction.

Median and min–max timings describe this machine and run. CPU frequency, caches,
JIT behavior, and other workloads can affect the numbers. The host is not
isolated. Identical machine instructions in the raw and assembly cases are
expected to perform similarly; minor differences are not evidence that manually
writing bytes is intrinsically faster.

The memory-loading method follows Apple's
[Apple silicon JIT guidance](https://developer.apple.com/documentation/apple-silicon/porting-just-in-time-compilers-to-apple-silicon).
