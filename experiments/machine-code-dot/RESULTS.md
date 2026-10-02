# Recorded result

Run: 2026-10-02T09:04:24.090093+00:00

Float64 dot product, 4,096 elements. Median of five samples; lower is faster.

| Implementation | Median per call |
| --- | ---: |
| Rust | 3.835 µs |
| C | 3.850 µs |
| C++ | 3.853 µs |
| Go | 3.904 µs |
| Hand-encoded machine code | 3.913 µs |
| Assembly (identical instructions) | 3.939 µs |
| F# | 3.948 µs |
| C# | 3.952 µs |
| Julia | 4.111 µs |
| TypeScript / Node | 6.574 µs |
| Python | 264.339 µs |

All 286 correctness checks passed, including the checksums of every timed sample.

The machine-code function and its assembly control are byte-for-byte identical.
Their small timing difference is within the observed run-to-run variation.
Most compiled implementations cluster together at this size. Clang is faster
on tiny vectors, where its unrolling and SIMD products help more.

This is an exploratory run on a non-isolated host using the installed toolchains.
It is not a ranking of the languages in general. The complete per-sample data,
four vector lengths, and environment versions are in results/results.json.
Open results/index.html for the dashboard.
