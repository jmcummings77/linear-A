# API compatibility: 0.1.0

This is the shared float64 contract for these release kits. A check mark is
a declared capability validated by the shared conformance suite in CI.
Package installation tests separately exercise solving and multiplication.

| Port | Required runtime/toolchain | arithmetic | determinants | eigenpairs | geometry | linear-solves | least-squares |
| --- | --- | --- | --- | --- | --- | --- | --- |
| assembly | ARM64 + C11 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| c | C11 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| cpp | C++17 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| csharp | .NET 10; legacy netstandard2.1 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| fsharp | .NET 10 | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| go | Go 1.22+ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| julia | Julia 1.10+ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| python | Python 3.9+ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| rust | Rust 1.69+ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| typescript | Node 22+ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| wasm | Emscripten; Node 22+ or browser | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

C# entries apply to `Matrix<double>` on .NET 10. The netstandard2.1
target contains the older compatibility classes and does not promise this
complete contract. ARM64 and WebAssembly use the C API/kernels; C++ shares
the C general-eigenvalue and factorization kernels.
