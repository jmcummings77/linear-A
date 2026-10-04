# Getting started

[Project overview](../README.md) · [Design and numerical scope](design.md) ·
[Reports and exploration](reports.md)

Choose the implementation you want to use. Each has a reusable matrix API, its
own tests, and a runner for the shared comparison protocol. You only need the
toolchain for the implementations you select.

## Get the source

```sh
git clone https://github.com/jmcummings77/linear-A.git
cd linear-A
```

Run the commands below from this repository root. For packaged artifacts and
isolated consumer examples, use [the release guide](../release/README.md).

## Choose an implementation

Each language guide explains its API, ownership and error rules, build steps,
and standalone tests.

| Implementation | Library location | Guide and toolchain |
| --- | --- | --- |
| ARM64 assembly | [ports/assembly](../ports/assembly) | [Assembly guide](../ports/assembly/README.md); ARM64 host and C11 compiler/assembler |
| C | [ports/c](../ports/c) | [C guide](../ports/c/README.md); C11 compiler |
| C# | [linear-A/linear-A](../linear-A/linear-A) | [C# guide](../linear-A/README.md); .NET SDK selected by [global.json](../global.json) |
| C++ | [ports/cpp](../ports/cpp) | [C++ guide](../ports/cpp/README.md); C++17 compiler |
| F# | [ports/fsharp](../ports/fsharp) | [F# guide](../ports/fsharp/README.md); .NET 10 SDK |
| Go | [ports/go](../ports/go) | [Go guide](../ports/go/README.md); Go 1.22+ |
| Julia | [ports/julia](../ports/julia) | [Julia guide](../ports/julia/README.md); Julia 1.10+ |
| Python | [ports/python](../ports/python) | [Python guide](../ports/python/README.md); Python 3.9+ |
| Rust | [ports/rust](../ports/rust) | [Rust guide](../ports/rust/README.md); Rust/Cargo 1.69+ |
| TypeScript | [ports/typescript](../ports/typescript) | [TypeScript guide](../ports/typescript/README.md); Node.js 22+ and npm |
| WebAssembly | [ports/wasm](../ports/wasm) | [WebAssembly guide](../ports/wasm/README.md); Emscripten for builds, Node.js or a browser to run |

The [shared API guide](../ports/README.md) explains the common float64 contract
and the differences between these interfaces. The [design guide](design.md)
describes which implementations share numerical kernels.

## Try the Python implementation

The [one-minute example](../README.md#one-minute-introduction) uses only Python's
standard library. Its shell command sets an import path for that invocation.
In another shell, or in an interactive Python session, run this from the
repository root instead:

```python
import sys
sys.path.insert(0, "ports/python")
from matrix import Matrix

a = Matrix(2, 3, [1, 2, 3, 4, 5, 6])
gram = a.multiply(a.transpose())
print(gram.row(0), gram.row(1))  # [14.0, 32.0] [32.0, 77.0]
print(gram.determinant())       # 54.0
```

Run its standalone tests:

```sh
python3 -m unittest discover -s ports/python -p 'test_*.py'
```

## Verify or compare implementations

The common harness checks selected implementations against independent
references before measuring them. A verification run needs Python 3.9+ and the
selected implementation's toolchain; this example needs only Python:

```sh
python3 benchmarks/run.py --verify-only --require-all --implementations python --output .build/verification-python
```

Add runner IDs such as `c`, `go`, or `rust` after `--implementations` to include
other installed toolchains. `--require-all` makes a missing selected toolchain
an error. Results go into the chosen output directory; use a different directory
for a new run to preserve earlier evidence.

For timing, profiling, executable overrides, specialized solver suites, and
publication captures, follow the [benchmark guide](../benchmarks/README.md).
Public captures require a clean committed source tree and the capture wrapper;
ordinary development runs can measure working changes. The
[runner protocol](../benchmarks/PROTOCOL.md) defines portable inputs, outputs,
and timing boundaries.

Before contributing, follow the [development checks](../CONTRIBUTING.md#checks-before-a-pull-request)
and the chosen language's formatting and testing instructions. Generated checks
and native fuzzing have their own [guide](../tests/bugfinding/README.md).
