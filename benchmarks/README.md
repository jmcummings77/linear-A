# Checking and comparing the matrix implementations

The harness builds eleven float64 implementations, checks their results against
independent references, and creates an offline performance report. Run commands
from the repository root. Python 3.9+ runs the harness; the individual toolchains
are listed in the [ports guide](../ports/README.md). ARM64 assembly requires a
native ARM64 host.

[View the live comparison report](https://jmcummings77.github.io/linear-A/latest/)
or [browse all reports](https://jmcummings77.github.io/linear-A/). The published
pages include the interactive charts, live WebAssembly controls, and geometry view.
Separate reports cover [determinants](https://jmcummings77.github.io/linear-A/determinants/),
[eigenvalues](https://jmcummings77.github.io/linear-A/eigen/), and
[vectors](https://jmcummings77.github.io/linear-A/vectors/).

The repository also includes the [saved HTML](reports/latest/index.html),
its [measurement data](reports/latest/results.json), and sanitized raw profiles.
Download the HTML file to use it offline. The **Publish reports** workflow copies
these saved files to GitHub Pages when they change on `main`, without rebuilding
the libraries or rerunning measurements. It can also be run manually from Actions.

```sh
python3 benchmarks/run.py --suite quick --profiles --require-all
```

All implementations are selected by default: `csharp`, `fsharp`, `rust`, `go`,
`typescript`, `python`, `cpp`, `c`, `assembly`, `julia`, and `wasm`. The initial
build may download compiler/package dependencies. Matrix arithmetic itself uses no
third-party numeric libraries or BLAS.

## Common commands

```sh
# Build and verify every implementation, without performance measurements.
python3 benchmarks/run.py --verify-only --require-all

# Work with only the toolchains installed on this machine.
python3 benchmarks/run.py --implementations csharp python typescript --suite quick

# Reuse artifacts from a previous successful build.
python3 benchmarks/run.py --no-build --suite quick --require-all

# Select operations and a reproducible input seed; save a separate report.
python3 benchmarks/run.py --implementations rust go cpp --operations multiply determinant --seed 42 --samples 5 --output benchmarks/reports/comparison

# Increase matrix sizes and the target duration per sample.
python3 benchmarks/run.py --suite full --require-all

# Save eigenvalue/eigenvector comparisons separately from earlier reports.
python3 benchmarks/run.py --operations eigen_symmetric --output benchmarks/reports/eigen --require-all
python3 benchmarks/run.py --operations eigen_general --output benchmarks/reports/eigen-general --require-all

# Compare cross products and rotation matrix construction in a separate report.
python3 benchmarks/run.py --operations cross rotation2d rotation3d --output benchmarks/reports/vectors --require-all
```

`--implementations` takes a space-separated list. `--require-all` makes missing
requested toolchains a failure; without it, missing toolchains are recorded as
unavailable. Build failures, correctness failures, and invalid timed results
still fail the run. Only implementations that pass correctness checks are timed.
`--no-build` skips building and dependency installation; it requires compatible
existing artifacts and the selected toolchains for discovery, version reporting,
and execution. `--verify-only`
also suppresses profiling even when `--profiles` is supplied.

## Selecting toolchains

Executables are found through `PATH`. Override an individual executable with an
environment variable; values must be executable paths, not shell command strings
with extra arguments.

| Variable | Default executable | Used for |
| --- | --- | --- |
| `DOTNET` | `dotnet` | C# and F# build/runtime |
| `CARGO`, `RUSTC` | `cargo`, `rustc` | Rust build and version reporting |
| `GO` | `go` | Go build/runtime tools |
| `NODE`, `NPM` | `node`, `npm` | TypeScript builds; TypeScript/WebAssembly runtimes |
| `EMCC` | `emcc` | Emscripten compiler for the WebAssembly target |
| `PYTHON` | `python3` | Python port runtime |
| `JULIA` | `julia` | Julia port runtime |
| `CC`, `CXX` | `clang`, `clang++` | C/assembly and C++ builds |

For example, `CC=gcc CXX=g++ python3 benchmarks/run.py --implementations c cpp`
uses GCC-compatible compilers. `PYTHON` selects the port runtime; the command
used to launch `run.py` selects the harness runtime. Keep Rust/Cargo from the same
toolchain. When `NODE` is set, the harness prepends that executable's directory
to the TypeScript build's `PATH` so npm subprocesses use the intended runtime.
The .NET SDK still follows `global.json`. Compiler/runtime ages can differ across
local installations; the report records actual versions rather than assuming
they are all current.

WebAssembly needs an activated [official Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html)
in addition to Node.js. The current SDK installer and compiler require Python
3.10 or newer; this is separate from the harness's Python 3.9 minimum. `EMCC` selects its compiler
and `NODE` selects the JavaScript runtime. To build and test the reusable port on
its own:

```sh
python3 ports/wasm/build.py --output .build/wasm
node --test ports/wasm/*.test.mjs
python3 benchmarks/run.py --implementations c wasm --suite quick --profiles
```

The build emits `.build/wasm/matrix.mjs` and `matrix.wasm`. The standalone runner
uses this module by default; `LINEAR_A_WASM_MODULE` selects another generated
module. The [WebAssembly guide](../ports/wasm/README.md) describes the Node/browser
API and memory ownership.

## What is checked and timed

The [protocol](PROTOCOL.md) specifies all input generation and output rules.
Correctness fixtures include rectangular and empty matrices, fractional values,
shape errors, singular matrices, pivot changes, triangular checks, and seeded
small determinants. Small determinant references use the independent permutation
formula with rational arithmetic. These checks supplement each port's unit tests.
There are 119 shared fixtures: 39 arithmetic, 23 symmetric eigenvalue, 20 general
eigenvalue, and 37
cross-product/rotation cases. Eigenvalue cases cover known and repeated spectra, indefinite and singular
matrices, extreme scales, empty results, and rejected unsupported inputs.
Eigenvalue checks compare the spectrum and verify `A Q = Q diag(values)` and
`Qᵀ Q = I`; they do not demand a particular eigenvector sign or basis within a
repeated eigenspace.
Vector checks cover handedness, mixed row/column orientations, finite extreme
magnitudes, invalid shapes, and zero axes. Cross references use rational
arithmetic; rotation references independently act on basis vectors with
quaternion multiplication. Quaternions are a test reference, not an exposed API.

Timed operations are addition, subtraction, scale, transpose, multiplication,
trace, determinant, symmetric eigendecomposition, cross products, and 2D/3D
rotation matrix construction. The quick suite uses matrix sizes
16 and 48, determinant sizes 8 and 16, and eigenvalue sizes 8 and 16. The full
suite uses 64, 128, and 256, determinant sizes 8, 24, and 48, and eigenvalue sizes
16, 32, and 48. In either suite, `cross` and `rotation3d` require size 3 and
`rotation2d` requires size 2. Input seed defaults to 17 and must be between 0 and
2147483646.

`eigen_symmetric` computes all eigenvalues and eigenvectors with cyclic Jacobi
rotations. Every implementation requires finite, exactly symmetric input and uses
relative Frobenius tolerance `1e-12` with a 50-sweep limit; nonconvergence fails
explicitly. Values are ascending and eigenvectors are orthonormal columns.
Use `eigen_general` for arbitrary finite real square matrices and complex results.
It uses balanced Hessenberg reduction and bounded real double-shift QR, followed
by complex back-substitution and normalization. Values are sorted by real then
imaginary part; right eigenvectors are unit columns, potentially dependent for
defective inputs. The shared analytic block-triangular workloads have known
real/complex spectra; every exact timed input is checked for spectrum,
normalization and complex residuals before measurement. C and C++ share their
numerical kernel; ARM64 and WASM use that C kernel too. Other ports implement it
in their own language. No external LAPACK/BLAS solver is used.
Accuracy is normwise, so tiny eigenvalues in mixed-scale inputs can have large
relative errors.

The eigenvalue workload is a symmetric tridiagonal Toeplitz matrix with diagonal
`d = 2 + (seed % 17)/16` and adjacent off-diagonals `-1`. Its ascending spectrum
is `d - 2 cos(k π/(n+1))`, for `k = 1..n`. Before calibration or timing, the
harness validates the complete eigenbasis for that exact input against this
analytic spectrum, eigenpair residuals, and orthogonality. The timed checksum
consumes every eigenvalue and vector entry: `sum((i+1)*values[i]) + sum(Q[i,j]^2)`.
Its independent reference is the weighted analytic spectrum plus `n`. This
checksum supplements the complete pre-timing validation; it does not replace it.

`cross` uses 3×1 vectors containing the first three usual generated entries from
`seed` and `seed+1`, with the right vector's third entry negated to avoid a
trivially zero sum checksum. `rotation2d` constructs a 2×2 matrix from an empty
input; `rotation3d` constructs a 3×3 matrix about axis `[1,2,3]`. Both rotation
benchmarks use 0.5 radians. They measure matrix construction; applying a rotation
to a vector is a separate multiplication. Before calibration, each exact timed
input receives a complete independent matrix-result check. Timed results use
the ordinary first/middle/last matrix-entry checksum.

Inputs are constructed before timing. Every batch warms up with
`max(5, min(iterations, 100))` operations. Two calibration batches choose iteration
counts and are excluded from reported samples. Each timed iteration consumes its
result in a checksum that the harness checks independently. Timing includes result
allocation, checksum consumption, deterministic disposal where applicable, and
garbage collection when the runtime schedules it. It excludes process startup,
input preparation, compilation warmup, and JSON serialization.

The harness runs implementations serially and shuffles their order with a seeded
generator. It reports the median, minimum, maximum, and median absolute deviation
across three quick or five full samples by default. `--samples` overrides this.
These measurements compare the implementations, their memory layouts, algorithms,
and toolchains. JIT behavior, GC, CPU scheduling, thermal state, and compiler age
also affect results. The quick suite is a small smoke test, not a universal
language ranking; shared CI machines are especially unsuitable for treating
small timing differences as stable conclusions.

The C# baseline uses the existing `Matrix<double>` API, with a copy before its
mutating scale operation. Its determinant now calls the generic type-aware Auto
algorithm directly; for double it uses scaled partial-pivot LU. Other ports use
Auto, with safe tiny or triangular shortcuts and partial-pivot LU otherwise.
Older saved reports retain their recorded implementation notes and measurements.
No fast-math flags are used. ARM64 results include the assembly kernels and their
C support code; they are not measurements of a standalone assembly runtime.

WebAssembly compiles the same float64 kernels as the C port. Its measured loop
runs in Node.js and includes JavaScript-to-WebAssembly calls, matrix allocation,
checksum readback, and explicit disposal. Input copying, module loading, and
initial compilation occur before timing. These results compare native C with the
WebAssembly runtime and its bridge costs; they are not browser benchmarks.

## Reports and profiles

The default output directory is `benchmarks/reports/latest`:

- `index.html` is a self-contained report with embedded data and no network assets.
  Open it directly in a browser, including offline. When a verified WebAssembly
  build is available, it also contains the executable live-testing bundle.
- `results.json` stores checks, raw timing samples, summary statistics, machine
  information, compiler/runtime versions, build commands, revision, and source hash.
- `profiles/` contains raw profiler output when `--profiles` is enabled.

Reports are prepared for public sharing by default. Exported metadata uses a
fixed set of fields: tool commands use executable names, project sources use
repository-relative paths, and external sources use portable labels. Home and
temporary paths, environment variables, credentials, process/thread identifiers,
machine names, and authorship annotations are omitted or redacted. Function
names, source lines, measurements, compiler versions, basic OS/CPU information,
timestamps, and Git provenance remain available. Raw profiles and the embedded
JavaScript/WebAssembly bundle receive the same treatment before publication.

Build artifacts under `.build/` and generated reports are ignored by Git. Keep
the JSON with any exported chart so its environment and measurements remain
reviewable. Use `--output` to preserve multiple runs.

Profiles come from separate, instrumented or sampled runs of 48×48 matrix
multiplication and do not contribute to comparative benchmark timings. Each
profile records its workload, iteration count, and input seed. Available backends include Python stack
sampling, V8's CPU profiler for JavaScript/WebAssembly, Go's CPU profiler, Julia's
`Profile`, and macOS native stack sampling. Native profiling is currently
unavailable on other operating systems; unavailable or failed profiles are shown
explicitly. `--require-all`
requires implementations, not a profiler for every implementation.

Profile weights represent sampled CPU time or sampling intervals, depending on the
backend. They are not exact wall-clock durations of individual function calls.
Some profiles include startup/warmup, and aggregate profiles do not preserve
chronology; each profile records these limits. A sample timeline is shown only
when the backend provides chronological samples. Sampling and symbol quality can
also differ between runtimes, so flame graphs explain where a run spent time
rather than supplying an independent language ranking.

Node/WebAssembly timelines follow V8's chronological sample order and include
JavaScript glue, WebAssembly execution, and Node startup/warmup in the separate
profiling run. WebAssembly frames may have incomplete names or appear as module
function identifiers. These samples do not describe browser execution or give
exact durations for each JavaScript-to-WebAssembly call.

## Live WebAssembly in the report

The live controls execute the C-backed WebAssembly API in your current browser.
Choose the math checks to run all 119 shared fixtures, or choose addition,
subtraction, scale, transpose, multiplication, trace, determinant, symmetric or general real
eigenvalues/eigenvectors, cross products, or rotation construction to measure
locally. Matrix sizes are limited to 16, 48,
64, 128, or 256; determinant uses 4, 8, 16, 24, or 48, and eigendecomposition uses
4, 8, 16, 32, or 48. Cross/3D-rotation size is fixed at 3 and 2D-rotation size at
2. Select three or five samples. A benchmark first passes all 119 fixtures, then
validates each batch's checksum. Eigenvalue workloads also validate the exact
timed input's complete eigenbasis before measurement; vector/rotation workloads
check the complete result matrix independently before starting the clock.

Each run uses a fresh Web Worker so the page stays responsive. Stop terminates
the worker, and a 45-second timeout bounds each run. The live section needs a
modern browser with WebAssembly and module-worker support. Browser restrictions
or initialization failures are reported in the page.

Live timings use the browser's `performance.now()` clock. Browser scheduling,
clock resolution, JIT compilation, and machine load affect them. They remain
separate from the historical comparison charts and Node measurements. The live
JSON export contains the current browser runs; it does not replace `results.json`
or recorded profiles. It includes only a recognized browser name/version, without
the full user-agent string or arbitrary runtime diagnostics. Live runs do not
collect a stack profile.

Report generation includes live controls automatically when the recorded results
contain a passed `wasm` implementation and `.build/wasm` contains a valid
`matrix.wasm` and its generated `matrix.mjs`. The HTML embeds the binary,
JavaScript glue, API wrapper, worker, and fixture data, so the generated file
needs no server or network connection. When the verified result or build is
missing, the report explains why live execution is unavailable.

Regenerate a public report from existing measurements and an existing WebAssembly
build, without running benchmarks:

```sh
python3 benchmarks/report.py benchmarks/reports/latest/results.json benchmarks/reports/latest/index.html --wasm-dir .build/wasm

# Keep the report as a viewer for recorded results only.
python3 benchmarks/report.py benchmarks/reports/latest/results.json benchmarks/reports/latest/index.html --without-live
```

`--wasm-dir` selects an alternate build directory; it defaults to `.build/wasm`.
Regeneration writes sanitized `results.json` and raw profiles alongside the HTML.
With the same output directory it updates those artifacts in place; with a new
directory it leaves the input artifacts untouched. Historical measurements and
their provenance are retained exactly. The embedded code can be newer than the recorded measurements, so the
live bundle has its own SHA-256 identifier displayed and exported separately
from the recorded source hash.

## Timing chart language filter

Use the **Languages** dropdown in **Operation timings** to check which languages
are shown. Unchecking Python, for example, removes its bars and rescales the chart
to the remaining medians. The sample table and comparison baseline use the same
selection; choices remain in effect when changing operation or matrix size.
**Select all** and **Clear all** let you reset or build a smaller comparison.
In determinant reports, each language checkbox includes all its algorithm variants.
Downloads retain the complete recorded run.

## Report appearance

The **Appearance** selector offers System, Light, and Dark themes. System is the
initial default and follows changes to the device's color preference. An explicit
choice is saved locally when browser storage is available; the controls also work
when storage is blocked. Theme changes preserve live results, animation progress,
and flamechart zoom. The theme script and palettes are embedded for offline use.

The report's emblem and decorative graphics adapt six patterns from
[Book of Shapes](https://bookofshapes.com/). They are embedded as static SVG, with
colors matched to each theme. See [artwork sources and reuse terms](ARTWORK.md).
These illustrations do not represent benchmark measurements or replace the live
WASM geometry view.

## Animated geometry and calculation steps

The report includes an interactive 3D view of 125 vectors evaluated by the live
WebAssembly library. Edit two 3×3 matrices or choose rotation, shear, reflection,
collapse, axis-rotation, and cross-direction presets. Multiplication, addition,
subtraction, transpose, scaling, rotation, and crossing with a direction produce
a new operator; the input and output vector fields are checked
against independent references before they are displayed.

**Randomize matrices** fills both M and N with values from −2 to 2 in quarter
steps and computes a new field, keeping the selected operation and its settings.
Use **Show overlays** to toggle the axes, bounding box, and eigenvector ellipse
independently. Visibility changes apply immediately during playback and keep the
current calculation progress. The ellipse appears for complex eigenvalue pairs.

Play a 10, 20, or 30 second transformation, pause, or scrub to any point. Drag the
plot to orbit and scroll to zoom; when the plot has keyboard focus, arrow keys
rotate it and `+`/`-` zoom. The animation starts only when requested. The wire
unit cube and RGB basis vectors show the same matrix acting on geometry; arrows
show `A*x` anchored at `x`, with a fixed display scale. Column vectors mean that
`(M*N)*x` applies `N` first. Ordinary operator interpolation is `(1-t)*M + t*R`,
so even invertible endpoints can pass through a singular matrix. The displayed determinant is the
computed signed volume scale, not an assertion of exact numerical rank.

The axis-rotation preset uses an editable direction `u` and an angle displayed
in degrees, converted to radians for the library. Playback follows
`R(u, tθ) M`: the worker computes and independently verifies 901 WASM matrices
along this angular path. These are sampled rotation states, not captured CPU
function frames. The vector-field and cube move through the actual rotation
path rather than a linear blend of the endpoint matrices.

The cross-direction preset crosses each column of `M` with a chosen direction
`u`, forming a linear operator whose field is `(M x) × u`. This applies the
library's vector cross product; it does not define a cross product of arbitrary
matrices `M` and `N`. A zero cross direction is valid and collapses the field;
a rotation axis must be nonzero.

For every real input `M`, the view also computes its eigenpairs in
WebAssembly and checks the residual and orthogonality independently. Choose
**Symmetric stretch · eigenvectors** to see eigenvalues `0.5`, `1`, and `3`, with
purple arrows along the corresponding unit eigenvectors. For nonsymmetric
inputs, real eigenvectors remain purple arrows; each complex conjugate pair
gets a gold ellipse in the plane spanned by its real and imaginary components.
The inspector lists the full complex components. Defective inputs may repeat
a direction; an independent basis is not promised. Expand **Eigenvalues &
eigenvectors of M** for their components. These arrows describe the fixed input
`M` throughout playback; they are not eigenvectors of each interpolated matrix
or captured Jacobi iteration frames.

A separate instrumented build enables **Calculation steps** for multiplication:

```sh
python3 ports/wasm/build.py --trace --output .build/wasm-trace
python3 benchmarks/report.py benchmarks/reports/latest/results.json benchmarks/reports/latest/index.html --trace-wasm-dir .build/wasm-trace
```

The trace records operands, products, running sums, indices, and source lines
inside the C multiplication kernel. The report reconstructs the partial matrix
from those captured accumulators, highlights the active inputs and C source
line, and displays the logical call path. There are 27 updates for a 3×3 product.
The trace build is capped at 8×8 inputs and 512 events; its source snapshot is
saved with the build and embedded with the report. The ordinary WASM build has
no trace hook or trace API. Only multiplication exposes these captured C function
steps; rotations and cross products have no captured instruction trace. Recorded
performance measurements and sampled profiles are unchanged.

**Locate multiplication in profile** selects the matching function in the saved
WASM flamegraph. That profile describes a separate 48×48 run: its CPU samples
cannot be aligned to the current arithmetic steps. Step playback has no timing
claims. If the optional trace build is absent, geometric animation still works.
All code and data are embedded, with the same publication sanitization as the
other report assets; no external rendering library or network asset is needed.

## Continuous integration

[The ports workflow](../.github/workflows/ports.yml) runs unit tests and all-runner
verification on pull requests using the documented
[`ubuntu-24.04-arm` GitHub-hosted runner](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).
This executes the ARM64 kernels natively on Linux. The existing C# workflow
continues to test its supported operating systems.

A manual workflow run can enable the `performance` input to collect a quick
benchmark report and available profiles. Reports are uploaded as workflow
artifacts; there is no timing threshold that fails a pull request. CI installs
the SDK recorded in `global.json`, current Rust/Go/Julia stable toolchains, a
current Python 3, Node LTS, the pinned TypeScript dependency, and Emscripten through
the [official setup action](https://github.com/emscripten-core/setup-emsdk).
The workflow builds the WebAssembly module, runs its JavaScript API and live-worker
tests, and includes `wasm` in verification of all eleven runners.

## Determinant algorithm comparisons

```sh
python3 benchmarks/determinants.py --suite quick --require-all
python3 benchmarks/determinants.py --suite full --samples 5 --require-all
```

This separate suite defaults to `benchmarks/reports/determinants/index.html` and
uses the same build, verification, calibration, timing and public sanitization
pipeline. It supports `--implementations`, `--no-build`, `--verify-only`, `--seed`,
and `--output`. It does not overwrite the ordinary saved comparison.

Choose a workload in the report:

- `determinant`: Auto versus explicit LU on identical general matrices.
- `determinant_spd`: LU versus Cholesky on identical symmetric positive-definite matrices.
- `determinant_small`: C# Auto, LU and the former cofactor algorithm at sizes 3, 5 and 7.

The positive-definite input is `(A + transpose(A))/2 + 4n I`, using the usual
unshifted generated A. Its entries are exact multiples of 1/32. Symmetry, positive
diagonals and strict diagonal dominance guarantee positive definiteness. The
reference computes an exact integer Bareiss determinant before converting to a
float; small correctness fixtures also use an independent permutation formula.
Compare speed ratios only within one workload and size. Cholesky's requirement
is checked explicitly, and it is never selected automatically for arbitrary input.

The live report offers `determinant · Auto`, `determinant · LU`, `SPD determinant · LU`,
and `SPD determinant · Cholesky`, using the same inputs and independent integer
oracle in the browser. These live samples remain separate from saved measurements.
