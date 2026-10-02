# WebAssembly float64 matrices

The WebAssembly build compiles the [C matrix implementation](../c/matrix.c)
with Emscripten. The reusable ES module works in Node.js and modern browsers.
Arithmetic runs in WebAssembly; JavaScript handles initialization, input/output
copies, and object ownership. It uses float64 throughout, sharing C's determinant
algorithms, finite-value checks, and floating-point limitations.
This is a compiled C backend, not a separate matrix algorithm.

Build from the repository root after activating Emscripten:

```sh
python3 ports/wasm/build.py --emcc emcc --output .build/wasm
node --test ports/wasm/*.test.mjs
node ports/wasm/runner.mjs bench multiply 48 1000 17
```

The build produces `.build/wasm/matrix.mjs` and `.build/wasm/matrix.wasm` for
Node.js, browser pages, and browser workers.
Emscripten 6.0.10 is tested; its compiler requires Python 3.10 or newer. When
necessary, set `EMSDK_PYTHON` to a supported interpreter. Node.js 22 or newer is
recommended. The build retains WASM function names for profiling, disables fused
multiply/add and fast-math, and permits memory growth. Allocation failures are
reported explicitly. Generated artifacts belong in the ignored `.build` folder.

Initialize one API instance and dispose every owned matrix:

```js
import { createMatrixAPI } from "./ports/wasm/matrix.mjs";

const { Matrix } = await createMatrixAPI();
const a = new Matrix(2, 3, [1, 2, 3, 4, 5, 6]);
const b = new Matrix(3, 2, [7, 8, 9, 10, 11, 12]);
let product;
try {
  product = a.multiply(b);
  console.log([...product.toArray()]); // [58, 64, 139, 154]
} finally {
  product?.dispose();
  a.dispose();
  b.dispose();
}
```

`Matrix` provides `rows`, `cols`, `get`, `set`, `copy`, `add`, `subtract`, `scale`,
`transpose`, `multiply`, `cross`, `row`, `column`, `trace`, `determinant`,
`eigenSymmetric`, `eigenGeneral`, `triangular`, and `toArray`;
`Matrix.identity(size)` creates an identity matrix. Arithmetic creates
independent matrices. `row`, `column`, and `toArray` return independent
`Float64Array` copies. Dimensions and indices are checked, and inputs must be
finite JavaScript numbers. Zero-size dimensions are supported. The wasm32 address
space limits dimensions and allocations; oversized requests fail.

`a.determinant()` or `a.determinant("auto")` uses guarded formulas for sizes 0
through 3 and an exact triangular shortcut, falling back to partial-pivot LU.
`a.determinant("lu")` always uses LU. `a.determinant("cholesky")` uses Cholesky
only when explicitly requested and requires exact symmetry and positive
computed pivots; otherwise it throws. Determinant factors are scaled to avoid
intermediate product overflow or underflow. Factorization can still lose range,
and a finite final determinant follows ordinary double rounding and underflow.
These operations preserve the matrix. Unknown algorithm names throw a
`RangeError`.

`a.eigenSymmetric()` returns `{ values, vectors }` for an exactly symmetric real
matrix: `values` is a copied `Float64Array` sorted ascending, and `vectors` is an
independently owned `Matrix` with corresponding unit eigenvectors as columns.
The result satisfies `A V = V diag(values)`. Dispose `vectors` after use; `values`
needs no disposal. Options can be supplied as
`a.eigenSymmetric({ tolerance: 1e-12, maxSweeps: 50 })`; tolerance must be finite
and between zero and one, and `maxSweeps` an integer from 1 through 10000.

The solver uses cyclic Jacobi rotations with power-of-two scaling. The default
relative tolerance is `1e-12`, measured by the off-diagonal Frobenius norm against
the original matrix norm, with at most 50 sweeps. It preserves the input. The
largest-magnitude component of each eigenvector is made nonnegative; repeated
eigenvalues can have any orthonormal basis of their eigenspace. Accuracy is
relative to the matrix norm, so very small eigenvalues beside much larger ones
need not have small relative error. Exact diagonal inputs retain their original
values even across extreme scales. Empty input returns empty eigenvalues and a
0-by-0 eigenvector matrix. Nonsymmetric matrices, invalid options, failure to
converge, and eigenvalues outside the finite double range are rejected.

`a.cross(b)` creates the vector cross product. Rotation factories are
`Matrix.rotation2D(radians)`, `Matrix.rotationX(radians)`,
`Matrix.rotationY(radians)`, `Matrix.rotationZ(radians)`, and
`Matrix.rotationAxisAngle(axis, radians)`. Each returns an owned matrix that
must be disposed. The axis must belong to the same initialized API instance.

The cross product is defined for three-dimensional vectors, represented as
3-by-1 columns or 1-by-3 rows. The operands may use different orientations; the
result keeps the left operand's shape and owns independent storage. Components
use the right-hand rule: `x × y = z`. This is a vector operation, not a matrix
product. Arithmetic rejects nonfinite results, including product overflow.

Rotations are right-handed **active transformations of column vectors**, with
angles in radians. In 2D, `R = [[cos θ, -sin θ], [sin θ, cos θ]]`; a positive
quarter-turn maps x to y. In 3D, positive X rotates y toward z, positive Y
rotates z toward x, and positive Z rotates x toward y. An arbitrary axis may be
a row or column three-vector and must be finite and nonzero. Axis normalization
uses its largest component before computing its norm, so huge, tiny, and
subnormal axes work without magnitude overflow. Rodrigues' formula uses a
stable `1 - cos θ` expression for small angles. Rotation matrices can be applied
with ordinary matrix multiplication. Both operands and axes stay unchanged.

Each matrix owns C heap memory. `dispose()` is idempotent; subsequent access
throws. Garbage collection does not dispose matrices automatically. Where
available, `Symbol.dispose` is supported. Matrices from different API instances
cannot be combined. The wrapper refreshes heap views after memory growth; it
does not expose views that would become detached.

For a browser, serve the files over HTTP(S) and import the same wrapper. The
default module URL resolves `.build/wasm/matrix.mjs` relative to the repository.
For a deployed layout, pass the generated module's URL and optionally locate its
binary separately:

```js
const { Matrix } = await createMatrixAPI({
  moduleUrl: new URL("./assets/matrix.mjs", import.meta.url),
  locateFile: name => new URL(`./assets/${name}`, import.meta.url).href,
});
```

Initialization also accepts `wasmBinary`, a `Uint8Array` containing the compiled
binary, to load embedded bytes without fetching a separate `.wasm` file. The
generated JavaScript module is still supplied through `moduleUrl`. The offline
report uses this option inside its worker.

Serve `.wasm` as `application/wasm`. Browser loading follows normal same-origin
and CORS rules. The wrapper itself does not import Node-specific modules.

The Node runner implements [the shared protocol](../../benchmarks/PROTOCOL.md).
Set `LINEAR_A_WASM_MODULE` to an alternate generated `.mjs` file when needed.
Benchmark timing includes the JavaScript operation loop, JS/WASM calls, C result
allocation, checksum consumption, result transfer, and explicit disposal.
Ordinary matrix operations use a three-element checksum computed in WASM.
Eigendecomposition transfers all eigenvalues and eigenvectors into JavaScript;
its checksum consumes every eigenvalue and vector entry. Module loading, WASM compilation, and input copies
are outside the timer. The comparison therefore measures the reusable API's
interop cost as well as its math kernels. Node's V8 CPU profiler can observe both
JavaScript and retained WASM function frames in a separate profiling run.

The [comparison report](../../benchmarks/README.md#live-webassembly-in-the-report)
can embed this API, the generated module and binary, and a worker in one offline
HTML file. Its live controls run the 39 shared correctness fixtures and bounded
benchmarks directly in the browser. Each run uses a fresh worker with a Stop
button and a timeout. Live browser results and their bundle identifier are
separate from recorded Node timings and profiles; no live stack profile is
collected. A report with a passed WebAssembly result and the current build
available includes this bundle automatically. Standalone browser applications
can continue to serve the module and binary as described above.

### Instrumented multiplication

`python3 ports/wasm/build.py --trace` creates a separate `.build/wasm-trace`
module for the report's calculation stepper. Only this build exposes
`matrix.multiplyWithTrace(other)`, returning `{ matrix, steps, source }`. Dispose
the returned matrix normally; steps are independent JavaScript snapshots of
actual C multiply-add events. Inputs are limited to dimensions of at most eight.
The ordinary build and benchmark API do not include instrumentation.


`a.eigenGeneral({ maxIterations: 1000 })` handles finite real square matrices,
including complex spectra. It returns `{ valuesReal, valuesImag, vectorsReal,
vectorsImag }`: the first two are copied `Float64Array`s; the latter two are
independently owned `Matrix` objects that both require `dispose()`. Complex
values sort by `(real, imag)` and vector columns have unit complex norm,
satisfying `A (Vr + i Vi) = (Vr + i Vi) diag(real + i imag)`. Results preserve the
input, and all temporary allocations are released on failure.

`maxIterations` must be an integer from 1 through 100000 and bounds QR steps
between root deflations. The C backend uses power-of-two balancing and scaling,
Hessenberg reduction, and implicit real double-shift QR with machine-epsilon
deflation, adapted from [public-domain JAMA](https://math.nist.gov/javanumerics/jama/).
Failures to converge or produce finite output throw. General eigenvectors need
not be orthogonal or independent for defective matrices; rounding and
conditioning still limit accuracy. The symmetric `eigenSymmetric()` API remains
available. The Node runner supports the `eigen_general` check and benchmark
operation with split real/imaginary eigenvalues and eigenvector matrices.

## Linear systems

See [the shared solver guide](../SOLVING.md) for this implementation’s reusable
LU, Cholesky, and column-pivoted QR APIs, multiple right-hand sides, least squares,
condition diagnostics, ownership, and numerical limits.
