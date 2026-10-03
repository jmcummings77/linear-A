# Numerical playground

[Open the playground](https://jmcummings77.github.io/linear-A/applications/) ·
[Portable HTML](index.html)

Three local browser applications use the public WebAssembly matrix API:

- **PCA:** drag/add/remove points or edit coordinates with keyboard controls.
  The sample covariance uses n−1; symmetric eigenvectors give principal axes,
  variances and projections. Zero variance and tied directions are explicit.
- **Least squares:** fit degree 1–3 polynomials to editable observations using
  column-pivoted QR. Coefficients are for t=x/5, not raw x. Residual segments,
  training RMSE and R² expose the fit; rank-deficient designs fail explicitly.
- **Image rank:** paint a bounded 24×24 grayscale matrix or choose a preset.
  The rank slider reconstructs it from its leading singular directions, showing
  error, energy and factor storage. No input leaves the page.

## Numerical scope

The image demo obtains U from the symmetric eigendecomposition of AAᵀ and
reconstructs Uₖ(UₖᵀA), avoiding division by small singular values. This is an
educational Gram-matrix route to an SVD reconstruction, not a general-purpose SVD
implementation. It squares the condition number and cannot accurately resolve
arbitrarily small singular values. Tiny negative computed Gram eigenvalues are
clamped to zero for the displayed singular values.

Errors use the unclipped reconstruction; only grayscale display values are
clipped to [0,1]. Factor storage counts k(2n+1) floating-point scalars versus n²
entries. It does not estimate a compressed image file, and high ranks can use
more storage. Images are generated or painted locally, not borrowed assets.

PCA directions can change sign without changing their meaning. For tied
variances there is no unique first direction; axis overlays are suppressed.
Curve fitting shows training error only, without confidence bands or prediction
claims. Fixed coordinate bounds and normalized x keep these examples small;
they are not guarantees for arbitrary real matrices.

## Build and verify

```sh
python3 ports/wasm/build.py --output .build/wasm
node --test ports/wasm/applications.test.mjs
python3 applications/build.py
python3 -m http.server 8000 --directory applications
```

The builder verifies all shared WASM fixtures before embedding the sanitized
runtime into a standalone `index.html`. Serve it over HTTP(S) to support module
workers consistently. CI tests the applications with the actual WASM kernels
and regenerates the page; Pages publishes the committed HTML without running
benchmarks. Application calculations run in workers with 300 ms debouncing,
stale-result rejection, bounded inputs and timeouts. Native allocations are
released in `finally` blocks. No external scripts, network computation, user
tracking or dependencies are needed at runtime.

Tests cover analytic PCA, translation invariance, tied/zero variance, known
polynomial coefficients, residual orthogonality, rank rejection, monotone
reconstruction error, full-rank recovery, rank-one and blank images, storage
counts and repeated resource use. The embedded bundle hash identifies its
runtime and application math. Existing benchmark measurements stay separate.
