# Contributing

This is a hobby project for learning linear algebra and comparing languages
through eleven matrix implementations. Contributions to any language, the shared
contract, examples, math checks, and reports are welcome.

Work on a branch and open a pull request against `main` with a short explanation
of the change and how it was checked. Start with the
[implementation guide](ports/README.md) and the README for the language you are
changing. Each implementation has its own build and test commands; working on one
does not require installing every toolchain.

## Shared behavior and math changes

The [runner protocol](benchmarks/PROTOCOL.md) defines common float64 inputs,
outputs, errors, and timing boundaries. Language-specific APIs can expose
additional types or operations; document those differences in their own guides.
Changes to the shared contract should describe their effect on all implementations
and update the relevant runners and fixtures together.

For math changes, include a small worked example with expected values calculated
independently. Add a regression test for a bug before changing the algorithm.
Use rectangular and nonsymmetric matrices where appropriate: square or identity
matrices can hide indexing mistakes. For floating-point assertions, choose a
tolerance suitable for the values under test and explain it when it is not
obvious. Keep overflow, rounding, shape errors, and solver nonconvergence distinct.

Write numerical kernels as named phases: input validation, preparation,
factorization or iteration, and result acceptance. Give nested loops and failure
branches their own blocks, and name intermediate quantities by their role.
Explain choices such as reorthogonalization, pivot acceptance, scaling, dropped
fill, and true-residual stopping where they occur. State the numerical reason for
a threshold and distinguish machine epsilon from the smallest representable
positive value. The C# [GMRES implementation](linear-A/linear-A/Gmres.cs) is one
example of this structure.

For readability refactors, preserve arithmetic order, tolerances, ownership,
exceptions, and stopping precedence. Verify residuals and failure cases against
the existing contract; do not silently change numerical policy during cleanup.

Shared correctness references live in `benchmarks/*reference.py`; standalone
unit tests live beside each implementation or its test project. Verification
compares results with those independent references. No implementation is the
numerical oracle for another. For eigenpairs, check spectra and residuals;
require orthogonality only for the symmetric solver, and allow sign changes and
valid bases for repeated eigenspaces.

## Checks before a pull request

Run the affected implementation's standalone tests and formatting checks as
documented in its guide. Then verify the relevant comparison runners from the
repository root, selecting their IDs explicitly:

```sh
# Example subset; choose the implementations affected by your change.
python3 benchmarks/run.py --verify-only --implementations go python rust --require-all --output benchmarks/reports/verification

# Full verification when all toolchains are available, on an ARM64 host.
python3 benchmarks/run.py --verify-only --require-all --output benchmarks/reports/verification
```

Use `--require-all` so a missing selected toolchain is reported as a failure.
For changes to the harness or report controls, also run:

```sh
python3 -m unittest discover -s tests/benchmark_harness -v
node --test tests/report-theme.test.mjs tests/timing-report.test.mjs
```

The [Matrix ports workflow](.github/workflows/ports.yml) builds and checks all
eleven implementations on ARM64 Linux. The [CI workflow](.github/workflows/ci.yml)
also checks the C# runtime targets on Linux, Windows, and macOS. Their platform
coverage differs; each implementation follows the same shared math contract.

## Dependencies and formatting

Keep numerical kernels dependency-free and use the repository's formatting
conventions. When updating a toolchain or dependency, update its manifest and
committed lockfile where applicable, then rebuild and run its tests using the
documented setup. Toolchain floors, ownership rules, and runtime-specific checks
belong in the corresponding language guide. The
[C# guide](linear-A/README.md#development-checks) covers its SDK, analyzers,
formatting, generic numeric types, and legacy .NET APIs.

## Release documentation

Treat [release/manifest.json](release/manifest.json) as the release-version and
capability source of truth; `python3 release/build.py --version` prints its
version. When preparing a release, update the package versions, current
`docs/api/<version>/README.md`, root README's API link, and affected shared
contracts together. Regenerate the compatibility table with
`python3 release/build.py --write-table`, then run
`python3 -m unittest discover -s tests/release -v`. Keep older API documents
versioned and distinguish capability introduction dates from the current release.
Use the manifest version in installation examples instead of copying a literal
version into the general release guide.

## Benchmarks and public reports

Use a separate output directory for new measurements and verification runs so
earlier results stay intact. Do not edit timings to fit expectations. Record
algorithm, compiler, and runtime differences, and distinguish instrumented
profiles from comparative timing runs. The [benchmark guide](benchmarks/README.md)
documents the methodology and report generation pipeline.

Generate reports through that pipeline so local paths and identifying metadata
are sanitized. When changing presentation only, regenerate HTML from the saved
JSON without rerunning benchmarks or changing measurements. Check the interactive
controls in a browser. Selected report snapshots are tracked in Git and publish
to Pages when their committed files change; build outputs, caches, and new run
directories remain ignored.

Capture new measurements for publication with `benchmarks/reproduce.py` after
committing the source. It builds in a fresh clone of the recorded commit and
records the Git tree, invocation, and artifact hashes in `provenance.json`.
Publish that file with the complete report and keep the source commit in the
published branch's history. Run `python3 benchmarks/check_provenance.py` before
submitting report changes; this also covers the separately published experiments.
See [the capture commands](benchmarks/README.md#publishing-reproducible-source-records).

Historical measurements whose dirty source was not archived cannot be made
recoverable by adding a hash or labeling them clean. Preserve their disclosed
limitations and fixed artifact digests until replacing them with a new capture.
Do not expand the historical exceptions to admit new measurements.
