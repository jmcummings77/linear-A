# Shared numerical bug finding

The generated checks use the same protocol for all eleven implementations.
They supplement the exact, analytic, and hand-written regression tests; they do
not replace them. A passing campaign establishes only that the tested inputs
passed. Agreement between implementations is not an independent mathematical
reference: several ports share native numerical kernels.

## Run, replay, and extend

```sh
# Build every port, then check two consecutive seeds and saved regressions.
python3 benchmarks/bugfinding.py

# Use existing runner builds for a larger, reproducible campaign.
python3 benchmarks/bugfinding.py --no-build --seed 20261002 --seeds 20

# Fast local subset; no implicit skipping of requested implementations.
python3 benchmarks/bugfinding.py --no-build --implementations c python --seeds 20

# Replay a reduced counterexample against any selected ports.
python3 benchmarks/bugfinding.py --no-build --implementations c python \
  --replay .build/bugfinding/failure-EXAMPLE.json
```

The default output is `.build/bugfinding/`, outside tracked content. A summary
records each property/implementation result, compiler versions, timestamp and
Git provenance. Failure artifacts include the generator version, seed, original
case, reduced case, failure category, reduction budget and repeat confirmation.
They contain matrix data and sanitized diagnostics, not environments or local
executable paths. CI uploads these JSON files even if a property fails.

To retain an actual counterexample, inspect the artifact and copy it to
`tests/bugfinding/regressions/` with a descriptive name. Every ordinary campaign
replays those files. Review the cause and expected behavior before committing;
do not change the oracle merely to make a failing implementation pass. The
initial zero-inner-dimension fixture is a contract regression, not a claimed
newly discovered defect. Deliberately faulty implementations in the harness
unit tests verify that detection and reduction actually work.

## Families and evidence

Each seed generates ten cases with dimensions at most six:

- Rectangular and empty products, checked using exact rational dot products.
- Transpose involution and `(AB)ᵀ = BᵀAᵀ`, also independently checked against
  reference results so two compensating bugs cannot satisfy the identity alone.
- Dense and duplicate-row determinants, with an independent Leibniz oracle,
  plus `det(Aᵀ) = det(A)`.
- Exact dyadic cancellation in addition.
- Strictly diagonally dominant systems with multiple right-hand sides.
- Gram matrices plus identity for Cholesky solves and symmetric eigenchecks.
- Full-column-rank tall matrices with noisy observations for least squares.

Solver answers use exact rational Gauss–Jordan elimination. The least-squares
reference uses normal equations **only in exact arithmetic**. The existing
shared comparison requires finite, correctly shaped outputs with `1e-9`
relative / `1e-10` absolute tolerance. These conservative limits are suitable for
these small, controlled families; they are not a universal error bound for
arbitrarily ill-conditioned systems. Eigenchecks use the existing residual and
orthogonality criteria, allowing different bases for repeated eigenspaces.
Cross-port comparisons cover scalar/matrix results, not eigenvector equality.

Generation version 1 uses Python's seeded `random.Random` and records actual
inputs on failure. Saved case data is the definitive replay input if a future
Python version or generator changes its sequence. Raw invalid protocol text and
unbounded dimensions are not passed through the replay interface.

## Reduction

The reducer tries dimension removal, then replacing entries with zero or unit
magnitude. Every accepted change reduces dimensions or value magnitude. It
preserves conformable shapes, exact nonsingularity for solver references, and
symmetry/positive definiteness for Cholesky inputs. It checks candidates against
the same failing implementation and failure category; differential failures
also rerun the comparison implementation.

Reduction defaults to 100 attempts and a 30-second search budget, plus any
in-flight runner request (each limited to 30 seconds). The result is a locally
simplified example, **not a globally minimal proof** or a guarantee of the same
root cause. Failure categories can include different defects; inspect the saved
original and reduced inputs together. A final repeat records whether the reduced
failure reproduced. Unstable failures still fail the campaign and remain saved.

## Native fuzzing and memory checks

```sh
# Clang with libFuzzer: coverage-guided mutations + ASan + UBSan.
python3 tests/fuzz/run.py --runs 10000
python3 tests/fuzz/run.py --assembly --runs 10000  # ARM64 hosts

# Replay a crash input, keeping the exact bytes.
python3 tests/fuzz/run.py --replay .build/fuzz/crash-EXAMPLE

# Apple Clang without libFuzzer: deterministic sanitizer smoke campaign.
python3 tests/fuzz/run.py --smoke --runs 10000
```

The native target decodes at most 128 bytes into bounded dimensions and dyadic
values, exercising multiplication, transpose, input preservation, invalid shape
and index errors, nonfinite input rejection, LU/Cholesky/QR factorization,
repeated solves and disposal. Products have an exact bounded dyadic oracle;
SPD systems have known solutions. A four-input starting corpus covers empty,
small and odd shapes and byte extremes; mutated inputs persist under `.build`.
The decoder intentionally keeps pointers valid; this is not arbitrary-pointer
fault injection. Assembly instructions themselves are not compiler-instrumented,
although the C allocation, validation, and shared solver paths are instrumented.

Linux CI runs real coverage-guided C and ARM64 campaigns. The deterministic
fallback is explicitly not coverage-guided and must not be reported as such.
Raw sanitizer logs stay local because runtimes can print identifying host paths.
Crash bytes contain only generated numerical inputs; retain/replay them when
investigating a failure. No bug is implied when the campaign finds none.
