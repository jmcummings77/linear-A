# Compare benchmark snapshots

[Public comparison](https://jmcummings77.github.io/linear-A/compare/) ·
[Portable HTML](index.html)

Compare saved runs without executing benchmarks or checking out old code. Select
baseline/candidate snapshots, or load two local `results.json` files in the
browser. Filter by implementation, operation and size, inspect individual sample
times, swap the direction, and download a numeric comparison JSON.

## Files and Git revisions

```sh
# Generate a standalone comparison from two saved runs.
python3 benchmarks/compare.py /path/to/baseline/results.json \
  /path/to/candidate/results.json --output /tmp/comparison.html

# Read the report file from two commits; no checkout or benchmark execution.
python3 benchmarks/compare.py --from-git HEAD~3 HEAD \
  --path benchmarks/reports/latest/results.json --output /tmp/revisions.html

# Refresh the public catalog from the existing committed snapshots.
python3 benchmarks/compare.py --catalog

node --test tests/comparison*.test.mjs
python3 -m unittest discover -s tests/benchmark_harness
```

The browser has no Git filesystem access: use the CLI for arbitrary Git refs or
load JSON files from those revisions. “Saved in” identifies the commit containing
the report file. The report's **recorded revision** and source fingerprint identify
the measured source, which can be older or dirty. The initial public catalog has
multiple runs recording the same revision; it is not a fabricated before/after
experiment. Verification-only reports without timing rows are omitted from the
catalog. Nothing here rewrites historical measurement JSON or raw profiles.

## Matching and interpretation

Rows join on exact implementation ID, operation, and matrix size. Missing or
failed rows remain visible and never receive invented zero timings. Each side's
median, median absolute deviation (MAD), min and max are recomputed from saved
`ns_per_op` samples. The loader checks positive finite timings, positive integer
iteration counts, consistency with elapsed time / iterations, duplicate row
identities and schema version. It does not trust a stale stored median.

The ratio is **candidate median / baseline median**: less than one is faster.
The percent change is `(ratio − 1) × 100`. The displayed envelope is
`[candidate min / baseline max, candidate max / baseline min]`. It describes the
observed batches; it is **not a confidence interval**, a probability of regression,
or a significance test. There are no performance gates or automatic regression
claims. Repeated batches on one machine do not establish general performance.

Different seeds, numeric/layout contracts, timing boundaries, warmup rules or
explicit workload versions withhold ratios. Toolchain and build-command changes
warn per row. Hardware/OS, suite, sampling and provenance differences warn at the
snapshot level. Basic hardware details cannot prove identical hosts. Legacy
reports lack explicit workload versions and process-unit metadata; the dashboard
exposes that uncertainty instead of backfilling historical evidence. Before
attributing a ratio to code, inspect the relevant protocol and source changes.

New standard and determinant reports record `sample_unit: fresh_process_per_batch`
and an explicit workload version. These identify the existing subprocess-per-
sample behavior without changing timings. Bump the workload version when changing
generated inputs or the operation contract; changes to implementation algorithms
alone do not require a new workload version. Profiles remain separate evidence.

Allocation is already included in the timing boundary, but the standard snapshots
do not record allocation counts. The comparison reports this limitation and does
not invent allocation deltas. The locality experiment uses a different schema and
is deliberately rejected by this loader.

## Privacy and export

The Python renderer reuses public-report sanitization and an explicit field list:
no runner paths, environment dictionaries, raw profiles or failure diagnostics
are embedded. Labels for arbitrary file inputs are generic. Browser-loaded files
stay local and are never uploaded; their metadata is displayed as text, not HTML.
A downloaded comparison includes numerical samples, validated revision/source
hashes and fixed diagnostic messages. It omits arbitrary toolchain strings,
command metadata and other imported free text. It is a comparison artifact, not
a lossless export of the original reports. Original inputs are never changed.

The page is standalone, supports light/dark themes, and uses no network requests
or external assets at runtime. Pages publishes the committed HTML as `/compare/`.
