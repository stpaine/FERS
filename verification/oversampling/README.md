# Oversampling Audit Harness

This directory contains a standalone C++ executable for auditing the oversampling implementation used by `libfers`.

## Practical Envelope

The copied resampler uses a single-stage Blackman-windowed sinc with `2 * filter_length` taps and cutoff `1 / oversample_ratio`.

- For modest ratios, this is adequate.
- For high ratios with short filters, the FIR sum falls well below `oversample_ratio`, so the harness now emits a runtime warning with the measured FIR sum and estimated round-trip gain.
- As a practical rule of thumb, settings around `filter_length >= 4 * oversample_ratio` keep the DC-gain error small enough for the current design to behave predictably.

Warnings are advisory only in this harness version. Out-of-envelope cases still run so the matrix can record how far performance degrades.

## Scope

- Source-faithful copies of the core oversampling path:
  - `packages/libfers/src/signal/dsp_filters.*`
  - `packages/libfers/src/signal/radar_signal.*`
  - `packages/libfers/src/interpolation/interpolation_filter.*`
  - `packages/libfers/src/serial/response.*`
  - `packages/libfers/src/processing/signal_processor.*`
  - `packages/libfers/src/processing/finalizer_pipeline.*`
- Formula-level probes for oversampling-dependent callers:
  - transmitter PRF quantization
  - receiver window PRF/skip quantization
  - CW timestep and CW buffer sizing
- Deferred executable coverage:
  - `packages/libfers/src/noise/falpha_branch.cpp` uses `signal::DecadeUpsampler`, but this harness only documents that path in manifests and source maps.

## Usage

Configure and build:

```bash
cmake -S verification/oversampling -B verification/oversampling/build
cmake --build verification/oversampling/build
```

Run a single suite:

```bash
./verification/oversampling/build/src/oversampling_audit suite render --ratio 4 --filter-length 33 --seed 42 --output-dir verification/oversampling/output
```

Run all suites for one parameter set:

```bash
./verification/oversampling/build/src/oversampling_audit all --ratio 4 --filter-length 33 --seed 42 --output-dir verification/oversampling/output
```

Run a matrix:

```bash
./verification/oversampling/build/src/oversampling_audit matrix --ratios 1,2,4,8 --filter-lengths 8,33 --seed 42 --output-dir verification/oversampling/output
```

Run the comprehensive Python driver and post-analysis:

```bash
python3 verification/oversampling/run_and_analyze.py --clean-output
```

## Output

Each suite writes:

- `manifest.json`
- `summary.csv`
- suite-specific CSV artifacts
- `README.txt` describing the generated files

The matrix and `all` commands create per-suite subdirectories so each process can initialize the interpolation-filter singleton with one filter length only.

## Python Analysis

`run_and_analyze.py` automates the end-to-end workflow:

- optionally rebuilds the standalone harness
- runs a matrix over ratios and filter lengths
- parses every generated `manifest.json`, `summary.csv`, and suite-specific CSV
- emits aggregated analysis artifacts under `output/.../analysis/`

Generated analysis artifacts:

- `case_analysis.csv`: one row per `(ratio, filter_length, suite)` case
- `suite_summary.csv`: pass/fail rollup by suite
- `global_summary.json`: machine-readable aggregate report
- `report.md`: human-readable Markdown report with key findings and failing cases

Analysis notes:

- Render expected-value checks include the upsample FIR DC gain, so they isolate interpolation behavior instead of conflating interpolation with resampler gain loss.
- Resampler cases can still fail by policy when the FIR sum deviates materially from the requested ratio; those failures should be read as out-of-envelope configuration warnings, not harness-oracle bugs.
- Pipeline downsampling RMS uses a `0.06` threshold, and resampler peak alignment allows a `+/-1` sample shift.

Default comprehensive sweep:

- ratios: `1,2,4,8,16,24,32,64`
- filter lengths: `8,16,33,64`
- seed: `42`
