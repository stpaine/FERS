# Oversampling audit report

## Implementation and recorded run

The standalone harness uses the source copies in `verification/oversampling/src/`.
The recorded sweep evaluated 128 cases across four suites, oversampling ratios
from 1 to 64, and filter lengths from 8 to 64. Its overall pass rate was **64.84%**.

Current `libfers` accepts oversampling ratios **1 through 8** and uses a default
render filter length of **33**. `params::validateOversampleRatio` rejects values
outside that range. The harness extends the ratio range to characterize filter
behaviour; its source copies and recorded results form a separate audit workflow.

## Recorded suite results

| Suite | Pass rate | Observation |
| --- | --- | --- |
| Clocking | 100% | PRF quantization, skip, and CW buffer-sizing cases passed. |
| Pipeline | 100% | Thermal-noise and ADC pipeline cases passed. |
| Resampler | 37.5% | Short filters reduced round-trip gain at larger ratios. |
| Render | 21.8% | Phase/amplitude ripple and fractional-delay errors exceeded the configured tolerances in the failing cases. |

## Filter response

The copied resampler uses a single-stage Blackman-windowed sinc filter with
`2 * filter_length` taps and cutoff `1 / oversample_ratio`. The audit uses a
coefficient sum equal to the ratio as its DC-gain reference. Its working envelope
is `filter_length >= 4 * oversample_ratio`.

For `ratio = 64` and `filter_length = 8`, the recorded FIR sum was approximately
**6.7**, and `estimated_roundtrip_gain` was **0.0109**. The harness computes this
quantity as `(fir_sum / oversample_ratio)^2`.

The recorded values at `filter_length = 33` were:

| Ratio | Reported round-trip gain |
| --- | ---: |
| 8 | 1.0006 |
| 16 | 0.966 |
| 24 | 0.751 |
| 32 | 0.535 |
| 64 | 0.171 |

At ratio 8, the constant-render check exceeded its `< 1e-5` tolerance even though
the reported round-trip gain was close to one. The higher-ratio cases characterize
the standalone harness beyond the current library's accepted range.

## Current verification workflow

Run the standalone commands in [README.md](README.md) to generate a fresh sweep
and its manifests, CSV results, and analysis report. The main library's Debug and
Release CTest suites exercise its supported configurations, including the
all-zero ADC regression and streaming sample-boundary checks.

The ratio cap proposed in the original audit is implemented in
[`parameters.h`](../../packages/libfers/src/core/parameters.h). Filter redesign
and support for larger ratios remain separate implementation work.
