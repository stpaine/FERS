# Using FERS

FERS runs a `.fersxml` scenario and writes receiver I/Q data to HDF5 files. A normal workflow is:

1. Create or choose a scenario. You can write XML directly or use [[FERS UI]] to build it visually.
2. Run it with `fers-cli` or from the UI Simulation Run view.
3. Open the generated HDF5 files in an analysis script.
4. Compare the simulated signals with the radar effect you expected.

## Run A Scenario

```bash
./build/release/packages/fers-cli/fers-cli path/to/scenario.fersxml --out-dir=./results
```

If `--out-dir` is omitted, FERS writes results beside the scenario file.

Generate KML instead of running the simulation:

```bash
./build/release/packages/fers-cli/fers-cli path/to/scenario.fersxml --kml=scene.kml
```

## What A Scenario Describes

A scenario describes the radar scene, not the analysis you want to run afterward.

At a high level it contains:

- Global simulation settings, such as start time, end time, sample rate, and coordinate system.
- Waveforms, such as pulsed waveform files, CW tones, or FMCW chirps.
- Timing sources, which model oscillator frequency, phase, and noise behavior.
- Antennas, which determine gain as a function of direction.
- Platforms, which move through the scene.
- Radar objects on platforms: transmitters, receivers, monostatic radars, and targets.
- Optional schedules, which turn transmitters and receivers on only during selected time intervals.

The complete XML reference is in [[XML Schema Reference]].

## Choose The Right Radar Mode

| Mode | Use when | Scenario pieces |
| --- | --- | --- |
| Pulsed | You have a pulse waveform and want range-gated receiver windows. | `<pulsed_from_file>`, `<pulsed_mode>`, `prf`, `window_skip`, `window_length`. |
| CW | You want continuous-wave Doppler-style output. | `<cw/>`, `<cw_mode/>`. |
| FMCW | You want chirp-based ranging or range-Doppler analysis. | `<fmcw_linear_chirp>` or `<fmcw_triangle>`, `<fmcw_mode>`, optional dechirp settings. |

The waveform type and radar mode must match. For example, a `<cw/>` waveform must be used with `<cw_mode/>`.

## Monostatic And Bistatic Layouts

Use `<monostatic>` when the transmitter and receiver are colocated and share the same antenna, waveform, and timing reference.

Use separate `<transmitter>` and `<receiver>` elements when you need a bistatic or multistatic scene. The transmitter and receiver can be on different platforms, use different antennas, and have separate schedules.

## Output Files

FERS writes one HDF5 file per receiver:

```text
<receiver-name>_results.h5
```

Monostatic radars also write through their receiver side, using the monostatic name.

### Pulsed Output

Pulsed receivers write one pair of datasets per receiver window:

```text
chunk_000000_I
chunk_000000_Q
chunk_000001_I
chunk_000001_Q
...
```

Each chunk represents one receiver window. Use `window_skip` and `window_length` from the scenario, plus the metadata stored in the HDF5 file, to convert samples to range.

### CW And FMCW Output

Streaming receivers write:

```text
I_data
Q_data
```

For FMCW with dechirping enabled, these datasets contain the dechirped IF output. Without dechirping, they contain the received complex baseband signal.

## Reconstruct Physical I/Q Values

FERS normalizes stored samples. To recover the physical sample values, multiply by the `fullscale` attribute.

```python
import h5py

with h5py.File("Rx_results.h5", "r") as h5:
    fullscale = h5.attrs["fullscale"]
    iq = (h5["I_data"][:] + 1j * h5["Q_data"][:]) * fullscale
```

For pulsed chunks, `fullscale` is stored on the chunk datasets:

```python
with h5py.File("PulsedRadar_results.h5", "r") as h5:
    i = h5["chunk_000000_I"][:]
    q = h5["chunk_000000_Q"][:]
    fullscale = h5["chunk_000000_I"].attrs["fullscale"]
    iq = (i + 1j * q) * fullscale
```

## Metadata

FERS stores metadata in the HDF5 output, including receiver mode, sample rates, sample counts, time span, fullscale values, and FMCW settings when applicable.

Use the metadata whenever possible instead of hard-coding assumptions in analysis scripts.

## Sample Rates

The most important rates are:

| Setting | Meaning |
| --- | --- |
| `<rate>` | Main output sample rate in Hz. |
| `<simSamplingRate>` | Rate used for geometry interpolation in pulsed response generation. |
| `<oversample>` | Internal oversampling multiplier. Final output is still usually written at `<rate>`. |
| `<if_sample_rate>` | Optional FMCW IF output rate after dechirping and resampling. |

Use a high enough `<rate>` to represent the signal bandwidth you expect. For FMCW, the waveform definition must also satisfy the anti-aliasing checks described in [[XML Schema Reference]].

## Validation

By default, FERS validates the XML structure before running:

```bash
./build/release/packages/fers-cli/fers-cli scenario.fersxml
```

You can skip XML validation:

```bash
./build/release/packages/fers-cli/fers-cli scenario.fersxml --no-validate
```

Use `--no-validate` only when diagnosing a file. A file that skips schema validation can still fail later if references, numbers, radar modes, or physics settings are invalid.

## KML Export

KML export helps inspect platform paths and scenario geometry in a map or geospatial viewer:

```bash
./build/release/packages/fers-cli/fers-cli scenario.fersxml --out-dir=./results --kml
```

Default output:

```text
<out-dir>/<scenario-name>.kml
```

KML uses the scenario's `<origin>` and `<coordinatesystem>` settings.

If KML generation fails, `fers-cli` exits with a nonzero status and logs the reason. Common causes are an invalid output path, missing parent directory, scenario loading failure, or invalid coordinate settings.

## Troubleshooting

| Symptom | Likely cause |
| --- | --- |
| XML validation fails | Elements are out of order, required attributes are missing, or a mode block does not match the schema. |
| Scenario loads but does not run | A referenced waveform, antenna, timing source, or platform object name is wrong. |
| Output is all zeros | Receiver schedule/window may miss the return, antenna pointing may be wrong, target may be outside the simulated time span, or path loss may make the signal very small. |
| Pulsed range looks offset | Check `window_skip`, `window_length`, waveform sample rate, and speed of propagation `c`. |
| FMCW run fails validation | Check chirp bandwidth, chirp duration, chirp period, `<rate>`, and `<oversample>`. |
| KML looks wrong geographically | Check `<origin>`, `<coordinatesystem>`, UTM zone, and hemisphere. |

## Old XML Files

Run the migration script for legacy files:

```bash
python3 scripts/migrate_fers_xml.py old.fersxml migrated.fersxml
```

Then inspect the result manually. Old `<export>` and `<multipath>` settings are not supported in the current scenario format.
