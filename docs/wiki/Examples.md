# Examples

The `examples/` directory contains working scenarios and analysis scripts. These are the best starting points for learning FERS because they show complete XML files, waveform generation, CLI usage, and post-processing.

Python analysis scripts generally need:

- `numpy`
- `h5py`
- `matplotlib`
- `scipy` for some signal-processing examples

## Mixed Pulsed And CW Scenario

Directory:

```text
examples/mixed_scenario
```

This scenario contains:

- One pulsed monostatic radar.
- One CW monostatic radar.
- One moving target.
- A generated pulsed waveform stored as HDF5.
- Python analysis for range and Doppler behavior.

Run it:

```bash
cd examples/mixed_scenario
python3 genpulse.py
../../build/release/packages/fers-cli/fers-cli example.fersxml --out-dir=. --log-level=INFO -n=2
python3 analysis.py
```

Expected generated files:

```text
pulse.h5
PulsedRadar_results.h5
CWRadar_results.h5
pulsed_analysis.png
cw_analysis.png
```

What to look at:

- `genpulse.py` shows how to create a pulse file with `/I/value` and `/Q/value`.
- `example.fersxml` shows two radar modes in one scenario.
- `analysis.py` shows how to reconstruct I/Q with `fullscale` and estimate range/Doppler.

Use this example when you want to learn:

- How pulsed and CW radars can coexist in one scenario.
- How to use a file-backed pulsed waveform.
- How receiver output differs between pulsed chunk output and streaming output.

## FMCW Monostatic Dechirp

Directory:

```text
examples/fmcw_monostatic_dechirp
```

This scenario contains:

- One monostatic FMCW radar.
- One target.
- A linear chirp waveform.
- Receiver-side dechirping inside FERS.
- IF output analysis.

Run it:

```bash
./build/release/packages/fers-cli/fers-cli examples/fmcw_monostatic_dechirp/example.fersxml --out-dir=./results/mono_fmcw --log-level=INFO -n=2
python3 examples/fmcw_monostatic_dechirp/analysis.py --results-dir ./results/mono_fmcw --output-dir ./results/mono_fmcw
```

Expected generated files:

```text
MonoFmcwRadar_results.h5
monostatic_dechirp_analysis.png
```

What to look at:

- The waveform uses `<fmcw_linear_chirp>`.
- The monostatic radar uses `<fmcw_mode dechirp_mode="physical">`.
- The dechirp reference is `source="attached"`, which is the normal monostatic setup.
- The analysis script reads dechirped IF output and compares it to the expected beat signal.

Use this example when you want to learn:

- Basic FMCW XML structure.
- How to request dechirped IF output.
- How `if_sample_rate` affects the generated result.

## FMCW Bistatic External Dechirp

Directory:

```text
examples/fmcw_bistatic_external_dechirp
```

This scenario contains:

- A separate FMCW transmitter.
- A separate receiver.
- One target.
- Raw received FMCW output from FERS.
- Python-side dechirping in the analysis script.

Run it:

```bash
./build/release/packages/fers-cli/fers-cli examples/fmcw_bistatic_external_dechirp/example.fersxml --out-dir=./results/bistatic_fmcw --log-level=INFO -n=2
python3 examples/fmcw_bistatic_external_dechirp/analysis.py --results-dir ./results/bistatic_fmcw --output-dir ./results/bistatic_fmcw
```

Expected generated files:

```text
BistaticFmcwRx_results.h5
bistatic_external_dechirp_analysis.png
```

What to look at:

- The transmitter and receiver are separate objects.
- The receiver uses `<fmcw_mode dechirp_mode="none">`.
- The HDF5 output contains raw received complex baseband.
- The analysis script performs the reference mixing outside FERS.

Use this example when you want to learn:

- Bistatic scenario structure.
- How to export raw FMCW receive data.
- How external dechirping compares with built-in dechirping.

## XML Antenna Pattern Example

File:

```text
examples/example_antenna_pattern.xml
```

This is a standalone antenna pattern asset. A scenario can reference it like this:

```xml
<antenna name="ExampleXmlAntenna" pattern="xml" filename="example_antenna_pattern.xml"/>
```

The pattern file defines azimuth and elevation gain samples. Use this when an analytic antenna model is not enough.

What to look at:

- `unit="deg"` or `unit="rad"` controls angle units.
- `format="dBi"` or `format="linear"` controls gain value format.
- `symmetry="mirrored"` can mirror positive-angle samples to negative angles.

## Starting A New Scenario From An Example

The fastest way to create a new scenario is to copy the closest example and change one thing at a time:

1. Copy the example directory.
2. Rename the scenario and output objects.
3. Change platform positions and motion.
4. Change waveform parameters.
5. Run the scenario.
6. Update the analysis script only after the simulation output looks reasonable.

For early debugging:

- Use `pattern="isotropic"` antennas.
- Use simple static platforms.
- Disable schedules until the geometry works.
- Start with `nopropagationloss="true"` only if you need to confirm signal timing without realistic link loss.
- Add antenna pointing, schedules, noise, and propagation-loss realism after the basic result is correct.

## Common Example Mistakes

| Problem | Fix |
| --- | --- |
| `pulse.h5` is missing | Run `python3 genpulse.py` in `examples/mixed_scenario`. |
| Analysis script cannot find results | Pass the same output directory to `fers-cli` and the analysis script. |
| HDF5 file exists but analysis gives zeros | Reconstruct I/Q with `fullscale` before calculating power or spectra. |
| FMCW output sample count is unexpected | Check `if_sample_rate`, chirp duration, chirp count, and schedule duration. |
| KML geometry looks misplaced | Set the KML/geospatial `<origin>` and `<coordinatesystem>` explicitly. |
