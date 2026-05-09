# FERS UI

`fers-ui` is the graphical desktop interface for FERS. It is designed for building, inspecting, running, and exporting radar scenarios without editing XML by hand.

Use the UI when you want a visual scenario workbench. Use `fers-cli` when you already have a `.fersxml` file and want a simple terminal run.

## Current Status

The UI is usable but still marked beta and under active development. Expect some incomplete workflows and rough edges. The CLI remains the most stable way to run production scenario batches.

Known practical limitations:

- The 3D Scenario view requires working WebGL2 support.
- Some macOS systems may experience WebGL context loss on launch.
- The UI can edit the current FERS scenario model, but advanced or unusual XML files may still need manual review after import/export.
- The UI does not author standalone antenna-pattern XML files directly. It can reference XML or HDF5 antenna-pattern files.
- The UI is not an HDF5 plotting tool. It shows output metadata after a run, but detailed signal analysis is still done in Python, MATLAB, Julia, or another analysis environment.

## Running The UI From Source

Install the normal FERS build requirements first. See [[Installing FERS]].

From the repository root:

```bash
bun install
bun ui:dev
```

Build a release bundle:

```bash
bun ui:build
```

The UI uses the same FERS simulation engine as the CLI. You still need vcpkg, CMake, a C++23 compiler, Rust, and the Tauri prerequisites for your platform.

## Main Views

The left application rail switches between the main UI workspaces:

| View | Purpose |
| --- | --- |
| Scenario Builder | Build and inspect the radar scene. |
| Asset Library | Save reusable waveforms, timings, antennas, and platforms. |
| Simulation Run | Choose output settings, run the simulation, generate KML, and inspect output metadata. |

The lower rail buttons open:

| Button | Purpose |
| --- | --- |
| Raw Logs | Show simulator log output and set the FERS log level. |
| Settings | Set application-level options such as preview playback duration and log level. |
| About | Version, license, and third-party license information. |

## Scenario Builder

Scenario Builder is the main editing workspace. It has three areas:

- Scenario Explorer on the left.
- 3D scenario viewport and timeline in the middle.
- Properties panel on the right.

### Scenario Explorer

The Scenario Explorer lists the major scenario objects:

- Global Parameters.
- Waveforms.
- Timings.
- Antennas.
- Platforms.

Use the add buttons beside each section to create new items. Select an item to edit it in the Properties panel.

Platforms can contain:

- Monostatic radars.
- Transmitters.
- Receivers.
- Targets.

The import and export buttons at the top of the Scenario Explorer load or save `.fersxml` files.

### 3D Viewport

The 3D viewport previews the scene at the current timeline time. It can show:

- Platform positions.
- Platform labels.
- Body axes.
- Motion paths.
- Velocity vectors.
- Antenna patterns.
- Boresight arrows.
- RF preview links.
- Link labels with approximate power or illumination values.

The viewport is a preview, not the simulation result. It helps you check geometry, pointing, schedules, and obvious link problems before running the full signal simulation.

### View Controls

The View Controls panel can:

- Frame the whole scene.
- Focus on the selected item.
- Follow the selected item as time changes.
- Toggle platform labels, body axes, velocity vectors, antenna patterns, boresight arrows, motion paths, and RF links.
- Filter RF link display by monostatic, bistatic illuminator, bistatic scattered, and direct interference paths.

### Timeline

The timeline previews platform motion and pointing over simulation time. It does not run the signal simulation.

Use it to:

- Scrub through time.
- Play/pause scene preview.
- Step backward or forward.
- Check whether moving platforms and pointing behavior look plausible.

The playback duration setting controls how quickly the preview plays in real time.

## Editing Scenario Properties

The Properties panel changes based on the selected item.

### Global Parameters

Global Parameters controls:

- Simulation name.
- Start and end time.
- Output sample rate.
- Internal simulation sampling rate.
- Propagation speed.
- Random seed.
- ADC bits.
- Oversampling ratio.
- Rotation angle unit.
- Origin latitude, longitude, and altitude.
- Coordinate system for KML/geospatial export: ENU, UTM, or ECEF.

The coordinate system setting affects KML/geospatial export. Platform `x`, `y`, and `altitude` remain meter-valued coordinates and are used directly by the preview and simulation. The UI does not use latitude/longitude as platform waypoints.

If you change the rotation angle unit after entering rotation values, the UI asks whether to convert existing values or keep the numbers unchanged.

### Waveforms

Waveform editing supports:

- Pulse file waveforms from `.csv` or `.h5`.
- CW waveforms.
- FMCW linear chirps.
- FMCW triangular chirps.

For FMCW waveforms, the UI warns or blocks when settings violate major FMCW constraints, such as chirp period shorter than chirp duration or sweep bandwidth too close to the effective sample-rate limit.

### Timings

Timing editing supports:

- Frequency.
- Frequency offset.
- Random frequency offset standard deviation.
- Phase offset.
- Random phase offset standard deviation.
- Noise entries with alpha and weight values.

For simple scenarios, a default timing source is often enough.

### Antennas

Antenna editing supports:

- Isotropic.
- Sinc.
- Gaussian.
- Square horn.
- Parabolic.
- XML pattern file.
- HDF5 pattern file.

The Properties panel includes `Mesh Scale Multiplier`, which affects only the size of the antenna-pattern preview mesh in the 3D viewport.

For square horn and parabolic antennas, `Design Frequency` helps the UI preview a frequency-dependent pattern when the antenna is not attached to a waveform. It is a preview aid, not a replacement for the waveform carrier frequency used in simulation.

### Platforms

Platform editing supports:

- Static, linear, or cubic motion paths.
- Position waypoints with `x`, `y`, `altitude`, and `time`.
- Fixed-rate rotation.
- Waypoint-based rotation.
- Rotation interpolation using static, linear, or cubic mode.

For static motion or rotation, the UI uses the first waypoint. Use linear or cubic interpolation when a platform moves or changes pointing through the simulation.

### Radar Components

Monostatic radars, transmitters, and receivers share common fields:

- Component name.
- Radar mode: pulsed, CW, or FMCW.
- Compatible waveform selection.
- Antenna selection.
- Timing source selection.
- Optional operating schedule.

Pulsed radar components also expose PRF. Pulsed receivers and monostatic radars expose receiver window skip and window length.

Receivers and monostatic radars also expose:

- Noise temperature.
- Ignore direct paths.
- Ignore propagation loss.

Use `Ignore Propagation Loss` only for debugging and early setup checks. It is not a realistic link-budget setting.

### FMCW Receiver Settings

For FMCW receivers and monostatic radars, the UI supports:

- Dechirp mode: none, physical, or ideal.
- Dechirp reference: attached, transmitter, or custom waveform.
- Optional IF sample rate.
- Optional IF filter bandwidth.
- Optional IF transition width.

The UI blocks simulations and KML generation when it detects invalid FMCW settings that would fail later.

### Targets

Target editing supports:

- Isotropic RCS with a constant value in square meters.
- File-backed RCS.
- RCS model: constant, chi-square, or gamma.
- `k` value for chi-square and gamma models.

The 3D viewport currently visualizes simple isotropic RCS as a wireframe sphere. It does not visualize every file-backed RCS pattern.

## Importing And Exporting Scenarios

Import supports `.xml` and `.fersxml` files. Importing a scenario replaces the current in-memory scenario. If the current scenario has unsaved changes, the UI asks for confirmation first.

Export writes the current scenario as `.fersxml`. The suggested filename is based on the simulation name.

Important notes:

- Imported scenarios can produce warnings. Review them before running.
- Exported XML should still be treated as the authoritative scenario file.
- Advanced XML features that are not exposed in the UI may be simplified, omitted, or require manual review.

## Asset Library

The Asset Library stores reusable scenario pieces:

- Waveforms.
- Timings.
- Antennas.
- Platforms.

To save an asset, select it in Scenario Builder and use `Save to Asset Library` in the Properties panel.

To load an asset, open Asset Library and press `Load` on the saved item.

Platform assets include their dependencies, such as referenced waveforms, timings, and antennas. When loaded into another scenario, the UI creates renamed copies as needed to avoid name conflicts.

Asset library files can be imported and exported as JSON. Exported files normally use:

```text
fers-assets.fersasset.json
```

The local catalog is stored by the app as `asset-library.json` in the application data location for your operating system. When desktop file access is not available, the UI falls back to browser local storage.

## Running Simulations In The UI

Open the Simulation Run view to run the current scenario.

1. Choose an output directory.
2. Press `Run Simulation`.
3. Watch progress for the main run and output export tasks.
4. Review the output metadata table after completion.
5. Open the generated HDF5 files in your analysis tools.

The UI writes the same HDF5 result files as `fers-cli`.

The output metadata table shows:

- Receiver name.
- Radar mode.
- Sample rate.
- Total samples.
- Sample range.
- Pulse count or streaming segment summary.
- Output file path.
- FMCW details when available.

Use `Export JSON` to save the output metadata beside the generated HDF5 files.

## Generating KML

The Simulation Run view can generate KML without running the full signal simulation.

Use this to inspect:

- Platform motion paths.
- Geospatial placement.
- Antenna pointing.

If KML geometry looks wrong, check Global Parameters, especially origin and coordinate system.

## Logs And Diagnostics

Open Raw Logs from the left rail to inspect simulator output. You can change the FERS log level from the log drawer or Settings.

Recommended log levels:

| Level | Use |
| --- | --- |
| `INFO` | Normal runs. |
| `DEBUG` | Scenario troubleshooting. |
| `TRACE` | Very detailed diagnostics. Can produce a lot of output. |
| `WARNING`, `ERROR`, `FATAL` | Quiet runs where only problems should be shown. |
| `OFF` | Disable FERS logs in the UI. |

If a scenario fails to run, open Raw Logs and check the latest errors before editing the XML by hand.

## UI Versus XML

The UI is a scenario editor and runner. The `.fersxml` file is still the portable scenario format.

Use this rule of thumb:

- Use the UI to build and inspect scenarios visually.
- Export `.fersxml` when you want to archive, share, review, or batch-run a scenario.
- Use `fers-cli` for scripted runs and repeatable command-line workflows.
- Use Python or another analysis tool for detailed HDF5 signal analysis.
