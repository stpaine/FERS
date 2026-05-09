# FERS Documentation

FERS is the Flexible Extensible Radar Simulator. It simulates radar signals at the I/Q sample level so you can build a radar scenario, run it, and analyze the generated receiver data.

FERS is useful when you need to study how radar waveforms, platform motion, timing, antenna gain, propagation loss, target RCS, and receiver configuration affect the signal that appears at a receiver.

## What FERS Can Simulate

FERS can model:

- Pulsed radar using complex baseband waveform files.
- Continuous-wave radar.
- FMCW linear chirps and triangular chirps.
- Monostatic radars, where the transmitter and receiver are attached to the same platform.
- Bistatic or multistatic layouts, where transmitters, receivers, and targets are on separate moving platforms.
- Multiple transmitters, receivers, radars, and targets in the same scenario.
- Platform motion with static, linear, or cubic interpolation.
- Platform pointing with fixed-rate rotation or waypoint-based rotation.
- Isotropic, analytic, XML, and file-backed antenna patterns.
- Constant, file-backed, and fluctuating target RCS.
- Receiver thermal noise.
- Timing frequency offset, phase offset, random offset, and phase-noise entries.
- Optional schedules that turn transmitters and receivers on only during selected time periods.
- KML scene export for geospatial inspection.
- HDF5 result output for Python, MATLAB, Julia, C++, or other analysis tools.

## What FERS Does Not Do

FERS is not a full electromagnetic solver. It does not calculate detailed scattering from CAD geometry, terrain, buildings, sea clutter, weather, or material models.

Current limitations to keep in mind:

- Multipath from terrain, ground bounce, buildings, and reflectors is not part of the current XML model.
- Clutter and environmental propagation effects are not modeled as first-class scenario objects.
- Target shape is represented through RCS values or RCS pattern files, not physical meshes.
- Built-in FMCW support is available but still marked alpha in the project README.
- CW support is available but still marked beta in the project README.
- The XML schema validates structure, but some physical limits are checked only when FERS loads or runs the scenario.

## Common Workflow

1. Install and build FERS.
2. Write a `.fersxml` scenario, or build one visually with `fers-ui`.
3. Run the scenario with `fers-cli` or the Simulation Run view in `fers-ui`.
4. Inspect the generated HDF5 result files.
5. Analyze the I/Q data with your own scripts.
6. Optionally export KML to inspect platform geometry.

```bash
./build/release/packages/fers-cli/fers-cli examples/fmcw_monostatic_dechirp/example.fersxml --out-dir=./results
```

## Main Documentation Pages

- [[Installing FERS]]: install prerequisites, build the CLI, build tests, and install locally.
- [[Using FERS]]: the normal end-to-end workflow from scenario file to HDF5 output.
- [[FERS CLI]]: command-line options and examples.
- [[FERS UI]]: visual scenario building, asset reuse, simulation runs, KML generation, and UI limitations.
- [[libfers]]: when and how to use the FERS library from another application.
- [[XML Schema Reference]]: every XML element, attribute, parameter, unit, and practical caveat.
- [[Examples]]: walkthroughs of the scenarios in the `examples/` directory.
- [[Simulation Pipelines]]: Mermaid diagrams showing how FERS runs pulsed, CW, and FMCW simulations.

## Migrating Old Scenarios

Older FERS XML used names such as `<pulse>`, `pulse="..."`, and inline radar mode attributes. Current FERS uses `<waveform>`, `waveform="..."`, and explicit mode blocks such as `<pulsed_mode>`, `<cw_mode>`, and `<fmcw_mode>`.

Use the migration script as a starting point:

```bash
python3 scripts/migrate_fers_xml.py old_scenario.fersxml new_scenario.fersxml
```

Always inspect migrated files before using them. The migrator removes deprecated `<export>` settings and unsupported legacy `<multipath>` blocks because they do not map to the current scenario model.
