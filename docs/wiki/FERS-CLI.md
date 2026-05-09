# FERS CLI

`fers-cli` is the normal way to run FERS scenarios from a terminal. It loads a `.fersxml` file, validates it unless told not to, and either runs the simulation or exports KML.

## Basic Command

```bash
./build/release/packages/fers-cli/fers-cli scenario.fersxml --out-dir=./results
```

Windows:

```powershell
.\build\release\packages\fers-cli\fers-cli.exe scenario.fersxml --out-dir=.\results
```

## Options

| Option | What it does |
| --- | --- |
| `--help`, `-h` | Show help. |
| `--version`, `-v` | Show the FERS version. |
| `--out-dir=<dir>` | Choose where HDF5 result files are written. |
| `--kml` | Export KML and do not run the simulation. |
| `--kml=<file>` | Export KML to a specific file and do not run the simulation. |
| `--no-validate` | Skip XML schema validation before loading. |
| `--log-level=<level>` | Set logging detail. Use `TRACE`, `DEBUG`, `INFO`, `WARNING`, `ERROR`, or `FATAL`. |
| `--log-file=<file>` | Write logs to a `.log` or `.txt` file as well as the terminal. |
| `-n=<threads>` | Choose how many worker threads to use. |

If `--out-dir` is not supplied, results are written beside the scenario file.

## Run Examples

Run the mixed pulsed/CW example:

```bash
cd examples/mixed_scenario
python3 genpulse.py
../../build/release/packages/fers-cli/fers-cli example.fersxml --out-dir=.
python3 analysis.py
```

Run an FMCW example into a separate output folder:

```bash
./build/release/packages/fers-cli/fers-cli examples/fmcw_monostatic_dechirp/example.fersxml --out-dir=./results --log-level=INFO -n=4
```

Export KML for a scenario:

```bash
./build/release/packages/fers-cli/fers-cli examples/fmcw_monostatic_dechirp/example.fersxml --out-dir=./results --kml
```

## Choosing Thread Count

Use `-n=<threads>` to control the worker pool:

```bash
./build/release/packages/fers-cli/fers-cli scenario.fersxml -n=8
```

Higher values can reduce runtime for larger scenarios, but they do not guarantee a speedup for small scenarios. Values above the machine's hardware concurrency are clamped.

## Logging

For normal use:

```bash
--log-level=INFO
```

For diagnosing a scenario:

```bash
--log-level=DEBUG --log-file=run.log
```

Use `TRACE` only when you need very detailed logs. It can produce a large amount of output.

## Exit Status

`fers-cli` returns `0` when help/version is shown or a simulation succeeds. It returns `1` for argument errors, scenario loading errors, output-directory errors, and simulation failures.

Practical note: always check the log output for KML export errors. The current CLI can log a KML failure even when the process exit status is still `0`.

## Output Directory Behavior

`--out-dir` controls:

- HDF5 result files from simulation runs.
- The default location for KML when `--kml` is used without a file path.

Examples:

Writes results under `scenes/`.

```bash
fers-cli scenes/demo.fersxml
```

Writes results under `results/`.

```bash
fers-cli scenes/demo.fersxml --out-dir=results
```

Writes `results/demo.kml`.

```bash
fers-cli scenes/demo.fersxml --out-dir=results --kml
```

Writes `maps/demo.kml`.

```bash
fers-cli scenes/demo.fersxml --out-dir=results --kml=maps/demo.kml
```
