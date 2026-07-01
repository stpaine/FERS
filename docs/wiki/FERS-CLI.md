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
| `--vita49 host:port` | Stream receiver output as the FERS VITA 49.2 UDP profile. |
| `--vita49-fullscale <positive-real>` | Required with `--vita49`; fixed ADC full-scale for int16 IQ scaling. |
| `--vita49-epoch <unix-nanoseconds>` | Optional deterministic VITA stream epoch for replay. |
| `--vita49-max-udp-payload <bytes>` | Optional VITA UDP payload cap, `64..65507` bytes. Default is `1400`. |
| `--vita49-queue-depth <packets>` | Optional VITA sender queue depth. Must be greater than zero. Default is `1024`. |
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

Run an SFCW example:

```bash
./build/release/packages/fers-cli/fers-cli examples/sfcw_monostatic/example.fersxml --out-dir=. --log-level=INFO
```

Export KML for a scenario:

```bash
./build/release/packages/fers-cli/fers-cli examples/fmcw_monostatic_dechirp/example.fersxml --out-dir=./results --kml
```

Stream receiver output as VITA 49.2 UDP:

```bash
./build/release/packages/fers-cli/fers-cli scenario.fersxml --out-dir=./results --vita49 127.0.0.1:4991 --vita49-fullscale 1.0
```

Use a deterministic replay epoch:

```bash
./build/release/packages/fers-cli/fers-cli scenario.fersxml --vita49 127.0.0.1:4991 --vita49-fullscale 1.0 --vita49-epoch 1700000000123456789
```

Tune VITA UDP packet sizing and queue depth:

```bash
./build/release/packages/fers-cli/fers-cli scenario.fersxml --vita49 127.0.0.1:4991 --vita49-fullscale 1.0 --vita49-max-udp-payload 1400 --vita49-queue-depth 1024
```

For the packet format, metadata fields, counters, and source-backed implementation contract, see [[VITA49 Streaming Implementation]].

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

`fers-cli` returns `0` when help/version is shown, a simulation succeeds, or KML export succeeds.

It returns `1` for argument errors, scenario loading errors, output-directory errors, simulation failures, and KML export failures. This means scripts can treat a nonzero exit status as a failed run or failed KML export.

## Output Directory Behavior

`--out-dir` controls:

- HDF5 result files from simulation runs.
- Run metadata, logs, and default KML paths for VITA 49.2 runs.
- The default location for KML when `--kml` is used without a file path.

`.fersxml` files remain scenario descriptions. They do not select VITA/network output; use the runtime switches above.

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
