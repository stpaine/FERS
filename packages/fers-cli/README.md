# FERS-CLI: Command-Line Interface for FERS

[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://github.com/stpaine/FERS/blob/master/LICENSE)

**FERS-CLI** is the command-line interface for the FERS simulation engine. It is a lightweight executable that acts as a
client to the core `libfers` library.

Its primary purpose is to provide a familiar, scriptable interface for users who prefer to run simulations from the
terminal and to maintain backward compatibility with the original FERS workflow.

## Features

- Full access to all core `libfers` simulation capabilities.
- Parses command-line arguments for simulation control.
- Loads scenarios from FERS XML files.
- Runs pulsed, CW, FMCW, and SFCW scenarios using the same XML workflow, including HDF5-backed CW/FMCW waveforms.
- Writes HDF5 receiver results by default or streams paced VITA 49.2 UDP output when selected at runtime.
- Displays real-time simulation progress in the console.
- Automatically outputs results to the scenario's directory by default.
- Generates KML files for scenario visualization.

## Building

The `fers-cli` executable is built as part of the main C++ build process for the monorepo. Please see the build
instructions in the [root `README.md`](https://github.com/stpaine/FERS/blob/master/README.md) file.

After a successful build, the executable can be found at:

- Linux/macOS: `build/release/packages/fers-cli/fers-cli`
- Windows: `build/release/packages/fers-cli/fers-cli.exe`

## Usage

Run the simulator from the command line, providing the path to a scenario XML file and any desired options.

### Linux / macOS:
```bash
./build/release/packages/fers-cli/fers-cli path/to/your/scenario.fersxml [options]
```

### Windows:
```powershell
.\build\release\packages\fers-cli\fers-cli.exe path\to\your\scenario.fersxml [options]
```

### Options

| Flag                                 | Description                                                                                                                          |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------ |
| `--help`, `-h`                       | Show the help message and exit.                                                                                                      |
| `--version`, `-v`                    | Show version information and exit.                                                                                                   |
| `--no-validate`                      | Disable the validation of the scenario file before running.                                                                          |
| `--kml[=<file>]`                     | Generate a KML visualization of the scenario and exit. If a filename is provided, it will be used. Otherwise, it defaults to the scenario name with a `.kml` extension in the output directory. |
| `--out-dir=<dir>`                    | Set the output directory for simulation results, VITA run metadata, and default KML output. Defaults to the directory containing the input scenario file. |
| `--vita49 <host:port>`               | Stream receiver output using the FERS VITA 49.2 UDP profile instead of writing HDF5 receiver files.                                  |
| `--vita49-fullscale <positive-real>` | Set the required fixed ADC full-scale used for VITA int16 IQ conversion.                                                             |
| `--vita49-epoch <unix-nanoseconds>`  | Set an optional fixed VITA UTC epoch; otherwise the epoch is selected when the first data batch is ready.                            |
| `--vita49-max-udp-payload <bytes>`   | Set the maximum VITA UDP datagram size. Valid range: `64..65507`; default: `1400`.                                                    |
| `--vita49-queue-depth <packets>`     | Set the paced sender's steady-state backpressure watermark. Must be greater than zero; default: `1024`.                              |
| `--log-level=<level>`                | Set the logging level (`TRACE`, `DEBUG`, `INFO`, `WARNING`, `ERROR`, `FATAL`).                                                       |
| `--log-file=<file>`                  | Log output to the specified `.log` or `.txt` file in addition to the console.                                                        |
| `-n=<threads>`                       | Set the number of threads to use for the simulation.                                                                                 |

### Example

```bash
./build/release/packages/fers-cli/fers-cli examples/mixed_scenario/example.fersxml --out-dir=./results --log-level=INFO -n=8
```

Stream receiver output as VITA 49.2 UDP:

```bash
./build/release/packages/fers-cli/fers-cli scenario.fersxml --out-dir=./results \
  --vita49 127.0.0.1:4991 --vita49-fullscale 1.0
```

VITA transport selection is a runtime option and is not stored in `.fersxml`. See the
[VITA 49.2 streaming implementation guide](../../docs/wiki/VITA49-Streaming-Implementation.md) for packet layout,
pacing, metadata, and counter semantics.
