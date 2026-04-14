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

| Flag                  | Description                                                                                                                          |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| `--help`, `-h`        | Show the help message and exit.                                                                                                      |
| `--version`, `-v`     | Show version information and exit.                                                                                                   |
| `--no-validate`       | Disable the validation of the scenario file before running.                                                                          |
| `--kml[=<file>]`      | Generate a KML visualization of the scenario and exit. If a filename is provided, it will be used. Otherwise, it defaults to the scenario name with a `.kml` extension in the output directory. |
| `--out-dir=<dir>`     | Set the output directory for simulation results and default KML output. Defaults to the directory containing the input scenario file. |
| `--log-level=<level>` | Set the logging level (`TRACE`, `DEBUG`, `INFO`, `WARNING`, `ERROR`, `FATAL`).                                                       |
| `--log-file=<file>`   | Log output to the specified `.log` or `.txt` file in addition to the console.                                                        |
| `-n=<threads>`        | Set the number of threads to use for the simulation.                                                                                 |

### Example

```bash
./build/release/packages/fers-cli/fers-cli examples/mixed_scenario/example.fersxml --out-dir=./results --log-level=INFO -n=8
```
