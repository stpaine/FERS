# FERS-CLI: Command-Line Interface for FERS

[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://github.com/davidbits/FERS/blob/master/LICENSE)

**FERS-CLI** is the command-line interface for the FERS simulation engine. It is a lightweight executable that acts as a
client to the core `libfers` library.

Its primary purpose is to provide a familiar, scriptable interface for users who prefer to run simulations from the
terminal and to maintain backward compatibility with the original FERS workflow.

## Features

- Full access to all core `libfers` simulation capabilities.
- Parses command-line arguments for simulation control.
- Loads scenarios from FERS XML files.
- Displays real-time simulation progress in the console.
- Generates KML files for scenario visualization.

## Building

The `fers-cli` executable is built as part of the main C++ build process for the monorepo. Please see the build
instructions in the [`packages/libfers/README.md`](https://github.com/davidbits/FERS/blob/master/packages/libfers/README.md) file.

After a successful build, the executable can be found at `build/packages/fers-cli/fers-cli`.

## Usage

Run the simulator from the command line, providing the path to a scenario XML file and any desired options.

```bash
# From the build/packages/ directory
./fers-cli/fers-cli path/to/your/scenario.fersxml [options]
```

### Options

| Flag                  | Description                                                                                                                          |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| `--help`, `-h`        | Show the help message and exit.                                                                                                      |
| `--version`, `-v`     | Show version information and exit.                                                                                                   |
| `--no-validate`       | Disable the validation of the scenario file before running.                                                                          |
| `--kml`               | Generate a KML visualization of the scenario and exit. The output will have the same name as the input file with a `.kml` extension. |
| `--log-level=<level>` | Set the logging level (`TRACE`, `DEBUG`, `INFO`, `WARNING`, `ERROR`, `FATAL`).                                                       |
| `--log-file=<file>`   | Log output to the specified `.log` or `.txt` file in addition to the console.                                                        |
| `-n=<threads>`        | Set the number of threads to use for the simulation.                                                                                 |

### Example

```bash
# From the build/packages directory
./fers-cli/fers-cli example.fersxml --log-level=DEBUG -n=8
```
