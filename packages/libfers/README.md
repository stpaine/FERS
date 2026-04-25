# libfers: The FERS Core Simulation Library

[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://github.com/stpaine/FERS/blob/master/LICENSE)

**libfers** is the core C++ simulation engine for the **Flexible Extensible Radar Simulator (FERS)**. It contains all
the logic for parsing scenarios, modeling physics, executing simulations, and processing output data.

This library is built with modern C++23 standards for optimal performance and maintainability. It exposes a stable
C-style Foreign Function Interface (FFI) in `include/libfers/api.h`, allowing it to be used by various frontends,
including the official `fers-cli` and `fers-ui` applications.

## Features

- **C-API:** A stable C-style API for creating and managing simulation contexts.
- **JSON & XML Serialization:** Robust serialization of the entire simulation state to and from JSON (for UI interoperability) and XML (for file-based workflows).
- **Signal-Level Modeling:** Creation of radar signal returns, including Doppler and phase modeling.
- **System Simulation:** Support for monostatic, multistatic, continuous wave (CW), pulsed, and native FMCW streaming
  radar systems.
- **Data Export:** Advanced data export in HDF5, CSV, and XML formats.
- **Geographic Visualization:** Generate KML files from scenarios.
- **Performance:** A unified event-driven architecture for efficient simulation of both pulsed and continuous-wave scenarios, with a global thread pool for parallelizing tasks.

## Building the Library

The C++ components are built using CMake from the repository root. For a complete list of prerequisites and environment setup instructions (including setting `VCPKG_ROOT`), please refer to the [root `README.md`](https://github.com/stpaine/FERS/blob/master/README.md).

### 1. Configure and Build

From the root of the repository, use the CMake presets to configure and build the project. This will automatically use `vcpkg` to manage dependencies.

```bash
cmake --preset=release
cmake --build --preset=release
```

The compiled artifacts will be located in the `build/release/` directory. Libraries (`.a`, `.lib`, `.so`, `.dll`) are in `build/release/packages/libfers/` or `build/release/bin/`.

**Windows Note:** Native Windows builds use MSVC. Ensure you are running these commands from the **Developer PowerShell for VS**.

### Running Tests

The `release` preset builds the library and CLI but does not enable unit tests. Use the `coverage` preset when you need to compile and run the Catch2 test suite:

```bash
cmake --fresh --preset=coverage
cmake --build --preset=coverage --parallel
ctest --preset=coverage --output-on-failure
```

### Build Options

You can customize the build by passing CMake options alongside the preset:

| Option                        | Description                                     | Default |
| ----------------------------- | ----------------------------------------------- | ------- |
| `-DFERS_BUILD_SHARED_LIBS=ON` | Build the shared library (`.so`/`.dll`).        | ON      |
| `-DFERS_BUILD_STATIC_LIBS=ON` | Build the static library (`.a`/`.lib`).         | ON      |
| `-DFERS_BUILD_DOCS=ON`        | Enable the `doc` target for Doxygen generation. | OFF     |

### 2. Install (Optional)

You can install the libraries, headers, and CLI executable to your system.

#### Linux / macOS:
```bash
sudo cmake --install build/release
sudo ldconfig # On Linux, to update the library cache
```

#### Windows:
Open an **Administrator** Developer PowerShell and run:
```powershell
cmake --install build/release
```

## Documentation

The source code documentation is automatically built and deployed to our [GitHub Pages site](https://davidbits.github.io/FERS/). If you wish to build it locally, you will need **Doxygen** and **Graphviz** installed.

### Configure with documentation enabled:

```bash
cmake --preset release -DFERS_BUILD_DOCS=ON
```

### Build the `doc` target:
```bash
cmake --build --preset release --target doc
```

The HTML output will be generated in the `build/release/docs/html/` directory. You can open `index.html` in your browser to view the documentation.
