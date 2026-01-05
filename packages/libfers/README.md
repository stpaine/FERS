# libfers: The FERS Core Simulation Library

[![Build Status](https://github.com/davidbits/FERS/actions/workflows/CMake.yml/badge.svg)](https://github.com/davidbits/FERS/actions/workflows/CMake.yml)
[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://github.com/davidbits/FERS/blob/master/LICENSE)

**libfers** is the core C++ simulation engine for the **Flexible Extensible Radar Simulator (FERS)**. It contains all
the logic for parsing scenarios, modeling physics, executing simulations, and processing output data.

This library is built with modern C++23 standards for optimal performance and maintainability. It exposes a stable
C-style Foreign Function Interface (FFI) in `include/libfers/api.h`, allowing it to be used by various frontends,
including the official `fers-cli` and `fers-ui` applications.

## Features

- **C-API:** A stable C-style API for creating and managing simulation contexts.
- **JSON & XML Serialization:** Robust serialization of the entire simulation state to and from JSON (for UI interoperability) and XML (for file-based workflows).
- **Signal-Level Modeling:** Creation of radar signal returns, including Doppler and phase modeling.
- **System Simulation:** Support for Monostatic, Multistatic, Continuous Wave (CW), and Pulsed radar systems.
- **Data Export:** Advanced data export in HDF5, CSV, and XML formats.
- **Geographic Visualization:** Generate KML files from scenarios.
- **Performance:** A unified event-driven architecture for efficient simulation of both pulsed and continuous-wave scenarios, with a global thread pool for parallelizing tasks.

## Dependencies

- A C++23 compatible compiler (e.g., GCC 11+, Clang 14+)
- **CMake** (3.22+)
- **libhdf5** (HDF5 support)
- **libxml2** (XML handling)
- **HighFive**, **GeographicLib**, & **nlohmann/json** (included as git submodules)
- **Doxygen** & **Graphviz** (for documentation generation, optional)
- **doxygen-awesome-css** (for enhanced Doxygen themes, optional)

## Building the Library and CLI

The C++ components are built using CMake from the repository root. For a complete development setup guide, please refer to the **[root `README.md`](../../README.md)**.

### 1. Install System Dependencies

**On Ubuntu/Debian:**

```bash
sudo apt-get update && sudo apt-get install build-essential cmake libhdf5-dev libxml2-dev xxd doxygen graphviz
```

**On macOS (using Homebrew):**
A modern LLVM toolchain is required, as the default system Clang may be outdated.

```bash
brew install cmake hdf5 libxml2 llvm doxygen graphviz
```

### 2. Configure and Build

From the root of the repository, create a build directory and run CMake/make.

**On Linux:**

```bash
# From the root FERS directory
mkdir build && cd build
cmake ..
make -j$(nproc)
```

**On macOS:**
You must point CMake to the Homebrew LLVM toolchain.

```bash
# From the FERS/ directory
mkdir build && cd build
CC=/usr/local/opt/llvm/bin/clang CXX=/usr/local/opt/llvm/bin/clang++ cmake ..
make -j$(sysctl -n hw.ncpu)
```

The compiled artifacts will be located in the `build/` directory. Libraries (`.a`, `.so`) are in `libfers/`,
and the `fers-cli` executable is in `fers-cli/`. DLLs on Windows will be placed in `bin/`.

### Build Options

You can customize the build using the following CMake options:

| Option                        | Description                                     | Default |
| ----------------------------- | ----------------------------------------------- | ------- |
| `-DCMAKE_BUILD_TYPE=Debug`    | Create a debug build with symbols.              | Release |
| `-DFERS_BUILD_SHARED_LIBS=ON` | Build the shared library (`.so`/`.dll`).        | ON      |
| `-DFERS_BUILD_STATIC_LIBS=ON` | Build the static library (`.a`).                | ON      |
| `-DFERS_BUILD_DOCS=ON`        | Enable the `doc` target for Doxygen generation. | OFF     |

Example of a debug build that only creates a static library:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Debug -DFERS_BUILD_SHARED_LIBS=OFF
```

### 3. Install (Optional)

You can install the libraries, headers, and CLI executable to your system.

```bash
# From the build directory
sudo make install
sudo ldconfig # On Linux, to update the cache
```

## Documentation

The source code documentation is automatically built and deployed to our [GitHub Pages site](https://davidbits.github.io/FERS/). If you wish to build it locally, you will need **Doxygen** and **Graphviz** installed. They are included in the dependency installation commands above.

1. **Configure with documentation enabled:**

    ```bash
    # From your build directory
    cmake .. -DFERS_BUILD_DOCS=ON
    ```

2. **Build the `doc` target:**
    ```bash
    make doc
    ```

The HTML output will be generated in the `build/docs/html/` directory. You can open `index.html` in your browser to view the documentation.
