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
- **System Simulation:** Support for Monostatic, Multistatic, Continuous Wave (CW), and Pulsed radar systems.
- **Data Export:** Advanced data export in HDF5, CSV, and XML formats.
- **Geographic Visualization:** Generate KML files from scenarios.
- **Performance:** A unified event-driven architecture for efficient simulation of both pulsed and continuous-wave scenarios, with a global thread pool for parallelizing tasks.

## Dependencies

- A C++23 compatible compiler (e.g., GCC 11+, Clang 14+)
- **CMake** (3.22+)
- **vcpkg** for C++ package management. All third-party libraries (like HDF5, libxml2, etc.) are managed through the `vcpkg.json` manifest at the repository root.
- **Doxygen** & **Graphviz** (for documentation generation, optional)

## Building the Library and CLI

The C++ components are built using CMake from the repository root. For a complete development setup guide, please refer to the **[root `README.md`](../../README.md)**.

### 1. Configure and Build

From the root of the repository, use the CMake presets to configure and build the project. This will automatically use `vcpkg` to manage dependencies. Ensure to install [vcpkg](https://learn.microsoft.com/en-us/vcpkg/get_started/get-started?pivots=shell-bash) before running the following commands.

```bash
# From the root FERS directory
cmake --preset=release
cmake --build --preset=release
```

The compiled artifacts will be located in the `build/release/` directory. Libraries (`.a`, `.so`) are in `build/release/packages/libfers/`,
and the `fers-cli` executable is in `build/release/packages/fers-cli/`. DLLs on Windows will be placed in `build/release/bin/`.

### Build Options

You can customize the build by passing CMake options alongside the preset:

| Option                        | Description                                     | Default |
| ----------------------------- | ----------------------------------------------- | ------- |
| `-DFERS_BUILD_SHARED_LIBS=ON` | Build the shared library (`.so`/`.dll`).        | ON      |
| `-DFERS_BUILD_STATIC_LIBS=ON` | Build the static library (`.a`).                | ON      |
| `-DFERS_BUILD_DOCS=ON`        | Enable the `doc` target for Doxygen generation. | OFF     |

Example of a debug build that only creates a static library:

```bash
# From the root FERS directory
cmake --preset debug -DFERS_BUILD_SHARED_LIBS=OFF
cmake --build --preset debug
```

### 2. Install (Optional)

You can install the libraries, headers, and CLI executable to your system.

```bash
# From the root FERS directory
sudo cmake --install build/release
sudo ldconfig # On Linux, to update the cache
```

## Documentation

The source code documentation is automatically built and deployed to our [GitHub Pages site](https://davidbits.github.io/FERS/). If you wish to build it locally, you will need **Doxygen** and **Graphviz** installed.

1. **Configure with documentation enabled:**

    ```bash
    # From the root FERS directory
    cmake --preset release -DFERS_BUILD_DOCS=ON
    ```

2. **Build the `doc` target:**
    ```bash
    cmake --build --preset release --target doc
    ```

The HTML output will be generated in the `build/release/docs/html/` directory. You can open `index.html` in your browser to view the documentation.
