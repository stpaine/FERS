# The Flexible Extensible Radar Simulator

<!-- TODO: Build status should be present for both ui and core -->

[![FERS Core CI](https://github.com/davidbits/FERS/actions/workflows/CMake.yml/badge.svg)](https://github.com/davidbits/FERS/actions/workflows/CMake.yml)
[![Documentation](https://github.com/davidbits/FERS/actions/workflows/docs.yml/badge.svg)](https://github.com/davidbits/FERS/actions/workflows/docs.yml)
[![GitHub issues](https://img.shields.io/github/issues/davidbits/FERS.svg)](https://github.com/davidbits/FERS/issues)
[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://github.com/davidbits/FERS/blob/master/LICENSE)

![Tech Stack](https://img.shields.io/badge/Core-C%2B%2B%2023-00599C?logo=cplusplus)
![Tech Stack](https://img.shields.io/badge/UI-Tauri%20%7C%20React-20232A?logo=react)

FERS is a comprehensive suite of tools for signal-level radar simulation. It consists of a high-performance C++
simulation engine and a modern, intuitive graphical user interface for building and visualizing complex scenarios.

This repository is structured as a monorepo, containing the core simulator, the UI, and the shared data schema as
semi-independent packages.

---

## Key Features

- **Unified Event-Driven Core:** A modernized C++23 engine featuring a unified event-driven architecture for efficient
  simulation of both pulsed and continuous-wave scenarios, with optimized multithreading.
- **Visual Scenario Builder:** An intuitive 3D interface to construct, configure, and visualize radar scenarios.
- **Flexible System Modeling:** Simulate a wide range of radar systems, including monostatic, multistatic, pulsed, and
  continuous wave (CW).
- **Advanced Data Export:** Output simulation data in HDF5, CSV, and XML formats for analysis.
- **Geographic Visualization:** Generate KML files from scenarios for accurate visualization in tools like Google Earth.
- **Modern Documentation:** A continuously updated and
  deployed [documentation site](https://davidbits.github.io/FERS/)
  with a searchable interface, generated directly from the source code.
- **Unified Schema:** A central XML schema ensures consistency and serves as the single source of truth for scenarios
  across the simulator and the UI.

## Repository Structure

This monorepo contains the following packages:

- `packages/libfers`: The core C++ radar simulation library. It contains all the core logic, physics, and file
  processing capabilities, exposed through a C-style API.
- `packages/fers-cli`: A lightweight command-line interface that acts as a wrapper around `libfers`, providing backward
  compatibility with the original FERS executable.
- `packages/fers-ui`: A modern desktop application built with Tauri and React. It provides a graphical interface for
  creating, editing, and visualizing FERS scenarios by linking against `libfers`.
- `packages/schema`: The XML Schema Definition (XSD) and Document Type Definition (DTD) that define the structure of
  FERS scenario files. This schema is the contract between the UI and the core simulator.

## Development Setup

Follow these steps to set up a development environment for building the C++ core and running the UI.

### 1. Prerequisites

Ensure you have the following tools installed on your system:

- A C++23 compatible compiler (e.g., GCC 11+, Clang 14+) and **CMake** (3.22+).
- [**Node.js**](https://nodejs.org/) (v18+) and [**pnpm**](https://pnpm.io/).
- The [**Rust toolchain**](https://www.rust-lang.org/tools/install).
- [**Tauri prerequisites**](https://tauri.app/start/prerequisites/) for your operating system.
- [**clang-format**](https://clang.llvm.org/docs/ClangFormat.html) (for code formatting).

### 2. Clone the Repository

Clone the repository and its submodules from the root of the monorepo.

```bash
git clone --recursive https://github.com/davidbits/FERS.git
cd FERS
```

### 3. Install System Dependencies (for C++ Core)

**On Ubuntu/Debian:**

```bash
sudo apt-get update && sudo apt-get install build-essential cmake libhdf5-dev libxml2-dev xxd
```

**On macOS (using Homebrew):**

```bash
brew install cmake hdf5 libxml2 llvm
```

### 4. Install Node.js Dependencies

From the **root of the repository**, install all JavaScript dependencies. This also sets up pre-commit hooks using
Husky.

```bash
pnpm install
```

### 5. Build the C++ Core

Create a build directory and compile the C++ libraries using CMake.

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

### 6. Run the UI

Navigate back to the root of the repository and start the development server.

```bash
# From the root FERS directory
pnpm ui:dev
```

## Contributing

We welcome contributions to the FERS project! Please read our [CONTRIBUTING.md](CONTRIBUTING.md) guide to get started.

Note that this repository uses **Husky** to enforce code quality with pre-commit hooks. When you commit, your staged
files will be automatically formatted and linted. Ensure you have `clang-format`, `prettier`, and the Rust toolchain
installed.

## License

- Copyright (C) 2006-2008 Marc Brooker and Michael Inggs.
- Copyright (C) 2008-present FERS contributors (see [AUTHORS.md](AUTHORS.md)).

This program is free software; you can redistribute it and/or modify it under the terms of the
[GNU General Public License](https://github.com/davidbits/FERS/blob/master/LICENSE) as published by the Free
Software Foundation; version 2 of the License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
the [GNU General Public License](https://github.com/davidbits/FERS/blob/master/LICENSE) for
more details.

You should have received a copy of
the [GNU General Public License](https://github.com/davidbits/FERS/blob/master/LICENSE) along with this program;
if not, write to
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

### User-Generated Files

Please note that this license only covers the source code, program binaries, and build system of FERS. Any input files
you create (such as simulation scenarios) and any results generated by the simulator are not covered by this license and
remain the copyright of their original author.

### Third-Party Libraries

FERS incorporates code from the following third-party libraries, which are provided under their own licenses. The full
text for these licenses can be found in the `THIRD_PARTY_LICENSES` directory.

- **libxml2:** Used for XML parsing. Licensed under the [MIT License](THIRD_PARTY_LICENSES/libxml2-LICENSE.txt).
- **HighFive:** A C++ header-only library for HDF5. Licensed under
  the [Boost Software License 1.0](THIRD_PARTY_LICENSES/HighFive-LICENSE.txt).
- **GeographicLib:** A library for geographic calculations. Licensed under
  the [MIT License](THIRD_PARTY_LICENSES/GeographicLib-LICENSE.txt).
- **libhdf5:** Used for HDF5 file handling. Licensed under
  the [BSD 3-Clause License](THIRD_PARTY_LICENSES/libhdf5-LICENSE.txt).
- **nlohmann/json:** A JSON library for C++. Licensed under the
  [MIT License](THIRD_PARTY_LICENSES/nlohmann-json-LICENSE.txt).
- **Tauri:** A framework for building desktop applications. Licensed under the
  [MIT License](THIRD_PARTY_LICENSES/tauri-LICENSE.txt).
- **React:** A JavaScript library for building user interfaces. Licensed under the
  [MIT License](THIRD_PARTY_LICENSES/react-LICENSE.txt).
- **MUI:** A React component library. Licensed under the [MIT License](THIRD_PARTY_LICENSES/mui-LICENSE.txt).
- **Three.js:** A 3D JavaScript library. Licensed under the [MIT License](THIRD_PARTY_LICENSES/threejs-LICENSE.txt).
- **Zustand:** A small, fast state-management library for React. Licensed under the
  [MIT License](THIRD_PARTY_LICENSES/zustand-LICENSE.txt).
- **Zod:** A TypeScript-first schema declaration and validation library. Licensed under the
  [MIT License](THIRD_PARTY_LICENSES/zod-LICENSE.txt).
- **React Three Fiber:** A React renderer for Three.js. Licensed under the
  [MIT License](THIRD_PARTY_LICENSES/react-three-fiber-LICENSE.txt).

### Historical Notice from Original Distribution

The following notice was part of the original FERS distribution:

> Should you wish to acquire a copy of FERS not covered by these terms, please contact the Department of
> Electrical Engineering at the University of Cape Town.
