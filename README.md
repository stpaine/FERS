# The Flexible Extensible Radar Simulator

[![Core CI](https://github.com/stpaine/FERS/actions/workflows/core.yml/badge.svg)](https://github.com/stpaine/FERS/actions/workflows/core.yml)
[![UI CI](https://github.com/stpaine/FERS/actions/workflows/ui.yml/badge.svg)](https://github.com/stpaine/FERS/actions/workflows/ui.yml)
[![Documentation](https://github.com/stpaine/FERS/actions/workflows/docs.yml/badge.svg)](https://github.com/stpaine/FERS/actions/workflows/docs.yml)
[![GitHub issues](https://img.shields.io/github/issues/stpaine/FERS.svg)](https://github.com/stpaine/FERS/issues)
[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://github.com/stpaine/FERS/blob/master/LICENSE)

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
- **Flexible System Modeling:** Simulate a wide range of radar systems, including monostatic, multistatic, pulsed,
  continuous wave (CW), and native FMCW streaming modes.
- **Advanced Data Export:** Output simulation data in HDF5 format for analysis.
- **Geographic Visualization:** Generate KML files from scenarios for accurate visualization in tools like Google Earth.
- **Modern Documentation:** A continuously updated and
  deployed [documentation site](https://davidbits.github.io/FERS/)
  with a searchable interface, generated directly from the source code.
- **Unified Schema:** A central XML schema ensures consistency and serves as the single source of truth for scenarios
  across the simulator and the UI.

## Feature Stability

Last reviewed: 2026-04-30.

These labels describe the current state of the `master` branch and should be read alongside
the [development status disclaimer](#disclaimer--development-status). "Verified" means the feature has been checked
against the current project expectations; "unverified" means it is available but still needs validation before users
should rely on it for important results.

| Feature | Stability | Verification |
| ------- | --------- | ------------ |
| fers-cli | Stable | N/A |
| fers-ui | Beta | N/A |
| Pulsed radar | Stable | Verified |
| CW radar | Beta | Unverified |
| FMCW radar | Alpha | Unverified |
| Multimodal | Alpha | Unverified |
| KML generation | Stable | Verified |

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

FERS supports native development on Linux, macOS, and Windows.

### Supported Development Targets

| Platform | Core build | UI build | Notes |
| -------- | ---------- | -------- | ----- |
| Linux x64 / ARM64 | Supported | Supported | Requires standard Tauri Linux system packages for the UI. |
| macOS Intel / Apple Silicon | Supported | Supported | `MACOSX_DEPLOYMENT_TARGET=14.0` is recommended locally. |
| Windows x64 | Supported | Supported | Use Visual Studio 2022 Build Tools (or later) and the MSVC Rust toolchain. |
| Windows ARM64 | Supported | Supported | Use native ARM64 MSVC tools and the `aarch64-pc-windows-msvc` Rust target. |

### 1. Prerequisites

Please install the following external dependencies using their official guides:

- A C++23 compatible compiler (GCC 11+, Clang 14+, or MSVC v143+)
- [**CMake**](https://cmake.org/download/) (3.22+) and [**Ninja**](https://ninja-build.org/)
- [**vcpkg**](https://vcpkg.io/en/getting-started.html) (C++ package manager)
- [**Bun**](https://bun.sh/) (JavaScript runtime and package manager)
- [**Rust**](https://www.rust-lang.org/tools/install) (Required for the Tauri UI)
- [**Tauri Prerequisites**](https://tauri.app/start/prerequisites/) (OS-specific system dependencies for the UI)

Then also review [Platform-Specific Notes](#platform-specific-notes) below.

### 2. Environment Configuration

**Important:** FERS relies on `vcpkg` to manage its C++ dependencies. You **must** set the `VCPKG_ROOT` environment variable to point to your vcpkg installation directory before building.

*Linux / macOS:*
```bash
export VCPKG_ROOT=/path/to/vcpkg
```

*Windows (PowerShell):*
```powershell
$env:VCPKG_ROOT = "C:\path\to\vcpkg"
```

> [!TIP]
> Installing Visual Studio should provide and link vcpkg if selected during installation, in which case the VCPKG_ROOT envvar is automatically set.

### 3. Clone and Install JS Dependencies

Clone the repository and install the frontend dependencies. This will also set up pre-commit hooks.

```bash
git clone https://github.com/stpaine/FERS.git
cd FERS
bun install
```

### 4. Build the Standalone C++ Core

You can configure and compile the C++ libraries and CLI directly using our CMake presets. This will automatically invoke `vcpkg` to fetch and build the required C++ dependencies (like HDF5 and libxml2).

```bash
cmake --preset=release
cmake --build --preset=release
```

Expected C++ artifacts:

- Linux/macOS CLI: `build/release/packages/fers-cli/fers-cli`
- Windows CLI: `build/release/packages/fers-cli/fers-cli.exe`

> [!TIP]
> You can install the built `fers-cli` release to your system using:
>
> **Linux / macOS:**
> ```bash
> sudo cmake --install build/release
> sudo ldconfig # For Linux only, to update the library cache
> ```
> This installs to `/usr/local` by default.
>
> **Windows:**
> Open an **Administrator** Developer PowerShell and run:
> ```powershell
> cmake --install build/release
> ```
> This installs to `C:\Program Files (x86)\FERS` by default.

### 5. Run the UI

The UI build process is completely self-contained. When you run the UI, Cargo will automatically invoke CMake to build the C++ backend in an isolated directory and link it to the Tauri application.

To start the development server:
```bash
bun ui:dev
```

To build a release UI bundle:

```bash
bun ui:build
```

> [!WARNING]
> On some Intel macOS systems, WebGL may be unavailable at startup due to a system WebKit limitation. FERS detects this and disables the 3D viewport gracefully. See issue [#181](https://github.com/davidbits/FERS/issues/181) for details.

### 6. Run C++ Unit Tests (optional)

Use the `coverage` preset to compile and run the Catch2 unit tests. This preset enables `FERS_BUILD_TESTS` across all platforms.

```bash
cmake --preset=coverage
cmake --build --preset=coverage --parallel
ctest --preset=coverage --output-on-failure
```

## Using Old XML Scenarios

The new FERS uses a different XML schema for scenarios than the original version. If you have existing scenarios in the old format, you can convert them to the new format using the included migration script:

```bash
python3 scripts/migrate_fers_xml.py old_scenario.fersxml new_scenario.fersxml
```

## Platform-Specific Notes

- **Windows:** FERS requires native MSVC (Visual Studio 2022 Build Tools, or later). MinGW and WSL are not officially supported. You should use the **Developer PowerShell for VS** when running build commands so that `cl.exe` and the Windows SDK are correctly prioritized in your `PATH`. Ensure you install the MSVC Rust toolchain.
- **macOS:** It is highly recommended to set `MACOSX_DEPLOYMENT_TARGET=14.0` in your environment to ensure modern C++ filesystem support.

> [!IMPORTANT]
> **Windows: Ensuring the x64 Toolchain**
>
> The default Developer PowerShell in Visual Studio may launch with the **x86** host toolchain, even on a 64-bit
> system. This causes CMake to detect a 32-bit compiler, which will fail when linking against the 64-bit vcpkg
> packages (you will see an error like `HighFiveConfig.cmake, version: 3.3.0 (64bit) ... not compatible`).
>
> To fix this, you must ensure the x64 toolchain is active. Before running any build commands, execute:
> ```powershell
> & "C:\Program Files\Microsoft Visual Studio\<VS_VERSION>\<EDITION>\Common7\Tools\Launch-VsDevShell.ps1" -Arch amd64 -HostArch amd64
> ```
> Replace `<VS_VERSION>` and `<EDITION>` with your installation values (e.g., `18` and `Community`).
>
> You can verify the correct compiler is active by running:
> ```powershell
> (Get-Command cl.exe).Source  # Should contain Hostx64\x64
> ```
>
> **Note:** If you get a script execution policy error, run
> `Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned` first.
>
> To make this permanent, edit the Developer PowerShell profile in **Tools → Options → Environment → Terminal
> Profiles** and set the arguments to:
> ```
> -NoExit -Command "& { & 'C:\Program Files\Microsoft Visual Studio\<VS_VERSION>\<EDITION>\Common7\Tools\Launch-VsDevShell.ps1' -Arch amd64 -HostArch amd64 }"
> ```

## Contributing

We welcome contributions to the FERS project! Please read our [CONTRIBUTING.md](CONTRIBUTING.md) guide to get started.

Note that this repository uses **Husky** to enforce code quality with pre-commit hooks. When you commit, your staged files will be automatically formatted and linted.

## Versioning

FERS uses a single repo-wide semantic version across the C++ core, CLI, UI, and release metadata.

- `version.txt` is the single source of truth for the managed release version.
- `CMakeLists.txt`, `vcpkg.json`, `packages/fers-ui/package.json`, and `packages/fers-ui/src-tauri/Cargo.toml` are kept aligned with that shared version.
- `packages/fers-ui/src-tauri/tauri.conf.json` reads the UI package version indirectly via `../package.json`.

## License

- Copyright (C) 2006-2008 Marc Brooker and Michael Inggs.
- Copyright (C) 2008-present FERS contributors (see [AUTHORS.md](AUTHORS.md)).

This program is free software; you can redistribute it and/or modify it under the terms of the
[GNU General Public License](https://github.com/stpaine/FERS/blob/master/LICENSE) as published by the Free
Software Foundation; version 2 of the License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
the [GNU General Public License](https://github.com/stpaine/FERS/blob/master/LICENSE) for
more details.

You should have received a copy of
the [GNU General Public License](https://github.com/stpaine/FERS/blob/master/LICENSE) along with this program;
if not, write to
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

### Third-Party Libraries

FERS incorporates code from several third-party libraries provided under their own licenses (MIT, BSD, Boost, MPL-2.0). License texts are in the `THIRD_PARTY_LICENSES` directory.

The JS and Rust license files are auto-generated and must be kept in sync with the dependency lock files:

```bash
bun run licenses:js    # regenerate THIRD_PARTY_LICENSES/js-licenses.txt
bun run licenses:rust  # regenerate THIRD_PARTY_LICENSES/rust-licenses.html
```

Run these after updating `packages/fers-ui/package.json` or `packages/fers-ui/src-tauri/Cargo.toml`. CI will fail if they are stale.

## Disclaimer & Development Status

> [!WARNING]
> Please be aware that FERS is currently undergoing a significant modernization and re-architecture. The `master` branch is under **heavy active development** and should be considered an **alpha-stage** project.

This means:

- **Stability:** Expect bugs, crashes, and incomplete features.
- **Breaking Changes:** The C-API, JSON/XML schemas, and internal architecture are subject to change without notice as the new foundation is stabilized.
- **Use Case:** This version is intended for development, testing, and community feedback. It is **not yet recommended for production or critical simulation work.**
