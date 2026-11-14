# Installing FERS

This page covers the normal build used by FERS users who want to run scenarios with `fers-cli`. The desktop UI has extra requirements, but the simulator and CLI are built with CMake, Ninja, and vcpkg.

## Requirements

Install:

- A C++23 compiler and standard library providing `std::format` and `std::expected`.
- CMake 3.22 or newer.
- Ninja.
- vcpkg.
- Bun, if you want to use the repository's JavaScript tooling or the desktop UI.
- Rust and the Tauri platform prerequisites, only if you want to build the desktop UI.

The current Linux Debug and Release builds have been tested with GCC 16.2.1 and CMake 4.5.0.

Set `VCPKG_ROOT` before configuring FERS:

```bash
export VCPKG_ROOT=/path/to/vcpkg
```

Windows PowerShell:

```powershell
$env:VCPKG_ROOT = "C:\path\to\vcpkg"
```

FERS can sometimes find vcpkg automatically, but setting `VCPKG_ROOT` avoids ambiguity.

## Get The Repository

```bash
git clone https://github.com/stpaine/FERS.git
cd FERS
```

If you plan to use the UI or repository tooling, install the JavaScript dependencies:

```bash
bun install
```

## Build The Simulator And CLI

From the repository root:

```bash
cmake --preset=release
cmake --build --preset=release
```

The CLI executable is created at:

| Platform | Path |
| --- | --- |
| Linux/macOS | `build/release/packages/fers-cli/fers-cli` |
| Windows | `build/release/packages/fers-cli/fers-cli.exe` |

Run a quick check:

```bash
./build/release/packages/fers-cli/fers-cli --version
```

On Windows:

```powershell
.\build\release\packages\fers-cli\fers-cli.exe --version
```

## Install FERS

You do not have to install FERS system-wide. You can run the CLI directly from the build directory.

If you do want to install it:

Linux:

```bash
sudo cmake --install build/release
sudo ldconfig
```

macOS:

```bash
sudo cmake --install build/release
```

Windows, from an Administrator Developer PowerShell:

```powershell
cmake --install build/release
```

The default install prefix is usually `/usr/local` on Linux/macOS and `C:\Program Files (x86)\FERS` on Windows.

## Run The Test Suite

Use this when you want to verify your local build:

```bash
cmake --preset=debug
cmake --build --preset=debug --parallel
ctest --preset=debug --output-on-failure

cmake --preset=release -DFERS_BUILD_TESTS=ON
cmake --build --preset=release --parallel
ctest --preset=release --output-on-failure
```

The `debug` preset enables tests. Release uses the explicit `FERS_BUILD_TESTS=ON` setting above. Use the `coverage` preset for instrumented builds and `scripts/run_tests_gen_coverage.sh` to generate an HTML coverage report.

## Build The Desktop UI

The UI is optional. It requires Bun, Rust, and Tauri prerequisites for your operating system. See [[FERS UI]] for the user workflow and current UI limitations.

Development mode:

```bash
bun ui:dev
```

Release bundle:

```bash
bun ui:build
```

Run these commands from the repository root.

## Platform Notes

### Windows

Use the MSVC toolchain. MinGW is not supported by the top-level build.

Use a Visual Studio Developer PowerShell configured for the architecture you are building:

| Native build target | Compiler path suffix | vcpkg target triplet | Rust target for the UI |
| --- | --- | --- | --- |
| Windows x64 | `Hostx64\x64\cl.exe` | `x64-windows-static-md` | `x86_64-pc-windows-msvc` |
| Windows ARM64 | `HostARM64\ARM64\cl.exe` | `arm64-windows-static-md` | `aarch64-pc-windows-msvc` |

For an x64 build on an x64 host, select the x64 tools explicitly if needed:

```powershell
& "C:\Program Files\Microsoft Visual Studio\<VS_VERSION>\<EDITION>\Common7\Tools\Launch-VsDevShell.ps1" -Arch amd64 -HostArch amd64
```

For a native ARM64 build, use Visual Studio for ARM64 with the **Desktop development with C++** workload and its native ARM64 compiler environment. Microsoft documents the toolset and compiler paths in [Visual Studio on Arm-powered devices](https://learn.microsoft.com/en-us/visualstudio/install/visual-studio-on-arm-devices).

Confirm the selected compiler before configuring:

```powershell
(Get-Command cl.exe).Source
```

The path should match the native target in the table. For ARM64, configure the matching dependencies explicitly:

```powershell
cmake --preset=release -DVCPKG_TARGET_TRIPLET=arm64-windows-static-md
cmake --build --preset=release
```

### macOS

A modern compiler is required. Homebrew LLVM is often the easiest way to get current C++23 support. The project recommends `MACOSX_DEPLOYMENT_TARGET=14.0` for local builds.

### Linux

Install the normal build toolchain plus CMake, Ninja, and vcpkg. UI builds also need the Linux packages required by Tauri.

## Docker

Older FERS documentation described a legacy Docker image. That image is not part of the current repository and should not be treated as the supported install path. Use the native CMake/vcpkg build above.
