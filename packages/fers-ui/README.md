# fers-ui: A Graphical Interface for the FERS Core Simulator

![Framework](https://img.shields.io/badge/Framework-React-61DAFB?logo=react)
![Language](https://img.shields.io/badge/Language-TypeScript-3178C6?logo=typescript)
![Desktop App](https://img.shields.io/badge/Tauri-v2-FFC336)
[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://github.com/stpaine/FERS/blob/master/LICENSE)

**fers-ui** is the official graphical user interface for the **Flexible Extensible Radar Simulator (FERS)**. It provides a professional-grade, visual workbench designed to streamline the entire simulation pipeline—from asset creation and scenario construction to simulation execution and results analysis. By offering an intuitive and powerful toolset, fers-ui dramatically improves the usability and accessibility of the core FERS engine.

## Key Features

- **Multi-Panel Scenario Builder:** A fully resizable workspace featuring an interactive 3D viewport for object placement, a hierarchical scene tree for organization, and a context-sensitive property inspector for detailed configuration.
- **Interactive Timeline:** A dedicated panel for visualizing and editing time-based events, such as platform motion paths and radar pulse schedules, with full playback controls.
- **Hierarchical Scene Tree:** An intuitive tree view for all simulation elements, enabling easy selection, parenting, and organization of complex scenarios.
- **Integrated Simulation Runner:** A focused workspace to configure global simulation parameters, trigger the FERS core engine, and monitor progress.
- **FERS XML Import/Export:** Generate a valid FERS XML configuration file directly from the visual scenario, or load an existing one to begin editing.

## Technology Stack

fers-ui is built with a modern technology stack to provide a powerful, cross-platform, and maintainable application.

- **Application Framework:** [**Tauri**](https://tauri.app/) (v2)
- **UI Library:** [**React**](https://react.dev/) with [**TypeScript**](https://www.typescriptlang.org/)
- **Component Library:** [**Material-UI (MUI)**](https://mui.com/material-ui/)
- **3D Rendering:** [**Three.js**](https://threejs.org/) via [**React Three Fiber**](https://docs.pmnd.rs/react-three-fiber)
- **State Management:** [**Zustand**](https://docs.pmnd.rs/zustand)
- **Schema Validation:** [**Zod**](https://zod.dev/)

## Getting Started

This application is part of a monorepo that includes the core C++ `libfers` library. The UI build process is entirely self-contained: Cargo will automatically invoke CMake to compile the C++ backend in an isolated directory during the build.

To set up the complete development environment, please follow the unified [Development Setup guide in the root README.md](https://github.com/stpaine/FERS/blob/master/README.md).

### Important Reminders:
1. You must have `vcpkg` installed and the `VCPKG_ROOT` environment variable set.
2. On Windows, you must use the **Developer PowerShell for VS** to ensure the MSVC toolchain (`cl.exe`, `link.exe`) is used instead of MinGW.

Once the environment is set up, you can run the UI from the **repository root** with:

```bash
bun ui:dev
```

To build a release bundle:

```bash
bun ui:build
```

> [!WARNING]
> The UI is currently in active development and may be unstable. Expect crashes and incomplete features.
> In particular, there is a known issue causing WebGL context loss on macOS on launch. See https://github.com/davidbits/FERS/issues/181 for details.
