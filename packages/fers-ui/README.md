# fers-ui: A Graphical Interface for the FERS Core Simulator

![Framework](https://img.shields.io/badge/Framework-React-61DAFB?logo=react)
![Language](https://img.shields.io/badge/Language-TypeScript-3178C6?logo=typescript)
![Desktop App](https://img.shields.io/badge/Tauri-v2-FFC336)
[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://github.com/davidbits/FERS/blob/master/LICENSE)

**fers-ui** is the official graphical user interface for the **Flexible Extensible Radar Simulator (FERS)**. It provides a professional-grade, visual workbench designed to streamline the entire simulation pipeline—from asset creation and scenario construction to simulation execution and results analysis. By offering an intuitive and powerful toolset, fers-ui dramatically improves the usability and accessibility of the core FERS engine.

## Key Features

- **Multi-Panel Scenario Builder:** A fully resizable workspace featuring an interactive 3D viewport for object placement, a hierarchical scene tree for organization, and a context-sensitive property inspector for detailed configuration.
- **Interactive Timeline:** A dedicated panel for visualizing and editing time-based events, such as platform motion paths and radar pulse schedules, with full playback controls.
- **Hierarchical Scene Tree:** An intuitive tree view for all simulation elements, enabling easy selection, parenting, and organization of complex scenarios.
- **Integrated Simulation Runner:** A focused workspace to configure global simulation parameters, trigger the FERS core engine, and monitor progress.
- **FERS XML Import/Export:** Generate a valid FERS XML configuration file directly from the visual scenario, or load an existing one to begin editing.

## Future Roadmap

The following major features are planned for future releases:

- **Centralized Asset Library:** A dedicated view for creating, managing, and reusing simulation components like radar pulses, antenna patterns, and timing sources across different scenarios.
- **Dynamic Results Analysis:** A post-simulation view for loading FERS output data to visualize target trajectories, platform movements, and analyze collected signals.

## Technology Stack

fers-ui is built with a modern technology stack to provide a powerful, cross-platform, and maintainable application.

- **Application Framework:** [**Tauri**](https://tauri.app/) (v2) - Provides a lightweight, secure, and performant way to build a native desktop application using a web front-end.
- **UI Library:** [**React**](https://react.dev/) with [**TypeScript**](https://www.typescriptlang.org/) - Offers a robust, component-based architecture for building complex UIs with the safety of static typing.
- **Component Library:** [**Material-UI (MUI)**](https://mui.com/material-ui/) & [**MUI X**](https://mui.com/x/) - A comprehensive library of UI components that accelerates development and ensures a consistent, modern design.
- **3D Rendering:** [**Three.js**](https://threejs.org/) via [**React Three Fiber**](https://docs.pmnd.rs/react-three-fiber) - Simplifies managing a 3D scene within a React application, powering the visual scenario builder.
- **State Management:** [**Zustand**](https://docs.pmnd.rs/zustand) - A minimal, fast, and scalable state management solution for a centralized application state.
- **Schema Validation:** [**Zod**](https://zod.dev/) - A TypeScript-first schema declaration and validation library, used for ensuring the integrity of scenario data on the front-end.

## Software Architecture

The application is architected as a multi-modal "Workbench" to provide a clean, context-focused user experience.

- **App Rail & Views:** The primary navigation is a vertical `AppRail` that allows switching between distinct workspaces called "Views" (e.g., `ScenarioView`, `AssetLibraryView`). This isolates different stages of the simulation workflow, preventing UI clutter.
- **Layouts:** The top-level `MainLayout` component orchestrates the active view and global UI structure.
- **Views:** These are high-level components that represent a complete user workspace. Each view is self-contained and manages its own panel arrangement (e.g., `ScenarioView` combines the 3D world, scene tree, and timeline).
- **Components:** The UI is built from small, reusable React components (`SceneTree`, `Timeline`, `PropertyInspector`, `ResizablePanel`) that encapsulate specific functionality and are composed within Views.

## Design Principles

- **Context-Focused Workflow:** The UI is designed to only present tools and information relevant to the user's current task. Switching from the "Scenario" to the "Asset Library" view changes the entire workspace context.
- **Component Modularity:** By breaking the UI into a clear hierarchy of `Layouts`, `Views`, and `Components`, we ensure the codebase is maintainable, testable, and easy to extend.
- **State-Driven UI:** The application relies on a centralized state managed by Zustand. UI components react to state changes, ensuring consistency across all panels and views.
- **Ergonomic Layout:** The use of resizable panels gives users full control over their workspace, allowing them to tailor the interface to their specific needs and monitor size.

## Getting Started

This application is part of a monorepo that includes the core C++ `libfers` library. To set up the complete development environment, please follow the unified **[Development Setup guide in the root README.md](https://github.com/davidbits/FERS/blob/master/README.md)**.

Once the environment is set up, you can run the UI from the **repository root** with:
`bash
    pnpm ui:dev
    `
