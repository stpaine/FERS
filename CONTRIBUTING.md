# Contributing to FERS

First off, thank you for considering contributing to FERS! We welcome contributions from the community to help improve
the C++ core, the UI, documentation, and more. Every contribution is appreciated.

To ensure a smooth and effective process, please read these guidelines before you start.

## Table of Contents

- [Ways to Contribute](#ways-to-contribute)
- [Reporting Bugs](#reporting-bugs)
- [Suggesting Enhancements](#suggesting-enhancements)
- [Your First Code Contribution](#your-first-code-contribution)
- [Development Setup](#development-setup)
    - [Core Simulator (`fers`) Setup](#core-simulator-fers-setup)
    - [User Interface (`fers-ui`) Setup](#user-interface-fers-ui-setup)
- [Pull Request Process](#pull-request-process)
- [Style Guides](#style-guides)
    - [C++ (`fers`) Style Guide](#c-fers-style-guide)
    - [TypeScript/React (`fers-ui`) Style Guide](#typescriptreact-fers-ui-style-guide)
- [Commit Message Guidelines](#commit-message-guidelines)

## Ways to Contribute

You can contribute in many ways:

- Reporting bugs and suggesting new features.
- Improving documentation.
- Adding or improving tests.
- Writing code to fix bugs or implement new features.

## Reporting Bugs

If you find a bug, please check the [existing issues](https://github.com/davidbits/FERS/issues) to see if it has
already been reported. If not,
please [open a new bug report](https://github.com/davidbits/FERS/issues/new/choose).

When filing a bug report, please include as many details as possible:

- A clear and descriptive title.
- A detailed description of the problem and steps to reproduce it.
- The expected behavior and what actually happened.
- Screenshots or logs if applicable.
- Information about your environment (OS, versions, etc.).

## Suggesting Enhancements

If you have an idea for a new feature or an improvement, we'd love to hear about it! Please check
the [existing issues and feature requests](https://github.com/davidbits/FERS/issues) first. If your idea is new,
please [open a new feature request](https://github.com/davidbits/FERS/issues/new/choose).

Provide a clear description of the proposed enhancement and explain the problem it solves or the value it adds.

## Your First Code Contribution

Unsure where to begin contributing to FERS? You can start by looking through `good first issue` and `help wanted`
issues:

- [Good First Issues](https://github.com/davidbits/FERS/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22) - issues
  which should only require a few lines of code, and a test or two.
- [Help Wanted Issues](https://github.com/davidbits/FERS/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22) - issues which
  should be a bit more involved than `good first issue` issues.

## Development Setup

This is a monorepo containing two main packages. Please follow the setup instructions for the package you intend to work
on.

### Core Simulator (`fers`) Setup

The core simulator is written in C++.

1. **Prerequisites**: Ensure you have a C++23 compiler, CMake, `libhdf5`, and `libxml2`.
2. **Clone the repo**: `git clone --recursive https://github.com/davidbits/FERS.git`
3. **Build**: Follow the detailed build instructions in the [`packages/fers/README.md`](packages/fers/README.md).

### User Interface (`fers-ui`) Setup

The UI is a Tauri desktop application built with React and TypeScript.

1. **Prerequisites**: Ensure you have Node.js, pnpm, and the Rust toolchain installed.
2. **Setup Tauri**: Follow the [Tauri prerequisites guide](https://tauri.app/start/prerequisites/) for your OS.
3. **Install & Run**:
    ```bash
    cd packages/fers-ui
    pnpm install
    pnpm tauri dev
    ```

For more details, see the [`packages/fers-ui/README.md`](packages/fers-ui/README.md).

## Pull Request Process

1. Fork the repository and create your branch from `main`.
2. If you've added code that should be tested, add tests.
3. Ensure the test suite passes (`run_sim_tests.py` for `fers`, and relevant tests for `fers-ui`).
4. Make sure your code lints and follows the style guides below.
5. Open a pull request with a clear title and a detailed description of your changes. Link to any relevant issues.
6. Be prepared to address feedback from the maintainers. The CI build must pass before a pull request can be merged.

## Style Guides

### C++ (`fers`) Style Guide

Please adhere to the existing code style. We use `.clang-format` to enforce formatting.

- **Standards**: Use modern C++23 features where appropriate (smart pointers, concepts, etc.).
- **Naming Conventions**:
    - Variables and functions: `snake_case`
    - Classes and structs: `UpperCamelCase`
    - Class member functions: `lowerCamelCase`
    - Constants: `ALL_UPPER_SNAKE_CASE`
- **Comments**: Use Doxygen-style comments for functions and classes.

### TypeScript/React (`fers-ui`) Style Guide

We use ESLint and Prettier to enforce a consistent code style. Please run `pnpm lint` and `pnpm format` before
committing.

- **Naming Conventions**:
    - Components and Types: `PascalCase`
    - Functions, variables, hooks: `camelCase`
- **Structure**: Follow the existing project structure (separating components, views, hooks, etc.).
- **Principles**: Adhere to principles like Single Responsibility and Don't Repeat Yourself (DRY).

## Commit Message Guidelines

We follow the [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) specification. This helps us
automate releases and makes the project history easier to read. Each commit message should consist of a header, a body,
and a footer.

```
<type>[optional scope]: <description>

[optional body]

[optional footer(s)]
```

**Types:** `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`, `build`, `ci`.

**Example:**

```
feat(ui): add property inspector for antennas

Implemented a new component in the left sidebar that displays
and allows editing of properties for the selected antenna element.

Fixes #42
```
