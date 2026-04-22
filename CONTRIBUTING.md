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
- [Continuous Integration](#continuous-integration)
- [Pull Request Process](#pull-request-process)
- [Style Guides](#style-guides)
- [Commit Message Guidelines](#commit-message-guidelines)

## Ways to Contribute

You can contribute in many ways:

- Reporting bugs and suggesting new features.
- Improving documentation.
- Adding or improving tests.
- Writing code to fix bugs or implement new features.

## Reporting Bugs

If you find a bug, please check the [existing issues](https://github.com/stpaine/FERS/issues) to see if it has
already been reported. If not,
please [open a new bug report](https://github.com/stpaine/FERS/issues/new/choose).

When filing a bug report, please include as many details as possible:

- A clear and descriptive title.
- A detailed description of the problem and steps to reproduce it.
- The expected behavior and what actually happened.
- Screenshots or logs if applicable.
- Information about your environment (OS, versions, etc.).

## Suggesting Enhancements

If you have an idea for a new feature or an improvement, we'd love to hear about it! Please check
the [existing issues and feature requests](https://github.com/stpaine/FERS/issues) first. If your idea is new,
please [open a new feature request](https://github.com/stpaine/FERS/issues/new/choose).

Provide a clear description of the proposed enhancement and explain the problem it solves or the value it adds.

## Your First Code Contribution

Unsure where to begin contributing to FERS? You can start by looking through `good first issue` and `help wanted`
issues:

- [Good First Issues](https://github.com/stpaine/FERS/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22) - issues
  which should only require a few lines of code, and a test or two.
- [Help Wanted Issues](https://github.com/stpaine/FERS/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22) - issues which
  should be a bit more involved than `good first issue` issues.

## Development Setup

This is a monorepo containing the C++ core simulator (`libfers` / `fers-cli`) and the Tauri desktop application (`fers-ui`).

For complete instructions on installing prerequisites (CMake, vcpkg, Bun, Rust) and building the project across Windows, macOS, and Linux, please refer to the **[Development Setup section in the root README.md](README.md#development-setup)**.

## Continuous Integration

Continuous Integration (CI) workflows run on every pull request and push to the `master` branch to verify the core build and tests, as well as the UI build. Please ensure your code passes locally before submitting a PR.

If you touch version metadata, run `python3 scripts/verify_versions.py` locally before pushing.
Release PRs additionally regenerate `bun.lock` and tracked license files automatically via `python3 scripts/refresh_release_pr_artifacts.py`.

## Pull Request Process

1. Fork the repository and create your branch from `master`.
2. If you've added code that should be tested, add tests.
3. Make sure your code lints and follows the style guides below.
4. Open a pull request with a clear title and a detailed description of your changes. Link to any relevant issues.
5. Be prepared to address feedback from the maintainers. The CI build must pass before a pull request can be merged.

## Style Guides

### C++ (`libfers`) Style Guide

Please adhere to the existing code style. We use `.clang-format` to enforce formatting.

- **Standards**: Use modern C++23 features where appropriate (smart pointers, concepts, etc.).
- **Naming Conventions**:
    - Variables and functions: `snake_case`
    - Classes and structs: `UpperCamelCase`
    - Class member functions: `lowerCamelCase`
    - Constants: `ALL_UPPER_SNAKE_CASE`
- **Comments**: Use Doxygen-style comments for functions and classes.

### TypeScript/React (`fers-ui`) Style Guide

We use Biome and Bun to enforce a consistent code style. Please run the UI lint task before committing:

```bash
bun lint:js
```

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

Closes #42
```
