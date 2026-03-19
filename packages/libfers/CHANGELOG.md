# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.4] - 2025-07-31

### Changed

- Updated `CMakeLists.txt` to support Python 3.11.13.
- Updated `python_extension.cpp` to support Python 3.11.13.

## [1.0.3] - 2025-05-09

### Changed

- Consolidated the `addDirect` and `simulateTarget` methods in `sim_threading.cpp`.
- Clean up some minor formatting issues in the codebase.
- Applied small static code analysis fixes.

### Removed

- Removed the unused `TaskThread` class.

## [1.0.2] - 2025-04-21

### Changed

- Removed `noexcept` from the FAlphaBranch constructor and flush() methods

## [1.0.1] - 2025-04-11

### Added

- .clang-format configuration file for consistent code formatting.

### Changed

- Updated CMake configuration to support Python 3.11.12
- Updated python_extension.cpp to support Python 3.11.12

## [1.0.0] - 2024-10-24

### Added

- Modernized the entire codebase to C++23.
- Added support for Python 3.7 to 3.11.10.
- Added HighFive for HDF5 support.
- Improved documentation.

### Changed

- Improved CMake setup.
- Removed TinyXML and replaced it with libxml2.
- Removed Boost and FFTW3 dependencies.
- Migrate to semantic versioning.

### Fixed

- Eliminated memory leaks.
- Fixed known multipath surface SIGSEGV bug.
- Fixed unknown h5 antenna gain pattern bug.

## [0.23] - 2013-12-10

### Changed

- Fixed for current versions of Boost.
- General improvements with the CMake setup.
- Build configuration defaults to Release.

## [0.22] - 2013-03-15

### Fixed

- Fixed "[BUG] Requested delay filter value out of range" problem caused by outbound arrays indexes. This could still be revised more elegantly.

### Changed

- Tidied up general program output.
- Added title banner to program.

## [0.21] - 2012-07-16

### Removed

- Removed included TinyXML—TinyXML is now in Debian and Red Hat.

## [0.20] - 2012-07-06

### Changed

- Various patches from RRSG @ UCT.
- Fix up CMakeLists to allow both in-source and out-of-source building.

### Removed

- Remove debian directory.
- Remove `src/FindPythonLibs.cmake` - use CMake's module instead.

## [0.3] - 2006-11-03

### Added

- Added the ability to export to HDF5 files.

### Fixed

- Fixed phase and doppler calculations.
- Fixed up some examples to the new XML format. More work needs to be done here.
