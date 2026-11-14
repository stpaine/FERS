# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

This changelog tracks monorepo releases starting with `v0.1.0`.
Legacy `libfers` history remains in `packages/libfers/CHANGELOG.md`.

## [Unreleased]

### Added

- CLI VITA run metadata files containing output configuration and final stream counters.

### Fixed

- Apply XML random seeds before constructing stochastic models and retain generated seeds for replay.
- Preserve the accepted world, parameters, and master generator after rejected JSON replacements.
- Correct parabolic antenna boresight gain and preserve its near-axis continuity.
- Retain zero-valued samples when ADC quantization is enabled.
- Join pulsed finalizer workers during normal completion and exception cleanup.
- Process every event at an admitted timestamp before evaluating termination.
- Clip FMCW IF output to the requested simulation extent.
- Snap exact streaming sample boundaries while retaining ceiling behaviour off the grid.
- Retain scheduled HDF5 receiver intervals with one final scale and one metadata entry per file.
- Send VITA closing context after pending data with the exclusive stream-end timestamp.
- Assign VITA packet counts in transmission order, including failed send attempts.
- Export moving KML paths as altitude-aware lines with attached-object descriptions and built-in symbols.

### Documentation

- Update XML/API contracts, output interpretation, build and test instructions, and feature verification guidance.
