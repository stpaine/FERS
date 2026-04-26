# FERS XML Schema

This directory contains the schema definitions for FERS XML inputs. These schemas cover both top-level simulation
scenario files and standalone XML asset formats that are validated by `libfers`.

## Files

- `fers-xml.xsd`: The **XML Schema Definition (XSD)** used for modern, robust validation of scenario files. This is the
  primary schema used by both the core simulator and the UI.
- `fers-xml.dtd`: The **Document Type Definition (DTD)**, an older schema format retained for legacy purposes.
- `antenna-pattern.xsd`: The XSD for standalone XML antenna-pattern assets loaded by `XmlAntenna`.
- `antenna-pattern.dtd`: The legacy DTD counterpart for standalone XML antenna-pattern assets.

## Role in the Monorepo

The scenario schema serves as the formal contract between the `libfers` core library and its clients (`fers-cli`,
`fers-ui`). The standalone antenna-pattern schema defines the contract for XML antenna asset files referenced by
scenario documents.

- **`fers-ui`** generates `.fersxml` files that conform to this schema.
- **`libfers`** validates incoming `.fersxml` files against this schema before running a simulation.
- **`libfers`** also validates standalone XML antenna-pattern assets against the embedded antenna-pattern schemas when
  those files are loaded.

During the build process of the `libfers` package, the schema files are converted into C header files and embedded
directly into the library binary. This ensures that the simulator always has access to the correct schema versions for
validation without relying on external files.

## Modifying the Schema

When making changes to a schema definition, update both the XSD and DTD for that format to maintain consistency. Any
changes must also be propagated throughout the codebase, particularly in the XML parsing logic within `libfers` and,
for scenario-schema changes, the XML generation logic in `fers-ui`.

The scenario XSD/DTD intentionally validate XML structure, not every numeric or cross-field physics rule. Scalar scenario
values are parsed by `libfers`, which performs the strict numeric checks that schemas cannot express portably across both
XSD and DTD, such as FMCW `chirp_period >= chirp_duration` and integer-positive `chirp_count`.
