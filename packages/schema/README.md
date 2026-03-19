# FERS XML Schema

This directory contains the schema definitions that govern the structure of FERS scenario files. These schemas are the
single source of truth for defining simulation scenarios.

## Files

- `fers-xml.xsd`: The **XML Schema Definition (XSD)** used for modern, robust validation of scenario files. This is the
  primary schema used by both the core simulator and the UI.
- `fers-xml.dtd`: The **Document Type Definition (DTD)**, an older schema format retained for legacy purposes.

## Role in the Monorepo

The schema serves as the formal contract between the `libfers` core library and its clients (`fers-cli`, `fers-ui`).

- **`fers-ui`** generates `.fersxml` files that conform to this schema.
- **`libfers`** validates incoming `.fersxml` files against this schema before running a simulation.

During the build process of the `libfers` package, the XSD and DTD files are converted into C header files and embedded
directly into the library binary. This ensures that the simulator always has access to the correct schema version for
validation without relying on external files.

## Modifying the Schema

When making changes to the scenario definition, both `fers-xml.xsd` and `fers-xml.dtd` should be updated to maintain
consistency. Any changes must be propagated throughout the codebase, particularly in the XML parsing logic within the
`libfers` package and the XML generation logic in the `fers-ui` package.
