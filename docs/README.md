# FERS documentation

User guides live in [wiki/](wiki/). The wiki workflow copies these tracked pages to the repository's GitHub wiki. API documentation is generated from source comments and the project READMEs using [Doxyfile](Doxyfile).

## Build API documentation

Install Doxygen and Graphviz, then run from the repository root:

```bash
cmake -S . -B build/docs -DFERS_DOCS_ONLY=ON -DFERS_BUILD_DOCS=ON
cmake --build build/docs --target doc
```

Open `build/docs/docs/html/index.html` to view the result. To add documentation to an existing core build, configure that build with `FERS_BUILD_DOCS=ON` and build its `doc` target.

## Theme

The Doxygen Awesome theme is vendored in `doxygen-awesome-css/`. Update its `doxygen-awesome*` files together from the [theme repository](https://github.com/jothepro/doxygen-awesome-css). Project-specific styling and extensions are configured in `custom.css`, `header.html.in`, and `CMakeLists.txt`.
