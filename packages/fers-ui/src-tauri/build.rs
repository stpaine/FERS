// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

//! # Tauri Application Build Script
//!
//! This build script configures the linking and code generation for the Tauri
//! desktop application that wraps the `libfers` C++ library. It performs three
//! key tasks:
//!
//! 1. **Library Linking**: Configures Cargo to link the `libfers` static library
//!    and its dependencies (both project-internal and system-provided).
//! 2. **FFI Binding Generation**: Uses `bindgen` to automatically generate Rust
//!    bindings from the C-style API header (`api.h`).
//! 3. **Tauri Integration**: Invokes the Tauri build process to prepare the
//!    application bundle.
//!
//! ## Build Dependencies
//!
//! This script assumes:
//! - The C++ libraries have been built and placed in the `build/` directory (via CMake).
//! - System dependencies (`libxml2`, `hdf5`) are available via `pkg-config`.
//! - The `bindgen` crate is available for FFI code generation.

use std::env;
use std::path::PathBuf;

fn main() {
    // --- 1. Link C++ libraries ---

    // -- Link libraries built within the project --

    // The directory containing the `Cargo.toml` file for this crate.
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());

    // The root build directory where CMake outputs are located.
    // This assumes a parallel build directory structure: `packages/fers-ui/src-tauri`
    // at the same level as `build/`.
    let build_dir = manifest_dir.join("../../../build");

    // Directory containing the compiled `libfers.a` static library.
    let libfers_lib_dir = build_dir.join("packages/libfers");

    // Directory containing the compiled `GeographicLib` static library.
    // CMake's FetchContent outputs build artifacts to the `_deps` directory.
    let geographiclib_lib_dir = build_dir.join("_deps/geographiclib-build/src");

    // Tell Cargo where to find our custom-built libraries.
    println!("cargo:rustc-link-search=native={}", libfers_lib_dir.display());
    println!("cargo:rustc-link-search=native={}", geographiclib_lib_dir.display());

    // Link the `libfers` static library (compiled from C++).
    println!("cargo:rustc-link-lib=static=fers");

    // Link the `GeographicLib` static library (third-party dependency).
    println!("cargo:rustc-link-lib=static=GeographicLib");

    // -- Find and link system dependencies using pkg-config --

    // Probe for `libxml2` and configure linking automatically.
    pkg_config::probe_library("libxml-2.0").unwrap();

    // Probe for `hdf5` and configure linking automatically.
    pkg_config::probe_library("hdf5").unwrap();

    // Link correct C++ standard library based on the target OS.
    // The C++ standard library is required because `libfers` is written in C++.
    if cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=c++");
    } else if cfg!(target_os = "linux") {
        println!("cargo:rustc-link-lib=stdc++");
    }
    // On Windows with MSVC, the C++ standard library is linked automatically,
    // so no explicit configuration is needed.

    // --- 2. Generate Rust bindings for the C++ API ---

    // Path to the C-style API header file that defines the FFI boundary.
    let header_path = manifest_dir.join("../../libfers/include/libfers/api.h");

    let libfers_path = libfers_lib_dir.join("libfers.a");
    if libfers_path.exists() {
        println!("cargo:rerun-if-changed={}", libfers_path.display());
    }

    // Also watch the GeographicLib library.
    let geographiclib_path = geographiclib_lib_dir.join(if cfg!(target_os = "windows") {
        "GeographicLib.lib"
    } else {
        "libGeographicLib.a"
    });

    if geographiclib_path.exists() {
        println!("cargo:rerun-if-changed={}", geographiclib_path.display());
    }

    // Generate Rust FFI bindings using `bindgen`.
    // - Only functions matching the `fers_*` pattern are exposed.
    // - Callbacks are used to emit Cargo configuration automatically.
    let bindings = bindgen::Builder::default()
        .header(header_path.to_str().unwrap())
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .allowlist_function("fers_.*")
        .generate()
        .expect("Unable to generate bindings");

    // The output directory where Cargo places build artifacts.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    // Write the generated bindings to a file that can be included at compile time.
    bindings.write_to_file(out_path.join("bindings.rs")).expect("Couldn't write bindings!");

    // --- 3. Tauri build step ---

    // Invoke the Tauri build process to prepare the desktop application.
    // This handles bundling, icon generation, and platform-specific configuration.
    tauri_build::build();
}
