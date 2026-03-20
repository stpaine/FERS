// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

use std::env;
use std::fs;
use std::path::PathBuf;

fn main() {
    // Load environment variables from .env file if it exists
    dotenv::dotenv().ok();

    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let repo_root = manifest_dir.join("../../..");

    // --- 1. Ensure VCPKG_ROOT is set ---
    let vcpkg_root = env::var("VCPKG_ROOT")
        .expect("VCPKG_ROOT must be set in your environment or in packages/fers-ui/src-tauri/.env");

    // --- 2. Invoke CMake to build libfers ---
    let mut config = cmake::Config::new(&repo_root);
    let vcpkg_toolchain = format!("{}/scripts/buildsystems/vcpkg.cmake", vcpkg_root);

    config
        .generator("Ninja")
        .define("CMAKE_TOOLCHAIN_FILE", &vcpkg_toolchain)
        .define("FERS_BUILD_TESTS", "OFF")
        .define("FERS_BUILD_SHARED_LIBS", "OFF")
        .define("FERS_BUILD_STATIC_LIBS", "ON")
        .build_target("fers_static")
        .profile("Release");

    let dst = config.build();
    let build_dir = dst.join("build");

    // --- 3. Locate vcpkg installed directory ---
    let vcpkg_installed_dir = build_dir.join("vcpkg_installed");

    let mut triplet_dir = None;
    if vcpkg_installed_dir.exists() {
        for entry in fs::read_dir(&vcpkg_installed_dir).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            if path.is_dir() && entry.file_name() != "vcpkg" {
                triplet_dir = Some(path);
                break;
            }
        }
    }

    let triplet_dir =
        triplet_dir.expect("Could not find vcpkg triplet directory. Did CMake run successfully?");
    let vcpkg_lib_dir = triplet_dir.join("lib");

    // --- 4. Link the libraries ---
    let libfers_lib_dir = build_dir.join("packages/libfers");
    println!("cargo:rustc-link-search=native={}", libfers_lib_dir.display());
    println!("cargo:rustc-link-lib=static=fers");

    println!("cargo:rustc-link-search=native={}", vcpkg_lib_dir.display());

    // We bypass pkg-config entirely for vcpkg dependencies because vcpkg's .pc files
    // for static libraries often omit transitive dependencies like szip and libaec.
    // Order is critical here: Dependents must appear BEFORE their dependencies.
    let vcpkg_libs = [
        "Geographic",
        "GeographicLib",
        "hdf5_hl_cpp",
        "hdf5_cpp",
        "hdf5_hl",
        "hdf5",
        "sz",
        "aec",
        "xml2",
        "lzma",
        "iconv",
        "z",
    ];

    for lib in vcpkg_libs {
        let lib_a = format!("lib{}.a", lib);
        let lib_lib = format!("{}.lib", lib);
        let lib_lib2 = format!("lib{}.lib", lib);

        // Only link if vcpkg actually built it
        if vcpkg_lib_dir.join(&lib_a).exists()
            || vcpkg_lib_dir.join(&lib_lib).exists()
            || vcpkg_lib_dir.join(&lib_lib2).exists()
        {
            println!("cargo:rustc-link-lib=static={}", lib);
        }
    }

    if cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=c++");
    } else if cfg!(target_os = "linux") {
        println!("cargo:rustc-link-lib=stdc++");
    }

    // --- 5. Generate Rust bindings for the C++ API ---
    let header_path = repo_root.join("packages/libfers/include/libfers/api.h");

    println!("cargo:rerun-if-changed={}", header_path.display());
    println!("cargo:rerun-if-env-changed=VCPKG_ROOT");
    println!("cargo:rerun-if-changed={}", repo_root.join("vcpkg.json").display());
    println!("cargo:rerun-if-changed={}", repo_root.join("packages/libfers/src").display());
    println!("cargo:rerun-if-changed={}", repo_root.join("packages/libfers/include").display());
    println!(
        "cargo:rerun-if-changed={}",
        repo_root.join("packages/libfers/CMakeLists.txt").display()
    );

    let bindings = bindgen::Builder::default()
        .header(header_path.to_str().unwrap())
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .allowlist_function("fers_.*")
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings.write_to_file(out_path.join("bindings.rs")).expect("Couldn't write bindings!");

    // --- 6. Tauri build step ---
    tauri_build::build();
}
