// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

use std::collections::HashSet;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

fn env_var(key: &str) -> Option<String> {
    env::var(key).ok().filter(|value| !value.is_empty())
}

fn apple_triplet(target: &str) -> Option<&'static str> {
    match target {
        "aarch64-apple-darwin" => Some("arm64-osx"),
        "x86_64-apple-darwin" => Some("x64-osx"),
        _ => None,
    }
}

fn apple_architecture(target: &str) -> Option<&'static str> {
    match target {
        "aarch64-apple-darwin" => Some("arm64"),
        "x86_64-apple-darwin" => Some("x86_64"),
        _ => None,
    }
}

fn linux_triplet(target: &str) -> Option<&'static str> {
    if target.starts_with("aarch64") {
        Some("arm64-linux")
    } else if target.starts_with("x86_64") {
        Some("x64-linux")
    } else {
        None
    }
}

fn windows_triplet(target: &str) -> Option<&'static str> {
    if target.starts_with("x86_64-pc-windows") {
        Some("x64-windows-static-md")
    } else if target.starts_with("aarch64-pc-windows") {
        Some("arm64-windows-static-md")
    } else {
        None
    }
}

fn static_library_name(path: &Path) -> Option<String> {
    let file_name = path.file_name()?.to_str()?;
    let name = file_name.strip_suffix(".a").or_else(|| file_name.strip_suffix(".lib"))?;

    let name = name.strip_prefix("lib").unwrap_or(name);
    if name.is_empty() {
        None
    } else {
        Some(name.to_owned())
    }
}

fn emit_static_links_from_file(links_file: &Path) {
    let contents = fs::read_to_string(links_file).unwrap_or_else(|err| {
        panic!("failed to read generated native link manifest {}: {err}", links_file.display())
    });

    let mut seen_dirs = HashSet::new();
    let mut seen_libs = HashSet::new();

    for line in contents.lines().map(str::trim).filter(|line| !line.is_empty()) {
        let lib_path = Path::new(line);

        if let Some(lib_dir) = lib_path.parent() {
            if !seen_dirs.contains(lib_dir) {
                println!("cargo:rustc-link-search=native={}", lib_dir.display());
                seen_dirs.insert(lib_dir.to_path_buf());
            }
        }

        let lib_name = static_library_name(lib_path).unwrap_or_else(|| {
            panic!("unsupported static library entry in {}: {}", links_file.display(), line)
        });

        if !seen_libs.contains(&lib_name) {
            println!("cargo:rustc-link-lib=static={}", lib_name);
            seen_libs.insert(lib_name);
        }
    }
}

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let repo_root = manifest_dir.join("../../..");

    // --- 1. Ensure VCPKG_ROOT is set ---
    // Resolution order:
    //   1. VCPKG_ROOT environment variable (explicit, always wins)
    //   2. repo-local vcpkg/ directory (convenient for contributors who clone vcpkg alongside)
    //   3. ~/vcpkg or %USERPROFILE%\vcpkg (common global install locations)
    // Note: a sibling ../vcpkg is intentionally not searched; contributors should set VCPKG_ROOT
    // or install vcpkg into one of the two locations above.
    let vcpkg_root = env::var("VCPKG_ROOT")
        .or_else(|_| {
            let mut candidates = vec![repo_root.join("vcpkg")];
            if let Some(home) = env_var("HOME") {
                candidates.push(PathBuf::from(home).join("vcpkg"));
            }
            if let Some(userprofile) = env_var("USERPROFILE") {
                candidates.push(PathBuf::from(userprofile).join("vcpkg"));
            }

            candidates
                .into_iter()
                .find(|p| p.exists())
                .map(|p| p.to_string_lossy().into_owned())
                .ok_or("VCPKG_ROOT not found")
        })
        .expect(
            "VCPKG_ROOT must be set in your environment, or vcpkg must exist at <repo>/vcpkg, ~/vcpkg, or %USERPROFILE%\\vcpkg",
        );

    // --- 2. Invoke CMake to build libfers ---
    // Note: libfers is always built in Release mode regardless of the Cargo profile.
    // This keeps the C++ build fast and avoids debug-symbol bloat in the native library.
    // Rust-side debug information is unaffected by this choice.
    let mut config = cmake::Config::new(&repo_root);
    let vcpkg_toolchain = format!("{}/scripts/buildsystems/vcpkg.cmake", vcpkg_root);
    let target = env::var("TARGET").unwrap();

    let mut vcpkg_triplet = if target.contains("apple") {
        apple_triplet(&target).map(str::to_owned)
    } else if target.contains("linux") {
        linux_triplet(&target).map(str::to_owned)
    } else if target.contains("windows") {
        windows_triplet(&target).map(str::to_owned)
    } else {
        None
    };

    if let Some(explicit_triplet) = env_var("VCPKG_TARGET_TRIPLET") {
        vcpkg_triplet = Some(explicit_triplet);
    }

    let vcpkg_triplet =
        vcpkg_triplet.unwrap_or_else(|| panic!("unsupported target for Tauri C++ build: {target}"));

    config
        .define("CMAKE_TOOLCHAIN_FILE", &vcpkg_toolchain)
        .define("FERS_BUILD_TESTS", "OFF")
        .define("FERS_BUILD_SHARED_LIBS", "OFF")
        .define("FERS_BUILD_STATIC_LIBS", "ON")
        .build_target("fers_static")
        .profile("Release");

    config.define("VCPKG_TARGET_TRIPLET", &vcpkg_triplet);

    if target.contains("apple") {
        if let Some(arch) = apple_architecture(&target) {
            config.define("CMAKE_OSX_ARCHITECTURES", arch);
        }

        // Cargo and Tauri inject MACOSX_DEPLOYMENT_TARGET=10.13 into the build environment.
        // We must explicitly override this for our C++ library to support std::filesystem (10.15+)
        // and other modern features. We enforce 14.0 here.
        let deployment_target = "14.0";
        config.define("CMAKE_OSX_DEPLOYMENT_TARGET", deployment_target);
        config.define("VCPKG_OSX_DEPLOYMENT_TARGET", deployment_target);
        config.env("MACOSX_DEPLOYMENT_TARGET", deployment_target);

        if let Some(sysroot) = env_var("SDKROOT") {
            config.define("CMAKE_OSX_SYSROOT", &sysroot);
        }
    }

    let dst = config.build();
    let build_dir = dst.join("build");

    // --- 3. Link the libraries ---
    let libfers_lib_dir = build_dir.join("packages/libfers");
    println!("cargo:rustc-link-search=native={}", libfers_lib_dir.display());
    println!("cargo:rustc-link-lib=static=fers");

    let vcpkg_links_file = libfers_lib_dir.join("vcpkg_links.txt");
    emit_static_links_from_file(&vcpkg_links_file);

    // Use the TARGET env var (not cfg!) so this is correct during cross-compilation.
    if target.contains("apple") {
        println!("cargo:rustc-link-lib=c++");
    } else if target.contains("linux") {
        println!("cargo:rustc-link-lib=stdc++");
    }

    // --- 4. Generate Rust bindings for the C++ API ---
    let header_path = repo_root.join("packages/libfers/include/libfers/api.h");

    println!("cargo:rerun-if-changed={}", header_path.display());
    println!("cargo:rerun-if-env-changed=VCPKG_ROOT");
    println!("cargo:rerun-if-env-changed=VCPKG_TARGET_TRIPLET");
    println!("cargo:rerun-if-env-changed=MACOSX_DEPLOYMENT_TARGET");
    println!("cargo:rerun-if-env-changed=CMAKE_OSX_DEPLOYMENT_TARGET");
    println!("cargo:rerun-if-env-changed=SDKROOT");
    println!("cargo:rerun-if-changed={}", repo_root.join("vcpkg.json").display());
    println!("cargo:rerun-if-changed={}", repo_root.join("CMakeLists.txt").display());
    println!("cargo:rerun-if-changed={}", repo_root.join("packages/CMakeLists.txt").display());
    println!("cargo:rerun-if-changed={}", repo_root.join("packages/libfers/src").display());
    println!("cargo:rerun-if-changed={}", repo_root.join("packages/libfers/include").display());
    println!(
        "cargo:rerun-if-changed={}",
        repo_root.join("packages/libfers/CMakeLists.txt").display()
    );
    println!("cargo:rerun-if-changed={}", vcpkg_links_file.display());

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
