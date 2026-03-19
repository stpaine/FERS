// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

//! # Tauri Desktop Application Entry Point
//!
//! This is the minimal `main.rs` entry point for the FERS Tauri desktop application.
//! Its sole responsibility is to delegate to the library's `run()` function, which
//! handles all application initialization and event loop management.
//!
//! ## Architecture Rationale
//!
//! By keeping `main.rs` minimal and placing all application logic in `lib.rs`, we:
//! * Enable the application code to be unit-tested (the `main` function itself cannot be tested).
//! * Allow the same code to be used in different contexts (e.g., integration tests, benchmarks).
//! * Follow Rust best practices for library-centric project structure.
//!
//! ## Platform-Specific Configuration
//!
//! The `#![cfg_attr]` attribute suppresses the console window on Windows release builds,
//! ensuring a native desktop application experience without a terminal appearing behind
//! the UI.

// Prevents additional console window on Windows in release, DO NOT REMOVE!!
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

/// Launches the FERS Tauri desktop application.
///
/// This function immediately delegates to `fers_ui_lib::run()`, which:
/// 1. Initializes the FFI connection to `libfers` (the C++ simulation engine).
/// 2. Sets up Tauri plugins for file dialogs, filesystem access, and shell operations.
/// 3. Registers all Tauri commands that the frontend can invoke.
/// 4. Starts the Tauri event loop, which handles UI rendering and IPC.
///
/// # Panics
///
/// This function will panic if:
/// * The `libfers` C++ library cannot be linked or initialized.
/// * The Tauri application fails to start due to configuration errors.
///
/// These panics are intentional, as the application cannot function without
/// a valid simulation context or UI framework.
///
/// # Example
///
/// This is the standard entry point and should not typically be modified:
///
/// ```no_run
/// fers_ui_lib::run();
/// ```
fn main() {
    fers_ui_lib::run()
}
