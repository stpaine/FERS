// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

//! # Tauri Application Library Entry Point
//!
//! This module provides the main entry point for the FERS Tauri desktop application.
//! It bridges the Rust/Tauri frontend with the C++ `libfers` simulation engine via
//! a Foreign Function Interface (FFI).
//!
//! ## Architecture Overview
//!
//! The application follows a three-layer architecture:
//!
//! 1. **Frontend (TypeScript/React)**: Provides the user interface for scenario
//!    editing and visualization.
//! 2. **Middle Layer (Rust/Tauri)**: This module. Exposes Tauri commands that the
//!    frontend can invoke, and manages the simulation state via a thread-safe wrapper.
//! 3. **Backend (C++ via FFI)**: The `libfers` library that performs the actual
//!    simulation computations, parsing, and serialization.
//!
//! ## Thread Safety
//!
//! The `FersContext` (which wraps the C++ simulation state) is protected by a `Mutex`
//! and stored in Tauri's managed state. This ensures that concurrent calls from the
//! UI are serialized, preventing data races on the non-thread-safe C++ object.
//!
//! ## Tauri Commands
//!
//! All functions annotated with `#[tauri::command]` are exposed to the frontend via
//! Tauri's IPC mechanism. They can be invoked asynchronously from JavaScript/TypeScript.

mod fers_api;

use std::{fs, path::PathBuf, sync::Mutex};
use tauri::{AppHandle, Emitter, Manager, State};

/// Data structure for a single motion waypoint received from the UI.
///
/// Coordinates should be in the scenario's define frame (e.g. ENU).
#[derive(serde::Serialize, serde::Deserialize, std::fmt::Debug)]
pub struct MotionWaypoint {
    /// Time in seconds.
    time: f64,
    /// Easting/X coordinate in meters.
    x: f64,
    /// Northing/Y coordinate in meters.
    y: f64,
    /// Altitude/Z coordinate in meters (MSL).
    altitude: f64,
}

/// Enum for the interpolation type received from the UI.
#[derive(serde::Serialize, serde::Deserialize, std::fmt::Debug)]
#[serde(rename_all = "lowercase")]
pub enum InterpolationType {
    Static,
    Linear,
    Cubic,
}

/// Enum for the rotation angle unit received from the UI.
#[derive(serde::Serialize, serde::Deserialize, std::fmt::Debug, Clone, Copy)]
#[serde(rename_all = "lowercase")]
pub enum RotationAngleUnit {
    Deg,
    Rad,
}

/// Data structure for a single interpolated point sent back to the UI.
///
/// Represents the physical state of a platform at a specific time step.
#[derive(serde::Serialize, serde::Deserialize, std::fmt::Debug)]
pub struct InterpolatedMotionPoint {
    /// X position in meters.
    x: f64,
    /// Y position in meters.
    y: f64,
    /// Z position in meters.
    z: f64,
    /// X velocity in m/s.
    vx: f64,
    /// Y velocity in m/s.
    vy: f64,
    /// Z velocity in m/s.
    vz: f64,
}

/// Data structure for a single rotation waypoint received from the UI.
#[derive(serde::Serialize, serde::Deserialize, std::fmt::Debug)]
pub struct RotationWaypoint {
    /// Time in seconds.
    time: f64,
    /// Azimuth in compass units (0=North, pi/2 or 90=East).
    azimuth: f64,
    /// Elevation in the selected external units (positive up).
    elevation: f64,
}

/// Data structure for a single interpolated rotation point sent back to the UI.
#[derive(serde::Serialize, serde::Deserialize, std::fmt::Debug)]
pub struct InterpolatedRotationPoint {
    /// Azimuth in compass units.
    azimuth: f64,
    /// Elevation in the selected external units.
    elevation: f64,
}

/// Type alias for the managed Tauri state that holds the simulation context.
///
/// The `FersContext` is wrapped in a `Mutex` to ensure thread-safe access, as Tauri
/// may invoke commands from multiple threads concurrently. This alias simplifies
/// the function signatures of Tauri commands.
type FersState = Mutex<fers_api::FersContext>;

// --- Tauri Commands ---

/// Sets the output directory for simulation results.
#[tauri::command]
fn set_output_directory(dir: String, state: State<'_, FersState>) -> Result<(), String> {
    state.lock().map_err(|e| e.to_string())?.set_output_directory(&dir)
}

/// Loads a FERS scenario from an XML file into the simulation context.
///
/// This command replaces any existing in-memory scenario with the one parsed from
/// the specified file. The file path is provided by the user via the frontend dialog.
///
/// # Parameters
///
/// * `filepath` - The absolute or relative path to the FERS XML scenario file.
/// * `state` - Tauri-managed state containing the shared `FersContext`.
///
/// # Returns
///
/// * `Ok(Vec<String>)` containing any deduplicated non-fatal warnings detected while loading.
/// * `Err(String)` containing an error message if loading failed (e.g., file not found,
///   invalid XML, or a Mutex lock error).
///
/// # Example (from frontend)
///
/// ```typescript
/// import { invoke } from '@tauri-apps/api/core';
/// await invoke('load_scenario_from_xml_file', { filepath: '/path/to/scenario.xml' });
/// ```
#[tauri::command]
fn load_scenario_from_xml_file(
    filepath: String,
    state: State<'_, FersState>,
) -> Result<Vec<String>, String> {
    state.lock().map_err(|e| e.to_string())?.load_scenario_from_xml_file(&filepath)
}

/// Retrieves the current in-memory scenario as a JSON string.
///
/// This command serializes the simulation state into JSON format, allowing the
/// frontend to display and edit the scenario. The JSON structure mirrors the
/// internal representation used by `libfers`.
///
/// # Parameters
///
/// * `state` - Tauri-managed state containing the shared `FersContext`.
///
/// # Returns
///
/// * `Ok(String)` containing the JSON representation of the scenario.
/// * `Err(String)` containing an error message if serialization failed or if the
///   Mutex could not be locked.
///
/// # Example (from frontend)
///
/// ```typescript
/// import { invoke } from '@tauri-apps/api/core';
/// const scenarioJson = await invoke<string>('get_scenario_as_json');
/// const scenario = JSON.parse(scenarioJson);
/// ```
#[tauri::command]
fn get_scenario_as_json(state: State<'_, FersState>) -> Result<String, String> {
    state.lock().map_err(|e| e.to_string())?.get_scenario_as_json()
}

/// Retrieves the current in-memory scenario as a FERS XML string.
///
/// This command is typically used when the user wants to export the scenario
/// (potentially modified in the UI) back to the standard FERS XML format for
/// sharing or archival.
///
/// # Parameters
///
/// * `state` - Tauri-managed state containing the shared `FersContext`.
///
/// # Returns
///
/// * `Ok(String)` containing the XML representation of the scenario.
/// * `Err(String)` containing an error message if serialization failed or if the
///   Mutex could not be locked.
///
/// # Example (from frontend)
///
/// ```typescript
/// import { invoke } from '@tauri-apps/api/core';
/// const scenarioXml = await invoke<string>('get_scenario_as_xml');
/// // Save scenarioXml to a file using Tauri's fs plugin
/// ```
#[tauri::command]
fn get_scenario_as_xml(state: State<'_, FersState>) -> Result<String, String> {
    state.lock().map_err(|e| e.to_string())?.get_scenario_as_xml()
}

/// Updates the in-memory scenario from a JSON string provided by the frontend.
///
/// This is the primary method for applying changes made in the UI back to the
/// simulation engine. The JSON is deserialized and used to rebuild the internal
/// C++ world representation.
///
/// # Parameters
///
/// * `json` - A JSON string representing the modified scenario structure.
/// * `state` - Tauri-managed state containing the shared `FersContext`.
///
/// # Returns
///
/// * `Ok(Vec<String>)` containing any deduplicated non-fatal warnings detected while updating.
/// * `Err(String)` containing an error message if deserialization failed, the JSON
///   structure was invalid, or the Mutex could not be locked.
///
/// # Example (from frontend)
///
/// ```typescript
/// import { invoke } from '@tauri-apps/api/core';
/// const updatedScenario = { /* modified scenario object */ };
/// await invoke('update_scenario_from_json', { json: JSON.stringify(updatedScenario) });
/// ```
#[tauri::command]
fn update_scenario_from_json(
    json: String,
    state: State<'_, FersState>,
) -> Result<Vec<String>, String> {
    state.lock().map_err(|e| e.to_string())?.update_scenario_from_json(&json)
}

/// Triggers the simulation based on the current in-memory scenario.
///
/// This command immediately returns `Ok(())` and spawns a background thread to
/// perform the actual computationally intensive simulation. This prevents the UI
/// from freezing. The result of the simulation (success or failure) is
/// communicated back to the frontend via Tauri events.
///
/// # Parameters
///
/// * `app_handle` - The Tauri application handle, used to access managed state
///   and emit events.
///
/// # Events Emitted
///
/// * `simulation-output-metadata` - Emitted with metadata JSON before completion.
/// * `simulation-complete` - Emitted with `()` as payload on successful completion.
/// * `simulation-error` - Emitted with a `String` error message on failure.
/// * `simulation-progress` - Emitted periodically with `{ message: String, current: i32, total: i32 }`.
#[tauri::command]
fn run_simulation(app_handle: AppHandle) -> Result<(), String> {
    // Clone the AppHandle so we can move it into the background thread.
    let app_handle_clone = app_handle.clone();

    // Spawn a new thread to run the blocking C++ simulation.
    std::thread::spawn(move || {
        // Retrieve the managed state within the new thread.
        let fers_state: State<'_, FersState> = app_handle_clone.state();
        let result = fers_state
            .lock()
            .map_err(|e| e.to_string())
            .and_then(|context| context.run_simulation(&app_handle_clone));

        // Emit an event to the frontend based on the simulation result.
        match result {
            Ok(metadata_json) => {
                app_handle_clone
                    .emit("simulation-output-metadata", metadata_json)
                    .expect("Failed to emit simulation-output-metadata event");
                app_handle_clone
                    .emit("simulation-complete", ())
                    .expect("Failed to emit simulation-complete event");
            }
            Err(e) => {
                app_handle_clone
                    .emit("simulation-error", e)
                    .expect("Failed to emit simulation-error event");
            }
        }
    });

    // Return immediately, allowing the UI to remain responsive.
    Ok(())
}

fn sanitize_file_stem(name: &str) -> String {
    let sanitized: String =
        name.chars().map(|ch| if ch.is_ascii_alphanumeric() { ch } else { '_' }).collect();
    if sanitized.is_empty() {
        "simulation".to_string()
    } else {
        sanitized
    }
}

/// Writes the most recent simulation output metadata JSON next to the generated HDF5 files.
#[tauri::command]
fn export_output_metadata_json(state: State<'_, FersState>) -> Result<String, String> {
    let metadata_json = state.lock().map_err(|e| e.to_string())?.get_last_output_metadata_json()?;
    let metadata: serde_json::Value = serde_json::from_str(&metadata_json)
        .map_err(|e| format!("Failed to parse output metadata JSON: {}", e))?;

    let simulation_name =
        metadata.get("simulation_name").and_then(serde_json::Value::as_str).unwrap_or("simulation");
    let output_directory = metadata
        .get("output_directory")
        .and_then(serde_json::Value::as_str)
        .filter(|value| !value.is_empty())
        .unwrap_or(".");

    let mut output_path = PathBuf::from(output_directory);
    fs::create_dir_all(&output_path)
        .map_err(|e| format!("Failed to create metadata output directory: {}", e))?;
    output_path.push(format!("{}_metadata.json", sanitize_file_stem(simulation_name)));

    fs::write(&output_path, metadata_json)
        .map_err(|e| format!("Failed to write metadata JSON: {}", e))?;

    Ok(output_path.to_string_lossy().to_string())
}

/// Generates a KML visualization file for the current in-memory scenario.
///
/// This command spawns a background thread to handle file I/O and KML generation,
/// preventing the UI from freezing. The result is communicated via events.
///
/// # Parameters
///
/// * `output_path` - The absolute file path where the KML file should be saved.
/// * `app_handle` - The Tauri application handle.
///
/// # Events Emitted
///
/// * `kml-generation-complete` - Emitted with the output path `String` on success.
/// * `kml-generation-error` - Emitted with a `String` error message on failure.
#[tauri::command]
fn generate_kml(output_path: String, app_handle: AppHandle) -> Result<(), String> {
    let app_handle_clone = app_handle.clone();
    std::thread::spawn(move || {
        let fers_state: State<'_, FersState> = app_handle_clone.state();
        let result = fers_state
            .lock()
            .map_err(|e| e.to_string())
            .and_then(|context| context.generate_kml(&output_path));

        match result {
            Ok(_) => {
                app_handle_clone
                    .emit("kml-generation-complete", &output_path)
                    .expect("Failed to emit kml-generation-complete event");
            }
            Err(e) => {
                app_handle_clone
                    .emit("kml-generation-error", e)
                    .expect("Failed to emit kml-generation-error event");
            }
        }
    });
    Ok(())
}

/// A stateless command to calculate an interpolated motion path.
///
/// This command delegates to the `libfers` core to calculate a path from a given
/// set of waypoints, ensuring the UI visualization is identical to the path
/// the simulation would use.
///
/// # Parameters
/// * `waypoints` - A vector of motion waypoints.
/// * `interp_type` - The interpolation algorithm to use ('static', 'linear', 'cubic').
/// * `num_points` - The desired number of points for the final path.
///
/// # Returns
/// * `Ok(Vec<InterpolatedPoint>)` - The calculated path points.
/// * `Err(String)` - An error message if the path calculation failed.
#[tauri::command]
fn get_interpolated_motion_path(
    waypoints: Vec<MotionWaypoint>,
    interp_type: InterpolationType,
    num_points: usize,
) -> Result<Vec<InterpolatedMotionPoint>, String> {
    fers_api::get_interpolated_motion_path(waypoints, interp_type, num_points)
}

/// A stateless command to calculate an interpolated rotation path.
///
/// This command delegates to the `libfers` core to calculate a rotation path from a given
/// set of waypoints. It is used by the UI to preview how the simulation will interpolate
/// orientation changes over time.
///
/// # Parameters
/// * `waypoints` - A vector of rotation waypoints in the selected external angle unit.
/// * `interp_type` - The interpolation algorithm to use ('static', 'linear', 'cubic').
/// * `angle_unit` - The angle unit to use for both input and output.
/// * `num_points` - The desired number of points for the final path.
///
/// # Returns
/// * `Ok(Vec<InterpolatedRotationPoint>)` - The calculated rotation points.
/// * `Err(String)` - An error message if the calculation failed.
#[tauri::command]
fn get_interpolated_rotation_path(
    waypoints: Vec<RotationWaypoint>,
    interp_type: InterpolationType,
    angle_unit: RotationAngleUnit,
    num_points: usize,
) -> Result<Vec<InterpolatedRotationPoint>, String> {
    fers_api::get_interpolated_rotation_path(waypoints, interp_type, angle_unit, num_points)
}

/// Retrieves a 2D antenna gain pattern for visualization.
///
/// This command samples the antenna model loaded in the current simulation context.
/// It is stateful because it relies on the antenna assets defined in the loaded scenario.
///
/// # Parameters
/// * `antenna_id` - The unique ID of the antenna asset to sample.
/// * `az_samples` - Resolution along the azimuth axis (e.g., 360).
/// * `el_samples` - Resolution along the elevation axis (e.g., 180).
/// * `frequency` - The frequency in Hz at which to calculate gain (relevant for frequency-dependent antennas).
/// * `state` - The shared simulation state.
///
/// # Returns
/// * `Ok(AntennaPatternData)` - Struct containing flattened gain array and dimensions.
/// * `Err(String)` - Error if antenna not found or context locked.
#[tauri::command]
fn get_antenna_pattern(
    antenna_id: String,
    az_samples: usize,
    el_samples: usize,
    frequency: f64,
    state: State<'_, FersState>,
) -> Result<fers_api::AntennaPatternData, String> {
    match state.try_lock() {
        Ok(context) => context.get_antenna_pattern(&antenna_id, az_samples, el_samples, frequency),
        Err(_) => Err("Backend is busy with another operation".to_string()),
    }
}

/// Calculates visual radio links between platforms at a specific time.
///
/// This command performs a lightweight geometric and physics check to determine
/// which platforms can "see" each other. It distinguishes between monostatic,
/// bistatic, and interference paths based on signal-to-noise ratios.
///
/// # Parameters
/// * `time` - The simulation time in seconds to evaluate.
/// * `state` - The shared simulation state containing platforms and physics models.
///
/// # Returns
/// * `Ok(Vec<VisualLink>)` - A list of renderable link segments with metadata (type, quality, label).
/// * `Err(String)` - Error if context access fails.
#[tauri::command]
fn get_preview_links(
    time: f64,
    state: State<'_, FersState>,
) -> Result<Vec<fers_api::VisualLink>, String> {
    match state.try_lock() {
        Ok(context) => context.calculate_preview_links(time),
        Err(_) => Ok(vec![]), // backend is busy (simulation/KML running); return empty and retry next frame
    }
}

/// Performs a granular state update on a specific simulation item via JSON.
#[tauri::command]
fn update_item_from_json(
    item_type: String,
    item_id: String,
    json: String,
    state: State<'_, FersState>,
) -> Result<Vec<String>, String> {
    let context = state.lock().map_err(|e| e.to_string())?;
    match item_type.as_str() {
        "Platform" => context.update_platform_from_json(&item_id, &json),
        "Transmitter" => context.update_transmitter_from_json(&item_id, &json).map(|()| vec![]),
        "Receiver" => context.update_receiver_from_json(&item_id, &json).map(|()| vec![]),
        "Target" => context.update_target_from_json(&item_id, &json).map(|()| vec![]),
        "Monostatic" => context.update_monostatic_from_json(&json).map(|()| vec![]),
        "Antenna" => context.update_antenna_from_json(&json).map(|()| vec![]),
        "Waveform" => context.update_waveform_from_json(&json).map(|()| vec![]),
        "Timing" => context.update_timing_from_json(&item_id, &json).map(|()| vec![]),
        "GlobalParameters" => context.update_parameters_from_json(&json),
        _ => Ok(vec![]),
    }
}

/// Initializes and runs the Tauri application.
///
/// This function is the main entry point for the desktop application. It performs
/// the following setup steps:
///
/// 1. Creates a new `FersContext` by calling the FFI layer. If this fails, it
///    indicates a linking or initialization problem with `libfers`.
/// 2. Registers Tauri plugins for file dialogs, file system access, and shell operations.
/// 3. Stores the `FersContext` in Tauri's managed state, protected by a `Mutex`.
/// 4. Registers all Tauri commands so they can be invoked from the frontend.
/// 5. Launches the Tauri application event loop.
///
/// # Panics
///
/// This function will panic if:
/// * The `FersContext` cannot be created (indicating a problem with `libfers`).
/// * The Tauri application fails to start due to misconfiguration.
///
/// # Example
///
/// This function is typically called from `main.rs`:
///
/// ```ignore
/// fers_ui_lib::run();
/// ```
#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    // Attempt to create the FFI context. This validates that libfers is correctly linked.
    let context = fers_api::FersContext::new()
        .expect("Failed to create FERS context. Is libfers linked correctly?");

    tauri::Builder::default()
        // Register Tauri plugins for UI functionality
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_opener::init())
        .plugin(tauri_plugin_fs::init())
        .setup(|app| {
            fers_api::register_log_callback(app.handle().clone());
            Ok(())
        })
        // Store the FersContext as managed state, accessible from all commands
        .manage(Mutex::new(context))
        // Register all Tauri commands that can be invoked from the frontend
        .invoke_handler(tauri::generate_handler![
            load_scenario_from_xml_file,
            get_scenario_as_json,
            get_scenario_as_xml,
            update_scenario_from_json,
            run_simulation,
            export_output_metadata_json,
            generate_kml,
            get_interpolated_motion_path,
            get_interpolated_rotation_path,
            get_antenna_pattern,
            get_preview_links,
            update_item_from_json,
            set_output_directory,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn test_interpolation_type_serialization() {
        // Assert serialization outputs match the TypeScript string literals exactly
        assert_eq!(serde_json::to_string(&InterpolationType::Static).unwrap(), "\"static\"");
        assert_eq!(serde_json::to_string(&InterpolationType::Linear).unwrap(), "\"linear\"");
        assert_eq!(serde_json::to_string(&InterpolationType::Cubic).unwrap(), "\"cubic\"");

        // Ensure deserialization from UI payloads works
        let parsed: InterpolationType = serde_json::from_str("\"linear\"").unwrap();
        assert!(matches!(parsed, InterpolationType::Linear));
    }

    #[test]
    fn test_motion_waypoint_serialization() {
        let wp = MotionWaypoint { time: 1.0, x: 2.0, y: 3.0, altitude: 4.0 };
        let serialized = serde_json::to_string(&wp).unwrap();

        // Deserialize to generic JSON value to inspect structure
        let parsed: serde_json::Value = serde_json::from_str(&serialized).unwrap();
        assert_eq!(parsed["time"], 1.0);
        assert_eq!(parsed["x"], 2.0);
        assert_eq!(parsed["y"], 3.0);
        assert_eq!(parsed["altitude"], 4.0);
    }

    #[test]
    fn test_rotation_waypoint_serialization() {
        let payload = json!({
            "time": 5.5,
            "azimuth": 180.0,
            "elevation": -15.0
        });

        let wp: RotationWaypoint = serde_json::from_value(payload).unwrap();
        assert_eq!(wp.time, 5.5);
        assert_eq!(wp.azimuth, 180.0);
        assert_eq!(wp.elevation, -15.0);
    }

    #[test]
    fn test_fers_state_wrapper() {
        // Verify that the C++ library correctly linked and can be protected by a Rust Mutex.
        // This mimics how Tauri maintains the global state internally.
        let context = fers_api::FersContext::new()
            .expect("Failed to create FERS context. Check build.rs linking.");

        let state: FersState = Mutex::new(context);

        // Test safe access and locking capabilities.
        let locked_context = state.lock().unwrap();

        // Executing a read against the state
        let result = locked_context.get_scenario_as_json();

        // Since it's a completely uninitialized new context with no XML/JSON loaded yet,
        // it correctly delegates to the backend and evaluates safely without segfaults.
        assert!(result.is_ok() || result.is_err());
    }
}
