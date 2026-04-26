#ifndef FERS_API_H
#define FERS_API_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief An opaque handle to an in-memory FERS simulation context.
 *
 * This pointer represents a live simulation instance. It is an "opaque" type,
 * meaning its internal structure is hidden from the API consumer. This design
 * is intentional to decouple the client from the C++ implementation details of
 * the simulator core. It ensures that the client code (e.g., in Rust) will not
 * break even if the internal C++ `FersContext` class is changed, maintaining ABI
 * stability across library versions.
 *
 * @note The client is responsible for managing the lifetime of this handle via
 *       `fers_context_create()` and `fers_context_destroy()`.
 */
struct fers_context;

typedef struct fers_context fers_context_t; // NOLINT(*-use-using)

/**
 * @brief A function pointer type for progress reporting callbacks.
 *
 * This callback can be implemented by the client to receive progress updates
 * during long-running operations like `fers_run_simulation`.
 *
 * @param message A descriptive message about the current operation.
 * @param current The current progress step.
 * @param total The total number of steps for the operation.
 * @param user_data An opaque pointer passed back to the caller, useful for
 *                  maintaining state (e.g., a class instance or application handle).
 */
typedef void (*fers_progress_callback_t)(const char* message, int current, int total,
										 void* user_data); // NOLINT(*-use-using)


// --- Context Lifecycle ---

/**
 * @brief Creates a new FERS simulation context.
 *
 * Allocates and initializes a new, empty simulation context in memory. This
 * context serves as the container for a scenario loaded via one of the
 * `fers_load_...` or `fers_update_...` functions.
 *
 * @note In a C API, RAII is not available. The caller is responsible for
 *       destroying the returned context using `fers_context_destroy()` to
 *       prevent resource leaks.
 *
 * @return A non-null opaque pointer (handle) to the simulation context on success.
 *         Returns NULL on failure (e.g., out of memory).
 */
fers_context_t* fers_context_create();

/**
 * @brief Destroys a FERS simulation context and releases all associated memory.
 *
 * This function must be called for every context created by `fers_context_create()`
 * to ensure proper cleanup of the underlying C++ objects. Accessing the context
 * handle after calling this function results in undefined behavior.
 *
 * @param context A valid pointer to a `fers_context_t` handle. If context is NULL,
 *                the function performs no action for safety and does not set an error.
 */
void fers_context_destroy(fers_context_t* context);

/**
 * @brief Sets the output directory for simulation results.
 * @param context A valid `fers_context_t` handle.
 * @param out_dir A null-terminated UTF-8 string for the output directory path.
 * @return 0 on success, non-zero on error.
 */
int fers_set_output_directory(fers_context_t* context, const char* out_dir);

/**
 * @brief Log levels for the FERS library.
 */
typedef enum // NOLINT(*-use-using)
{
	FERS_LOG_TRACE,
	FERS_LOG_DEBUG,
	FERS_LOG_INFO,
	FERS_LOG_WARNING,
	FERS_LOG_ERROR,
	FERS_LOG_FATAL,
	FERS_LOG_OFF
} fers_log_level_t;

/**
 * @brief A function pointer type for receiving formatted log lines.
 *
 * @param level The severity level of the log line.
 * @param line The full formatted log line, without a trailing newline.
 * @param user_data An opaque pointer passed back to the caller.
 */
typedef void (*fers_log_callback_t)(fers_log_level_t level, const char* line, void* user_data); // NOLINT(*-use-using)

/**
 * @brief Configures the internal logger.
 * @param level The minimum severity level to log.
 * @param log_file_path Optional path to a log file. Pass NULL to disable file logging.
 * @return 0 on success, non-zero on error.
 */
int fers_configure_logging(fers_log_level_t level, const char* log_file_path);

/**
 * @brief Returns the library version string.
 *
 * The returned pointer remains valid for the lifetime of the process and must
 * not be freed by the caller.
 */
const char* fers_get_version(void);

/**
 * @brief Returns the current internal logger level.
 */
fers_log_level_t fers_get_log_level();

/**
 * @brief Registers a callback for formatted log lines.
 * Pass NULL as callback to disable log callbacks.
 */
void fers_set_log_callback(fers_log_callback_t callback, void* user_data);

/**
 * @brief Submits a log message to the library's unified logging system.
 * This ensures CLI messages match the format (timestamps, alignment) of library messages.
 */
void fers_log(fers_log_level_t level, const char* message);

/**
 * @brief Sets the number of worker threads for the simulation.
 * @param num_threads The number of threads to use.
 * @return 0 on success, non-zero on error. This function clears any previous
 *         thread-local error message at entry, like the other fallible API calls.
 */
int fers_set_thread_count(unsigned num_threads);

// --- Scenario Loading & Serialization ---

/**
 * @brief Loads a scenario into the context from a FERS XML file.
 *
 * This is the standard method for initializing a simulation context from a file on
 * disk. It is essential for interoperability with the CLI and legacy workflows that
 * rely on the FERS XML format.
 *
 * @param context A valid `fers_context_t` handle.
 * @param xml_filepath A null-terminated UTF-8 string for the input XML file path.
 * @param validate A boolean (0 or 1) indicating whether to validate the XML
 *                 against the embedded FERS schema. Validation is recommended to
 *                 ensure scenario correctness.
 * @return 0 on success, a non-zero error code on failure. Use
 *         `fers_get_last_error_message()` to retrieve error details.
 */
int fers_load_scenario_from_xml_file(fers_context_t* context, const char* xml_filepath, int validate);

/**
 * @brief Loads a scenario into the context from a FERS XML string.
 *
 * This function provides a way to load a scenario from an in-memory string,
 * avoiding file I/O. It is useful for test harnesses or for UIs that manage
 * scenarios as text content before parsing.
 *
 * @param context A valid `fers_context_t` handle.
 * @param xml_content A null-terminated UTF-8 string containing the FERS scenario in XML format.
 * @param validate A boolean (0 or 1) indicating whether to validate the XML
 *                 against the embedded FERS schema.
 * @return 0 on success, a non-zero error code on failure. Use
 *         `fers_get_last_error_message()` to retrieve error details.
 */
int fers_load_scenario_from_xml_string(fers_context_t* context, const char* xml_content, int validate);

/**
 * @brief Serializes the current simulation scenario into a JSON string.
 *
 * This function is the primary method for the UI to retrieve the full state of
 * the simulation. JSON is used as the interchange format because it is lightweight,
 * human-readable, and natively supported by web technologies, making it trivial
 * to parse and use in the React/TypeScript frontend.
 *
 * @note Memory Management: The returned string is allocated by this library and
 *       its ownership is transferred to the caller. It is crucial to free this
 *       string using `fers_free_string()` to prevent memory leaks.
 *
 * @param context A valid `fers_context_t` handle.
 * @return A dynamically allocated, null-terminated C-string containing the
 *         JSON representation of the scenario. Returns NULL on failure.
 */
char* fers_get_scenario_as_json(fers_context_t* context);

/**
 * @brief Serializes the current simulation scenario into a FERS XML string.
 *
 * This function enables exporting the in-memory state back into the standard FERS
 * XML file format. This is essential for interoperability with legacy tools and
 * for allowing a user to save a scenario that was created or modified in the UI.
 *
 * @note Memory Management: The returned string is dynamically allocated and
 *       its ownership is transferred to the caller. It must be freed using
 *       `fers_free_string()` to prevent memory leaks.
 *
 * @param context A valid `fers_context_t` handle.
 * @return A dynamically allocated, null-terminated C-string containing the
 *         XML representation of the scenario. Returns NULL on failure.
 */
char* fers_get_scenario_as_xml(fers_context_t* context);

/**
 * @brief Returns JSON metadata for the most recent simulation output files.
 *
 * The returned JSON describes generated HDF5 file structure and sample ranges.
 * The caller owns the returned string and must free it with `fers_free_string`.
 * If no simulation has completed yet, the returned JSON contains an empty `files`
 * array.
 *
 * @param context A valid `fers_context_t` handle.
 * @return A heap-allocated JSON string, or NULL on error.
 */
char* fers_get_last_output_metadata_json(fers_context_t* context);

/**
 * @brief Returns a JSON projection of simulation startup memory and HDF5 payload size.
 *
 * The projection includes phase-noise lookup memory, streaming I/Q buffer memory,
 * rendered HDF5 payload size, current resident memory not attributed to streaming
 * I/Q buffers where available, and an aggregate projected total. The caller owns
 * the returned string and must free it with `fers_free_string`.
 *
 * @param context A valid `fers_context_t` handle.
 * @return A heap-allocated JSON string, or NULL on error.
 */
char* fers_get_memory_projection_json(fers_context_t* context);

/**
 * @brief Updates the simulation scenario from a JSON string.
 *
 * This is the primary method for the UI to push its state back to the C++
 * core. It performs a full replacement of the existing scenario.
 *
 * @param context A valid `fers_context_t` handle.
 * @param scenario_json A null-terminated UTF-8 string containing the FERS scenario in JSON format.
 * @return 0 on success.
 *         1 on generic logic error.
 *         2 on JSON parsing/schema validation error.
 *         Use `fers_get_last_error_message()` to retrieve error details.
 */
int fers_update_scenario_from_json(fers_context_t* context, const char* scenario_json);

/**
 * @brief Updates a single platform's paths and name from JSON without full context recreation.
 * @param context A valid `fers_context_t` handle.
 * @param id The unique ID of the platform.
 * @param json The JSON string for the platform.
 * @return 0 on success, non-zero on failure.
 */
int fers_update_platform_from_json(fers_context_t* context, uint64_t id, const char* json);

/**
 * @brief Updates the global simulation parameters from JSON without full context recreation.
 * @param context A valid `fers_context_t` handle.
 * @param json The JSON string for the parameters.
 * @return 0 on success, non-zero on failure.
 */
int fers_update_parameters_from_json(fers_context_t* context, const char* json);

/**
 * @brief Updates a single antenna from JSON without full context recreation.
 * @param context A valid `fers_context_t` handle.
 * @param json The JSON string for the antenna.
 * @return 0 on success, non-zero on failure.
 */
int fers_update_antenna_from_json(fers_context_t* context, const char* json);

/**
 * @brief Updates a single waveform from JSON without full context recreation.
 * @param context A valid `fers_context_t` handle.
 * @param json The JSON string for the waveform.
 * @return 0 on success, non-zero on failure.
 */
int fers_update_waveform_from_json(fers_context_t* context, const char* json);

/**
 * @brief Updates a single transmitter from JSON without full context recreation.
 * @param context A valid `fers_context_t` handle.
 * @param id The unique ID of the transmitter.
 * @param json The JSON string for the transmitter.
 * @return 0 on success, non-zero on failure.
 */
int fers_update_transmitter_from_json(fers_context_t* context, uint64_t id, const char* json);

/**
 * @brief Updates a single receiver from JSON without full context recreation.
 * @param context A valid `fers_context_t` handle.
 * @param id The unique ID of the receiver.
 * @param json The JSON string for the receiver.
 * @return 0 on success, non-zero on failure.
 */
int fers_update_receiver_from_json(fers_context_t* context, uint64_t id, const char* json);

/**
 * @brief Updates a single target from JSON without full context recreation.
 * @param context A valid `fers_context_t` handle.
 * @param id The unique ID of the target.
 * @param json The JSON string for the target.
 * @return 0 on success, non-zero on failure.
 */
int fers_update_target_from_json(fers_context_t* context, uint64_t id, const char* json);

/**
 * @brief Updates a monostatic radar from JSON without full context recreation.
 * @param context A valid `fers_context_t` handle.
 * @param json The JSON string for the monostatic component.
 * @return 0 on success, non-zero on failure.
 */
int fers_update_monostatic_from_json(fers_context_t* context, const char* json);

/**
 * @brief Updates a single timing source from JSON without full context recreation.
 * @param context A valid `fers_context_t` handle.
 * @param id The unique ID of the timing source.
 * @param json The JSON string for the timing source.
 * @return 0 on success, non-zero on failure.
 */
int fers_update_timing_from_json(fers_context_t* context, uint64_t id, const char* json);

// --- Error Handling ---

/**
 * @brief Retrieves the last error message that occurred on the current thread.
 *
 * Because C++ exceptions cannot safely propagate across the FFI boundary into
 * other languages, this function provides the standard C-style error reporting
 * mechanism. The error state is stored in a thread-local variable to ensure that
 * concurrent API calls from different threads do not overwrite each other's error
 * messages. The error is cleared at the start of each fallible API call.
 *
 * @note Memory Management: The returned string's ownership is transferred to the
 *       caller. It MUST be freed using `fers_free_string()` to prevent memory leaks.
 *
 * @return A dynamically allocated, null-terminated C-string containing the last
 *         error message, or NULL if no error has occurred.
 */
char* fers_get_last_error_message();

/**
 * @brief Returns the last deduplicated rotation-unit warning list for the calling thread as JSON.
 *
 * The returned value is a JSON array of strings. It is populated by successful XML/JSON
 * load and update calls that detect suspicious rotation values. The caller owns the string
 * and must free it with `fers_free_string()`.
 *
 * @return A dynamically allocated JSON array string, or NULL if no warnings are available.
 */
char* fers_get_last_warning_messages_json();

/**
 * @brief Frees a string that was allocated and returned by the libfers API.
 *
 * This function must be used to release memory for any string returned by
 * functions like `fers_get_scenario_as_json` or `fers_get_last_error_message`.
 * It exists to ensure that the memory deallocation mechanism (`free`) matches the
 * allocation mechanism (`strdup`/`malloc`) used within the C++ library, preventing
 * potential crashes from mismatched allocators across language boundaries.
 *
 * @param str A pointer to the string to be freed. If str is NULL, no action is taken.
 */
void fers_free_string(char* str);


// --- Simulation Execution ---

/**
 * @brief Runs the simulation defined in the provided context.
 *
 * This function is synchronous and will block the calling thread until the
 * simulation is complete. This design keeps the API simple. For use in a
 * responsive UI, it is the responsibility of the caller (e.g., the Tauri
 * backend) to invoke this function on a separate worker thread to avoid
 * freezing the user interface.
 *
 * @param context A valid `fers_context_t` handle containing a loaded scenario.
 * @param callback A function pointer to a progress callback. Can be NULL.
 * @param user_data An opaque pointer passed to the callback function.
 * @return 0 on success, a non-zero error code on failure. Use
 *         `fers_get_last_error_message()` to retrieve error details.
 */
int fers_run_simulation(fers_context_t* context, fers_progress_callback_t callback, void* user_data);


// --- Utility Functions ---

/**
 * @brief Generates a KML file for visualizing the scenario in the context.
 *
 * This utility exists to provide a simple, out-of-the-box method for users to
 * validate and visualize the geographic layout and motion paths of their
 * scenarios in common external tools like Google Earth.
 *
 * @param context A valid `fers_context_t` handle containing a loaded scenario.
 * @param output_kml_filepath A null-terminated UTF-8 string for the output KML file path.
 * @return 0 on success, a non-zero error code on failure. Use
 *         `fers_get_last_error_message()` to retrieve error details.
 */
int fers_generate_kml(const fers_context_t* context, const char* output_kml_filepath);

// --- Antenna Pattern Utilities ---

/**
 * @brief Represents a sampled 2D antenna gain pattern.
 * Gains are in linear scale (not dB), normalized to the antenna's peak gain.
 * The data is structured as a flat array in row-major order (elevation rows, then azimuth columns).
 * @note The `gains` array must be freed using `fers_free_antenna_pattern_data`.
 */
typedef struct // NOLINT(*-use-using)
{
	double* gains; // Flat array of gain values [el_count * az_count]
	size_t az_count; // Number of samples along azimuth (-180 to +180 deg)
	size_t el_count; // Number of samples along elevation (-90 to +90 deg)
	double max_gain; // The peak gain found in the pattern (linear scale)
} fers_antenna_pattern_data_t;


/**
 * @brief Samples the gain pattern of a specified antenna and provides the data.
 *
 * This function calculates the antenna's far-field gain at a specified resolution
 * over the full sphere of directions (azimuth and elevation). The resulting gain
 * values are linear (not in dB) and normalized relative to the pattern's peak gain.
 * This is a stateless utility useful for UI previews and analysis.
 *
 * @param context A valid `fers_context_t` handle containing a loaded scenario with the antenna.
 * @param antenna_id The unique ID of the antenna asset to sample.
 * @param az_samples The desired number of sample points along the azimuth axis.
 *                   Must be at least 2 to span the full azimuth range.
 * @param el_samples The desired number of sample points along the elevation axis.
 *                   Must be at least 2 to span the full elevation range.
 * @param frequency_hz The frequency in Hz to use for gain calculation (affects aperture antennas).
 * @return A pointer to a `fers_antenna_pattern_data_t` struct containing the results.
 *         Returns NULL on failure (e.g., antenna not found). The caller owns the
 *         returned struct and must free it with `fers_free_antenna_pattern_data`.
 */
fers_antenna_pattern_data_t* fers_get_antenna_pattern(const fers_context_t* context, uint64_t antenna_id,
													  size_t az_samples, size_t el_samples, double frequency_hz);

/**
 * @brief Frees the memory allocated for an antenna pattern data structure.
 * @param data A pointer to the `fers_antenna_pattern_data_t` struct to free.
 */
void fers_free_antenna_pattern_data(fers_antenna_pattern_data_t* data);


// --- Path Interpolation Utilities ---

/**
 * @brief Defines the interpolation methods available for path generation.
 * This enum provides a language-agnostic way to specify the desired
 * interpolation algorithm when calling the path generation functions.
 */
typedef enum // NOLINT(*-use-using)
{
	FERS_INTERP_STATIC,
	FERS_INTERP_LINEAR,
	FERS_INTERP_CUBIC
} fers_interp_type_t;

/**
 * @brief Units used for external rotation angles and rates.
 */
typedef enum // NOLINT(*-use-using)
{
	FERS_ANGLE_UNIT_DEG,
	FERS_ANGLE_UNIT_RAD
} fers_angle_unit_t;

/**
 * @brief Represents a single waypoint for a motion path.
 * Coordinates are in the scenario's defined coordinate system (e.g., ENU meters).
 */
typedef struct // NOLINT(*-use-using)
{
	double time; /**< Time in seconds. */
	double x; /**< X coordinate in meters (East in ENU). */
	double y; /**< Y coordinate in meters (North in ENU). */
	double z; /**< Z coordinate in meters (Up/Altitude in ENU). */
} fers_motion_waypoint_t;

/**
 * @brief Represents a single waypoint for a rotation path.
 * Angles are in compass units (CW from North) selected by the caller.
 */
typedef struct // NOLINT(*-use-using)
{
	double time; /**< Time in seconds. */
	double azimuth; /**< Azimuth angle in compass units (0=North, pi/2 or 90=East). */
	double elevation; /**< Elevation angle in compass units (positive up). */
} fers_rotation_waypoint_t;

/**
 * @brief Represents a single interpolated point on a motion path.
 * Includes position and velocity in the scenario's coordinate frame.
 */
typedef struct // NOLINT(*-use-using)
{
	double x; /**< X position in meters. */
	double y; /**< Y position in meters. */
	double z; /**< Z position in meters. */
	double vx; /**< X velocity in m/s. */
	double vy; /**< Y velocity in m/s. */
	double vz; /**< Z velocity in m/s. */
} fers_interpolated_point_t;

/**
 * @brief Represents a single interpolated point on a rotation path.
 * Angles are in compass units (CW from North) selected by the caller.
 */
typedef struct // NOLINT(*-use-using)
{
	double azimuth; /**< Azimuth angle in compass units. */
	double elevation; /**< Elevation angle in compass units. */
} fers_interpolated_rotation_point_t;


/**
 * @brief A container for an array of interpolated motion path points.
 * @note The `points` array must be freed using `fers_free_interpolated_motion_path`.
 */
typedef struct // NOLINT(*-use-using)
{
	fers_interpolated_point_t* points;
	size_t count;
} fers_interpolated_path_t;

/**
 * @brief A container for an array of interpolated rotation path points.
 * @note The `points` array must be freed using `fers_free_interpolated_rotation_path`.
 */
typedef struct // NOLINT(*-use-using)
{
	fers_interpolated_rotation_point_t* points;
	size_t count;
} fers_interpolated_rotation_path_t;


/**
 * @brief Calculates an interpolated motion path from a set of waypoints.
 * This function is a stateless utility that computes the path without needing a
 * full simulation context. It is useful for UI previews.
 *
 * @param waypoints An array of `fers_motion_waypoint_t` structs.
 * @param waypoint_count The number of waypoints in the array.
 * @param interp_type The interpolation algorithm to use.
 * @param num_points The desired number of points in the output interpolated path.
 * @return A pointer to a `fers_interpolated_path_t` struct containing the results.
 *         Returns NULL on failure. The caller owns the returned struct and must
 *         free it with `fers_free_interpolated_motion_path`.
 */
fers_interpolated_path_t* fers_get_interpolated_motion_path(const fers_motion_waypoint_t* waypoints,
															size_t waypoint_count, fers_interp_type_t interp_type,
															size_t num_points);

/**
 * @brief Frees the memory allocated for an interpolated motion path.
 * @param path A pointer to the `fers_interpolated_path_t` struct to free.
 */
void fers_free_interpolated_motion_path(fers_interpolated_path_t* path);

/**
 * @brief Calculates an interpolated rotation path from a set of waypoints.
 * This function is a stateless utility for UI previews.
 *
 * @param waypoints An array of `fers_rotation_waypoint_t` structs.
 * @param waypoint_count The number of waypoints in the array.
 * @param interp_type The interpolation algorithm to use (STATIC, LINEAR, CUBIC).
 * @param angle_unit The unit used by the waypoint angles and desired output angles.
 * @param num_points The desired number of points in the output interpolated path.
 * @return A pointer to a `fers_interpolated_rotation_path_t` struct containing the results.
 *         Returns NULL on failure. The caller owns the returned struct and must
 *         free it with `fers_free_interpolated_rotation_path`.
 */
fers_interpolated_rotation_path_t* fers_get_interpolated_rotation_path(const fers_rotation_waypoint_t* waypoints,
																	   size_t waypoint_count,
																	   fers_interp_type_t interp_type,
																	   fers_angle_unit_t angle_unit, size_t num_points);

/**
 * @brief Frees the memory allocated for an interpolated rotation path.
 * @param path A pointer to the `fers_interpolated_rotation_path_t` struct to free.
 */
void fers_free_interpolated_rotation_path(fers_interpolated_rotation_path_t* path);

// --- Preview Link Calculation ---

/**
 * @brief Quality of the radio link based on SNR.
 */
typedef enum // NOLINT(*-use-using)
{
	FERS_LINK_STRONG, // SNR > 0 dB
	FERS_LINK_WEAK // SNR < 0 dB (Geometric possibility but sub-noise)
} fers_link_quality_t;

/**
 * @brief Type of visual link to render.
 */
typedef enum // NOLINT(*-use-using)
{
	FERS_LINK_MONOSTATIC, // Combined Tx/Rx path
	FERS_LINK_BISTATIC_TX_TGT, // Illuminator path
	FERS_LINK_BISTATIC_TGT_RX, // Scattered path
	FERS_LINK_DIRECT_TX_RX // Interference path
} fers_link_type_t;

/**
 * @brief Represents a single renderable line segment metadata.
 * Geometry is resolved client-side.
 */
typedef struct // NOLINT(*-use-using)
{
	fers_link_type_t type; /**< Type of the link (Monostatic, Bistatic, etc.). */
	fers_link_quality_t quality; /**< Signal quality (Strong/Weak). */
	char label[128]; /**< Pre-formatted label (e.g., "-95 dBm"). */
	uint64_t source_id; /**< ID of the source component (e.g. Transmitter). */
	uint64_t dest_id; /**< ID of the destination component (e.g. Receiver/Target). */
	uint64_t origin_id; /**< ID of the originating Transmitter (for scattered paths). */
	double rcs; /**< RCS in m^2 for this path. Negative if not applicable (e.g., non-monostatic links). */
	double actual_power_dbm; /**< Received power in dBm with actual RCS applied. -999 if not applicable. */
} fers_visual_link_t;

/**
 * @brief A container for a list of visual links.
 * @note The `links` array is owned by this struct and must be freed using `fers_free_preview_links`.
 */
typedef struct // NOLINT(*-use-using)
{
	fers_visual_link_t* links;
	size_t count;
} fers_visual_link_list_t;

/**
 * @brief Calculates visual links for a specific simulation time.
 * @param context The simulation context.
 * @param time The simulation time in seconds.
 * @return A pointer to a link list. Caller must free with fers_free_preview_links.
 */
fers_visual_link_list_t* fers_calculate_preview_links(const fers_context_t* context, double time);

/**
 * @brief Frees the memory allocated for a preview link list.
 * @param list The list to free.
 */
void fers_free_preview_links(fers_visual_link_list_t* list);

#ifdef __cplusplus
}
#endif

#endif
