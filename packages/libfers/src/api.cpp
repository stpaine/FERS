// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

/**
 * @file api.cpp
 * @brief Implementation of the C-style FFI for the libfers core library.
 *
 * This file provides the C implementations for the functions declared in `api.h`.
 * It acts as the bridge between the C ABI and the C++ core, handling object
 * creation/destruction, exception catching, error reporting, and type casting.
 */

#include <algorithm>
#include <core/logging.h>
#include <core/parameters.h>
#include <core/sim_id.h>
#include <cstring>
#include <filesystem>
#include <format>
#include <functional>
#include <libfers/api.h>
#include <math/path.h>
#include <math/rotation_path.h>
#include <mutex>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/fers_context.h"
#include "core/memory_projection.h"
#include "core/sim_threading.h"
#include "core/thread_pool.h"
#include "fers_version.h"
#include "serial/json_serializer.h"
#include "serial/kml_generator.h"
#include "serial/rotation_angle_utils.h"
#include "serial/rotation_warning_utils.h"
#include "serial/xml_parser.h"
#include "serial/xml_serializer.h"
#include "signal/radar_signal.h"
#include "simulation/channel_model.h"

// The fers_context struct is defined here as an alias for our C++ class.
// This allows the C-API to return an opaque pointer, hiding the C++ implementation.
struct fers_context : public FersContext
{
};

// A thread-local error message string ensures that error details from one
// thread's API call do not interfere with another's. This is crucial for a
// thread-safe FFI layer.
thread_local std::string last_error_message;
thread_local std::vector<std::string> last_warning_messages;

/**
 * @brief Centralized exception handler for the C-API boundary.
 *
 * This function catches standard C++ exceptions, records their `what()` message
 * into the thread-local error storage, and logs the error. This prevents C++
 * exceptions from propagating across the FFI boundary, which would be undefined behavior.
 * @param e The exception that was caught.
 * @param function_name The name of the API function where the error occurred.
 */
static void handle_api_exception(const std::exception& e, const std::string& function_name)
{
	last_error_message = e.what();
	LOG(logging::Level::ERROR, "API Error in {}: {}", function_name, last_error_message);
}

static void begin_warning_capture() noexcept
{
	last_warning_messages.clear();
	serial::rotation_warning_utils::clear_captured_warnings();
}

static void complete_warning_capture()
{
	last_warning_messages = serial::rotation_warning_utils::take_captured_warnings();
}

static void discard_warning_capture() noexcept
{
	last_warning_messages.clear();
	serial::rotation_warning_utils::clear_captured_warnings();
}

extern "C" {

fers_context_t* fers_context_create()
{
	last_error_message.clear();
	discard_warning_capture();
	try
	{
		return new fers_context_t();
	}
	catch (const std::bad_alloc& e)
	{
		handle_api_exception(e, "fers_context_create");
		return nullptr;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_context_create");
		return nullptr;
	}
}

void fers_context_destroy(fers_context_t* context)
{
	if (context == nullptr)
	{
		return;
	}
	delete context;
}

// Helper to map C enum to internal C++ enum
static logging::Level map_api_log_level(fers_log_level_t level)
{
	switch (level)
	{
	case FERS_LOG_TRACE:
		return logging::Level::TRACE;
	case FERS_LOG_DEBUG:
		return logging::Level::DEBUG;
	case FERS_LOG_INFO:
		return logging::Level::INFO;
	case FERS_LOG_WARNING:
		return logging::Level::WARNING;
	case FERS_LOG_ERROR:
		return logging::Level::ERROR;
	case FERS_LOG_FATAL:
		return logging::Level::FATAL;
	case FERS_LOG_OFF:
		return logging::Level::OFF;
	default:
		return logging::Level::INFO;
	}
}

static fers_log_level_t map_internal_log_level(logging::Level level)
{
	switch (level)
	{
	case logging::Level::TRACE:
		return FERS_LOG_TRACE;
	case logging::Level::DEBUG:
		return FERS_LOG_DEBUG;
	case logging::Level::INFO:
		return FERS_LOG_INFO;
	case logging::Level::WARNING:
		return FERS_LOG_WARNING;
	case logging::Level::ERROR:
		return FERS_LOG_ERROR;
	case logging::Level::FATAL:
		return FERS_LOG_FATAL;
	case logging::Level::OFF:
		return FERS_LOG_OFF;
	default:
		return FERS_LOG_INFO;
	}
}

namespace
{
	std::mutex log_callback_mutex; ///< Guards C API log callback state.
	fers_log_callback_t log_callback = nullptr; ///< Registered C API log callback, if any.
	void* log_callback_user_data = nullptr; ///< Opaque user data passed to the registered log callback.

	/// Forwards an internal formatted log line to the registered C API callback.
	void forward_log_callback(const logging::Level level, const std::string& line, void* /*user_data*/)
	{
		fers_log_callback_t callback = nullptr;
		void* user_data = nullptr;

		{
			std::scoped_lock lock(log_callback_mutex);
			callback = log_callback;
			user_data = log_callback_user_data;
		}

		if (callback != nullptr)
		{
			callback(map_internal_log_level(level), line.c_str(), user_data);
		}
	}
}

int fers_configure_logging(fers_log_level_t level, const char* log_file_path)
{
	last_error_message.clear();
	try
	{
		logging::logger.setLevel(map_api_log_level(level));
		if ((log_file_path != nullptr) && ((*log_file_path) != 0))
		{
			auto result = logging::logger.logToFile(log_file_path);
			if (!result)
			{
				last_error_message = result.error();
				return 1;
			}
		}
		return 0;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_configure_logging");
		return 1;
	}
}

const char* fers_get_version(void) { return FERS_VERSION_STRING; }

fers_log_level_t fers_get_log_level() { return map_internal_log_level(logging::logger.getLevel()); }

void fers_set_log_callback(fers_log_callback_t callback, void* user_data)
{
	{
		std::scoped_lock lock(log_callback_mutex);
		log_callback = callback;
		log_callback_user_data = user_data;
	}

	logging::logger.setCallback(callback == nullptr ? nullptr : forward_log_callback, nullptr);
}

void fers_log(fers_log_level_t level, const char* message)
{
	if (message == nullptr)
		return;
	// We pass a default source_location because C-API calls don't provide C++ source info
	logging::logger.log(map_api_log_level(level), message, std::source_location::current());
}

int fers_set_thread_count(unsigned num_threads)
{
	last_error_message.clear();
	try
	{
		if (auto res = params::setThreads(num_threads); !res)
		{
			last_error_message = res.error();
			return 1;
		}
		return 0;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_set_thread_count");
		return 1;
	}
}

int fers_set_output_directory(fers_context_t* context, const char* out_dir)
{
	last_error_message.clear();
	if ((context == nullptr) || (out_dir == nullptr))
	{
		last_error_message = "Invalid arguments: context or out_dir is NULL.";
		LOG(logging::Level::ERROR, last_error_message);
		return -1;
	}
	auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		ctx->setOutputDir(out_dir);
		return 0;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_set_output_directory");
		return 1;
	}
}

int fers_load_scenario_from_xml_file(fers_context_t* context, const char* xml_filepath, const int validate)
{
	last_error_message.clear();
	begin_warning_capture();
	if ((context == nullptr) || (xml_filepath == nullptr))
	{
		last_error_message = "Invalid arguments: context or xml_filepath is NULL.";
		LOG(logging::Level::ERROR, last_error_message);
		discard_warning_capture();
		return -1;
	}

	auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		// Set default output directory to the scenario file's directory
		std::filesystem::path p(xml_filepath);
		auto parent = p.parent_path();
		if (parent.empty())
			parent = ".";
		ctx->setOutputDir(parent.string());

		serial::parseSimulation(xml_filepath, ctx->getWorld(), static_cast<bool>(validate), ctx->getMasterSeeder());

		// After parsing, seed the master random number generator. This is done
		// to ensure simulation reproducibility. If the scenario specifies a seed,
		// it is used; otherwise, a non-deterministic seed is generated so that
		// subsequent runs are unique by default.
		if (params::params.random_seed)
		{
			LOG(logging::Level::INFO, "Using master seed from scenario file: {}", *params::params.random_seed);
			ctx->getMasterSeeder().seed(*params::params.random_seed);
		}
		else
		{
			const auto seed = std::random_device{}();
			LOG(logging::Level::INFO, "No master seed provided in scenario. Using random_device seed: {}", seed);
			params::params.random_seed = seed;
			ctx->getMasterSeeder().seed(seed);
		}
		complete_warning_capture();
		return 0; // Success
	}
	catch (const std::exception& e)
	{
		discard_warning_capture();
		handle_api_exception(e, "fers_load_scenario_from_xml_file");
		return 1; // Error
	}
}

int fers_load_scenario_from_xml_string(fers_context_t* context, const char* xml_content, const int validate)
{
	last_error_message.clear();
	begin_warning_capture();
	if ((context == nullptr) || (xml_content == nullptr))
	{
		last_error_message = "Invalid arguments: context or xml_content is NULL.";
		LOG(logging::Level::ERROR, last_error_message);
		discard_warning_capture();
		return -1;
	}

	auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		serial::parseSimulationFromString(xml_content, ctx->getWorld(), static_cast<bool>(validate),
										  ctx->getMasterSeeder());

		// After parsing, seed the master random number generator. This ensures
		// that if the scenario provides a seed, the simulation will be
		// reproducible. If not, a random seed is used to ensure unique runs.
		if (params::params.random_seed)
		{
			LOG(logging::Level::INFO, "Using master seed from scenario string: {}", *params::params.random_seed);
			ctx->getMasterSeeder().seed(*params::params.random_seed);
		}
		else
		{
			const auto seed = std::random_device{}();
			LOG(logging::Level::INFO, "No master seed provided in scenario. Using random_device seed: {}", seed);
			params::params.random_seed = seed;
			ctx->getMasterSeeder().seed(seed);
		}

		complete_warning_capture();
		return 0; // Success
	}
	catch (const std::exception& e)
	{
		discard_warning_capture();
		handle_api_exception(e, "fers_load_scenario_from_xml_string");
		return 1; // Parsing or logic error
	}
}

char* fers_get_scenario_as_json(fers_context_t* context)
{
	last_error_message.clear();
	if (context == nullptr)
	{
		last_error_message = "Invalid context provided to fers_get_scenario_as_json.";
		LOG(logging::Level::ERROR, last_error_message);
		return nullptr;
	}

	const auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		const nlohmann::json j = serial::world_to_json(*ctx->getWorld());
		const std::string json_str = j.dump(2);
		// A heap-allocated copy of the string is returned. This is necessary
		// to transfer ownership of the memory across the FFI boundary to a
		// client that will free it using `fers_free_string`.
		return strdup(json_str.c_str());
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_get_scenario_as_json");
		return nullptr;
	}
}

char* fers_get_scenario_as_xml(fers_context_t* context)
{
	last_error_message.clear();
	if (context == nullptr)
	{
		last_error_message = "Invalid context provided to fers_get_scenario_as_xml.";
		LOG(logging::Level::ERROR, last_error_message);
		return nullptr;
	}

	const auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		const std::string xml_str = serial::world_to_xml_string(*ctx->getWorld());
		if (xml_str.empty())
		{
			throw std::runtime_error("XML serialization resulted in an empty string.");
		}
		// `strdup` is used to create a heap-allocated string that can be safely
		// passed across the FFI boundary. The client is responsible for freeing
		// this memory with `fers_free_string`.
		return strdup(xml_str.c_str());
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_get_scenario_as_xml");
		return nullptr;
	}
}

char* fers_get_last_output_metadata_json(fers_context_t* context)
{
	last_error_message.clear();
	if (context == nullptr)
	{
		last_error_message = "Invalid context provided to fers_get_last_output_metadata_json.";
		LOG(logging::Level::ERROR, last_error_message);
		return nullptr;
	}

	const auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		const std::string json_str = ctx->getLastOutputMetadataJson();
		return strdup(json_str.c_str());
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_get_last_output_metadata_json");
		return nullptr;
	}
}

char* fers_get_memory_projection_json(fers_context_t* context)
{
	last_error_message.clear();
	if (context == nullptr)
	{
		last_error_message = "Invalid context provided to fers_get_memory_projection_json.";
		LOG(logging::Level::ERROR, last_error_message);
		return nullptr;
	}

	auto* ctx = reinterpret_cast<FersContext*>(context);

	try
	{
		const auto projection = core::projectSimulationMemory(*ctx->getWorld());
		const std::string json_str = core::memoryProjectionToJsonString(projection);
		return strdup(json_str.c_str());
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_get_memory_projection_json");
		return nullptr;
	}
}

int fers_update_platform_from_json(fers_context_t* context, uint64_t id, const char* json)
{
	last_error_message.clear();
	begin_warning_capture();
	if ((context == nullptr) || (json == nullptr))
	{
		discard_warning_capture();
		return -1;
	}
	auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		auto* p = ctx->getWorld()->findPlatform(id);
		if (p == nullptr)
		{
			last_error_message = "Platform not found";
			discard_warning_capture();
			return 1;
		}
		auto j = nlohmann::json::parse(json);
		serial::update_platform_paths_from_json(j, p);
		if (j.contains("name"))
		{
			p->setName(j.at("name").get<std::string>());
		}
		complete_warning_capture();
		return 0;
	}
	catch (const std::exception& e)
	{
		discard_warning_capture();
		handle_api_exception(e, "fers_update_platform_from_json");
		return 1;
	}
}

int fers_update_parameters_from_json(fers_context_t* context, const char* json)
{
	last_error_message.clear();
	begin_warning_capture();
	if ((context == nullptr) || (json == nullptr))
	{
		discard_warning_capture();
		return -1;
	}
	auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		auto j = nlohmann::json::parse(json);
		serial::update_parameters_from_json(j, ctx->getMasterSeeder());
		complete_warning_capture();
		return 0;
	}
	catch (const std::exception& e)
	{
		discard_warning_capture();
		handle_api_exception(e, "fers_update_parameters_from_json");
		return 1;
	}
}

int fers_update_antenna_from_json(fers_context_t* context, const char* json)
{
	last_error_message.clear();
	if ((context == nullptr) || (json == nullptr))
		return -1;
	auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		auto j = nlohmann::json::parse(json);
		auto id = j.at("id").is_string() ? std::stoull(j.at("id").get<std::string>()) : j.at("id").get<uint64_t>();
		auto* ant = ctx->getWorld()->findAntenna(id);
		if (ant == nullptr)
		{
			last_error_message = "Antenna not found";
			return 1;
		}
		serial::update_antenna_from_json(j, ant, *ctx->getWorld());
		return 0;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_update_antenna_from_json");
		return 1;
	}
}

int fers_update_waveform_from_json(fers_context_t* context, const char* json)
{
	last_error_message.clear();
	if ((context == nullptr) || (json == nullptr))
		return -1;
	auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		auto j = nlohmann::json::parse(json);
		auto wf = serial::parse_waveform_from_json(j);
		if (wf)
		{
			ctx->getWorld()->replace(std::move(wf));
		}
		return 0;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_update_waveform_from_json");
		return 1;
	}
}

int fers_update_transmitter_from_json(fers_context_t* context, uint64_t id, const char* json)
{
	last_error_message.clear();
	if ((context == nullptr) || (json == nullptr))
		return -1;
	auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		auto* tx = ctx->getWorld()->findTransmitter(id);
		if (tx == nullptr)
		{
			last_error_message = "Transmitter not found";
			return 1;
		}
		auto j = nlohmann::json::parse(json);
		serial::update_transmitter_from_json(j, tx, *ctx->getWorld(), ctx->getMasterSeeder());
		return 0;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_update_transmitter_from_json");
		return 1;
	}
}

int fers_update_receiver_from_json(fers_context_t* context, uint64_t id, const char* json)
{
	last_error_message.clear();
	if ((context == nullptr) || (json == nullptr))
		return -1;
	auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		auto* rx = ctx->getWorld()->findReceiver(id);
		if (rx == nullptr)
		{
			last_error_message = "Receiver not found";
			return 1;
		}
		auto j = nlohmann::json::parse(json);
		serial::update_receiver_from_json(j, rx, *ctx->getWorld(), ctx->getMasterSeeder());
		return 0;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_update_receiver_from_json");
		return 1;
	}
}

int fers_update_target_from_json(fers_context_t* context, uint64_t id, const char* json)
{
	last_error_message.clear();
	if ((context == nullptr) || (json == nullptr))
		return -1;
	auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		auto* tgt = ctx->getWorld()->findTarget(id);
		if (tgt == nullptr)
		{
			last_error_message = "Target not found";
			return 1;
		}
		auto j = nlohmann::json::parse(json);
		serial::update_target_from_json(j, tgt, *ctx->getWorld(), ctx->getMasterSeeder());
		return 0;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_update_target_from_json");
		return 1;
	}
}

int fers_update_monostatic_from_json(fers_context_t* context, const char* json)
{
	last_error_message.clear();
	if ((context == nullptr) || (json == nullptr))
		return -1;
	auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		auto j = nlohmann::json::parse(json);
		uint64_t tx_id =
			j.at("tx_id").is_string() ? std::stoull(j.at("tx_id").get<std::string>()) : j.at("tx_id").get<uint64_t>();
		uint64_t rx_id =
			j.at("rx_id").is_string() ? std::stoull(j.at("rx_id").get<std::string>()) : j.at("rx_id").get<uint64_t>();
		auto* tx = ctx->getWorld()->findTransmitter(tx_id);
		auto* rx = ctx->getWorld()->findReceiver(rx_id);
		if ((tx == nullptr) || (rx == nullptr))
		{
			last_error_message = "Monostatic components not found";
			return 1;
		}
		serial::update_monostatic_from_json(j, tx, rx, *ctx->getWorld(), ctx->getMasterSeeder());
		return 0;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_update_monostatic_from_json");
		return 1;
	}
}

int fers_update_timing_from_json(fers_context_t* context, uint64_t id, const char* json)
{
	last_error_message.clear();
	if ((context == nullptr) || (json == nullptr))
		return -1;
	auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		if (ctx->getWorld()->findTiming(id) == nullptr)
		{
			last_error_message = "Timing not found";
			return 1;
		}
		auto j = nlohmann::json::parse(json);
		serial::update_timing_from_json(j, *ctx->getWorld(), id);
		return 0;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_update_timing_from_json");
		return 1;
	}
}

int fers_update_scenario_from_json(fers_context_t* context, const char* scenario_json)
{
	last_error_message.clear();
	begin_warning_capture();
	if ((context == nullptr) || (scenario_json == nullptr))
	{
		last_error_message = "Invalid arguments: context or scenario_json is NULL.";
		LOG(logging::Level::ERROR, last_error_message);
		discard_warning_capture();
		return -1;
	}

	auto* ctx = reinterpret_cast<FersContext*>(context);
	try
	{
		const nlohmann::json j = nlohmann::json::parse(scenario_json);
		serial::json_to_world(j, *ctx->getWorld(), ctx->getMasterSeeder());
		complete_warning_capture();

		return 0; // Success
	}
	catch (const nlohmann::json::exception& e)
	{
		// A specific catch block for JSON errors is used to provide more
		// detailed feedback to the client (e.g., the UI), which can help
		// developers diagnose schema or data format issues more easily.
		last_error_message = "JSON parsing/deserialization error: " + std::string(e.what());
		LOG(logging::Level::ERROR, "API Error in {}: {}", "fers_update_scenario_from_json", last_error_message);
		discard_warning_capture();
		return 2; // JSON error
	}
	catch (const std::exception& e)
	{
		discard_warning_capture();
		handle_api_exception(e, "fers_update_scenario_from_json");
		return 1; // Generic error
	}
}

char* fers_get_last_error_message()
{
	if (last_error_message.empty())
	{
		return nullptr; // No error to report
	}
	// `strdup` allocates with `malloc`, which is part of the C standard ABI,
	// making it safe to transfer ownership across the FFI boundary. The caller
	// must then free this memory using `fers_free_string`.
	return strdup(last_error_message.c_str());
}

char* fers_get_last_warning_messages_json()
{
	if (last_warning_messages.empty())
	{
		return nullptr;
	}

	const std::string warning_json = nlohmann::json(last_warning_messages).dump();
	last_warning_messages.clear();
	return strdup(warning_json.c_str());
}

void fers_free_string(char* str)
{
	if (str != nullptr)
	{
		free(str);
	}
}

int fers_run_simulation(fers_context_t* context, fers_progress_callback_t callback, void* user_data)
{
	last_error_message.clear();
	if (context == nullptr)
	{
		last_error_message = "Invalid context provided to fers_run_simulation.";
		LOG(logging::Level::ERROR, last_error_message);
		return -1;
	}

	auto* ctx = reinterpret_cast<FersContext*>(context);

	// Wrap the C-style callback in a std::function for easier use in C++.
	// This also handles the case where the callback is null.
	std::function<void(const std::string&, int, int)> progress_fn;
	if (callback != nullptr)
	{
		progress_fn = [callback, user_data](const std::string& msg, const int current, const int total)
		{ callback(msg.c_str(), current, total, user_data); };
	}

	try
	{
		pool::ThreadPool pool(params::renderThreads());

		ctx->clearLastOutputMetadata();
		const auto output_metadata = core::runEventDrivenSim(ctx->getWorld(), pool, progress_fn, ctx->getOutputDir());
		ctx->setLastOutputMetadata(output_metadata);

		return 0;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_run_simulation");
		return 1;
	}
}

int fers_generate_kml(const fers_context_t* context, const char* output_kml_filepath)
{
	last_error_message.clear();
	if ((context == nullptr) || (output_kml_filepath == nullptr))
	{
		last_error_message = "Invalid arguments: context or output_kml_filepath is NULL.";
		LOG(logging::Level::ERROR, last_error_message);
		return -1;
	}

	const auto* ctx = reinterpret_cast<const FersContext*>(context);

	try
	{
		if (serial::KmlGenerator::generateKml(*ctx->getWorld(), output_kml_filepath))
		{
			return 0; // Success
		}

		last_error_message = "KML generation failed for an unknown reason.";
		LOG(logging::Level::ERROR, last_error_message);
		return 2; // Generation failed
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_generate_kml");
		return 1; // Exception thrown
	}
}

// --- Helper to convert C-API enum to C++ enum ---
math::Path::InterpType to_cpp_interp_type(const fers_interp_type_t type)
{
	switch (type)
	{
	case FERS_INTERP_LINEAR:
		return math::Path::InterpType::INTERP_LINEAR;
	case FERS_INTERP_CUBIC:
		return math::Path::InterpType::INTERP_CUBIC;
	case FERS_INTERP_STATIC:
	default:
		return math::Path::InterpType::INTERP_STATIC;
	}
}

math::RotationPath::InterpType to_cpp_rot_interp_type(const fers_interp_type_t type)
{
	switch (type)
	{
	case FERS_INTERP_LINEAR:
		return math::RotationPath::InterpType::INTERP_LINEAR;
	case FERS_INTERP_CUBIC:
		return math::RotationPath::InterpType::INTERP_CUBIC;
	case FERS_INTERP_STATIC:
	default:
		return math::RotationPath::InterpType::INTERP_STATIC;
	}
}


fers_interpolated_path_t* fers_get_interpolated_motion_path(const fers_motion_waypoint_t* waypoints,
															const size_t waypoint_count,
															const fers_interp_type_t interp_type,
															const size_t num_points)
{
	last_error_message.clear();
	if ((waypoints == nullptr) || waypoint_count == 0 || num_points == 0)
	{
		last_error_message = "Invalid arguments: waypoints cannot be null and counts must be > 0.";
		LOG(logging::Level::ERROR, last_error_message);
		return nullptr;
	}
	if (interp_type == FERS_INTERP_CUBIC && waypoint_count < 2)
	{
		last_error_message = "Cubic interpolation requires at least 2 waypoints.";
		LOG(logging::Level::ERROR, last_error_message);
		return nullptr;
	}

	try
	{
		math::Path path;
		path.setInterp(to_cpp_interp_type(interp_type));

		for (size_t i = 0; i < waypoint_count; ++i)
		{
			math::Coord c;
			c.t = waypoints[i].time;
			c.pos.x = waypoints[i].x;
			c.pos.y = waypoints[i].y;
			c.pos.z = waypoints[i].z;
			path.addCoord(c);
		}

		path.finalize();

		auto* result_path = new fers_interpolated_path_t();
		result_path->points = new fers_interpolated_point_t[num_points];
		result_path->count = num_points;

		const double start_time = waypoints[0].time;
		const double end_time = waypoints[waypoint_count - 1].time;
		const double duration = end_time - start_time;

		// Handle static case separately
		if (waypoint_count < 2 || duration <= 0)
		{
			const math::Vec3 pos = path.getPosition(start_time);
			for (size_t i = 0; i < num_points; ++i)
			{
				result_path->points[i] = {pos.x, pos.y, pos.z, 0.0, 0.0, 0.0};
			}
			return result_path;
		}

		const double time_step =
			duration / static_cast<double>(num_points > 1 ? num_points - 1 : static_cast<size_t>(1));

		for (size_t i = 0; i < num_points; ++i)
		{
			const double t = start_time + static_cast<double>(i) * time_step;
			const math::Vec3 pos = path.getPosition(t);
			const math::Vec3 vel = path.getVelocity(t);
			result_path->points[i] = {pos.x, pos.y, pos.z, vel.x, vel.y, vel.z};
		}

		return result_path;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_get_interpolated_motion_path");
		return nullptr;
	}
}

void fers_free_interpolated_motion_path(fers_interpolated_path_t* path)
{
	if (path != nullptr)
	{
		delete[] path->points;
		delete path;
	}
}

fers_interpolated_rotation_path_t* fers_get_interpolated_rotation_path(const fers_rotation_waypoint_t* waypoints,
																	   const size_t waypoint_count,
																	   const fers_interp_type_t interp_type,
																	   const fers_angle_unit_t angle_unit,
																	   const size_t num_points)
{
	last_error_message.clear();
	last_warning_messages.clear();
	serial::rotation_warning_utils::clear_captured_warnings();
	if ((waypoints == nullptr) || waypoint_count == 0 || num_points == 0)
	{
		last_error_message = "Invalid arguments: waypoints cannot be null and counts must be > 0.";
		LOG(logging::Level::ERROR, last_error_message);
		return nullptr;
	}
	if (interp_type == FERS_INTERP_CUBIC && waypoint_count < 2)
	{
		last_error_message = "Cubic interpolation requires at least 2 waypoints.";
		LOG(logging::Level::ERROR, last_error_message);
		return nullptr;
	}

	try
	{
		const auto unit =
			angle_unit == FERS_ANGLE_UNIT_RAD ? params::RotationAngleUnit::Radians : params::RotationAngleUnit::Degrees;
		math::RotationPath path;
		path.setInterp(to_cpp_rot_interp_type(interp_type));

		for (size_t i = 0; i < waypoint_count; ++i)
		{
			serial::rotation_warning_utils::maybe_warn_about_rotation_value(
				waypoints[i].azimuth, unit, serial::rotation_warning_utils::ValueKind::Angle, "C-API",
				std::format("rotation waypoint {}", i), "azimuth");
			serial::rotation_warning_utils::maybe_warn_about_rotation_value(
				waypoints[i].elevation, unit, serial::rotation_warning_utils::ValueKind::Angle, "C-API",
				std::format("rotation waypoint {}", i), "elevation");
			path.addCoord(serial::rotation_angle_utils::external_rotation_to_internal(
				waypoints[i].azimuth, waypoints[i].elevation, waypoints[i].time, unit));
		}

		path.finalize();

		auto* result_path = new fers_interpolated_rotation_path_t();
		result_path->points = new fers_interpolated_rotation_point_t[num_points];
		result_path->count = num_points;

		const double start_time = waypoints[0].time;
		const double end_time = waypoints[waypoint_count - 1].time;
		const double duration = end_time - start_time;

		// Handle static case separately
		if (waypoint_count < 2 || duration <= 0)
		{
			const math::SVec3 rot = path.getPosition(start_time);
			for (size_t i = 0; i < num_points; ++i)
			{
				result_path->points[i] = fers_interpolated_rotation_point_t{
					serial::rotation_angle_utils::internal_azimuth_to_external(rot.azimuth, unit),
					serial::rotation_angle_utils::internal_elevation_to_external(rot.elevation, unit)};
			}
			return result_path;
		}

		const double time_step =
			duration / static_cast<double>(num_points > 1 ? num_points - 1 : static_cast<size_t>(1));

		for (size_t i = 0; i < num_points; ++i)
		{
			const double t = start_time + static_cast<double>(i) * time_step;
			const math::SVec3 rot = path.getPosition(t);

			result_path->points[i] = fers_interpolated_rotation_point_t{
				serial::rotation_angle_utils::internal_azimuth_to_external(rot.azimuth, unit),
				serial::rotation_angle_utils::internal_elevation_to_external(rot.elevation, unit)};
		}

		return result_path;
	}
	catch (const std::exception& e)
	{
		serial::rotation_warning_utils::clear_captured_warnings();
		handle_api_exception(e, "fers_get_interpolated_rotation_path");
		return nullptr;
	}
}

void fers_free_interpolated_rotation_path(fers_interpolated_rotation_path_t* path)
{
	if (path != nullptr)
	{
		delete[] path->points;
		delete path;
	}
}

// --- Antenna Pattern Implementation ---

fers_antenna_pattern_data_t* fers_get_antenna_pattern(const fers_context_t* context, const uint64_t antenna_id,
													  const size_t az_samples, const size_t el_samples,
													  const double frequency_hz)
{
	last_error_message.clear();
	if ((context == nullptr) || az_samples < 2 || el_samples < 2)
	{
		last_error_message = "Invalid arguments: context must be non-null and sample counts must be >= 2.";
		LOG(logging::Level::ERROR, last_error_message);
		return nullptr;
	}

	try
	{
		const auto* ctx = reinterpret_cast<const FersContext*>(context);
		antenna::Antenna* ant = ctx->getWorld()->findAntenna(static_cast<SimId>(antenna_id));

		if (ant == nullptr)
		{
			last_error_message = "Antenna ID '" + std::to_string(antenna_id) + "' not found in the world.";
			LOG(logging::Level::ERROR, last_error_message);
			return nullptr;
		}

		// TODO: Currently only using the first-found waveform. This is incorrect but also difficult to represent
		//		 correctly in scenarios with multiple waveforms as the gain for squarehorn and parabolic antennas
		// depends on 		 the wavelength. Hence a decision needs to be made about whether to return multiple patterns
		// per waveform or 		 have the user specify a representative wavelength in the UI per antenna.
		// Calculate wavelength from the provided frequency.
		// Default to 1GHz (0.3m) if frequency is invalid/zero, though the UI should prevent this
		// for antennas that strictly require it (Horn/Parabolic).
		RealType wavelength = 0.3;
		if (frequency_hz > 0.0)
		{
			wavelength = params::c() / frequency_hz;
		}

		auto* data = new fers_antenna_pattern_data_t();
		data->az_count = az_samples;
		data->el_count = el_samples;
		const size_t total_samples = az_samples * el_samples;
		data->gains = new double[total_samples];

		// The reference angle (boresight) is implicitly the local X-axis in the FERS engine.
		// We pass a zero rotation to get the gain relative to this boresight.
		const math::SVec3 ref_angle(1.0, 0.0, 0.0);
		double max_gain = 0.0;

		const RealType az_denominator = static_cast<RealType>(az_samples - 1);
		const RealType el_denominator = static_cast<RealType>(el_samples - 1);

		for (size_t i = 0; i < el_samples; ++i)
		{
			// Elevation from -PI/2 to PI/2
			const RealType elevation = (static_cast<RealType>(i) / el_denominator) * PI - (PI / 2.0);
			for (size_t j = 0; j < az_samples; ++j)
			{
				// Azimuth from -PI to PI
				const RealType azimuth = (static_cast<RealType>(j) / az_denominator) * 2.0 * PI - PI;
				const math::SVec3 sample_angle(1.0, azimuth, elevation);
				const RealType gain = ant->getGain(sample_angle, ref_angle, wavelength);
				data->gains[i * az_samples + j] = gain;
				max_gain = std::max(gain, max_gain);
			}
		}

		data->max_gain = max_gain;

		// Normalize the gains
		if (max_gain > 0)
		{
			for (size_t i = 0; i < total_samples; ++i)
			{
				data->gains[i] /= max_gain;
			}
		}

		return data;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_get_antenna_pattern");
		return nullptr;
	}
}

void fers_free_antenna_pattern_data(fers_antenna_pattern_data_t* data)
{
	if (data != nullptr)
	{
		delete[] data->gains;
		delete data;
	}
}

// --- Preview Link Calculation Implementation ---

fers_visual_link_list_t* fers_calculate_preview_links(const fers_context_t* context, const double time)
{
	last_error_message.clear();
	if (context == nullptr)
	{
		last_error_message = "Invalid context passed to fers_calculate_preview_links";
		LOG(logging::Level::ERROR, last_error_message);
		return nullptr;
	}

	try
	{
		const auto* ctx = reinterpret_cast<const FersContext*>(context);
		// Call the core physics logic in channel_model.cpp
		const auto cpp_links = simulation::calculatePreviewLinks(*ctx->getWorld(), time);

		// Convert C++ vector to C-API struct
		auto* result = new fers_visual_link_list_t();
		result->count = cpp_links.size();

		if (!cpp_links.empty())
		{
			result->links = new fers_visual_link_t[result->count];
			for (size_t i = 0; i < result->count; ++i)
			{
				const auto& src = cpp_links[i];
				auto& dst = result->links[i];

				// Map enums
				switch (src.type)
				{
				case simulation::LinkType::Monostatic:
					dst.type = FERS_LINK_MONOSTATIC;
					break;
				case simulation::LinkType::BistaticTxTgt:
					dst.type = FERS_LINK_BISTATIC_TX_TGT;
					break;
				case simulation::LinkType::BistaticTgtRx:
					dst.type = FERS_LINK_BISTATIC_TGT_RX;
					break;
				case simulation::LinkType::DirectTxRx:
					dst.type = FERS_LINK_DIRECT_TX_RX;
					break;
				}

				dst.quality = (src.quality == simulation::LinkQuality::Strong) ? FERS_LINK_STRONG : FERS_LINK_WEAK;

				// Safe string copy
				std::strncpy(dst.label, src.label.c_str(), sizeof(dst.label) - 1);
				dst.label[sizeof(dst.label) - 1] = '\0';

				dst.source_id = static_cast<uint64_t>(src.source_id);
				dst.dest_id = static_cast<uint64_t>(src.dest_id);
				dst.origin_id = static_cast<uint64_t>(src.origin_id);
				dst.rcs = src.rcs;
				dst.actual_power_dbm = src.actual_power_dbm;
				dst.display_value = src.display_value;
			}
		}
		else
		{
			result->links = nullptr;
		}
		return result;
	}
	catch (const std::exception& e)
	{
		handle_api_exception(e, "fers_calculate_preview_links");
		return nullptr;
	}
}

void fers_free_preview_links(fers_visual_link_list_t* list)
{
	if (list != nullptr)
	{
		delete[] list->links;
		delete list;
	}
}
}
