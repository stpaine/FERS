// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

/**
 * @file main.cpp
 * @brief Entry point for the FERS command-line interface (CLI).
 *
 * This executable acts as a wrapper around the libfers core library. It parses
 * command-line arguments, uses the libfers C-API to load and run a simulation,
 * and reports progress to the console.
 */

#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <libfers/api.h>
#include <string>

#include "arg_parser.h"

// --- Shim to restore LOG() functionality via C-API ---
namespace logging
{
	inline std::string getLevelString(fers_log_level_t l)
	{
		switch (l)
		{
		case FERS_LOG_TRACE:
			return "TRACE";
		case FERS_LOG_DEBUG:
			return "DEBUG";
		case FERS_LOG_INFO:
			return "INFO";
		case FERS_LOG_WARNING:
			return "WARNING";
		case FERS_LOG_ERROR:
			return "ERROR";
		case FERS_LOG_FATAL:
			return "FATAL";
		default:
			return "UNKNOWN";
		}
	}
}

// Define LOG macro to format arguments and pass them to the library's logger
#define LOG(level, ...)                                                                                                \
	do                                                                                                                 \
	{                                                                                                                  \
		std::string _msg = std::format(__VA_ARGS__);                                                                   \
		fers_log(level, _msg.c_str());                                                                                 \
	}                                                                                                                  \
	while (0)
// -----------------------------------------------------

int main(const int argc, char* argv[])
{
	// Parse command-line arguments using the local arg parser
	const auto config_result = core::parseArguments(argc, argv);
	if (!config_result)
	{
		if (config_result.error() != "Help requested." && config_result.error() != "Version requested." &&
			config_result.error() != "No arguments provided.")
		{
			// Use basic stderr here because logging isn't configured yet
			std::cerr << "[ERROR] Argument parsing error: " << config_result.error() << std::endl;
			return 1;
		}
		return 0;
	}

	const auto& [script_file, log_level, num_threads, validate, log_file, generate_kml] = config_result.value();

	// Configure logging via the C API
	const char* log_file_ptr = log_file ? log_file->c_str() : nullptr;
	if (fers_configure_logging(log_level, log_file_ptr) != 0)
	{
		// If we can't configure logging, we must print to stderr manually
		char* err = fers_get_last_error_message();
		std::cerr << "[ERROR] Failed to configure logging: " << (err ? err : "Unknown error") << std::endl;
		fers_free_string(err);
		return 1;
	}

	LOG(FERS_LOG_INFO, "FERS CLI started. Using libfers backend.");

	LOG(FERS_LOG_DEBUG,
		"Running FERS with arguments: script_file={}, log_level={}, num_threads={}, validate={}, log_file={}",
		script_file, logging::getLevelString(log_level), num_threads, validate, log_file.value_or("None"));

	// Create a simulation context using the C-API
	fers_context_t* context = fers_context_create();
	if (!context)
	{
		LOG(FERS_LOG_FATAL, "Failed to create FERS simulation context.");
		return 1;
	}

	// Load the scenario from file via the C-API
	LOG(FERS_LOG_INFO, "Loading scenario from '{}'...", script_file);
	if (fers_load_scenario_from_xml_file(context, script_file.c_str(), validate ? 1 : 0) != 0)
	{
		char* err = fers_get_last_error_message();
		LOG(FERS_LOG_FATAL, "Failed to load scenario: {}", err ? err : "Unknown error");
		fers_free_string(err);
		fers_context_destroy(context);
		return 1;
	}

	if (generate_kml)
	{
		std::filesystem::path kml_output_path = script_file;
		kml_output_path.replace_extension(".kml");
		const std::string kml_output_file = kml_output_path.string();

		LOG(FERS_LOG_INFO, "Generating KML file for scenario: {}", kml_output_file);
		if (fers_generate_kml(context, kml_output_file.c_str()) == 0)
		{
			LOG(FERS_LOG_INFO, "KML file generated successfully: {}", kml_output_file);
		}
		else
		{
			char* err = fers_get_last_error_message();
			LOG(FERS_LOG_FATAL, "Failed to generate KML file: {}", err ? err : "Unknown error");
			fers_free_string(err);
		}

		fers_context_destroy(context);
		return 0; // Exit after generating KML
	}

	// Set thread count via the C-API
	if (fers_set_thread_count(num_threads) != 0)
	{
		char* err = fers_get_last_error_message();
		LOG(FERS_LOG_ERROR, "Failed to set number of threads: {}", err ? err : "Unknown error");
		fers_free_string(err);
	}

	// Run the simulation via the C-API
	LOG(FERS_LOG_INFO, "Starting simulation...");
	if (fers_run_simulation(context, nullptr, nullptr) != 0)
	{
		char* err = fers_get_last_error_message();
		LOG(FERS_LOG_FATAL, "Simulation run failed: {}", err ? err : "Unknown error");
		fers_free_string(err);
		fers_context_destroy(context);
		return 1;
	}
	LOG(FERS_LOG_INFO, "Simulation completed successfully.");

	fers_context_destroy(context);

	return 0;
}
