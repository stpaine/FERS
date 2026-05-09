// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

/**
 * @file cli_runner.cpp
 * @brief Testable implementation of the FERS command-line interface.
 */

#include "cli_runner.h"

#include <filesystem>
#include <format>
#include <iostream>
#include <libfers/api.h>
#include <string>
#include <utility>

#include "arg_parser.h"
#include "cli_paths.h"

namespace
{
	std::string getLevelString(fers_log_level_t level)
	{
		switch (level)
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

	template <typename... Args>
	void log(fers_log_level_t level, std::format_string<Args...> fmt, Args&&... args)
	{
		const std::string message = std::format(fmt, std::forward<Args>(args)...);
		fers_log(level, message.c_str());
	}
}

namespace core
{
	int runCli(const int argc, char* argv[])
	{
		const auto config_result = parseArguments(argc, argv);
		if (!config_result)
		{
			if (config_result.error() != "Help requested." && config_result.error() != "Version requested." &&
				config_result.error() != "No arguments provided.")
			{
				std::cerr << "[ERROR] Argument parsing error: " << config_result.error() << '\n';
				return 1;
			}
			return 0;
		}

		const auto& [script_file, log_level, num_threads, validate, log_file, generate_kml, kml_file, output_dir] =
			config_result.value();

		const char* log_file_ptr = log_file ? log_file->c_str() : nullptr;
		if (fers_configure_logging(log_level, log_file_ptr) != 0)
		{
			char* err = fers_get_last_error_message();
			std::cerr << "[ERROR] Failed to configure logging: " << ((err != nullptr) ? err : "Unknown error") << '\n';
			fers_free_string(err);
			return 1;
		}

		log(FERS_LOG_INFO, "FERS CLI started. Using libfers backend.");

		log(FERS_LOG_DEBUG,
			"Running FERS with arguments: script_file={}, log_level={}, num_threads={}, validate={}, log_file={}",
			script_file, getLevelString(log_level), num_threads, validate, log_file.value_or("None"));

		const std::filesystem::path final_out_dir = resolveOutputDir(script_file, output_dir);

		fers_context_t* context = fers_context_create();
		if (context == nullptr)
		{
			log(FERS_LOG_FATAL, "Failed to create FERS simulation context.");
			return 1;
		}

		log(FERS_LOG_INFO, "Loading scenario from '{}'...", script_file);
		if (fers_load_scenario_from_xml_file(context, script_file.c_str(), validate ? 1 : 0) != 0)
		{
			char* err = fers_get_last_error_message();
			log(FERS_LOG_FATAL, "Failed to load scenario: {}", err ? err : "Unknown error");
			fers_free_string(err);
			fers_context_destroy(context);
			return 1;
		}

		if (fers_set_output_directory(context, final_out_dir.string().c_str()) != 0)
		{
			char* err = fers_get_last_error_message();
			log(FERS_LOG_FATAL, "Failed to set output directory: {}", err ? err : "Unknown error");
			fers_free_string(err);
			fers_context_destroy(context);
			return 1;
		}

		if (generate_kml)
		{
			const std::filesystem::path kml_output_path = resolveKmlOutputPath(script_file, final_out_dir, kml_file);
			const std::string kml_output_file = kml_output_path.string();

			log(FERS_LOG_INFO, "Generating KML file for scenario: {}", kml_output_file);
			const int kml_result = fers_generate_kml(context, kml_output_file.c_str());
			if (kml_result == 0)
			{
				log(FERS_LOG_INFO, "KML file generated successfully: {}", kml_output_file);
			}
			else
			{
				char* err = fers_get_last_error_message();
				log(FERS_LOG_FATAL, "Failed to generate KML file: {}", err ? err : "Unknown error");
				fers_free_string(err);
			}

			fers_context_destroy(context);
			return kml_result == 0 ? 0 : 1;
		}

		if (fers_set_thread_count(num_threads) != 0)
		{
			char* err = fers_get_last_error_message();
			log(FERS_LOG_ERROR, "Failed to set number of threads: {}", err ? err : "Unknown error");
			fers_free_string(err);
		}

		log(FERS_LOG_INFO, "Starting simulation...");
		if (fers_run_simulation(context, nullptr, nullptr) != 0)
		{
			char* err = fers_get_last_error_message();
			log(FERS_LOG_FATAL, "Simulation run failed: {}", err ? err : "Unknown error");
			fers_free_string(err);
			fers_context_destroy(context);
			return 1;
		}
		log(FERS_LOG_INFO, "Simulation completed successfully.");

		fers_context_destroy(context);

		return 0;
	}
}
