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

	class ContextHandle
	{
	public:
		explicit ContextHandle(fers_context_t* context) noexcept : _context(context) {}
		ContextHandle(const ContextHandle&) = delete;
		ContextHandle& operator=(const ContextHandle&) = delete;
		ContextHandle(ContextHandle&&) = delete;
		ContextHandle& operator=(ContextHandle&&) = delete;
		~ContextHandle()
		{
			if (_context != nullptr)
			{
				fers_context_destroy(_context);
			}
		}

		[[nodiscard]] fers_context_t* get() const noexcept { return _context; }

	private:
		fers_context_t* _context;
	};

	bool isParseExit(const std::string& error)
	{
		return error == "Help requested." || error == "Version requested." || error == "No arguments provided.";
	}

	int logApiFailure(const fers_log_level_t level, const char* message, const int result)
	{
		char* err = fers_get_last_error_message();
		log(level, "{}: {}", message, (err != nullptr) ? err : "Unknown error");
		fers_free_string(err);
		return result;
	}

	int configureLogging(const core::Config& config)
	{
		const char* log_file_ptr = config.log_file ? config.log_file->c_str() : nullptr;
		if (fers_configure_logging(config.log_level, log_file_ptr) == 0)
		{
			return 0;
		}

		char* err = fers_get_last_error_message();
		std::cerr << "[ERROR] Failed to configure logging: " << ((err != nullptr) ? err : "Unknown error") << '\n';
		fers_free_string(err);
		return 1;
	}

	int loadScenario(fers_context_t* context, const core::Config& config)
	{
		log(FERS_LOG_INFO, "Loading scenario from '{}'...", config.script_file);
		if (fers_load_scenario_from_xml_file(context, config.script_file.c_str(), config.validate ? 1 : 0) == 0)
		{
			return 0;
		}
		return logApiFailure(FERS_LOG_FATAL, "Failed to load scenario", 1);
	}

	int configureVita49Output(fers_context_t* context, const core::Config& config)
	{
		if (!config.vita49_enabled)
		{
			return 0;
		}
		if (!config.vita49_fullscale.has_value())
		{
			log(FERS_LOG_FATAL, "Failed to configure VITA49 fullscale: missing required fullscale.");
			return 1;
		}
		if (fers_enable_vita49_udp_output(context, config.vita49_host.c_str(), config.vita49_port) != 0)
		{
			return logApiFailure(FERS_LOG_FATAL, "Failed to configure VITA49 endpoint", 1);
		}
		if (fers_set_vita49_fullscale(context, config.vita49_fullscale.value_or(0.0)) != 0)
		{
			return logApiFailure(FERS_LOG_FATAL, "Failed to configure VITA49 fullscale", 1);
		}
		if (config.vita49_epoch_unix_nanoseconds.has_value() &&
			fers_set_vita49_epoch_unix_nanoseconds(context, *config.vita49_epoch_unix_nanoseconds) != 0)
		{
			return logApiFailure(FERS_LOG_FATAL, "Failed to configure VITA49 epoch", 1);
		}
		if (config.vita49_max_udp_payload.has_value() &&
			fers_set_vita49_max_udp_payload(context, *config.vita49_max_udp_payload) != 0)
		{
			return logApiFailure(FERS_LOG_FATAL, "Failed to configure VITA49 max UDP payload", 1);
		}
		if (config.vita49_queue_depth.has_value() &&
			fers_set_vita49_queue_depth(context, *config.vita49_queue_depth) != 0)
		{
			return logApiFailure(FERS_LOG_FATAL, "Failed to configure VITA49 queue depth", 1);
		}
		return 0;
	}

	int setOutputDirectory(fers_context_t* context, const std::filesystem::path& final_out_dir)
	{
		if (fers_set_output_directory(context, final_out_dir.string().c_str()) == 0)
		{
			return 0;
		}
		return logApiFailure(FERS_LOG_FATAL, "Failed to set output directory", 1);
	}

	int generateKml(fers_context_t* context, const core::Config& config, const std::filesystem::path& final_out_dir)
	{
		const std::filesystem::path kml_output_path =
			core::resolveKmlOutputPath(config.script_file, final_out_dir, config.kml_file);
		const std::string kml_output_file = kml_output_path.string();

		log(FERS_LOG_INFO, "Generating KML file for scenario: {}", kml_output_file);
		if (fers_generate_kml(context, kml_output_file.c_str()) == 0)
		{
			log(FERS_LOG_INFO, "KML file generated successfully: {}", kml_output_file);
			return 0;
		}
		return logApiFailure(FERS_LOG_FATAL, "Failed to generate KML file", 1);
	}

	void configureThreadCount(const core::Config& config)
	{
		if (fers_set_thread_count(config.num_threads) == 0)
		{
			return;
		}
		char* err = fers_get_last_error_message();
		log(FERS_LOG_ERROR, "Failed to set number of threads: {}", (err != nullptr) ? err : "Unknown error");
		fers_free_string(err);
	}

	int runSimulation(fers_context_t* context)
	{
		log(FERS_LOG_INFO, "Starting simulation...");
		if (fers_run_simulation(context, nullptr, nullptr) != 0)
		{
			return logApiFailure(FERS_LOG_FATAL, "Simulation run failed", 1);
		}
		log(FERS_LOG_INFO, "Simulation completed successfully.");
		return 0;
	}
}

namespace core
{
	int runCli(const int argc, char* argv[])
	{
		const auto config_result = parseArguments(argc, argv);
		if (!config_result)
		{
			if (!isParseExit(config_result.error()))
			{
				std::cerr << "[ERROR] Argument parsing error: " << config_result.error() << '\n';
				return 1;
			}
			return 0;
		}

		const auto& config = config_result.value();

		if (configureLogging(config) != 0)
		{
			return 1;
		}

		log(FERS_LOG_INFO, "FERS CLI started. Using libfers backend.");

		log(FERS_LOG_DEBUG,
			"Running FERS with arguments: script_file={}, log_level={}, num_threads={}, validate={}, log_file={}",
			config.script_file, getLevelString(config.log_level), config.num_threads, config.validate,
			config.log_file.value_or("None"));

		const std::filesystem::path final_out_dir = resolveOutputDir(config.script_file, config.output_dir);

		ContextHandle context(fers_context_create());
		if (context.get() == nullptr)
		{
			log(FERS_LOG_FATAL, "Failed to create FERS simulation context.");
			return 1;
		}

		if (loadScenario(context.get(), config) != 0)
		{
			return 1;
		}

		if (configureVita49Output(context.get(), config) != 0)
		{
			return 1;
		}

		if (setOutputDirectory(context.get(), final_out_dir) != 0)
		{
			return 1;
		}

		if (config.generate_kml)
		{
			return generateKml(context.get(), config, final_out_dir);
		}

		configureThreadCount(config);

		return runSimulation(context.get());
	}
}
