// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2024-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file arg_parser.cpp
 * @brief Implementation of the command-line argument parser for the application.
 */

#include "arg_parser.h"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace
{
	/**
	 * @brief Checks if the given file has a valid log file extension.
	 *
	 * @param filePath The path of the log file.
	 * @return true if the file has a valid extension, false otherwise.
	 */
	bool isValidLogFileExtension(const std::string& filePath) noexcept
	{
		static const std::vector<std::string> VALID_EXTENSIONS = {".log", ".txt"};
		std::string extension = std::filesystem::path(filePath).extension().string();
		std::ranges::transform(extension, extension.begin(), tolower);
		return std::ranges::find(VALID_EXTENSIONS, extension) != VALID_EXTENSIONS.end();
	}

	/**
	 * @brief Parses the logging level from a string representation.
	 *
	 * @param level The string representation of the logging level.
	 * @return std::optional<Level> The corresponding logging level, or `std::nullopt` if invalid.
	 */
	std::optional<fers_log_level_t> parseLogLevel(const std::string& level) noexcept
	{
		static const std::unordered_map<std::string, fers_log_level_t> LEVEL_MAP = {
			{"TRACE", FERS_LOG_TRACE},	   {"DEBUG", FERS_LOG_DEBUG}, {"INFO", FERS_LOG_INFO},
			{"WARNING", FERS_LOG_WARNING}, {"ERROR", FERS_LOG_ERROR}, {"FATAL", FERS_LOG_FATAL}};

		if (const auto it = LEVEL_MAP.find(level); it != LEVEL_MAP.end())
		{
			return it->second;
		}
		return std::nullopt;
	}

	/**
	 * @brief Handles the log-level argument and sets the logging level.
	 *
	 * @param arg The log-level argument string.
	 * @param config The configuration object to update.
	 * @return std::expected<void, std::string> An expected object with an error message if the log level is invalid.
	 */
	std::expected<void, std::string> handleLogLevel(const std::string& arg, core::Config& config) noexcept
	{
		const std::string level_str = arg.substr(12);
		if (const auto level = parseLogLevel(level_str))
		{
			config.log_level = *level;
			return {};
		}

		std::cerr << "[ERROR] Invalid log level '" << level_str << "'\n";
		return std::unexpected("Invalid log level: " + level_str);
	}

	/**
	 * @brief Handles the log-file argument and sets the log file path.
	 *
	 * @param arg The log-file argument string.
	 * @param config The configuration object to update.
	 * @return std::expected<void, std::string> An expected object with an error message if the log file path is
	 * invalid.
	 */
	std::expected<void, std::string> handleLogFile(const std::string& arg, core::Config& config) noexcept
	{
		if (std::string log_file_path = arg.substr(11); isValidLogFileExtension(log_file_path))
		{
			config.log_file = log_file_path;
			return {};
		}
		else
		{
			std::cerr << "[ERROR] Invalid log file extension. Must be .log or .txt\n";
			return std::unexpected("Invalid log file extension: " + log_file_path);
		}
	}

	/**
	 * @brief Handles the number of threads argument and sets the number of threads.
	 *
	 * @param arg The number of threads argument string.
	 * @param config The configuration object to update.
	 * @return std::expected<void, std::string> An expected object with an error message if the number of threads is
	 * invalid.
	 */
	std::expected<void, std::string> handleNumThreads(const std::string& arg, core::Config& config) noexcept
	{
		try
		{
			config.num_threads = std::stoi(arg.substr(3));
			if (config.num_threads == 0)
			{
				return std::unexpected("Number of threads must be greater than 0");
			}

			if (const unsigned max_threads = std::thread::hardware_concurrency();
				max_threads > 0 && config.num_threads > max_threads)
			{
				std::cerr << "[WARNING] Thread count exceeds available processors. Clamping.\n";
				config.num_threads = max_threads;
			}
			return {};
		}
		catch (const std::exception&)
		{
			return std::unexpected("Invalid number of threads specified.");
		}
	}

	/**
	 * @brief Handles the command-line argument and updates the configuration.
	 *
	 * @param arg The command-line argument string.
	 * @param config The configuration object to update.
	 * @param scriptFileSet A flag indicating if the script file has been set.
	 * @param programName The name of the program executable.
	 * @return std::expected<void, std::string> An expected object with an error message if the argument is invalid.
	 */
	std::expected<void, std::string> handleArgument(const std::string& arg, core::Config& config, bool& scriptFileSet,
													const char* programName) noexcept
	{
		if (arg == "--help" || arg == "-h")
		{
			core::showHelp(programName);
			return std::unexpected("Help requested.");
		}
		if (arg == "--version" || arg == "-v")
		{
			core::showVersion();
			return std::unexpected("Version requested.");
		}
		if (arg.rfind("--log-level=", 0) == 0)
		{
			return handleLogLevel(arg, config);
		}
		if (arg.rfind("--log-file=", 0) == 0)
		{
			return handleLogFile(arg, config);
		}
		if (arg.rfind("-n=", 0) == 0)
		{
			return handleNumThreads(arg, config);
		}
		if (arg == "--no-validate")
		{
			config.validate = false;
			return {};
		}
		if (arg == "--kml")
		{
			config.generate_kml = true;
			return {};
		}
		if (arg[0] != '-' && !scriptFileSet)
		{
			config.script_file = arg;
			scriptFileSet = true;
			return {};
		}

		std::cerr << "[ERROR] Unrecognized option: '" << arg << "'\n";
		return std::unexpected("Unrecognized argument: " + arg);
	}
}

namespace core
{
	void showHelp(const char* programName) noexcept
	{
		std::cout << R"(/------------------------------------------------\
| FERS - The Flexible Extensible Radar Simulator |
| Version 1.0.0                                  |
\------------------------------------------------/
Usage: )" << programName
				  << R"( <scriptfile> [options]

Options:
  --help, -h              Show this help message and exit
  --version, -v           Show version information and exit
  --no-validate           Disable XML schema validation before running.
  --kml                   Generate a KML visualization of the scenario and exit. The output file
                          will have the same name as the input file with a .kml extension.
  --log-level=<level>     Set the logging level (TRACE, DEBUG, INFO, WARNING, ERROR, FATAL)
  --log-file=<file>       Log output to the specified .log or .txt file as well as the console.
  -n=<threads>            Number of threads to use

Arguments:
  <scriptfile>            Path to the simulation script file (XML)

Example:
  )" << programName
				  << R"( simulation.fersxml --log-level=DEBUG --log-file=output.log -n=4

This program runs radar simulations based on an XML script file.
Make sure the script file follows the correct format to avoid errors.
)";
	}

	void showVersion() noexcept
	{
		std::cout << R"(
/------------------------------------------------\
| FERS - The Flexible Extensible Radar Simulator |
| Version 1.0.0                                  |
| Author: Marc Brooker                           |
\------------------------------------------------/
)" << std::endl;
	}

	std::expected<Config, std::string> parseArguments(const int argc, char* argv[]) noexcept
	{
		Config config;
		bool script_file_set = false;

		if (argc < 2)
		{
			showHelp(argv[0]);
			return std::unexpected("No arguments provided.");
		}

		for (int i = 1; i < argc; ++i)
		{
			if (const auto result = handleArgument(argv[i], config, script_file_set, argv[0]); !result)
			{
				return std::unexpected(result.error());
			}
		}

		if (!script_file_set)
		{
			return std::unexpected("No script file provided.");
		}
		return config;
	}
}
