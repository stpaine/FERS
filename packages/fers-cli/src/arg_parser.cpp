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
#include <cctype>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <format>
#include <iostream>
#include <limits>
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

	bool isUnsignedDecimal(const std::string& value) noexcept
	{
		return !value.empty() &&
			std::ranges::all_of(value, [](const unsigned char ch) { return std::isdigit(ch) != 0; });
	}

	std::expected<std::uint64_t, std::string> parseUnsigned64(const std::string& value,
															  const std::string& field_name) noexcept
	{
		if (!isUnsignedDecimal(value))
		{
			return std::unexpected(field_name + " must be an unsigned decimal integer");
		}
		try
		{
			std::size_t consumed = 0;
			const auto parsed = std::stoull(value, &consumed, 10);
			if (consumed != value.size())
			{
				return std::unexpected(field_name + " must be an unsigned decimal integer");
			}
			return static_cast<std::uint64_t>(parsed);
		}
		catch (const std::exception&)
		{
			return std::unexpected(field_name + " is out of range");
		}
	}

	std::expected<double, std::string> parsePositiveReal(const std::string& value,
														 const std::string& field_name) noexcept
	{
		try
		{
			std::size_t consumed = 0;
			const double parsed = std::stod(value, &consumed);
			if (consumed != value.size() || !std::isfinite(parsed) || parsed <= 0.0)
			{
				return std::unexpected(field_name + " must be a positive real number");
			}
			return parsed;
		}
		catch (const std::exception&)
		{
			return std::unexpected(field_name + " must be a positive real number");
		}
	}

	std::expected<void, std::string> handleVita49Endpoint(const std::string& value, core::Config& config) noexcept
	{
		const std::size_t separator = value.rfind(':');
		if (separator == std::string::npos || separator == 0 || separator + 1 == value.size())
		{
			return std::unexpected("Invalid VITA49 endpoint: expected host:port");
		}

		const std::string host = value.substr(0, separator);
		const std::string port_text = value.substr(separator + 1);
		if (host.find(':') != std::string::npos)
		{
			return std::unexpected("Invalid VITA49 endpoint: expected host:port");
		}

		const auto parsed_port = parseUnsigned64(port_text, "VITA49 port");
		if (!parsed_port)
		{
			return std::unexpected(parsed_port.error());
		}
		if (*parsed_port == 0 || *parsed_port > std::numeric_limits<std::uint16_t>::max())
		{
			return std::unexpected("VITA49 port must be in the range 1..65535");
		}

		config.vita49_enabled = true;
		config.vita49_host = host;
		config.vita49_port = static_cast<std::uint16_t>(*parsed_port);
		return {};
	}

	std::expected<void, std::string> handleVita49Fullscale(const std::string& value, core::Config& config) noexcept
	{
		const auto parsed = parsePositiveReal(value, "VITA49 fullscale");
		if (!parsed)
		{
			return std::unexpected(parsed.error());
		}
		config.vita49_fullscale = *parsed;
		return {};
	}

	std::expected<void, std::string> handleVita49Epoch(const std::string& value, core::Config& config) noexcept
	{
		constexpr std::uint64_t nanoseconds_per_second = 1'000'000'000ULL;
		constexpr std::uint64_t max_vrt_utc_epoch_ns =
			static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max()) * nanoseconds_per_second +
			(nanoseconds_per_second - 1ULL);
		const auto parsed = parseUnsigned64(value, "VITA49 epoch");
		if (!parsed)
		{
			return std::unexpected(parsed.error());
		}
		if (*parsed > max_vrt_utc_epoch_ns)
		{
			return std::unexpected("VITA49 epoch must fit the VRT 32-bit UTC seconds timestamp field");
		}
		config.vita49_epoch_unix_nanoseconds = *parsed;
		return {};
	}

	std::expected<void, std::string> handleVita49MaxUdpPayload(const std::string& value, core::Config& config) noexcept
	{
		const auto parsed = parseUnsigned64(value, "VITA49 max UDP payload");
		if (!parsed)
		{
			return std::unexpected(parsed.error());
		}
		if (*parsed < 64u || *parsed > 65'507u)
		{
			return std::unexpected("VITA49 max UDP payload must be between 64 and 65507 bytes");
		}
		config.vita49_max_udp_payload = static_cast<std::uint16_t>(*parsed);
		return {};
	}

	std::expected<void, std::string> handleVita49QueueDepth(const std::string& value, core::Config& config) noexcept
	{
		const auto parsed = parseUnsigned64(value, "VITA49 queue depth");
		if (!parsed)
		{
			return std::unexpected(parsed.error());
		}
		if (*parsed == 0u || *parsed > std::numeric_limits<std::uint32_t>::max())
		{
			return std::unexpected("VITA49 queue depth must be in the range 1..4294967295");
		}
		config.vita49_queue_depth = static_cast<std::uint32_t>(*parsed);
		return {};
	}

	using ConfigValueHandler = std::expected<void, std::string> (*)(const std::string&, core::Config&) noexcept;

	std::optional<ConfigValueHandler> valueOptionHandler(const std::string& arg) noexcept
	{
		if (arg == "--vita49")
		{
			return handleVita49Endpoint;
		}
		if (arg == "--vita49-fullscale")
		{
			return handleVita49Fullscale;
		}
		if (arg == "--vita49-epoch")
		{
			return handleVita49Epoch;
		}
		if (arg == "--vita49-max-udp-payload")
		{
			return handleVita49MaxUdpPayload;
		}
		if (arg == "--vita49-queue-depth")
		{
			return handleVita49QueueDepth;
		}
		return std::nullopt;
	}

	std::expected<void, std::string> validateVita49Options(const core::Config& config) noexcept
	{
		if (!config.vita49_enabled)
		{
			if (config.vita49_fullscale.has_value())
			{
				return std::unexpected("--vita49-fullscale requires --vita49");
			}
			if (config.vita49_epoch_unix_nanoseconds.has_value())
			{
				return std::unexpected("--vita49-epoch requires --vita49");
			}
			if (config.vita49_max_udp_payload.has_value())
			{
				return std::unexpected("--vita49-max-udp-payload requires --vita49");
			}
			if (config.vita49_queue_depth.has_value())
			{
				return std::unexpected("--vita49-queue-depth requires --vita49");
			}
			return {};
		}
		if (!config.vita49_fullscale.has_value())
		{
			return std::unexpected("--vita49-fullscale is required when --vita49 is used");
		}
		return {};
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
		std::string const log_file_path = arg.substr(11);
		if (isValidLogFileExtension(log_file_path))
		{
			config.log_file = log_file_path;
			return {};
		}

		std::cerr << "[ERROR] Invalid log file extension. Must be .log or .txt\n";
		return std::unexpected("Invalid log file extension: " + log_file_path);
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
			const int requested_threads = std::stoi(arg.substr(3));
			if (requested_threads <= 0)
			{
				return std::unexpected("Number of threads must be greater than 0");
			}

			config.num_threads = static_cast<unsigned>(requested_threads);
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
		if (arg.rfind("--out-dir=", 0) == 0)
		{
			config.output_dir = arg.substr(10);
			return {};
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
		if (arg.rfind("--kml=", 0) == 0)
		{
			config.generate_kml = true;
			config.kml_file = arg.substr(6);
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
		const char* version = fers_get_version();

		std::cout << "/------------------------------------------------\\\n"
				  << "| FERS - The Flexible Extensible Radar Simulator |\n"
				  << std::format("| Version {:<40}|\n", version)
				  << "\\------------------------------------------------/\n"
				  << "Usage: " << programName << R"( <scriptfile> [options]

Options:
  --help, -h              Show this help message and exit
  --version, -v           Show version information and exit
  --no-validate           Disable XML schema validation before running.
  --kml[=<file>]          Generate a KML visualization of the scenario and exit. If a filename
                          is provided, it will be used. Otherwise, it defaults to the scenario
                          name with a .kml extension in the output directory.
  --out-dir=<dir>         Set the output directory for simulation results and default KML output.
                          Defaults to the directory containing the script file.
  --vita49 host:port      Stream receiver output as the FERS VITA 49.2 UDP profile.
  --vita49-fullscale <x>  Set required positive fixed ADC full-scale for VITA int16 IQ output.
  --vita49-epoch <ns>     Set optional deterministic VITA epoch as Unix nanoseconds.
  --vita49-max-udp-payload <bytes>
                          Set optional VITA UDP payload cap, 64..65507 bytes.
  --vita49-queue-depth <packets>
                          Set optional VITA sender queue depth, greater than zero.
  --log-level=<level>     Set the logging level (TRACE, DEBUG, INFO, WARNING, ERROR, FATAL)
  --log-file=<file>       Log output to the specified .log or .txt file as well as the console.
  -n=<threads>            Number of threads to use

Arguments:
  <scriptfile>            Path to the simulation script file (XML)

Example:
  )" << programName
				  << R"( simulation.fersxml --out-dir=./results --log-level=DEBUG -n=4

This program runs radar simulations based on an XML script file.
Make sure the script file follows the correct format to avoid errors.
)";
	}

	void showVersion() noexcept
	{
		const char* version = fers_get_version();

		std::cout << '\n'
				  << "/------------------------------------------------\\\n"
				  << "| FERS - The Flexible Extensible Radar Simulator |\n"
				  << std::format("| Version {:<40}|\n", version)
				  << "| Authors: Marc Brooker, Michael Inggs,         |\n"
				  << "|          and FERS Contributors                |\n"
				  << "\\------------------------------------------------/\n\n";
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
			const std::string arg = argv[i];
			const auto require_value = [&](const std::string& option) -> std::expected<std::string, std::string>
			{
				if (i + 1 >= argc)
				{
					return std::unexpected(option + " requires a value");
				}
				++i;
				return std::string{argv[i]};
			};

			if (const auto handler = valueOptionHandler(arg))
			{
				const auto value = require_value(arg);
				if (!value)
				{
					return std::unexpected(value.error());
				}
				if (const auto result = (*handler)(*value, config); !result)
				{
					return std::unexpected(result.error());
				}
				continue;
			}

			if (const auto result = handleArgument(arg, config, script_file_set, argv[0]); !result)
			{
				return std::unexpected(result.error());
			}
		}

		if (!script_file_set)
		{
			return std::unexpected("No script file provided.");
		}
		if (const auto result = validateVita49Options(config); !result)
		{
			return std::unexpected(result.error());
		}
		return config;
	}
}
