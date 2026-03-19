// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2024-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file logging.h
 * @brief Header file for the logging system.
 *
 * This file defines a logging system with multiple log levels, supporting thread-safe logging
 * to a file or standard output, and allowing formatted log messages with source location details.
 * The logging system is customizable with log levels and file output for effective debugging
 * and tracking of events in an application.
 */

#pragma once

#define LOG(level, ...) log(level, std::source_location::current(), __VA_ARGS__)

#include <expected>
#include <format>
#include <fstream>
#include <mutex>
#include <source_location>
#include <string>
#include <utility>

namespace logging
{
	/**
	 * @class Level
	 * @brief Enum class representing the log levels.
	 */
	enum class Level
	{
		TRACE, ///< Trace level for detailed debugging information.
		DEBUG, ///< Debug level for general debugging information.
		INFO, ///< Info level for informational messages.
		WARNING, ///< Warning level for potentially harmful situations.
		ERROR, ///< Error level for error events.
		FATAL ///< Fatal level for severe error events.
	};

	/**
	 * @class Logger
	 * @brief Thread-safe logger class for handling logging operations.
	 */
	class Logger
	{
	public:
		/**
		 * @brief Sets the logging level.
		 *
		 * @param level The logging level to set.
		 */
		void setLevel(const Level level) noexcept { _log_level = level; }

		/**
		 * @brief Logs a message with a specific log level and source location.
		 *
		 * @param level The log level.
		 * @param message The message to log.
		 * @param location The source location of the log call.
		 */
		void log(Level level, const std::string& message,
				 const std::source_location& location = std::source_location::current()) noexcept;

		/**
		 * @brief Logs a formatted message with a specific log level and source location.
		 *
		 * @tparam Args Variadic template for format arguments.
		 * @param level The log level.
		 * @param location The source location of the log call.
		 * @param formatStr The format string.
		 * @param args The format arguments.
		 */
		template <typename... Args>
		void log(const Level level, const std::source_location& location, const std::string& formatStr,
				 Args&&... args) noexcept
		{
			if (level >= _log_level)
			{
				const std::string message = std::vformat(formatStr, std::make_format_args(args...));
				log(level, message, location);
			}
		}

		/**
		 * @brief Sets the log file path to log messages to a file.
		 *
		 * @param filePath The path to the log file.
		 * @return std::expected indicating success or error message on failure.
		 */
		std::expected<void, std::string> logToFile(const std::string& filePath) noexcept;

	private:
		Level _log_level = Level::INFO; ///< Current log level.
		std::optional<std::ofstream> _log_file; ///< Output file stream for logging to a file.
		std::mutex _log_mutex; ///< Mutex for thread-safe logging.

		/**
		 * @brief Gets the current timestamp as a string.
		 *
		 * @return The current timestamp.
		 */
		static std::string getCurrentTimestamp() noexcept;
	};

	/**
	 * @brief Externally available logger object.
	 */
	extern Logger logger;

	/**
	 * @brief Converts a log level enum value to its string representation.
	 *
	 * @param level The log level.
	 * @return A string representing the log level.
	 */
	inline std::string getLevelString(const Level level) noexcept
	{
		switch (level)
		{
		case Level::TRACE:
			return "TRACE";
		case Level::DEBUG:
			return "DEBUG";
		case Level::INFO:
			return "INFO";
		case Level::WARNING:
			return "WARNING";
		case Level::ERROR:
			return "ERROR";
		case Level::FATAL:
			return "FATAL";
		default:
			return "UNKNOWN";
		}
	}

	/**
	 * @brief Logs a formatted message with a specific log level and source location.
	 *
	 * @tparam Args Variadic template for format arguments.
	 * @param level The log level.
	 * @param location The source location of the log call.
	 * @param formatStr The format string.
	 * @param args The format arguments.
	 */
	template <typename... Args>
	void log(Level level, const std::source_location& location, const std::string& formatStr, Args&&... args) noexcept
	{
		logger.log(level, location, formatStr, std::forward<Args>(args)...);
	}
}
