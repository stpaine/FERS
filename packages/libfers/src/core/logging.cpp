// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2024-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file logging.cpp
 * @brief Implementation of the logging system.
 */

#include "logging.h"

#include <chrono>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace logging
{
	Logger logger;

	std::string Logger::getCurrentTimestamp() noexcept
	{
		const auto now = std::chrono::system_clock::now();
		const std::time_t time = std::chrono::system_clock::to_time_t(now);
		std::tm tm{};
		localtime_r(&time, &tm);

		std::ostringstream oss;
		oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
		return oss.str();
	}

	void Logger::log(const Level level, const std::string& message, const std::source_location& location) noexcept
	{
		if (level >= _log_level)
		{
			Callback callback = nullptr;
			void* callback_user_data = nullptr;
			std::string line;

			{
				std::scoped_lock lock(_log_mutex);

				const std::string filename = std::filesystem::path(location.file_name()).filename().string();
				const std::string file_line = filename + ":" + std::to_string(location.line());

				std::ostringstream oss;
				oss << "[" << getCurrentTimestamp() << "] " << "[" << std::setw(7) << std::left << getLevelString(level)
					<< "] " << "[" << std::setw(30) << std::left << file_line << "] " << message;
				line = oss.str();

				std::cerr << line << '\n';

				if (_log_file && _log_file->is_open())
				{
					*_log_file << line << '\n';
					_log_file->flush();
				}

				callback = _callback;
				callback_user_data = _callback_user_data;
			}

			if (callback != nullptr)
			{
				callback(level, line, callback_user_data);
			}
		}
	}

	void Logger::setCallback(Callback callback, void* user_data) noexcept
	{
		std::scoped_lock lock(_log_mutex);
		_callback = callback;
		_callback_user_data = user_data;
	}

	std::expected<void, std::string> Logger::logToFile(const std::string& filePath) noexcept
	{
		std::scoped_lock lock(_log_mutex);

		std::ofstream file(filePath, std::ios::out | std::ios::trunc);
		if (!file)
		{
			return std::unexpected("Unable to open log file: " + filePath);
		}

		_log_file = std::move(file);
		return {};
	}
}
