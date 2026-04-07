// SPDX-License-Identifier: GPL-2.0-only
//
// Local harness support code.

#pragma once

#include <iostream>
#include <mutex>
#include <sstream>
#include <string>

namespace logging
{
	enum class Level
	{
		TRACE,
		DEBUG,
		INFO,
		WARNING,
		ERROR,
		FATAL,
		OFF
	};

	class Logger
	{
	public:
		void setLevel(const Level level) noexcept { _level = level; }

		template <typename... Args>
		void log(const Level level, Args&&... args) noexcept
		{
			if (level < _level || _level == Level::OFF)
			{
				return;
			}

			std::ostringstream stream;
			(stream << ... << std::forward<Args>(args));

			std::scoped_lock lock(_mutex);
			std::cerr << "[" << levelString(level) << "] " << stream.str() << '\n';
		}

	private:
		static const char* levelString(Level level) noexcept;

		Level _level = Level::INFO;
		std::mutex _mutex;
	};

	extern Logger logger;
}

#define LOG(level, ...) ::logging::logger.log(level, __VA_ARGS__)

