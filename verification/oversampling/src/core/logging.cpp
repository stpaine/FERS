// SPDX-License-Identifier: GPL-2.0-only
//
// Local harness support code.

#include "logging.h"

namespace logging
{
	Logger logger;

	const char* Logger::levelString(const Level level) noexcept
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
		case Level::OFF:
			return "OFF";
		}
		return "UNKNOWN";
	}
}

