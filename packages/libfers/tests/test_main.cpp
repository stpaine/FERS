// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>
#include <cstdlib>
#include <iostream>
#include <string_view>

#include "core/logging.h"

// Helper to parse log level from an environment variable string
logging::Level levelFromString(std::string_view levelStr)
{
	if (levelStr == "TRACE")
		return logging::Level::TRACE;
	if (levelStr == "DEBUG")
		return logging::Level::DEBUG;
	if (levelStr == "INFO")
		return logging::Level::INFO;
	if (levelStr == "WARNING")
		return logging::Level::WARNING;
	if (levelStr == "ERROR")
		return logging::Level::ERROR;
	if (levelStr == "FATAL")
		return logging::Level::FATAL;
	if (levelStr == "OFF")
		return logging::Level::OFF;
	return logging::Level::OFF; // Default to OFF if invalid
}

int main(int argc, char* argv[])
{
	// By default, suppress all logging during tests to keep output clean.
	logging::Level testLogLevel = logging::Level::OFF;

	// Allow overriding the log level with an environment variable for debugging.
	// Example: FERS_TEST_LOG_LEVEL=DEBUG ctest -R my_failing_test
	if (const char* env_level = std::getenv("FERS_TEST_LOG_LEVEL"))
	{
		testLogLevel = levelFromString(env_level);
		if (testLogLevel != logging::Level::OFF)
		{
			std::cout << "[Test Runner] Overriding log level to: " << logging::getLevelString(testLogLevel)
					  << std::endl;
		}
	}

	logging::logger.setLevel(testLogLevel);

	// Run the Catch2 test session
	int result = Catch::Session().run(argc, argv);

	return result;
}
