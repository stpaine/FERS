#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "core/logging.h"

namespace
{
	struct CerrCapture
	{
		std::ostringstream buffer;
		std::streambuf* old{nullptr};
		CerrCapture() { old = std::cerr.rdbuf(buffer.rdbuf()); }
		~CerrCapture() { std::cerr.rdbuf(old); }
		[[nodiscard]] std::string str() const { return buffer.str(); }
	};

	struct LogLevelGuard
	{
		explicit LogLevelGuard(logging::Level level) { logging::logger.setLevel(level); }
		~LogLevelGuard() { logging::logger.setLevel(logging::Level::INFO); }
	};
}

TEST_CASE("getLevelString maps levels", "[core][logging]")
{
	REQUIRE(logging::getLevelString(logging::Level::TRACE) == "TRACE");
	REQUIRE(logging::getLevelString(logging::Level::DEBUG) == "DEBUG");
	REQUIRE(logging::getLevelString(logging::Level::INFO) == "INFO");
	REQUIRE(logging::getLevelString(logging::Level::WARNING) == "WARNING");
	REQUIRE(logging::getLevelString(logging::Level::ERROR) == "ERROR");
	REQUIRE(logging::getLevelString(logging::Level::FATAL) == "FATAL");

	const auto unknown = static_cast<logging::Level>(-1);
	REQUIRE(logging::getLevelString(unknown) == "UNKNOWN");
}

TEST_CASE("Logger respects log level filtering", "[core][logging]")
{
	LogLevelGuard guard(logging::Level::ERROR);
	CerrCapture capture;

	logging::logger.log(logging::Level::INFO, "ignored message");
	REQUIRE(capture.str().empty());

	logging::logger.log(logging::Level::ERROR, "logged message");
	const std::string output = capture.str();
	REQUIRE(output.find("logged message") != std::string::npos);
	REQUIRE(output.find("ERROR") != std::string::npos);
}

TEST_CASE("Logger includes source location and formats messages", "[core][logging]")
{
	LogLevelGuard guard(logging::Level::TRACE);
	CerrCapture capture;

	const std::source_location location = std::source_location::current();
	logging::logger.log(logging::Level::INFO, location, "value {}", 7);

	const std::string output = capture.str();
	const std::string filename = std::filesystem::path(location.file_name()).filename().string();
	const std::string file_line = filename + ":" + std::to_string(location.line());

	REQUIRE(output.find("value 7") != std::string::npos);
	REQUIRE(output.find(file_line) != std::string::npos);
	REQUIRE(output.find("INFO") != std::string::npos);
}

TEST_CASE("Logger can write to a file", "[core][logging]")
{
	LogLevelGuard guard(logging::Level::INFO);

	const auto log_path = std::filesystem::temp_directory_path() / "fers_logging_test.log";
	const auto result = logging::logger.logToFile(log_path.string());
	REQUIRE(result.has_value());

	logging::logger.log(logging::Level::INFO, "file message");

	std::ifstream file(log_path);
	REQUIRE(file.good());
	std::stringstream contents;
	contents << file.rdbuf();
	const std::string text = contents.str();
	REQUIRE(text.find("file message") != std::string::npos);
	REQUIRE(text.find("INFO") != std::string::npos);
}

TEST_CASE("Logger reports file open errors", "[core][logging]")
{
	const auto bad_path = std::filesystem::temp_directory_path() / "fers_invalid" / "log.txt";
	const auto result = logging::logger.logToFile(bad_path.string());

	REQUIRE_FALSE(result.has_value());
}
