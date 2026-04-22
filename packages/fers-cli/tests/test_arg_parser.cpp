// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

#include <catch2/catch_test_macros.hpp>
#include <functional>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "arg_parser.h"

namespace
{
	std::expected<core::Config, std::string> parseArgs(const std::initializer_list<std::string>& args) noexcept
	{
		std::vector<std::string> storage;
		storage.reserve(args.size() + 1);
		storage.emplace_back("fers-cli");
		storage.insert(storage.end(), args.begin(), args.end());

		std::vector<char*> argv;
		argv.reserve(storage.size());
		for (std::string& arg : storage)
		{
			argv.push_back(arg.data());
		}

		return core::parseArguments(static_cast<int>(argv.size()), argv.data());
	}

	template <typename Fn>
	std::string captureStdout(Fn&& fn)
	{
		std::ostringstream buffer;
		auto* const original = std::cout.rdbuf(buffer.rdbuf());
		fn();
		std::cout.rdbuf(original);
		return buffer.str();
	}
}

TEST_CASE("parseArguments accepts a script file and preserves defaults", "[fers-cli][arg-parser]")
{
	const auto result = parseArgs({"scenario.fersxml"});

	REQUIRE(result);
	CHECK(result->script_file == "scenario.fersxml");
	CHECK(result->log_level == FERS_LOG_INFO);
	CHECK(result->num_threads == std::thread::hardware_concurrency());
	CHECK(result->validate);
	CHECK_FALSE(result->log_file.has_value());
	CHECK_FALSE(result->generate_kml);
	CHECK_FALSE(result->kml_file.has_value());
	CHECK_FALSE(result->output_dir.has_value());
}

TEST_CASE("parseArguments applies recognized options", "[fers-cli][arg-parser]")
{
	const unsigned hw_threads = std::thread::hardware_concurrency();
	const unsigned requested_threads = hw_threads > 0 ? hw_threads : 1;

	const auto n_arg = std::string{"-n="} + std::to_string(requested_threads);

	const auto result = parseArgs({"scenario.fersxml", "--log-level=DEBUG", "--log-file=runner.log", n_arg,
								   "--out-dir=results", "--no-validate", "--kml=preview.kml"});

	REQUIRE(result);
	CHECK(result->script_file == "scenario.fersxml");
	CHECK(result->log_level == FERS_LOG_DEBUG);
	CHECK(result->log_file == std::optional<std::string>{"runner.log"});
	CHECK(result->num_threads == requested_threads);
	CHECK_FALSE(result->validate);
	CHECK(result->generate_kml);
	CHECK(result->kml_file == std::optional<std::string>{"preview.kml"});
	CHECK(result->output_dir == std::optional<std::string>{"results"});
}

TEST_CASE("parseArguments accepts flags before and after the script", "[fers-cli][arg-parser]")
{
	const auto result = parseArgs({"--out-dir=results", "scenario.fersxml", "--kml"});

	REQUIRE(result);
	CHECK(result->script_file == "scenario.fersxml");
	CHECK(result->output_dir == std::optional<std::string>{"results"});
	CHECK(result->generate_kml);
	CHECK_FALSE(result->kml_file.has_value());
}

TEST_CASE("parseArguments reports help and prints usage text", "[fers-cli][arg-parser]")
{
	std::expected<core::Config, std::string> result = std::unexpected(std::string{});
	const std::string output = captureStdout([&]() { result = parseArgs({"--help"}); });

	REQUIRE_FALSE(result);
	CHECK(result.error() == "Help requested.");
	CHECK(output.contains("Usage: fers-cli"));
	CHECK(output.contains("--kml[=<file>]"));
}

TEST_CASE("parseArguments reports version and prints version text", "[fers-cli][arg-parser]")
{
	std::expected<core::Config, std::string> result = std::unexpected(std::string{});
	const std::string output = captureStdout([&]() { result = parseArgs({"--version"}); });

	REQUIRE_FALSE(result);
	CHECK(result.error() == "Version requested.");
	CHECK(output.contains(std::string{"Version "} + fers_get_version()));
	CHECK(output.contains("FERS - The Flexible Extensible Radar Simulator"));
}

TEST_CASE("parseArguments reports missing arguments and prints help text", "[fers-cli][arg-parser]")
{
	std::expected<core::Config, std::string> result = std::unexpected(std::string{});
	const std::string output = captureStdout([&]() { result = parseArgs({}); });

	REQUIRE_FALSE(result);
	CHECK(result.error() == "No arguments provided.");
	CHECK(output.contains("Usage: fers-cli"));
}

TEST_CASE("parseArguments rejects an invalid log level", "[fers-cli][arg-parser]")
{
	const auto result = parseArgs({"scenario.fersxml", "--log-level=VERBOSE"});

	REQUIRE_FALSE(result);
	CHECK(result.error() == "Invalid log level: VERBOSE");
}

TEST_CASE("parseArguments rejects an invalid log file extension", "[fers-cli][arg-parser]")
{
	const auto result = parseArgs({"scenario.fersxml", "--log-file=runner.json"});

	REQUIRE_FALSE(result);
	CHECK(result.error() == "Invalid log file extension: runner.json");
}

TEST_CASE("parseArguments rejects a non-numeric thread count", "[fers-cli][arg-parser]")
{
	const auto result = parseArgs({"scenario.fersxml", "-n=many"});

	REQUIRE_FALSE(result);
	CHECK(result.error() == "Invalid number of threads specified.");
}

TEST_CASE("parseArguments rejects a zero thread count", "[fers-cli][arg-parser]")
{
	const auto result = parseArgs({"scenario.fersxml", "-n=0"});

	REQUIRE_FALSE(result);
	CHECK(result.error() == "Number of threads must be greater than 0");
}

TEST_CASE("parseArguments rejects a negative thread count", "[fers-cli][arg-parser]")
{
	const auto result = parseArgs({"scenario.fersxml", "-n=-4"});

	REQUIRE_FALSE(result);
	CHECK(result.error() == "Number of threads must be greater than 0");
}

TEST_CASE("parseArguments requires a script file", "[fers-cli][arg-parser]")
{
	const auto result = parseArgs({"--out-dir=results"});

	REQUIRE_FALSE(result);
	CHECK(result.error() == "No script file provided.");
}

TEST_CASE("parseArguments rejects unknown options", "[fers-cli][arg-parser]")
{
	const auto result = parseArgs({"scenario.fersxml", "--bogus"});

	REQUIRE_FALSE(result);
	CHECK(result.error() == "Unrecognized argument: --bogus");
}

TEST_CASE("parseArguments rejects extra positional arguments", "[fers-cli][arg-parser]")
{
	const auto result = parseArgs({"scenario.fersxml", "other.fersxml"});

	REQUIRE_FALSE(result);
	CHECK(result.error() == "Unrecognized argument: other.fersxml");
}

TEST_CASE("parseArguments clamps the thread count to the hardware limit", "[fers-cli][arg-parser]")
{
	const unsigned max_threads = std::thread::hardware_concurrency();
	if (max_threads == 0)
	{
		SUCCEED("hardware_concurrency() returned 0; clamp behavior is not observable on this platform.");
		return;
	}

	const auto result = parseArgs({"scenario.fersxml", "-n=999999"});

	REQUIRE(result);
	CHECK(result->num_threads == max_threads);
}
