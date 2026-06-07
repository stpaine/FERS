// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <functional>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
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
		std::ostringstream const buffer;
		auto* const original = std::cout.rdbuf(buffer.rdbuf());
		std::invoke(std::forward<Fn>(fn));
		std::cout.rdbuf(original);
		return buffer.str();
	}

	void requireRecognizedOptions(const core::Config& config, const unsigned requested_threads)
	{
		CHECK(config.script_file == "scenario.fersxml");
		CHECK(config.log_level == FERS_LOG_DEBUG);
		CHECK(config.log_file == std::optional<std::string>{"runner.log"});
		CHECK(config.num_threads == requested_threads);
		CHECK_FALSE(config.validate);
		CHECK(config.generate_kml);
		CHECK(config.kml_file == std::optional<std::string>{"preview.kml"});
		CHECK(config.output_dir == std::optional<std::string>{"results"});
		CHECK(config.vita49_enabled);
		CHECK(config.vita49_host == "localhost");
		CHECK(config.vita49_port == 4991);
		REQUIRE(config.vita49_fullscale.has_value());
		CHECK(config.vita49_fullscale.value_or(0.0) == 2.5);
		CHECK(config.vita49_epoch_unix_nanoseconds == std::optional<std::uint64_t>{1700000000123456789ULL});
		CHECK(config.vita49_max_udp_payload == std::optional<std::uint16_t>{900});
		CHECK(config.vita49_queue_depth == std::optional<std::uint32_t>{17});
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
	CHECK_FALSE(result->vita49_enabled);
	CHECK(result->vita49_host.empty());
	CHECK(result->vita49_port == 0);
	CHECK_FALSE(result->vita49_fullscale.has_value());
	CHECK_FALSE(result->vita49_epoch_unix_nanoseconds.has_value());
	CHECK_FALSE(result->vita49_max_udp_payload.has_value());
	CHECK_FALSE(result->vita49_queue_depth.has_value());
}

TEST_CASE("parseArguments applies recognized options", "[fers-cli][arg-parser]")
{
	const unsigned hw_threads = std::thread::hardware_concurrency();
	const unsigned requested_threads = hw_threads > 0 ? hw_threads : 1;

	const auto n_arg = std::string{"-n="} + std::to_string(requested_threads);

	const auto result = parseArgs(
		{"scenario.fersxml", "--log-level=DEBUG", "--log-file=runner.log", n_arg, "--out-dir=results", "--no-validate",
		 "--kml=preview.kml", "--vita49", "localhost:4991", "--vita49-fullscale", "2.5", "--vita49-epoch",
		 "1700000000123456789", "--vita49-max-udp-payload", "900", "--vita49-queue-depth", "17"});

	REQUIRE(result);
	requireRecognizedOptions(*result, requested_threads);
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
	CHECK(output.contains("--vita49 host:port"));
	CHECK(output.contains("--vita49-max-udp-payload <bytes>"));
	CHECK(output.contains("--vita49-queue-depth <packets>"));
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

TEST_CASE("parseArguments rejects malformed VITA49 endpoints", "[fers-cli][arg-parser][vita49]")
{
	const auto missing_port = parseArgs({"scenario.fersxml", "--vita49", "localhost", "--vita49-fullscale", "1.0"});
	REQUIRE_FALSE(missing_port);
	CHECK(missing_port.error() == "Invalid VITA49 endpoint: expected host:port");

	const auto empty_host = parseArgs({"scenario.fersxml", "--vita49", ":4991", "--vita49-fullscale", "1.0"});
	REQUIRE_FALSE(empty_host);
	CHECK(empty_host.error() == "Invalid VITA49 endpoint: expected host:port");

	const auto zero_port = parseArgs({"scenario.fersxml", "--vita49", "localhost:0", "--vita49-fullscale", "1.0"});
	REQUIRE_FALSE(zero_port);
	CHECK(zero_port.error() == "VITA49 port must be in the range 1..65535");

	const auto high_port = parseArgs({"scenario.fersxml", "--vita49", "localhost:65536", "--vita49-fullscale", "1.0"});
	REQUIRE_FALSE(high_port);
	CHECK(high_port.error() == "VITA49 port must be in the range 1..65535");
}

TEST_CASE("parseArguments rejects invalid VITA49 fullscale usage", "[fers-cli][arg-parser][vita49]")
{
	const auto missing = parseArgs({"scenario.fersxml", "--vita49", "localhost:4991"});
	REQUIRE_FALSE(missing);
	CHECK(missing.error() == "--vita49-fullscale is required when --vita49 is used");

	const auto zero = parseArgs({"scenario.fersxml", "--vita49", "localhost:4991", "--vita49-fullscale", "0"});
	REQUIRE_FALSE(zero);
	CHECK(zero.error() == "VITA49 fullscale must be a positive real number");

	const auto without_endpoint = parseArgs({"scenario.fersxml", "--vita49-fullscale", "1.0"});
	REQUIRE_FALSE(without_endpoint);
	CHECK(without_endpoint.error() == "--vita49-fullscale requires --vita49");
}

TEST_CASE("parseArguments rejects invalid VITA49 epoch usage", "[fers-cli][arg-parser][vita49]")
{
	const auto negative = parseArgs(
		{"scenario.fersxml", "--vita49", "localhost:4991", "--vita49-fullscale", "1.0", "--vita49-epoch", "-1"});
	REQUIRE_FALSE(negative);
	CHECK(negative.error() == "VITA49 epoch must be an unsigned decimal integer");

	const auto too_large = parseArgs({"scenario.fersxml", "--vita49", "localhost:4991", "--vita49-fullscale", "1.0",
									  "--vita49-epoch", "4294967296000000000"});
	REQUIRE_FALSE(too_large);
	CHECK(too_large.error() == "VITA49 epoch must fit the VRT 32-bit UTC seconds timestamp field");

	const auto without_endpoint = parseArgs({"scenario.fersxml", "--vita49-epoch", "1700000000000000000"});
	REQUIRE_FALSE(without_endpoint);
	CHECK(without_endpoint.error() == "--vita49-epoch requires --vita49");
}

TEST_CASE("parseArguments rejects invalid VITA49 payload and queue usage", "[fers-cli][arg-parser][vita49]")
{
	const auto payload_too_small = parseArgs({"scenario.fersxml", "--vita49", "localhost:4991", "--vita49-fullscale",
											  "1.0", "--vita49-max-udp-payload", "63"});
	REQUIRE_FALSE(payload_too_small);
	CHECK(payload_too_small.error() == "VITA49 max UDP payload must be between 64 and 65507 bytes");

	const auto payload_too_large = parseArgs({"scenario.fersxml", "--vita49", "localhost:4991", "--vita49-fullscale",
											  "1.0", "--vita49-max-udp-payload", "65508"});
	REQUIRE_FALSE(payload_too_large);
	CHECK(payload_too_large.error() == "VITA49 max UDP payload must be between 64 and 65507 bytes");

	const auto payload_without_endpoint = parseArgs({"scenario.fersxml", "--vita49-max-udp-payload", "900"});
	REQUIRE_FALSE(payload_without_endpoint);
	CHECK(payload_without_endpoint.error() == "--vita49-max-udp-payload requires --vita49");

	const auto zero_queue = parseArgs(
		{"scenario.fersxml", "--vita49", "localhost:4991", "--vita49-fullscale", "1.0", "--vita49-queue-depth", "0"});
	REQUIRE_FALSE(zero_queue);
	CHECK(zero_queue.error() == "VITA49 queue depth must be in the range 1..4294967295");

	const auto queue_without_endpoint = parseArgs({"scenario.fersxml", "--vita49-queue-depth", "1"});
	REQUIRE_FALSE(queue_without_endpoint);
	CHECK(queue_without_endpoint.error() == "--vita49-queue-depth requires --vita49");
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
