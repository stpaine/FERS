// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "cli_runner.h"

namespace
{
	std::filesystem::path uniqueTempPath(const std::string& prefix, const std::string& extension = "")
	{
		return std::filesystem::temp_directory_path() /
			(prefix + "_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + extension);
	}

	struct ScopedPath
	{
		std::filesystem::path path;

		explicit ScopedPath(std::filesystem::path value) : path(std::move(value)) {}

		~ScopedPath()
		{
			std::error_code ec;
			std::filesystem::remove_all(path, ec);
		}
	};

	void writeTextFile(const std::filesystem::path& path, const std::string& content)
	{
		std::ofstream out(path, std::ios::binary);
		REQUIRE(out.is_open());
		out << content;
	}

	std::string minimalScenarioXml(const std::string& simulationName = "CLI KML Scenario")
	{
		return "<simulation name=\"" + simulationName +
			"\">\n"
			"  <parameters>\n"
			"    <starttime>0</starttime>\n"
			"    <endtime>1</endtime>\n"
			"    <rate>1000</rate>\n"
			"    <origin latitude=\"-33.957652\" longitude=\"18.4611991\" altitude=\"111.01\"/>\n"
			"    <coordinatesystem frame=\"ENU\"/>\n"
			"  </parameters>\n"
			"  <waveform name=\"cli_wave\">\n"
			"    <power>1000</power>\n"
			"    <carrier_frequency>1000000000</carrier_frequency>\n"
			"    <cw/>\n"
			"  </waveform>\n"
			"  <timing name=\"cli_clock\">\n"
			"    <frequency>1000000</frequency>\n"
			"  </timing>\n"
			"  <antenna name=\"cli_iso\" pattern=\"isotropic\"/>\n"
			"  <platform name=\"cli_sensor\">\n"
			"    <motionpath interpolation=\"static\">\n"
			"      <positionwaypoint>\n"
			"        <x>0</x>\n"
			"        <y>0</y>\n"
			"        <altitude>100</altitude>\n"
			"        <time>0</time>\n"
			"      </positionwaypoint>\n"
			"    </motionpath>\n"
			"    <rotationpath interpolation=\"static\">\n"
			"      <rotationwaypoint>\n"
			"        <azimuth>0</azimuth>\n"
			"        <elevation>0</elevation>\n"
			"        <time>0</time>\n"
			"      </rotationwaypoint>\n"
			"    </rotationpath>\n"
			"    <transmitter name=\"cli_tx\" waveform=\"cli_wave\" antenna=\"cli_iso\" timing=\"cli_clock\">\n"
			"      <cw_mode/>\n"
			"    </transmitter>\n"
			"  </platform>\n"
			"</simulation>\n";
	}

	int runCliWithArgs(const std::initializer_list<std::string>& args)
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

		return core::runCli(static_cast<int>(argv.size()), argv.data());
	}
}

TEST_CASE("runCli returns nonzero when KML generation fails", "[fers-cli][runner][kml]")
{
	const auto scenario_path = uniqueTempPath("fers_cli_kml_exit", ".fersxml");
	ScopedPath scenario_guard(scenario_path);
	writeTextFile(scenario_path, minimalScenarioXml());

	const auto missing_parent = uniqueTempPath("fers_cli_missing_kml_dir");
	ScopedPath missing_parent_guard(missing_parent);
	std::error_code ec;
	std::filesystem::remove_all(missing_parent, ec);
	const auto kml_path = missing_parent / "out.kml";

	const int exit_code = runCliWithArgs({scenario_path.string(), "--log-level=FATAL", "--kml=" + kml_path.string()});

	CHECK(exit_code == 1);
	CHECK_FALSE(std::filesystem::exists(kml_path));
}
