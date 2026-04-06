// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <optional>
#include <string>

#include "cli_paths.h"

TEST_CASE("resolveOutputDir uses the script parent directory by default", "[fers-cli][paths]")
{
	CHECK(core::resolveOutputDir("inputs/scenario.fersxml", std::nullopt) == std::filesystem::path("inputs"));
}

TEST_CASE("resolveOutputDir falls back to the current directory for bare filenames", "[fers-cli][paths]")
{
	CHECK(core::resolveOutputDir("scenario.fersxml", std::nullopt) == std::filesystem::path("."));
}

TEST_CASE("resolveOutputDir respects an explicit override", "[fers-cli][paths]")
{
	CHECK(core::resolveOutputDir("inputs/scenario.fersxml", std::optional<std::string>{"results"}) ==
		  std::filesystem::path("results"));
}

TEST_CASE("resolveKmlOutputPath defaults to the scenario name in the output directory", "[fers-cli][paths]")
{
	CHECK(core::resolveKmlOutputPath("inputs/scenario.fersxml", std::filesystem::path("results"), std::nullopt) ==
		  std::filesystem::path("results/scenario.kml"));
}

TEST_CASE("resolveKmlOutputPath places filename-only overrides in the output directory", "[fers-cli][paths]")
{
	CHECK(core::resolveKmlOutputPath("inputs/scenario.fersxml", std::filesystem::path("results"),
									 std::optional<std::string>{"preview.kml"}) ==
		  std::filesystem::path("results/preview.kml"));
}

TEST_CASE("resolveKmlOutputPath preserves relative paths with parent directories", "[fers-cli][paths]")
{
	CHECK(core::resolveKmlOutputPath("inputs/scenario.fersxml", std::filesystem::path("results"),
									 std::optional<std::string>{"kml/preview.kml"}) ==
		  std::filesystem::path("kml/preview.kml"));
}

TEST_CASE("resolveKmlOutputPath preserves absolute overrides", "[fers-cli][paths]")
{
	CHECK(core::resolveKmlOutputPath("inputs/scenario.fersxml", std::filesystem::path("results"),
									 std::optional<std::string>{"/tmp/preview.kml"}) ==
		  std::filesystem::path("/tmp/preview.kml"));
}
