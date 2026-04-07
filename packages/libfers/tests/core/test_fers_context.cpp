// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

#include <catch2/catch_test_macros.hpp>
#include <random>

#include "core/fers_context.h"

TEST_CASE("FersContext constructs world and seeder", "[core][context]")
{
	FersContext context;

	core::World* world_ptr = context.getWorld();
	REQUIRE(world_ptr != nullptr);
	REQUIRE(context.getWorld() == world_ptr);
}

TEST_CASE("FersContext master seeder is deterministic", "[core][context]")
{
	FersContext context;
	std::mt19937& seeder = context.getMasterSeeder();

	seeder.seed(12345);
	const auto value1 = seeder();
	const auto value2 = seeder();

	seeder.seed(12345);
	REQUIRE(seeder() == value1);
	REQUIRE(seeder() == value2);
}

TEST_CASE("FersContext manages output directory", "[core][context]")
{
	FersContext context;

	// Should default to current directory
	REQUIRE(context.getOutputDir() == ".");

	context.setOutputDir("/custom/output/path");
	REQUIRE(context.getOutputDir() == "/custom/output/path");
}
