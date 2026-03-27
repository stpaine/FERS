// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <fstream>
#include <string>

#include "core/parameters.h"
#include "core/world.h"
#include "math/path.h"
#include "radar/platform.h"
#include "radar/target.h"
#include "serial/kml_generator.h"

using Catch::Matchers::ContainsSubstring;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
	};

	std::string getTempKmlPath(const std::string& name)
	{
		return (std::filesystem::temp_directory_path() / name).string();
	}

	std::string readFileContent(const std::string& path)
	{
		std::ifstream in(path);
		if (!in.is_open())
			return "";
		return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
	}

	void populateMinimalWorld(core::World& world, const math::Vec3& pos)
	{
		auto plat = std::make_unique<radar::Platform>("TestPlat");

		plat->getMotionPath()->addCoord({pos, 0.0});
		plat->getMotionPath()->finalize();
		plat->getRotationPath()->addCoord({0.0, 0.0, 0.0});
		plat->getRotationPath()->finalize();

		// Attach a target so the platform is picked up by the KML generator
		auto tgt = radar::createIsoTarget(plat.get(), "TestTarget", 1.0, 42);

		world.add(std::move(tgt));
		world.add(std::move(plat));
	}
}

TEST_CASE("KmlGenerator: Generates valid KML file in ENU frame", "[serial][kml][facade]")
{
	ParamGuard guard;
	params::params.reset();
	params::setCoordinateSystem(params::CoordinateFrame::ENU, 0, true);
	params::setOrigin(-33.957652, 18.4611991, 111.01); // UCT

	core::World world;
	// ENU coordinates are local offsets in meters
	populateMinimalWorld(world, {100.0, 200.0, 300.0});

	std::string out_path = getTempKmlPath("test_enu.kml");
	std::filesystem::remove(out_path);

	REQUIRE(serial::KmlGenerator::generateKml(world, out_path) == true);

	std::string content = readFileContent(out_path);
	REQUIRE_THAT(content, ContainsSubstring("<?xml version=\"1.0\""));
	REQUIRE_THAT(content, ContainsSubstring("<name>TestPlat</name>"));

	std::filesystem::remove(out_path);
}

TEST_CASE("KmlGenerator: Generates valid KML file in UTM frame", "[serial][kml][facade]")
{
	ParamGuard guard;
	params::params.reset();
	// Zone 34S (South Africa)
	params::setCoordinateSystem(params::CoordinateFrame::UTM, 34, false);

	core::World world;
	// Valid UTM coordinates: Easting = 500000 (center), Northing = 6234000 (~34 deg South)
	populateMinimalWorld(world, {500000.0, 6234000.0, 300.0});

	std::string out_path = getTempKmlPath("test_utm.kml");
	std::filesystem::remove(out_path);

	REQUIRE(serial::KmlGenerator::generateKml(world, out_path) == true);

	std::string content = readFileContent(out_path);
	REQUIRE_THAT(content, ContainsSubstring("<?xml version=\"1.0\""));
	REQUIRE_THAT(content, ContainsSubstring("<name>TestPlat</name>"));

	std::filesystem::remove(out_path);
}

TEST_CASE("KmlGenerator: Generates valid KML file in ECEF frame", "[serial][kml][facade]")
{
	ParamGuard guard;
	params::params.reset();
	params::setCoordinateSystem(params::CoordinateFrame::ECEF, 0, true);

	core::World world;
	// Valid ECEF coordinates: Surface of the Earth at the equator/prime meridian
	populateMinimalWorld(world, {6378137.0, 0.0, 0.0});

	std::string out_path = getTempKmlPath("test_ecef.kml");
	std::filesystem::remove(out_path);

	REQUIRE(serial::KmlGenerator::generateKml(world, out_path) == true);

	std::string content = readFileContent(out_path);
	REQUIRE_THAT(content, ContainsSubstring("<?xml version=\"1.0\""));
	REQUIRE_THAT(content, ContainsSubstring("<name>TestPlat</name>"));

	std::filesystem::remove(out_path);
}

TEST_CASE("KmlGenerator: Fails gracefully on invalid file path", "[serial][kml][facade]")
{
	ParamGuard guard;
	core::World world;

	std::string bad_path = "/root/invalid_dir_name_12345/test.kml";

	REQUIRE(serial::KmlGenerator::generateKml(world, bad_path) == false);
}

TEST_CASE("KmlGenerator: Catches exceptions from invalid UTM zones", "[serial][kml][facade]")
{
	ParamGuard guard;
	params::params.reset();
	// Zone 99 is invalid and will cause GeographicLib to throw an exception
	params::setCoordinateSystem(params::CoordinateFrame::UTM, 99, true);

	core::World world;
	populateMinimalWorld(world, {500000.0, 6234000.0, 300.0});

	std::string out_path = getTempKmlPath("test_utm_fail.kml");

	// The KmlGenerator should catch the GeographicLib exception and return false
	REQUIRE(serial::KmlGenerator::generateKml(world, out_path) == false);

	std::filesystem::remove(out_path);
}
