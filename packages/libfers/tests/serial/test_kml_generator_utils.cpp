// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/parameters.h"
#include "core/world.h"
#include "math/path.h"
#include "math/rotation_path.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/kml_generator_utils.h"
#include "signal/radar_signal.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace
{
	serial::kml_generator_utils::KmlContext createTestContext()
	{
		serial::kml_generator_utils::KmlContext ctx;
		ctx.parameters.simulation_name = "TestSim";
		ctx.parameters.start = 0.0;
		// Dummy converter: treats x, y, z directly as lat, lon, alt for easy string matching
		ctx.converter = [](const math::Vec3& pos, double& lat, double& lon, double& alt)
		{
			lat = pos.x;
			lon = pos.y;
			alt = pos.z;
		};
		return ctx;
	}
}

TEST_CASE("KML Math: Sinc Antenna Gain", "[serial][kml][math]")
{
	// At theta = 0, gain should be alpha
	REQUIRE_THAT(serial::kml_generator_utils::sincAntennaGain(0.0, 2.0, 1.0, 1.0), WithinAbs(2.0, 1e-6));

	// At theta = pi, sin(pi)/pi = 0, gain should be 0
	REQUIRE_THAT(serial::kml_generator_utils::sincAntennaGain(PI, 1.0, 1.0, 1.0), WithinAbs(0.0, 1e-6));
}

TEST_CASE("KML Math: 3dB Drop Angles", "[serial][kml][math]")
{
	SECTION("Sinc Antenna")
	{
		// For alpha=1, beta=1, gamma=2, G(theta) = sinc(theta)^2
		// The 3dB drop (half power) occurs when sinc(theta)^2 = 0.5 -> sinc(theta) = 0.707
		// This happens at approx 1.391 rad, which is ~79.7 degrees.
		// The function uses a 1000-point discretization, so we allow a 1.0 degree tolerance.
		double angle = serial::kml_generator_utils::find3DbDropAngle(1.0, 1.0, 2.0);
		REQUIRE_THAT(angle, WithinAbs(79.7, 1.0));
	}

	SECTION("Gaussian Antenna")
	{
		antenna::Gaussian valid_gaussian("gauss", 1.0, 1.0);
		// theta = sqrt(ln(2) / scale) * 180 / PI
		double expected_deg = std::sqrt(std::log(2.0) / 1.0) * 180.0 / PI;
		REQUIRE_THAT(serial::kml_generator_utils::findGaussian3DbDropAngle(&valid_gaussian),
					 WithinAbs(expected_deg, 1e-4));

		antenna::Gaussian invalid_gaussian("gauss_bad", -1.0, 1.0);
		REQUIRE_THAT(serial::kml_generator_utils::findGaussian3DbDropAngle(&invalid_gaussian), WithinAbs(0.0, 1e-6));
	}

	SECTION("Parabolic Antenna")
	{
		antenna::Parabolic valid_parabolic("dish", 1.0);
		// arg = 1.6 * lambda / (PI * D). For lambda=0.1, D=1.0 -> arg = 0.16 / PI
		double arg = 0.16 / PI;
		double expected_deg = std::asin(arg) * 180.0 / PI;
		REQUIRE_THAT(serial::kml_generator_utils::findParabolic3DbDropAngle(&valid_parabolic, 0.1),
					 WithinAbs(expected_deg, 1e-4));

		// Too wide beam (arg > 1.0) caps at 90 degrees
		REQUIRE_THAT(serial::kml_generator_utils::findParabolic3DbDropAngle(&valid_parabolic, 10.0),
					 WithinAbs(90.0, 1e-6));

		// Invalid diameter
		antenna::Parabolic invalid_parabolic("dish_bad", 0.0);
		REQUIRE_THAT(serial::kml_generator_utils::findParabolic3DbDropAngle(&invalid_parabolic, 0.1),
					 WithinAbs(0.0, 1e-6));
	}

	SECTION("SquareHorn Antenna")
	{
		antenna::SquareHorn valid_horn("horn", 1.0);
		// arg = 1.39155 * lambda / (PI * D). For lambda=0.1, D=1.0 -> arg = 0.139155 / PI
		double arg = 0.139155 / PI;
		double expected_deg = std::asin(arg) * 180.0 / PI;
		REQUIRE_THAT(serial::kml_generator_utils::findSquareHorn3DbDropAngle(&valid_horn, 0.1),
					 WithinAbs(expected_deg, 1e-4));

		// Too wide beam (arg > 1.0) caps at 90 degrees
		REQUIRE_THAT(serial::kml_generator_utils::findSquareHorn3DbDropAngle(&valid_horn, 10.0), WithinAbs(90.0, 1e-6));

		// Invalid dimension
		antenna::SquareHorn invalid_horn("horn_bad", -0.5);
		REQUIRE_THAT(serial::kml_generator_utils::findSquareHorn3DbDropAngle(&invalid_horn, 0.1), WithinAbs(0.0, 1e-6));
	}
}

TEST_CASE("KML Geodesic: Coordinate Calculations", "[serial][kml][geodesic]")
{
	SECTION("Destination Coordinate")
	{
		double dest_lat, dest_lon;
		// Start at equator/prime meridian, move East (90 deg) by ~1 degree of longitude (111319.49 meters)
		serial::kml_generator_utils::calculateDestinationCoordinate(0.0, 0.0, 90.0, 111319.49, dest_lat, dest_lon);
		REQUIRE_THAT(dest_lat, WithinAbs(0.0, 1e-4));
		REQUIRE_THAT(dest_lon, WithinAbs(1.0, 1e-4));
	}

	SECTION("Circle Coordinates")
	{
		auto circle = serial::kml_generator_utils::generateCircleCoordinates(0.0, 0.0, 20.0);
		REQUIRE(circle.size() == 100);
		// First point is bearing 0 (North) -> Longitude should be 0, Latitude should be positive
		REQUIRE_THAT(circle[0].second, WithinAbs(0.0, 1e-6));
		REQUIRE(circle[0].first > 0.0);
	}
}

TEST_CASE("KML Formatting: Strings and XML", "[serial][kml][formatting]")
{
	SECTION("formatCoordinates")
	{
		std::string formatted = serial::kml_generator_utils::formatCoordinates(1.2345678, 2.3456789, 3.4567891);
		REQUIRE(formatted == "1.234568,2.345679,3.456789");
	}

	SECTION("writeKmlHeaderAndStyles")
	{
		std::ostringstream oss;
		auto ctx = createTestContext();
		serial::kml_generator_utils::writeKmlHeaderAndStyles(oss, ctx);
		std::string out = oss.str();
		REQUIRE_THAT(out, ContainsSubstring("<?xml version=\"1.0\""));
		REQUIRE_THAT(out, ContainsSubstring("<name>TestSim</name>"));
		REQUIRE_THAT(out, ContainsSubstring("<Style id=\"target\">"));
	}

	SECTION("writePoint")
	{
		std::ostringstream oss;
		serial::kml_generator_utils::writePoint(oss, "  ", "MyPoint", "#style", "1,2,3", 100.0, 0.0);
		std::string out = oss.str();
		REQUIRE_THAT(out, ContainsSubstring("<name>MyPoint</name>"));
		REQUIRE_THAT(out, ContainsSubstring("<styleUrl>#style</styleUrl>"));
		REQUIRE_THAT(out, ContainsSubstring("<coordinates>1,2,3</coordinates>"));
		REQUIRE_THAT(out, ContainsSubstring("<extrude>1</extrude>")); // Because 100.0 > 0.0
	}
}

TEST_CASE("KML Logic: Platform Styles and Primary Radar", "[serial][kml][logic]")
{
	radar::Platform plat("plat");
	radar::Transmitter tx(&plat, "tx", radar::OperationMode::CW_MODE);
	radar::Receiver rx(&plat, "rx", 42, radar::OperationMode::CW_MODE);
	auto tgt = radar::createIsoTarget(&plat, "tgt", 1.0, 42);

	SECTION("getPlacemarkStyleForPlatform")
	{
		REQUIRE(serial::kml_generator_utils::getPlacemarkStyleForPlatform({&rx}) == "#receiver");
		REQUIRE(serial::kml_generator_utils::getPlacemarkStyleForPlatform({&tx}) == "#transmitter");
		REQUIRE(serial::kml_generator_utils::getPlacemarkStyleForPlatform({tgt.get()}) == "#target");
		REQUIRE(serial::kml_generator_utils::getPlacemarkStyleForPlatform({&tx, tgt.get()}) == "#transmitter");
		REQUIRE(serial::kml_generator_utils::getPlacemarkStyleForPlatform({&rx, &tx}) == "#receiver");
	}

	SECTION("getPrimaryRadar")
	{
		REQUIRE(serial::kml_generator_utils::getPrimaryRadar({tgt.get()}) == nullptr);
		REQUIRE(serial::kml_generator_utils::getPrimaryRadar({tgt.get(), &tx}) == &tx);
	}
}

TEST_CASE("KML Generation: Antenna Visualization", "[serial][kml][generation]")
{
	auto ctx = createTestContext();
	radar::Platform plat("plat");
	plat.getMotionPath()->addCoord({math::Vec3(0, 0, 100), 0.0});
	plat.getMotionPath()->finalize();
	plat.getRotationPath()->addCoord({0.0, 0.0, 0.0});
	plat.getRotationPath()->finalize();

	radar::Transmitter tx(&plat, "tx", radar::OperationMode::CW_MODE);

	SECTION("Isotropic Antenna")
	{
		antenna::Isotropic iso("iso");
		tx.setAntenna(&iso);

		std::ostringstream oss;
		serial::kml_generator_utils::generateAntennaKml(oss, &plat, &tx, ctx, "");
		std::string out = oss.str();

		REQUIRE_THAT(out, ContainsSubstring("<name>Isotropic pattern range</name>"));
		REQUIRE_THAT(out, ContainsSubstring("<Polygon>"));
	}

	SECTION("Directional Antenna (Gaussian)")
	{
		antenna::Gaussian gauss("gauss", 1.0, 1.0);
		tx.setAntenna(&gauss);

		std::ostringstream oss;
		serial::kml_generator_utils::generateAntennaKml(oss, &plat, &tx, ctx, "");
		std::string out = oss.str();

		REQUIRE_THAT(out, ContainsSubstring("<name>Antenna Boresight</name>"));
		REQUIRE_THAT(out, ContainsSubstring("<name>Antenna 3dB Beamwidth</name>"));
		REQUIRE_THAT(out, ContainsSubstring("<name>Antenna Arrow</name>"));
	}
}

TEST_CASE("KML Generation: Path Visualization", "[serial][kml][generation]")
{
	auto ctx = createTestContext();
	radar::Platform plat("plat");

	SECTION("Static Path")
	{
		plat.getMotionPath()->setInterp(math::Path::InterpType::INTERP_STATIC);
		plat.getMotionPath()->addCoord({math::Vec3(10, 20, 30), 0.0});
		plat.getMotionPath()->finalize();

		std::ostringstream oss;
		serial::kml_generator_utils::generatePlatformPathKml(oss, &plat, "#style", 0.0, ctx, "");
		std::string out = oss.str();

		REQUIRE_THAT(out, ContainsSubstring("<Point>"));
		REQUIRE_THAT(out, ContainsSubstring("<coordinates>20.000000,10.000000,30.000000</coordinates>"));
	}

	SECTION("Dynamic Path")
	{
		plat.getMotionPath()->setInterp(math::Path::InterpType::INTERP_LINEAR);
		plat.getMotionPath()->addCoord({math::Vec3(0, 0, 100), 0.0});
		plat.getMotionPath()->addCoord({math::Vec3(10, 10, 100), 10.0});
		plat.getMotionPath()->finalize();

		std::ostringstream oss;
		serial::kml_generator_utils::generatePlatformPathKml(oss, &plat, "#style", 0.0, ctx, "");
		std::string out = oss.str();

		REQUIRE_THAT(out, ContainsSubstring("<gx:Track>"));
		REQUIRE_THAT(out, ContainsSubstring("<when>0</when>"));
		REQUIRE_THAT(out, ContainsSubstring("<when>10</when>"));
		REQUIRE_THAT(out, ContainsSubstring("Start: plat"));
		REQUIRE_THAT(out, ContainsSubstring("End: plat"));
	}
}

TEST_CASE("KML Generation: Full Stream Facade", "[serial][kml][generation]")
{
	auto ctx = createTestContext();
	core::World world;

	auto plat = std::make_unique<radar::Platform>("UAV");
	plat->getMotionPath()->addCoord({math::Vec3(0, 0, 500), 0.0});
	plat->getMotionPath()->finalize();
	plat->getRotationPath()->addCoord({0.0, 0.0, 0.0});
	plat->getRotationPath()->finalize();

	auto tx = std::make_unique<radar::Transmitter>(plat.get(), "Tx", radar::OperationMode::CW_MODE);
	auto ant = std::make_unique<antenna::Isotropic>("Iso");
	tx->setAntenna(ant.get());

	world.add(std::move(tx));
	world.add(std::move(ant));
	world.add(std::move(plat));

	std::ostringstream oss;
	serial::kml_generator_utils::generateKmlToStream(oss, world, ctx);
	std::string out = oss.str();

	REQUIRE_THAT(out, ContainsSubstring("<?xml version=\"1.0\""));
	REQUIRE_THAT(out, ContainsSubstring("<name>Reference Coordinate</name>"));
	REQUIRE_THAT(out, ContainsSubstring("<name>UAV</name>"));
	REQUIRE_THAT(out, ContainsSubstring("Isotropic pattern range"));
	REQUIRE_THAT(out, ContainsSubstring("</kml>"));
}

TEST_CASE("KML Generation: Negative Azimuth Normalization", "[serial][kml][generation]")
{
	auto ctx = createTestContext();
	radar::Platform plat("plat");
	plat.getMotionPath()->addCoord({math::Vec3(0, 0, 100), 0.0});
	plat.getMotionPath()->finalize();

	// FERS azimuth = PI (180 deg).
	// KML start azimuth = 90 - 180 = -90. Normalized to 270.
	// Arrow heading = (270 + 180) % 360 = 450 % 360 = 90.
	plat.getRotationPath()->addCoord({PI, 0.0, 0.0});
	plat.getRotationPath()->finalize();

	radar::Transmitter tx(&plat, "tx", radar::OperationMode::CW_MODE);
	antenna::Gaussian gauss("gauss", 1.0, 1.0); // Directional to trigger arrow rendering
	tx.setAntenna(&gauss);

	std::ostringstream oss;
	serial::kml_generator_utils::generateAntennaKml(oss, &plat, &tx, ctx, "");
	REQUIRE_THAT(oss.str(), ContainsSubstring("<heading>90</heading>"));
}

TEST_CASE("KML Generation: Antenna Early Returns", "[serial][kml][generation]")
{
	auto ctx = createTestContext();
	radar::Platform plat("plat");
	radar::Transmitter tx(&plat, "tx", radar::OperationMode::CW_MODE);

	std::ostringstream oss;

	// 1. Empty motion path -> Early return
	serial::kml_generator_utils::generateAntennaKml(oss, &plat, &tx, ctx, "");
	REQUIRE(oss.str().empty());

	// 2. Motion path exists, but no antenna -> Early return
	plat.getMotionPath()->addCoord({math::Vec3(0, 0, 100), 0.0});
	plat.getMotionPath()->finalize();
	serial::kml_generator_utils::generateAntennaKml(oss, &plat, &tx, ctx, "");
	REQUIRE(oss.str().empty());
}

TEST_CASE("KML Generation: Antenna Wavelength and Dispatch", "[serial][kml][generation]")
{
	auto ctx = createTestContext();
	ctx.parameters.c = 3e8; // Speed of light

	radar::Platform plat("plat");
	plat.getMotionPath()->addCoord({math::Vec3(0, 0, 100), 0.0});
	plat.getMotionPath()->finalize();
	plat.getRotationPath()->addCoord({0.0, 0.0, 0.0});
	plat.getRotationPath()->finalize();

	auto sig = std::make_unique<fers_signal::CwSignal>();
	// Carrier = 3e8 Hz -> Wavelength = 1.0 m
	fers_signal::RadarSignal wave("wave", 1.0, 3e8, 1.0, std::move(sig), 1);

	radar::Transmitter tx(&plat, "tx", radar::OperationMode::CW_MODE);
	tx.setSignal(&wave);

	SECTION("Parabolic with Transmitter")
	{
		antenna::Parabolic dish("dish", 1.0); // D=1.0, lambda=1.0
		tx.setAntenna(&dish);
		std::ostringstream oss;
		serial::kml_generator_utils::generateAntennaKml(oss, &plat, &tx, ctx, "");
		REQUIRE_THAT(oss.str(), ContainsSubstring("Antenna 3dB Beamwidth"));
	}

	SECTION("SquareHorn with Receiver attached to Transmitter")
	{
		radar::Receiver rx(&plat, "rx", 42, radar::OperationMode::CW_MODE);
		rx.setAttached(&tx);
		antenna::SquareHorn horn("horn", 1.0); // D=1.0, lambda=1.0
		rx.setAntenna(&horn);
		std::ostringstream oss;
		serial::kml_generator_utils::generateAntennaKml(oss, &plat, &rx, ctx, "");
		REQUIRE_THAT(oss.str(), ContainsSubstring("Antenna 3dB Beamwidth"));
	}

	SECTION("Sinc Antenna")
	{
		antenna::Sinc sinc("sinc", 1.0, 1.0, 2.0);
		tx.setAntenna(&sinc);
		std::ostringstream oss;
		serial::kml_generator_utils::generateAntennaKml(oss, &plat, &tx, ctx, "");
		REQUIRE_THAT(oss.str(), ContainsSubstring("Antenna 3dB Beamwidth"));
	}

	SECTION("XmlAntenna (Symbolic)")
	{
		// Create a minimal valid XML file to satisfy the XmlAntenna constructor
		std::string temp_xml = (std::filesystem::temp_directory_path() / "dummy_ant.xml").string();
		std::ofstream out(temp_xml);
		out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
			<< "<antenna>\n"
			<< "  <azimuth>\n"
			<< "    <gainsample><angle>-180</angle><gain>0</gain></gainsample>\n"
			<< "    <gainsample><angle>0</angle><gain>1</gain></gainsample>\n"
			<< "    <gainsample><angle>180</angle><gain>0</gain></gainsample>\n"
			<< "  </azimuth>\n"
			<< "  <elevation>\n"
			<< "    <gainsample><angle>-90</angle><gain>0</gain></gainsample>\n"
			<< "    <gainsample><angle>0</angle><gain>1</gain></gainsample>\n"
			<< "    <gainsample><angle>90</angle><gain>0</gain></gainsample>\n"
			<< "  </elevation>\n"
			<< "</antenna>";
		out.close();

		antenna::XmlAntenna xml_ant("xml_ant", temp_xml);
		tx.setAntenna(&xml_ant);

		std::ostringstream oss;
		serial::kml_generator_utils::generateAntennaKml(oss, &plat, &tx, ctx, "");

		// Symbolic antennas should have the Boresight line but NO 3dB Beamwidth lines
		REQUIRE_THAT(oss.str(), ContainsSubstring("Antenna Boresight"));
		REQUIRE_THAT(oss.str(), !ContainsSubstring("Antenna 3dB Beamwidth"));

		std::filesystem::remove(temp_xml);
	}
}

TEST_CASE("KML Generation: Dynamic Path Edge Cases", "[serial][kml][generation]")
{
	auto ctx = createTestContext();
	radar::Platform plat("plat");
	plat.getMotionPath()->setInterp(math::Path::InterpType::INTERP_LINEAR);

	SECTION("Single coordinate (zero duration)")
	{
		plat.getMotionPath()->addCoord({math::Vec3(10, 20, 30), 5.0});
		plat.getMotionPath()->finalize();

		std::ostringstream oss;
		serial::kml_generator_utils::generatePlatformPathKml(oss, &plat, "#style", 0.0, ctx, "");
		std::string out = oss.str();

		// Should have gx:Track
		REQUIRE_THAT(out, ContainsSubstring("<gx:Track>"));
		// Should have exactly one <when>
		REQUIRE_THAT(out, ContainsSubstring("<when>5</when>"));
		// Should NOT have Start/End placemarks (endpoints early return due to size <= 1)
		REQUIRE_THAT(out, !ContainsSubstring("Start: plat"));
	}

	SECTION("Two coordinates, same time (zero duration)")
	{
		plat.getMotionPath()->addCoord({math::Vec3(10, 20, 30), 5.0});
		plat.getMotionPath()->addCoord({math::Vec3(10, 20, 30), 5.0});
		plat.getMotionPath()->finalize();

		std::ostringstream oss;
		serial::kml_generator_utils::generatePlatformPathKml(oss, &plat, "#style", 0.0, ctx, "");
		std::string out = oss.str();

		// Should have gx:Track
		REQUIRE_THAT(out, ContainsSubstring("<gx:Track>"));
		// Should have exactly one <when> because time_diff <= 0.0
		REQUIRE_THAT(out, ContainsSubstring("<when>5</when>"));
		// Should have Start/End placemarks because size > 1
		REQUIRE_THAT(out, ContainsSubstring("Start: plat"));
	}
}

TEST_CASE("KML Generation: processPlatform Early Return", "[serial][kml][generation]")
{
	auto ctx = createTestContext();
	radar::Platform plat("plat");
	std::vector<const radar::Object*> objects;

	std::ostringstream oss;
	// Platform has an empty motion path, so it should return immediately
	serial::kml_generator_utils::processPlatform(oss, &plat, objects, ctx, 0.0, "");
	REQUIRE(oss.str().empty());
}
