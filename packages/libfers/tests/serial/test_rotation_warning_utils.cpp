#include <catch2/catch_test_macros.hpp>
#include <string>

#include "core/parameters.h"
#include "serial/rotation_warning_utils.h"

TEST_CASE("rotation warning heuristic distinguishes obvious unit mismatches", "[serial][rotation_warning]")
{
	using serial::rotation_warning_utils::Confidence;
	using serial::rotation_warning_utils::ValueKind;

	SECTION("degree-like value declared in radians")
	{
		const auto result = serial::rotation_warning_utils::infer_unit_from_value(90.0, ValueKind::Angle);
		REQUIRE(result.inferred_unit == params::RotationAngleUnit::Degrees);
		REQUIRE(result.confidence == Confidence::High);
	}

	SECTION("radian-like value declared in degrees")
	{
		const auto result = serial::rotation_warning_utils::infer_unit_from_value(PI / 2.0, ValueKind::Angle);
		REQUIRE(result.inferred_unit == params::RotationAngleUnit::Radians);
		REQUIRE(result.confidence == Confidence::High);
	}

	SECTION("simple ambiguous value stays below warning threshold")
	{
		const auto result = serial::rotation_warning_utils::infer_unit_from_value(1.0, ValueKind::Angle);
		REQUIRE_FALSE(serial::rotation_warning_utils::should_warn(result.confidence,
																  serial::rotation_warning_utils::kWarningSensitivity));
	}
}

TEST_CASE("rotation warning collector deduplicates repeated messages", "[serial][rotation_warning]")
{
	serial::rotation_warning_utils::clear_captured_warnings();

	serial::rotation_warning_utils::maybe_warn_about_rotation_value(
		90.0, params::RotationAngleUnit::Radians, serial::rotation_warning_utils::ValueKind::Angle, "JSON",
		"platform 'collector' rotation waypoint 0", "azimuth");
	serial::rotation_warning_utils::maybe_warn_about_rotation_value(
		90.0, params::RotationAngleUnit::Radians, serial::rotation_warning_utils::ValueKind::Angle, "JSON",
		"platform 'collector' rotation waypoint 0", "azimuth");

	const auto warnings = serial::rotation_warning_utils::take_captured_warnings();
	REQUIRE(warnings.size() == 1);
	REQUIRE(warnings[0].find("platform 'collector' rotation waypoint 0") != std::string::npos);

	REQUIRE(serial::rotation_warning_utils::take_captured_warnings().empty());
}
