#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <libfers/api.h>

#include "api_test_helpers.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

TEST_CASE("API motion path interpolation validates arguments", "[api][paths]")
{
	api_test::clearLastError();
	const fers_motion_waypoint_t point{0.0, 1.0, 2.0, 3.0};

	SECTION("null waypoints")
	{
		api_test::MotionPath path(fers_get_interpolated_motion_path(nullptr, 1, FERS_INTERP_LINEAR, 2));
		REQUIRE(path.get() == nullptr);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("waypoints cannot be null"));
	}

	SECTION("zero waypoint count")
	{
		api_test::MotionPath path(fers_get_interpolated_motion_path(&point, 0, FERS_INTERP_LINEAR, 2));
		REQUIRE(path.get() == nullptr);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("counts must be > 0"));
	}

	SECTION("zero output count")
	{
		api_test::MotionPath path(fers_get_interpolated_motion_path(&point, 1, FERS_INTERP_LINEAR, 0));
		REQUIRE(path.get() == nullptr);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("counts must be > 0"));
	}

	SECTION("cubic needs more waypoints")
	{
		api_test::MotionPath path(fers_get_interpolated_motion_path(&point, 1, FERS_INTERP_CUBIC, 2));
		REQUIRE(path.get() == nullptr);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("Cubic interpolation requires at least 2 waypoints"));
	}
}

TEST_CASE("API motion path interpolation returns constant positions for static input", "[api][paths]")
{
	api_test::clearLastError();
	const fers_motion_waypoint_t waypoint{5.0, 1.5, -2.0, 8.25};
	api_test::MotionPath path(fers_get_interpolated_motion_path(&waypoint, 1, FERS_INTERP_STATIC, 4));

	REQUIRE(path.get() != nullptr);
	REQUIRE(path.get()->count == 4u);

	for (size_t i = 0; i < path.get()->count; ++i)
	{
		const auto& point = path.get()->points[i];
		REQUIRE_THAT(point.x, WithinAbs(1.5, 1e-12));
		REQUIRE_THAT(point.y, WithinAbs(-2.0, 1e-12));
		REQUIRE_THAT(point.z, WithinAbs(8.25, 1e-12));
		REQUIRE_THAT(point.vx, WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(point.vy, WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(point.vz, WithinAbs(0.0, 1e-12));
	}
}

TEST_CASE("API motion path interpolation falls back to repeated positions for non-increasing time spans",
		  "[api][paths]")
{
	api_test::clearLastError();
	const fers_motion_waypoint_t waypoints[] = {
		{1.0, 2.0, 3.0, 4.0},
		{1.0, 2.0, 3.0, 4.0},
	};
	api_test::MotionPath path(fers_get_interpolated_motion_path(waypoints, 2, FERS_INTERP_LINEAR, 3));

	REQUIRE(path.get() != nullptr);
	REQUIRE(path.get()->count == 3u);

	for (size_t i = 0; i < path.get()->count; ++i)
	{
		const auto& point = path.get()->points[i];
		REQUIRE_THAT(point.x, WithinAbs(2.0, 1e-12));
		REQUIRE_THAT(point.y, WithinAbs(3.0, 1e-12));
		REQUIRE_THAT(point.z, WithinAbs(4.0, 1e-12));
	}
}

TEST_CASE("API motion path interpolation returns expected linear positions and velocities", "[api][paths]")
{
	api_test::clearLastError();
	const fers_motion_waypoint_t waypoints[] = {
		{0.0, 0.0, 0.0, 0.0},
		{10.0, 10.0, 20.0, 30.0},
	};
	api_test::MotionPath path(fers_get_interpolated_motion_path(waypoints, 2, FERS_INTERP_LINEAR, 3));

	REQUIRE(path.get() != nullptr);
	REQUIRE(path.get()->count == 3u);

	const auto& start = path.get()->points[0];
	const auto& middle = path.get()->points[1];
	const auto& finish = path.get()->points[2];

	REQUIRE_THAT(start.x, WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(start.y, WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(start.z, WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(middle.x, WithinAbs(5.0, 1e-12));
	REQUIRE_THAT(middle.y, WithinAbs(10.0, 1e-12));
	REQUIRE_THAT(middle.z, WithinAbs(15.0, 1e-12));
	REQUIRE_THAT(finish.x, WithinAbs(10.0, 1e-12));
	REQUIRE_THAT(finish.y, WithinAbs(20.0, 1e-12));
	REQUIRE_THAT(finish.z, WithinAbs(30.0, 1e-12));

	for (const auto* point : {&start, &middle, &finish})
	{
		REQUIRE_THAT(point->vx, WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(point->vy, WithinAbs(2.0, 1e-12));
		REQUIRE_THAT(point->vz, WithinAbs(3.0, 1e-12));
	}
}

TEST_CASE("API motion path interpolation maps cubic enum to cubic behavior", "[api][paths]")
{
	api_test::clearLastError();
	const fers_motion_waypoint_t waypoints[] = {
		{0.0, 0.0, 0.0, 0.0},
		{1.0, 1.0, 0.0, 0.0},
		{2.0, 2.0, 0.0, 0.0},
	};
	api_test::MotionPath path(fers_get_interpolated_motion_path(waypoints, 3, FERS_INTERP_CUBIC, 3));

	REQUIRE(path.get() != nullptr);
	REQUIRE(path.get()->count == 3u);
	REQUIRE_THAT(path.get()->points[1].x, WithinAbs(1.0, 1e-9));
	REQUIRE_THAT(path.get()->points[1].y, WithinAbs(0.0, 1e-9));
}

TEST_CASE("API rotation path interpolation validates arguments", "[api][paths]")
{
	api_test::clearLastError();
	const fers_rotation_waypoint_t point{0.0, 0.0, 0.0};

	SECTION("null waypoints")
	{
		api_test::RotationPath path(fers_get_interpolated_rotation_path(nullptr, 1, FERS_INTERP_LINEAR, 2));
		REQUIRE(path.get() == nullptr);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("waypoints cannot be null"));
	}

	SECTION("zero waypoint count")
	{
		api_test::RotationPath path(fers_get_interpolated_rotation_path(&point, 0, FERS_INTERP_LINEAR, 2));
		REQUIRE(path.get() == nullptr);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("counts must be > 0"));
	}

	SECTION("zero output count")
	{
		api_test::RotationPath path(fers_get_interpolated_rotation_path(&point, 1, FERS_INTERP_LINEAR, 0));
		REQUIRE(path.get() == nullptr);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("counts must be > 0"));
	}

	SECTION("cubic needs more waypoints")
	{
		api_test::RotationPath path(fers_get_interpolated_rotation_path(&point, 1, FERS_INTERP_CUBIC, 2));
		REQUIRE(path.get() == nullptr);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("Cubic interpolation requires at least 2 waypoints"));
	}
}

TEST_CASE("API rotation path interpolation returns constant compass angles for static input", "[api][paths]")
{
	api_test::clearLastError();
	const fers_rotation_waypoint_t waypoint{0.0, 15.0, -10.0};
	api_test::RotationPath path(fers_get_interpolated_rotation_path(&waypoint, 1, FERS_INTERP_STATIC, 4));

	REQUIRE(path.get() != nullptr);
	REQUIRE(path.get()->count == 4u);

	for (size_t i = 0; i < path.get()->count; ++i)
	{
		REQUIRE_THAT(path.get()->points[i].azimuth_deg, WithinAbs(15.0, 1e-12));
		REQUIRE_THAT(path.get()->points[i].elevation_deg, WithinAbs(-10.0, 1e-12));
	}
}

TEST_CASE("API rotation path interpolation preserves compass-angle semantics", "[api][paths]")
{
	api_test::clearLastError();
	const fers_rotation_waypoint_t waypoints[] = {
		{0.0, 0.0, -10.0},
		{10.0, 180.0, 20.0},
	};
	api_test::RotationPath path(fers_get_interpolated_rotation_path(waypoints, 2, FERS_INTERP_LINEAR, 3));

	REQUIRE(path.get() != nullptr);
	REQUIRE(path.get()->count == 3u);

	REQUIRE_THAT(path.get()->points[0].azimuth_deg, WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(path.get()->points[1].azimuth_deg, WithinAbs(90.0, 1e-12));
	REQUIRE_THAT(path.get()->points[2].azimuth_deg, WithinAbs(180.0, 1e-12));
	REQUIRE_THAT(path.get()->points[0].elevation_deg, WithinAbs(-10.0, 1e-12));
	REQUIRE_THAT(path.get()->points[1].elevation_deg, WithinAbs(5.0, 1e-12));
	REQUIRE_THAT(path.get()->points[2].elevation_deg, WithinAbs(20.0, 1e-12));
}

TEST_CASE("API rotation path interpolation maps cubic enum to cubic behavior", "[api][paths]")
{
	api_test::clearLastError();
	const fers_rotation_waypoint_t waypoints[] = {
		{0.0, 90.0, 0.0},
		{1.0, 32.7042204869, 57.2957795131},
		{2.0, -24.5915590262, 114.5915590262},
	};
	api_test::RotationPath path(fers_get_interpolated_rotation_path(waypoints, 3, FERS_INTERP_CUBIC, 3));

	REQUIRE(path.get() != nullptr);
	REQUIRE(path.get()->count == 3u);
	REQUIRE_THAT(path.get()->points[1].azimuth_deg, WithinAbs(32.7042204869, 1e-6));
	REQUIRE_THAT(path.get()->points[1].elevation_deg, WithinAbs(57.2957795131, 1e-6));
}

TEST_CASE("API rotation path interpolation preserves unnormalized winding angles", "[api][paths]")
{
	api_test::clearLastError();
	const fers_rotation_waypoint_t waypoints[] = {
		{0.0, -30.0, 0.0},
		{10.0, 390.0, 0.0},
	};
	api_test::RotationPath path(fers_get_interpolated_rotation_path(waypoints, 2, FERS_INTERP_LINEAR, 3));

	REQUIRE(path.get() != nullptr);
	REQUIRE(path.get()->count == 3u);

	REQUIRE_THAT(path.get()->points[0].azimuth_deg, WithinAbs(-30.0, 1e-12));
	REQUIRE_THAT(path.get()->points[1].azimuth_deg, WithinAbs(180.0, 1e-12));
	REQUIRE_THAT(path.get()->points[2].azimuth_deg, WithinAbs(390.0, 1e-12));
}
