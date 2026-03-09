#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "math/path_utils.h"
#include "math/rotation_path.h"

using namespace math;
using Catch::Matchers::WithinAbs;

constexpr RealType PI_VAL = std::numbers::pi_v<RealType>;

TEST_CASE("RotationPath: State Management and Errors", "[math][rotation_path]")
{
	RotationPath rp;

	REQUIRE(rp.getCoords().empty());

	rp.addCoord(RotationCoord(0.0, 0.0, 0.0));
	REQUIRE_THROWS_AS(rp.getPosition(0.0), PathException);

	RotationPath rp_invalid(static_cast<RotationPath::InterpType>(999));
	rp_invalid.addCoord(RotationCoord(0.0, 0.0, 0.0));
	rp_invalid.finalize();
	REQUIRE_THROWS_AS(rp_invalid.getPosition(0.0), PathException);
}

TEST_CASE("RotationPath: Accessors and State", "[math][rotation_path]")
{
	RotationPath rp;

	SECTION("Type Accessor")
	{
		// Default is Static
		REQUIRE(rp.getType() == RotationPath::InterpType::INTERP_STATIC);

		rp.setInterp(RotationPath::InterpType::INTERP_LINEAR);
		REQUIRE(rp.getType() == RotationPath::InterpType::INTERP_LINEAR);
	}

	SECTION("Start Coordinate Accessor")
	{
		RotationCoord start_val(1.5, 2.5, 0.0);
		rp.setStart(start_val);

		RotationCoord s = rp.getStart();
		REQUIRE_THAT(s.azimuth, WithinAbs(1.5, 1e-9));
		REQUIRE_THAT(s.elevation, WithinAbs(2.5, 1e-9));
	}

	SECTION("Rate Accessor")
	{
		RotationCoord rate_val(0.1, 0.2, 0.0);
		rp.setRate(rate_val);

		RotationCoord r = rp.getRate();
		REQUIRE_THAT(r.azimuth, WithinAbs(0.1, 1e-9));
		REQUIRE_THAT(r.elevation, WithinAbs(0.2, 1e-9));
	}
}

TEST_CASE("RotationPath: Static Interpolation", "[math][rotation_path]")
{
	RotationPath rp(RotationPath::InterpType::INTERP_STATIC);
	rp.addCoord(RotationCoord(0.0, 0.0, 0.0));
	rp.addCoord(RotationCoord(PI_VAL, PI_VAL / 2, 1.0));
	rp.finalize();

	SVec3 pos = rp.getPosition(0.5);
	REQUIRE_THAT(pos.azimuth, WithinAbs(0.0, 1e-9));
}

TEST_CASE("RotationPath: Cubic Interpolation", "[math][rotation_path]")
{
	RotationPath rp(RotationPath::InterpType::INTERP_CUBIC);
	rp.addCoord(RotationCoord(0.0, 0.0, 0.0));
	rp.addCoord(RotationCoord(1.0, 1.0, 1.0));
	rp.addCoord(RotationCoord(2.0, 2.0, 2.0));
	rp.finalize();

	SVec3 pos = rp.getPosition(0.5);
	REQUIRE_THAT(pos.azimuth, WithinAbs(0.5, 1e-9));
}

TEST_CASE("RotationPath: Constant Rate", "[math][rotation_path]")
{
	RotationPath rp;
	RotationCoord start(0.0, 0.0, 0.0);
	RotationCoord rate(PI_VAL, 0.0, 0.0);

	rp.setConstantRate(start, rate);

	SECTION("Time 0.5 (Half Rotation)")
	{
		SVec3 pos = rp.getPosition(0.5);
		REQUIRE_THAT(pos.azimuth, WithinAbs(PI_VAL * 0.5, 1e-9));
	}

	SECTION("Time 2.0 (Full Rotation Wrap)")
	{
		SVec3 pos = rp.getPosition(2.0);
		REQUIRE_THAT(pos.azimuth, WithinAbs(0.0, 1e-9));
	}
}

TEST_CASE("RotationPath: Linear Interpolation", "[math][rotation_path]")
{
	RotationPath rp(RotationPath::InterpType::INTERP_LINEAR);
	rp.addCoord(RotationCoord(0.0, 0.0, 0.0));
	rp.addCoord(RotationCoord(PI_VAL, PI_VAL / 2, 1.0));
	rp.finalize();

	SECTION("Midpoint")
	{
		SVec3 pos = rp.getPosition(0.5);
		REQUIRE_THAT(pos.azimuth, WithinAbs(PI_VAL * 0.5, 1e-9));
		REQUIRE_THAT(pos.elevation, WithinAbs(PI_VAL * 0.25, 1e-9));
	}
}
