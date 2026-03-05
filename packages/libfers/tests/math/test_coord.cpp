#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "math/coord.h"

using namespace math;
using Catch::Matchers::WithinAbs;

TEST_CASE("Coord operations", "[math][coord]")
{
	SECTION("Initialization and assignment")
	{
		Coord c;
		c = 5.0;
		REQUIRE_THAT(c.t, WithinAbs(5.0, 1e-9));
		REQUIRE_THAT(c.pos.x, WithinAbs(5.0, 1e-9));
		REQUIRE_THAT(c.pos.y, WithinAbs(5.0, 1e-9));
		REQUIRE_THAT(c.pos.z, WithinAbs(5.0, 1e-9));
	}

	SECTION("Comparison")
	{
		Coord c1{Vec3(1, 2, 3), 1.0};
		Coord c2{Vec3(4, 5, 6), 2.0};
		REQUIRE(c1 < c2);
		REQUIRE_FALSE(c2 < c1);
	}

	SECTION("Arithmetic operations")
	{
		Coord c1{Vec3(1, 2, 3), 1.0};
		Coord c2{Vec3(4, 5, 6), 1.0};

		Coord add = c1 + c2;
		REQUIRE_THAT(add.pos.x, WithinAbs(5.0, 1e-9));
		REQUIRE_THAT(add.pos.y, WithinAbs(7.0, 1e-9));
		REQUIRE_THAT(add.pos.z, WithinAbs(9.0, 1e-9));
		REQUIRE_THAT(add.t, WithinAbs(1.0, 1e-9));

		Coord sub = c2 - c1;
		REQUIRE_THAT(sub.pos.x, WithinAbs(3.0, 1e-9));
		REQUIRE_THAT(sub.pos.y, WithinAbs(3.0, 1e-9));
		REQUIRE_THAT(sub.pos.z, WithinAbs(3.0, 1e-9));

		Coord mul = c1 * 2.0;
		REQUIRE_THAT(mul.pos.x, WithinAbs(2.0, 1e-9));
		REQUIRE_THAT(mul.pos.y, WithinAbs(4.0, 1e-9));
		REQUIRE_THAT(mul.pos.z, WithinAbs(6.0, 1e-9));
	}
}

TEST_CASE("RotationCoord operations", "[math][coord]")
{
	SECTION("Initialization and assignment")
	{
		RotationCoord rc(2.0);
		REQUIRE_THAT(rc.azimuth, WithinAbs(2.0, 1e-9));
		REQUIRE_THAT(rc.elevation, WithinAbs(2.0, 1e-9));
		REQUIRE_THAT(rc.t, WithinAbs(2.0, 1e-9));

		rc = 3.0;
		REQUIRE_THAT(rc.azimuth, WithinAbs(3.0, 1e-9));
	}

	SECTION("Arithmetic operations")
	{
		RotationCoord rc1(1.0, 2.0, 1.0);
		RotationCoord rc2(3.0, 4.0, 1.0);

		RotationCoord add = rc1 + rc2;
		REQUIRE_THAT(add.azimuth, WithinAbs(4.0, 1e-9));
		REQUIRE_THAT(add.elevation, WithinAbs(6.0, 1e-9));

		RotationCoord mul = rc1 * 2.0;
		REQUIRE_THAT(mul.azimuth, WithinAbs(2.0, 1e-9));
		REQUIRE_THAT(mul.elevation, WithinAbs(4.0, 1e-9));
	}
}
