#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "math/coord.h"

using namespace math;
using Catch::Matchers::WithinAbs;

TEST_CASE("Coord Struct Operations", "[math][coord]")
{
	SECTION("Assignment from Scalar")
	{
		Coord c;
		c = 5.0;
		REQUIRE_THAT(c.t, WithinAbs(5.0, 1e-9));
		REQUIRE_THAT(c.pos.x, WithinAbs(5.0, 1e-9));
		REQUIRE_THAT(c.pos.y, WithinAbs(5.0, 1e-9));
		REQUIRE_THAT(c.pos.z, WithinAbs(5.0, 1e-9));
	}

	SECTION("Comparison (Time based)")
	{
		Coord c1{Vec3(0, 0, 0), 1.0};
		Coord c2{Vec3(0, 0, 0), 2.0};
		REQUIRE(c1 < c2);
		REQUIRE_FALSE(c2 < c1);
	}

	SECTION("Coord Arithmetic")
	{
		Coord c1{Vec3(1, 2, 3), 1.0};
		Coord c2{Vec3(4, 5, 6), 1.0};

		Coord sum = c1 + c2;
		REQUIRE_THAT(sum.pos.x, WithinAbs(5.0, 1e-9));
		REQUIRE_THAT(sum.t, WithinAbs(1.0, 1e-9));

		Coord diff = c2 - c1;
		REQUIRE_THAT(diff.pos.x, WithinAbs(3.0, 1e-9));

		Coord prod = c1 * c2;
		REQUIRE_THAT(prod.pos.x, WithinAbs(4.0, 1e-9));

		Coord div = c2 / c1;
		REQUIRE_THAT(div.pos.x, WithinAbs(4.0, 1e-9));
	}

	SECTION("Coord Scalar Arithmetic")
	{
		Coord c{Vec3(2, 4, 8), 10.0};

		Coord res = c * 2.0;
		REQUIRE_THAT(res.pos.x, WithinAbs(4.0, 1e-9));
		REQUIRE_THAT(res.t, WithinAbs(10.0, 1e-9));

		res = c + 1.0;
		REQUIRE_THAT(res.pos.x, WithinAbs(3.0, 1e-9));

		res = c / 2.0;
		REQUIRE_THAT(res.pos.x, WithinAbs(1.0, 1e-9));

		res = 16.0 / c;
		REQUIRE_THAT(res.pos.z, WithinAbs(2.0, 1e-9));
	}
}

TEST_CASE("RotationCoord Struct Operations", "[math][coord]")
{
	SECTION("Constructors and Assignment")
	{
		RotationCoord r1;
		REQUIRE(r1.azimuth == 0.0);

		RotationCoord r2(5.0);
		REQUIRE_THAT(r2.azimuth, WithinAbs(5.0, 1e-9));
		REQUIRE_THAT(r2.t, WithinAbs(5.0, 1e-9));

		RotationCoord r3(1.0, 2.0, 3.0);
		REQUIRE_THAT(r3.azimuth, WithinAbs(1.0, 1e-9));
		REQUIRE_THAT(r3.elevation, WithinAbs(2.0, 1e-9));
		REQUIRE_THAT(r3.t, WithinAbs(3.0, 1e-9));

		RotationCoord r4;
		r4 = 7.5;
		REQUIRE_THAT(r4.azimuth, WithinAbs(7.5, 1e-9));
		REQUIRE_THAT(r4.elevation, WithinAbs(7.5, 1e-9));
		REQUIRE_THAT(r4.t, WithinAbs(7.5, 1e-9));
	}

	SECTION("Arithmetic")
	{
		RotationCoord a(1.0, 2.0, 10.0);
		RotationCoord b(2.0, 3.0, 20.0);

		RotationCoord c = a + b;
		REQUIRE_THAT(c.azimuth, WithinAbs(3.0, 1e-9));
		REQUIRE_THAT(c.t, WithinAbs(10.0, 1e-9));

		c = b - a;
		REQUIRE_THAT(c.azimuth, WithinAbs(1.0, 1e-9));

		c = a * b;
		REQUIRE_THAT(c.azimuth, WithinAbs(2.0, 1e-9));

		c = b / a;
		REQUIRE_THAT(c.azimuth, WithinAbs(2.0, 1e-9));
	}

	SECTION("Scalar Arithmetic")
	{
		RotationCoord a(2.0, 4.0, 10.0);

		RotationCoord c = a * 2.0;
		REQUIRE_THAT(c.elevation, WithinAbs(8.0, 1e-9));

		c = a + 1.0;
		REQUIRE_THAT(c.azimuth, WithinAbs(3.0, 1e-9));

		c = 10.0 / a;
		REQUIRE_THAT(c.azimuth, WithinAbs(5.0, 1e-9));

		c = a / 2.0;
		REQUIRE_THAT(c.azimuth, WithinAbs(1.0, 1e-9));
		REQUIRE_THAT(c.elevation, WithinAbs(2.0, 1e-9));
		REQUIRE_THAT(c.t, WithinAbs(10.0, 1e-9)); // Time should be copied, not divided
	}
}
