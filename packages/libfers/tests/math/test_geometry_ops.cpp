#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "math/geometry_ops.h"

using namespace math;
using Catch::Matchers::WithinAbs;

TEST_CASE("Vec3 operations", "[math][geometry]")
{
	SECTION("Initialization")
	{
		Vec3 v(1.0, 2.0, 3.0);
		REQUIRE_THAT(v.x, WithinAbs(1.0, 1e-9));
		REQUIRE_THAT(v.y, WithinAbs(2.0, 1e-9));
		REQUIRE_THAT(v.z, WithinAbs(3.0, 1e-9));
	}

	SECTION("Arithmetic")
	{
		Vec3 v1(1.0, 2.0, 3.0);
		Vec3 v2(4.0, 5.0, 6.0);

		v1 += v2;
		REQUIRE_THAT(v1.x, WithinAbs(5.0, 1e-9));
		REQUIRE_THAT(v1.y, WithinAbs(7.0, 1e-9));
		REQUIRE_THAT(v1.z, WithinAbs(9.0, 1e-9));

		Vec3 v3 = v1 - v2;
		REQUIRE_THAT(v3.x, WithinAbs(1.0, 1e-9));

		Vec3 v4 = v1 * 2.0;
		REQUIRE_THAT(v4.x, WithinAbs(10.0, 1e-9));
	}

	SECTION("Length and Dot Product")
	{
		Vec3 v(3.0, 4.0, 0.0);
		REQUIRE_THAT(v.length(), WithinAbs(5.0, 1e-9));

		Vec3 v2(1.0, 0.0, 0.0);
		REQUIRE_THAT(dotProduct(v, v2), WithinAbs(3.0, 1e-9));
	}
}

TEST_CASE("SVec3 operations", "[math][geometry]")
{
	SECTION("Conversion from Vec3")
	{
		Vec3 v(0.0, 1.0, 0.0); // Pointing along Y axis
		SVec3 sv(v);
		REQUIRE_THAT(sv.length, WithinAbs(1.0, 1e-9));
		REQUIRE_THAT(sv.azimuth, WithinAbs(PI / 2.0, 1e-9));
		REQUIRE_THAT(sv.elevation, WithinAbs(0.0, 1e-9));
	}

	SECTION("Conversion to Vec3")
	{
		SVec3 sv(2.0, PI / 2.0, 0.0);
		Vec3 v(sv);
		REQUIRE_THAT(v.x, WithinAbs(0.0, 1e-9));
		REQUIRE_THAT(v.y, WithinAbs(2.0, 1e-9));
		REQUIRE_THAT(v.z, WithinAbs(0.0, 1e-9));
	}

	SECTION("Subtraction with wrapping")
	{
		SVec3 sv1(1.0, PI - 0.1, 0.0);
		SVec3 sv2(1.0, -PI + 0.1, 0.0);
		SVec3 diff = sv1 - sv2;
		// Azimuth diff is 2*PI - 0.2, which wraps to -0.2
		REQUIRE_THAT(diff.azimuth, WithinAbs(-0.2, 1e-9));
	}
}
