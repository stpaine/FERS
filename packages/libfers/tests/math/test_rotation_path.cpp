#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "math/path_utils.h"
#include "math/rotation_path.h"

using namespace math;
using Catch::Matchers::WithinAbs;

TEST_CASE("RotationPath operations", "[math][rotation_path]")
{
	SECTION("Static Interpolation")
	{
		RotationPath rp(RotationPath::InterpType::INTERP_STATIC);
		rp.addCoord(RotationCoord(0.0, 0.0, 0.0));
		rp.addCoord(RotationCoord(PI, PI / 2, 1.0));

		REQUIRE_THROWS_AS(rp.getPosition(0.5), PathException);

		rp.finalize();

		SVec3 pos = rp.getPosition(0.5);
		REQUIRE_THAT(pos.azimuth, WithinAbs(0.0, 1e-9));
	}

	SECTION("Linear Interpolation")
	{
		RotationPath rp(RotationPath::InterpType::INTERP_LINEAR);
		rp.addCoord(RotationCoord(0.0, 0.0, 0.0));
		rp.addCoord(RotationCoord(PI, PI / 2, 1.0));
		rp.finalize();

		SVec3 pos = rp.getPosition(0.5);
		REQUIRE_THAT(pos.azimuth, WithinAbs(PI / 2, 1e-9));
		REQUIRE_THAT(pos.elevation, WithinAbs(PI / 4, 1e-9));
	}

	SECTION("Constant Rate Interpolation")
	{
		RotationPath rp;
		rp.setConstantRate(RotationCoord(0.0, 0.0, 0.0), RotationCoord(PI, 0.0, 0.0));
		// finalize is called automatically by setConstantRate

		SVec3 pos = rp.getPosition(0.5);
		REQUIRE_THAT(pos.azimuth, WithinAbs(PI / 2, 1e-9));
	}
}
