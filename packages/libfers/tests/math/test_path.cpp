#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "math/path.h"
#include "math/path_utils.h"

using namespace math;
using Catch::Matchers::WithinAbs;

TEST_CASE("Path operations", "[math][path]")
{
	SECTION("Static Interpolation")
	{
		Path p(Path::InterpType::INTERP_STATIC);
		p.addCoord(Coord{Vec3(1, 2, 3), 0.0});
		p.addCoord(Coord{Vec3(4, 5, 6), 1.0});

		REQUIRE_THROWS_AS(p.getPosition(0.5), PathException);

		p.finalize();

		Vec3 pos = p.getPosition(0.5);
		REQUIRE_THAT(pos.x, WithinAbs(1.0, 1e-9));

		Vec3 vel = p.getVelocity(0.5);
		REQUIRE_THAT(vel.x, WithinAbs(0.0, 1e-9));
	}

	SECTION("Linear Interpolation")
	{
		Path p(Path::InterpType::INTERP_LINEAR);
		p.addCoord(Coord{Vec3(0, 0, 0), 0.0});
		p.addCoord(Coord{Vec3(10, 0, 0), 1.0});
		p.finalize();

		Vec3 pos = p.getPosition(0.5);
		REQUIRE_THAT(pos.x, WithinAbs(5.0, 1e-9));

		Vec3 vel = p.getVelocity(0.5);
		REQUIRE_THAT(vel.x, WithinAbs(10.0, 1e-9));
	}

	SECTION("Cubic Interpolation")
	{
		Path p(Path::InterpType::INTERP_CUBIC);
		p.addCoord(Coord{Vec3(0, 0, 0), 0.0});
		p.addCoord(Coord{Vec3(10, 10, 0), 1.0});
		p.addCoord(Coord{Vec3(20, 0, 0), 2.0});
		p.finalize();

		Vec3 pos = p.getPosition(0.5);
		REQUIRE(pos.x > 0.0);
		REQUIRE(pos.x < 10.0);

		Vec3 vel = p.getVelocity(0.5);
		REQUIRE(vel.x > 0.0);
	}
}
