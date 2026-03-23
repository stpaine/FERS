#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "math/path.h"
#include "math/path_utils.h"

using namespace math;
using Catch::Matchers::WithinAbs;

TEST_CASE("Path: Error Handling and Edge Cases", "[math][path]")
{
	Path p;

	REQUIRE(p.getType() == Path::InterpType::INTERP_STATIC);

	REQUIRE(p.getCoords().empty());

	SECTION("GetPosition before Finalize")
	{
		p.addCoord({Vec3(0, 0, 0), 0});
		REQUIRE_THROWS_AS(p.getPosition(0), PathException);
	}

	SECTION("GetVelocity before Finalize")
	{
		p.addCoord({Vec3(0, 0, 0), 0});
		REQUIRE_THROWS_AS(p.getVelocity(0), PathException);
	}

	SECTION("Velocity on Empty Path")
	{
		p.finalize();
		Vec3 v = p.getVelocity(0);
		REQUIRE(v.length() == 0.0);
	}

	SECTION("Cubic with insufficient points")
	{
		p.setInterp(Path::InterpType::INTERP_CUBIC);
		p.addCoord({Vec3(0, 0, 0), 0});
		REQUIRE_THROWS_AS(p.finalize(), PathException);
	}

	SECTION("Invalid Interpolation Type Fallback")
	{
		Path p_invalid(static_cast<Path::InterpType>(999));
		p_invalid.addCoord({Vec3(0, 0, 0), 0.0});
		p_invalid.finalize();
		Vec3 v = p_invalid.getVelocity(0.0);
		REQUIRE(v.length() == 0.0);
	}
}

TEST_CASE("Path: Static Interpolation", "[math][path]")
{
	Path p(Path::InterpType::INTERP_STATIC);
	p.addCoord({Vec3(1, 2, 3), 0.0});
	p.addCoord({Vec3(4, 5, 6), 10.0});
	p.finalize();

	Vec3 pos = p.getPosition(5.0);
	REQUIRE_THAT(pos.x, WithinAbs(1.0, 1e-9));

	Vec3 vel = p.getVelocity(5.0);
	REQUIRE(vel.length() == 0.0);
}

TEST_CASE("Path: Linear Interpolation Velocity Edge Cases", "[math][path]")
{
	Path p(Path::InterpType::INTERP_LINEAR);

	SECTION("Single Coordinate Velocity")
	{
		p.addCoord({Vec3(10, 0, 0), 1.0});
		p.finalize();
		Vec3 vel = p.getVelocity(1.0);
		REQUIRE(vel.length() == 0.0); // No delta T, velocity is 0
	}

	SECTION("Duplicate Timestamps (dt <= EPSILON)")
	{
		p.addCoord({Vec3(0, 0, 0), 1.0});
		p.addCoord({Vec3(10, 0, 0), 1.0}); // Physically impossible, but mathematically handled
		p.finalize();
		Vec3 vel = p.getVelocity(1.0);
		REQUIRE(vel.length() == 0.0); // Prevent divide by zero
	}

	SECTION("Out of Bounds Clamping")
	{
		p.addCoord({Vec3(0, 0, 0), 1.0});
		p.addCoord({Vec3(10, 0, 0), 2.0});
		p.finalize();

		Vec3 vel_before = p.getVelocity(0.0);
		REQUIRE_THAT(vel_before.x, WithinAbs(10.0, 1e-9)); // Continues initial velocity backwards

		Vec3 vel_after = p.getVelocity(3.0);
		REQUIRE_THAT(vel_after.x, WithinAbs(10.0, 1e-9)); // Continues final velocity forwards
	}
}

TEST_CASE("Path: Cubic Interpolation Velocity Edge Cases", "[math][path]")
{
	Path p(Path::InterpType::INTERP_CUBIC);

	SECTION("Out of Bounds Clamping")
	{
		p.addCoord({Vec3(0, 0, 0), 1.0});
		p.addCoord({Vec3(10, 0, 0), 2.0});
		p.addCoord({Vec3(20, 0, 0), 3.0});
		p.finalize();

		Vec3 vel_before = p.getVelocity(0.0);
		REQUIRE(vel_before.x > 0.0); // Velocity should be correctly extrapolated

		Vec3 vel_after = p.getVelocity(4.0);
		REQUIRE(vel_after.x > 0.0);
	}

	SECTION("Duplicate Timestamps (h <= EPSILON)")
	{
		p.addCoord({Vec3(0, 0, 0), 1.0});
		p.addCoord({Vec3(10, 0, 0), 1.0}); // Duplicate time
		p.addCoord({Vec3(20, 0, 0), 2.0});
		p.finalize();

		Vec3 vel = p.getVelocity(0.5);
		REQUIRE(vel.length() == 0.0); // Prevent divide by zero
	}
}

TEST_CASE("Path Utils: Direct Template Calls", "[math][path]")
{
	std::vector<Coord> empty_coords;
	Coord out;
	std::vector<Coord> dd;

	REQUIRE_THROWS_AS(getPositionLinear(0.0, out, empty_coords), PathException);

	REQUIRE_THROWS_AS(getPositionCubic(0.0, out, empty_coords, dd), PathException);

	std::vector<Coord> coords = {Coord{Vec3(0, 0, 0), 1.0}, Coord{Vec3(10, 0, 0), 2.0}, Coord{Vec3(20, 0, 0), 3.0}};
	finalizeCubic(coords, dd);

	getPositionCubic(0.0, out, coords, dd);
	REQUIRE_THAT(out.pos.x, WithinAbs(0.0, 1e-9)); // Clamped to first point physically

	getPositionCubic(4.0, out, coords, dd);
	REQUIRE_THAT(out.pos.x, WithinAbs(20.0, 1e-9)); // Clamped to last point physically
}
