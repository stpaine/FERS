#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <memory>

#include "core/config.h"
#include "math/path.h"
#include "math/rotation_path.h"
#include "radar/platform.h"

using Catch::Matchers::WithinAbs;

namespace
{
	std::unique_ptr<math::Path> makeLinearPath(const math::Vec3& start_pos, RealType start_time,
											   const math::Vec3& end_pos, RealType end_time)
	{
		auto path = std::make_unique<math::Path>(math::Path::InterpType::INTERP_LINEAR);
		path->addCoord({start_pos, start_time});
		path->addCoord({end_pos, end_time});
		path->finalize();
		return path;
	}

	std::unique_ptr<math::RotationPath> makeLinearRotation(const math::RotationCoord& start,
														   const math::RotationCoord& end)
	{
		auto path = std::make_unique<math::RotationPath>(math::RotationPath::InterpType::INTERP_LINEAR);
		path->addCoord(start);
		path->addCoord(end);
		path->finalize();
		return path;
	}
}

TEST_CASE("Platform exposes default paths and identity", "[radar][platform]")
{
	radar::Platform platform("TestPlatform", 4242);

	REQUIRE(platform.getMotionPath() != nullptr);
	REQUIRE(platform.getRotationPath() != nullptr);
	REQUIRE(platform.getMotionPath()->getType() == math::Path::InterpType::INTERP_STATIC);
	REQUIRE(platform.getRotationPath()->getType() == math::RotationPath::InterpType::INTERP_STATIC);
	REQUIRE(platform.getName() == "TestPlatform");
	REQUIRE(platform.getId() == 4242);
}

TEST_CASE("Platform returns physically correct position and rotation", "[radar][platform]")
{
	radar::Platform platform("KinematicPlatform");

	platform.setMotionPath(makeLinearPath({0.0, 0.0, 0.0}, 0.0, {10.0, 0.0, 0.0}, 10.0));
	platform.setRotationPath(makeLinearRotation({0.0, 0.0, 0.0}, {0.5, 0.25, 2.0}));

	const math::Vec3 position = platform.getPosition(5.0);
	REQUIRE_THAT(position.x, WithinAbs(5.0, 1e-9));
	REQUIRE_THAT(position.y, WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(position.z, WithinAbs(0.0, 1e-12));

	const math::SVec3 rotation = platform.getRotation(1.0);
	REQUIRE_THAT(rotation.azimuth, WithinAbs(0.25, 1e-9));
	REQUIRE_THAT(rotation.elevation, WithinAbs(0.125, 1e-9));
}

TEST_CASE("Platform swaps motion paths correctly", "[radar][platform]")
{
	radar::Platform platform("PathSwap");
	platform.setMotionPath(makeLinearPath({0.0, 0.0, 0.0}, 0.0, {4.0, 0.0, 0.0}, 4.0));
	REQUIRE_THAT(platform.getPosition(2.0).x, WithinAbs(2.0, 1e-9));

	platform.setMotionPath(makeLinearPath({0.0, 0.0, 0.0}, 0.0, {6.0, 0.0, 0.0}, 6.0));
	REQUIRE_THAT(platform.getPosition(2.0).x, WithinAbs(2.0, 1e-9));
	REQUIRE_THAT(platform.getPosition(3.0).x, WithinAbs(3.0, 1e-9));
}
