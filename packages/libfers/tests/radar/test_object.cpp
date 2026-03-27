#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <memory>

#include "core/config.h"
#include "core/sim_id.h"
#include "math/path.h"
#include "math/rotation_path.h"
#include "radar/object.h"

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

TEST_CASE("Object exposes platform linkage and identity", "[radar][object]")
{
	radar::Platform platform("ObjectPlatform");
	radar::Object object(&platform, "TestObject", ObjectType::Transmitter, 9001);

	REQUIRE(object.getPlatform() == &platform);
	REQUIRE(object.getName() == "TestObject");
	REQUIRE(object.getId() == 9001);
}

TEST_CASE("Object allows updating its name", "[radar][object]")
{
	radar::Platform platform("ObjectPlatform");
	radar::Object object(&platform, "InitialName", ObjectType::Transmitter, 9001);

	REQUIRE(object.getName() == "InitialName");
	object.setName("UpdatedName");
	REQUIRE(object.getName() == "UpdatedName");
}

TEST_CASE("Object delegates position and rotation to platform", "[radar][object]")
{
	radar::Platform platform("MovingPlatform");
	platform.setMotionPath(makeLinearPath({0.0, 0.0, 0.0}, 0.0, {100.0, 0.0, 0.0}, 10.0));
	platform.setRotationPath(makeLinearRotation({0.0, 0.0, 0.0}, {1.0, 0.5, 4.0}));

	radar::Object object(&platform, "MovingObject", ObjectType::Target);

	const math::Vec3 position = object.getPosition(2.5);
	REQUIRE_THAT(position.x, WithinAbs(25.0, 1e-9));
	REQUIRE_THAT(position.y, WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(position.z, WithinAbs(0.0, 1e-12));

	const math::SVec3 rotation = object.getRotation(2.0);
	REQUIRE_THAT(rotation.azimuth, WithinAbs(0.5, 1e-9));
	REQUIRE_THAT(rotation.elevation, WithinAbs(0.25, 1e-9));

	REQUIRE(SimIdGenerator::getType(object.getId()) == ObjectType::Target);
}
