#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <memory>
#include <string>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "radar/platform.h"
#include "radar/radar_obj.h"
#include "timing/timing.h"

using Catch::Matchers::WithinAbs;

TEST_CASE("Radar enforces timing and antenna invariants", "[radar][radar_obj]")
{
	radar::Platform platform("RadarPlatform");
	radar::Radar radar(&platform, "RadarA");

	REQUIRE_THROWS_AS(radar.getTiming(), std::runtime_error);
	REQUIRE_THROWS_AS(radar.setTiming(nullptr), std::runtime_error);
	REQUIRE_THROWS_AS(radar.setAntenna(nullptr), std::logic_error);
}

TEST_CASE("Radar stores timing and antenna references", "[radar][radar_obj]")
{
	radar::Platform platform("RadarPlatform");
	radar::Radar radar(&platform, "RadarA", 7777);

	auto timing = std::make_shared<timing::Timing>("ClockA", 1234);
	radar.setTiming(timing);

	antenna::Isotropic antenna("IsoA");
	antenna.setEfficiencyFactor(0.85);
	radar.setAntenna(&antenna);

	REQUIRE(radar.getTiming() == timing);
	REQUIRE(radar.getAntenna() == &antenna);
	REQUIRE(radar.getId() == 7777);
}

TEST_CASE("Radar delegates gain and noise calculations", "[radar][radar_obj]")
{
	radar::Platform platform("RadarPlatform");
	radar::Radar radar(&platform, "RadarA");

	antenna::Isotropic antenna("IsoA");
	antenna.setEfficiencyFactor(0.7);
	radar.setAntenna(&antenna);

	const math::SVec3 angle(1.0, 0.1, 0.2);
	const math::SVec3 refangle(1.0, 0.1, 0.2);

	const RealType gain = radar.getGain(angle, refangle, 0.03);
	REQUIRE_THAT(gain, WithinAbs(0.7, 1e-12));

	const RealType noise_temp = radar.getNoiseTemperature(angle);
	REQUIRE_THAT(noise_temp, WithinAbs(0.0, 1e-12));
}

TEST_CASE("Radar attachment enforces single attached object", "[radar][radar_obj]")
{
	radar::Platform platform("RadarPlatform");
	radar::Radar radar_a(&platform, "RadarA");
	radar::Radar radar_b(&platform, "RadarB");
	radar::Radar radar_c(&platform, "RadarC");

	REQUIRE(radar_a.getAttached() == nullptr);

	radar_a.setAttached(&radar_b);
	REQUIRE(radar_a.getAttached() == &radar_b);

	REQUIRE_THROWS_AS(radar_a.setAttached(&radar_c), std::runtime_error);
}
