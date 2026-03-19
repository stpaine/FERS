#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "timing/prototype_timing.h"

using Catch::Matchers::WithinAbs;

TEST_CASE("PrototypeTiming stores basic metadata", "[timing][prototype]")
{
	const SimId id = 12345;
	timing::PrototypeTiming proto("reference", id);

	REQUIRE(proto.getName() == "reference");
	REQUIRE(proto.getId() == id);
	REQUIRE_FALSE(proto.getSyncOnPulse());
	REQUIRE_THAT(proto.getFrequency(), WithinAbs(0.0, 1e-12));
	REQUIRE_FALSE(proto.getPhaseOffset().has_value());
	REQUIRE_FALSE(proto.getFreqOffset().has_value());
}

TEST_CASE("PrototypeTiming manages offsets and noise parameters", "[timing][prototype]")
{
	timing::PrototypeTiming proto("clock");

	proto.setFrequency(10.0);
	proto.setSyncOnPulse();
	proto.setAlpha(2.0, 1.0);
	proto.setAlpha(-2.0, 0.5);
	proto.setFreqOffset(0.25);
	proto.setPhaseOffset(-0.5);
	proto.setRandomFreqOffsetStdev(0.0);
	proto.setRandomPhaseOffsetStdev(0.0);

	std::vector<RealType> alphas;
	std::vector<RealType> weights;
	proto.copyAlphas(alphas, weights);

	REQUIRE_THAT(proto.getFrequency(), WithinAbs(10.0, 1e-12));
	REQUIRE(proto.getSyncOnPulse());
	REQUIRE(alphas.size() == 2);
	REQUIRE(weights.size() == 2);
	REQUIRE_THAT(alphas[0], WithinAbs(2.0, 1e-12));
	REQUIRE_THAT(alphas[1], WithinAbs(-2.0, 1e-12));
	REQUIRE_THAT(weights[0], WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(weights[1], WithinAbs(0.5, 1e-12));
	REQUIRE(proto.getFreqOffset().has_value());
	REQUIRE(proto.getPhaseOffset().has_value());
	REQUIRE(proto.getRandomFreqOffsetStdev().has_value());
	REQUIRE(proto.getRandomPhaseOffsetStdev().has_value());
	REQUIRE_THAT(proto.getFreqOffset().value(), WithinAbs(0.25, 1e-12));
	REQUIRE_THAT(proto.getPhaseOffset().value(), WithinAbs(-0.5, 1e-12));
}
