#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <random>

#include "core/config.h"
#include "core/parameters.h"
#include "timing/prototype_timing.h"
#include "timing/timing.h"

using Catch::Matchers::WithinAbs;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
	};

	RealType expectedSample(const RealType phase_offset, const RealType freq_offset, const unsigned long count)
	{
		return phase_offset + 2.0 * PI * freq_offset * static_cast<RealType>(count) / params::rate();
	}
}

TEST_CASE("Timing defaults to disabled", "[timing]")
{
	timing::Timing timing_source("clk", 1234, 9);

	REQUIRE(timing_source.getName() == "clk");
	REQUIRE(timing_source.getId() == 9);
	REQUIRE_FALSE(timing_source.isEnabled());
	REQUIRE_THAT(timing_source.getNextSample(), WithinAbs(0.0, 1e-12));

	timing_source.skipSamples(10);
	timing_source.reset();
}

TEST_CASE("Timing initializes model and produces deterministic offsets", "[timing]")
{
	ParamGuard guard;
	params::setRate(1000.0);

	timing::PrototypeTiming proto("clock", 77);
	proto.setFrequency(10.0);
	proto.setFreqOffset(2.0);
	proto.setPhaseOffset(0.1);
	proto.setSyncOnPulse();

	timing::Timing timing_source("clock", 42, 77);
	timing_source.initializeModel(&proto);

	REQUIRE(timing_source.isEnabled());
	REQUIRE(timing_source.getSyncOnPulse());
	REQUIRE_THAT(timing_source.getFrequency(), WithinAbs(10.0, 1e-12));
	REQUIRE_THAT(timing_source.getFreqOffset(), WithinAbs(2.0, 1e-12));
	REQUIRE_THAT(timing_source.getPhaseOffset(), WithinAbs(0.1, 1e-12));

	const auto sample0 = timing_source.getNextSample();
	REQUIRE_THAT(sample0, WithinAbs(expectedSample(0.1, 2.0, 0), 1e-12));

	timing_source.skipSamples(9);
	const auto sample10 = timing_source.getNextSample();
	REQUIRE_THAT(sample10, WithinAbs(expectedSample(0.1, 2.0, 10), 1e-12));
}

TEST_CASE("Timing applies random frequency and phase offsets", "[timing]")
{
	ParamGuard guard;
	params::setRate(1000.0);

	const unsigned seed = 4242;
	const RealType base_freq_offset = 1.0;
	const RealType base_phase_offset = -0.2;
	const RealType freq_stdev = 0.5;
	const RealType phase_stdev = 0.25;

	std::mt19937 reference_rng(seed);
	std::normal_distribution<RealType> normal_dist(0.0, 1.0);
	const RealType expected_freq = base_freq_offset + normal_dist(reference_rng) * freq_stdev;
	const RealType expected_phase = base_phase_offset + normal_dist(reference_rng) * phase_stdev;

	timing::PrototypeTiming proto("random", 101);
	proto.setFrequency(5.0);
	proto.setFreqOffset(base_freq_offset);
	proto.setPhaseOffset(base_phase_offset);
	proto.setRandomFreqOffsetStdev(freq_stdev);
	proto.setRandomPhaseOffsetStdev(phase_stdev);

	timing::Timing timing_source("random", seed, 101);
	timing_source.initializeModel(&proto);

	REQUIRE(timing_source.isEnabled());
	REQUIRE_THAT(timing_source.getFreqOffset(), WithinAbs(expected_freq, 1e-12));
	REQUIRE_THAT(timing_source.getPhaseOffset(), WithinAbs(expected_phase, 1e-12));
	REQUIRE_THAT(timing_source.getNextSample(), WithinAbs(expectedSample(expected_phase, expected_freq, 0), 1e-12));
}

TEST_CASE("Timing handles zero frequency initialization", "[timing]")
{
	ParamGuard guard;
	params::setRate(1000.0);

	timing::PrototypeTiming proto("zero");
	proto.setFrequency(0.0);
	proto.setPhaseOffset(0.125);

	timing::Timing timing_source("zero", 7, 12);
	timing_source.initializeModel(&proto);

	REQUIRE_THAT(timing_source.getFrequency(), WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(timing_source.getPhaseOffset(), WithinAbs(0.125, 1e-12));
	REQUIRE(timing_source.isEnabled());
	REQUIRE_THAT(timing_source.getNextSample(), WithinAbs(expectedSample(0.125, 0.0, 0), 1e-12));
}

TEST_CASE("Timing reports disabled when no offsets or noise", "[timing]")
{
	ParamGuard guard;
	params::setRate(1000.0);

	timing::PrototypeTiming proto("quiet");
	proto.setFrequency(1.0);

	timing::Timing timing_source("quiet", 7, 88);
	timing_source.initializeModel(&proto);

	REQUIRE_FALSE(timing_source.isEnabled());
	REQUIRE_THAT(timing_source.getNextSample(), WithinAbs(0.0, 1e-12));
}

TEST_CASE("Timing clone requires initialized prototype", "[timing]")
{
	timing::Timing timing_source("clock", 11, 5);

	REQUIRE_THROWS_AS(timing_source.clone(), std::logic_error);
}

TEST_CASE("Timing clone reproduces prototype behavior", "[timing]")
{
	ParamGuard guard;
	params::setRate(1000.0);

	timing::PrototypeTiming proto("cloneable", 19);
	proto.setFrequency(10.0);
	proto.setFreqOffset(1.0);
	proto.setPhaseOffset(-0.25);

	timing::Timing timing_source("cloneable", 2024, 19);
	timing_source.initializeModel(&proto);

	auto clone = timing_source.clone();
	REQUIRE(clone->getId() == timing_source.getId());
	REQUIRE(clone->getName() == timing_source.getName());

	const auto sample0 = timing_source.getNextSample();
	const auto clone_sample0 = clone->getNextSample();
	REQUIRE_THAT(sample0, WithinAbs(clone_sample0, 1e-12));

	const auto sample1 = timing_source.getNextSample();
	const auto clone_sample1 = clone->getNextSample();
	REQUIRE_THAT(sample1, WithinAbs(clone_sample1, 1e-12));
}

TEST_CASE("Timing initializeModel is idempotent", "[timing]")
{
	ParamGuard guard;
	params::setRate(1000.0);

	timing::PrototypeTiming proto_a("clock", 3);
	proto_a.setFrequency(8.0);
	proto_a.setFreqOffset(1.0);
	proto_a.setPhaseOffset(0.05);

	timing::PrototypeTiming proto_b("clock", 4);
	proto_b.setFrequency(15.0);
	proto_b.setFreqOffset(7.0);
	proto_b.setPhaseOffset(0.9);

	timing::Timing timing_source("clock", 55, 3);
	timing_source.initializeModel(&proto_a);
	const auto sample0 = timing_source.getNextSample();

	timing_source.initializeModel(&proto_b);
	const auto sample1 = timing_source.getNextSample();

	REQUIRE_THAT(sample0, WithinAbs(expectedSample(0.05, 1.0, 0), 1e-12));
	REQUIRE_THAT(sample1, WithinAbs(expectedSample(0.05, 1.0, 1), 1e-12));
}
