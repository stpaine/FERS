#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <memory>
#include <vector>

#include "core/config.h"
#include "core/parameters.h"
#include "processing/finalizer_pipeline.h"
#include "radar/platform.h"
#include "radar/receiver.h"
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

	RealType expectedSample(RealType phase_offset, RealType freq_offset, unsigned long count)
	{
		return phase_offset + 2.0 * PI * freq_offset * static_cast<RealType>(count) / params::rate();
	}
}

TEST_CASE("advanceTimingModel safely ignores null and disabled timing models", "[processing][finalizer][timing]")
{
	processing::pipeline::advanceTimingModel(nullptr, nullptr, 1000.0);

	radar::Platform platform("RxPlatform");
	radar::Receiver receiver(&platform, "RxA", 1, radar::OperationMode::PULSED_MODE);
	auto timing_model = std::make_unique<timing::Timing>("clk", 123);

	processing::pipeline::advanceTimingModel(timing_model.get(), &receiver, 1000.0);

	REQUIRE_THAT(timing_model->getNextSample(), WithinAbs(0.0, 0.0));
}

TEST_CASE("advanceTimingModel resets sync-on-pulse timing before skipping receiver dead time",
		  "[processing][finalizer][timing]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);

	radar::Platform platform("RxPlatform");
	radar::Receiver receiver(&platform, "RxA", 2, radar::OperationMode::PULSED_MODE);
	receiver.setWindowProperties(0.002, 100.0, 0.003);

	timing::PrototypeTiming prototype("clk");
	prototype.setFrequency(10.0);
	prototype.setFreqOffset(2.0);
	prototype.setPhaseOffset(0.1);
	prototype.setSyncOnPulse();

	auto timing_model = std::make_unique<timing::Timing>("clk", 77);
	timing_model->initializeModel(&prototype);

	(void)timing_model->getNextSample();
	(void)timing_model->getNextSample();

	processing::pipeline::advanceTimingModel(timing_model.get(), &receiver, params::rate());

	REQUIRE_THAT(timing_model->getNextSample(), WithinAbs(expectedSample(0.1, 2.0, 3), 1e-12));
}

TEST_CASE("advanceTimingModel free-running mode skips only the inter-pulse interval", "[processing][finalizer][timing]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);

	radar::Platform platform("RxPlatform");
	radar::Receiver receiver(&platform, "RxA", 3, radar::OperationMode::PULSED_MODE);
	receiver.setWindowProperties(0.002, 100.0, 0.003);

	timing::PrototypeTiming prototype("clk");
	prototype.setFrequency(10.0);
	prototype.setFreqOffset(2.0);
	prototype.setPhaseOffset(0.1);

	auto timing_model = std::make_unique<timing::Timing>("clk", 88);
	timing_model->initializeModel(&prototype);

	(void)timing_model->getNextSample();
	(void)timing_model->getNextSample();

	processing::pipeline::advanceTimingModel(timing_model.get(), &receiver, params::rate());

	REQUIRE_THAT(timing_model->getNextSample(), WithinAbs(expectedSample(0.1, 2.0, 10), 1e-12));
}

TEST_CASE("calculateJitteredStart converts phase noise into rounded start time and fractional delay",
		  "[processing][finalizer][timing]")
{
	SECTION("positive phase noise delays the sample clock")
	{
		const auto [rounded_start, fractional_delay] = processing::pipeline::calculateJitteredStart(1.0, PI, 10.0, 8.0);

		REQUIRE_THAT(rounded_start, WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(fractional_delay, WithinAbs(0.4, 1e-12));
	}

	SECTION("negative phase noise advances the sample clock")
	{
		const auto [rounded_start, fractional_delay] =
			processing::pipeline::calculateJitteredStart(1.0, -PI, 10.0, 8.0);

		REQUIRE_THAT(rounded_start, WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(fractional_delay, WithinAbs(-0.4, 1e-12));
	}
}

TEST_CASE("addPhaseNoiseToWindow rotates each covered sample by the expected complex phase",
		  "[processing][finalizer][timing]")
{
	const std::vector<RealType> noise = {PI / 2.0, -PI / 2.0};
	std::vector<ComplexType> window = {
		ComplexType{1.0, 0.0},
		ComplexType{0.0, 1.0},
		ComplexType{3.0, -4.0},
	};

	processing::pipeline::addPhaseNoiseToWindow(noise, window);

	REQUIRE_THAT(window[0].real(), WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(window[0].imag(), WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(window[1].real(), WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(window[1].imag(), WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(window[2].real(), WithinAbs(3.0, 1e-12));
	REQUIRE_THAT(window[2].imag(), WithinAbs(-4.0, 1e-12));
}
