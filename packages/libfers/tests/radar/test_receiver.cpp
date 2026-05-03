#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <memory>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/parameters.h"
#include "core/rendering_job.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/transmitter.h"
#include "serial/response.h"
#include "signal/radar_signal.h"
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

	std::unique_ptr<serial::Response> makeResponse(const radar::Transmitter* tx,
												   std::vector<std::unique_ptr<fers_signal::RadarSignal>>& wave_store)
	{
		auto signal = std::make_unique<fers_signal::CwSignal>();
		auto wave = std::make_unique<fers_signal::RadarSignal>("Wave", 1.0, 1.0e9, 1.0, std::move(signal));
		const auto* wave_ptr = wave.get();
		wave_store.push_back(std::move(wave));
		return std::make_unique<serial::Response>(wave_ptr, tx);
	}

	void requireComplexVectorsNear(const std::vector<ComplexType>& actual, const std::vector<ComplexType>& expected,
								   const RealType tolerance)
	{
		REQUIRE(actual.size() == expected.size());
		for (std::size_t i = 0; i < actual.size(); ++i)
		{
			REQUIRE_THAT(actual[i].real(), WithinAbs(expected[i].real(), tolerance));
			REQUIRE_THAT(actual[i].imag(), WithinAbs(expected[i].imag(), tolerance));
		}
	}
}

TEST_CASE("Receiver basic accessors and flags", "[radar][receiver]")
{
	radar::Platform platform("RxPlatform");
	radar::Receiver rx(&platform, "RxA", 42, radar::OperationMode::PULSED_MODE, 6001);

	REQUIRE(rx.getId() == 6001);
	REQUIRE(rx.getMode() == radar::OperationMode::PULSED_MODE);
	REQUIRE_FALSE(rx.isActive());

	rx.setActive(true);
	REQUIRE(rx.isActive());

	rx.setMode(radar::OperationMode::CW_MODE);
	REQUIRE(rx.getMode() == radar::OperationMode::CW_MODE);

	rx.setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
	REQUIRE(rx.checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT));
	REQUIRE_FALSE(rx.checkFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS));

	rx.clearFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
	REQUIRE_FALSE(rx.checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT));
}

TEST_CASE("Receiver noise temperature guards and aggregation", "[radar][receiver]")
{
	radar::Platform platform("RxPlatform");
	radar::Receiver rx(&platform, "RxA", 1, radar::OperationMode::CW_MODE);
	antenna::Isotropic antenna("RxAnt");
	rx.setAntenna(&antenna);

	rx.setNoiseTemperature(50.0);
	REQUIRE_THAT(rx.getNoiseTemperature(), WithinAbs(50.0, 1e-12));

	const math::SVec3 angle(1.0, 0.1, 0.2);
	REQUIRE_THAT(rx.getNoiseTemperature(angle), WithinAbs(50.0, 1e-12));

	REQUIRE_THROWS_AS(rx.setNoiseTemperature(-1.0), std::runtime_error);
}

TEST_CASE("Receiver exposes RNG and mutable streaming data", "[radar][receiver]")
{
	radar::Platform platform("RxPlatform");
	radar::Receiver rx(&platform, "RxA", 99, radar::OperationMode::CW_MODE);

	const auto first_sample = rx.getRngEngine()();
	REQUIRE(first_sample != 0U);

	rx.prepareStreamingData(2);
	auto& mutable_data = rx.getMutableStreamingData();
	mutable_data[0] = ComplexType(2.0, -1.0);
	REQUIRE_THAT(rx.getStreamingData()[0].real(), WithinAbs(2.0, 1e-12));
	REQUIRE_THAT(rx.getStreamingData()[0].imag(), WithinAbs(-1.0, 1e-12));
}

TEST_CASE("Receiver windows quantize to sampling rate", "[radar][receiver]")
{
	ParamGuard guard;
	params::setTime(0.0, 10.0);
	params::setRate(1000.0);
	params::setOversampleRatio(1);

	radar::Platform platform("RxPlatform");
	radar::Receiver rx(&platform, "RxA", 7, radar::OperationMode::PULSED_MODE);

	rx.setWindowProperties(2.0, 333.0, 0.0015);

	const RealType expected_prf = 1.0 / (std::floor(1000.0 / 333.0) / 1000.0);
	const RealType expected_skip = std::floor(1000.0 * 0.0015) / 1000.0;

	REQUIRE_THAT(rx.getWindowPrf(), WithinAbs(expected_prf, 1e-12));
	REQUIRE_THAT(rx.getWindowSkip(), WithinAbs(expected_skip, 1e-12));
	REQUIRE_THAT(rx.getWindowLength(), WithinAbs(2.0, 1e-12));

	const unsigned expected_count = static_cast<unsigned>(std::ceil((10.0 - 0.0) * expected_prf));
	REQUIRE(rx.getWindowCount() == expected_count);

	REQUIRE_THROWS_AS(rx.getWindowStart(0), std::logic_error);

	auto timing = std::make_shared<timing::Timing>("RxClock", 99);
	rx.setTiming(timing);

	const RealType expected_start = 0.0 / expected_prf + expected_skip;
	REQUIRE_THAT(rx.getWindowStart(0), WithinAbs(expected_start, 1e-12));
}

TEST_CASE("Receiver inbox and interference log", "[radar][receiver]")
{
	radar::Platform platform("RxPlatform");
	radar::Receiver rx(&platform, "RxA", 5, radar::OperationMode::CW_MODE);
	radar::Transmitter tx(&platform, "TxA", radar::OperationMode::CW_MODE, 9001);

	std::vector<std::unique_ptr<fers_signal::RadarSignal>> waves;
	rx.addResponseToInbox(makeResponse(&tx, waves));
	rx.addResponseToInbox(makeResponse(&tx, waves));

	const auto drained = rx.drainInbox();
	REQUIRE(drained.size() == 2);
	REQUIRE(rx.drainInbox().empty());

	rx.addInterferenceToLog(makeResponse(&tx, waves));
	REQUIRE(rx.getPulsedInterferenceLog().size() == 1);
}

TEST_CASE("Receiver streaming data accumulation", "[radar][receiver]")
{
	radar::Platform platform("RxPlatform");
	radar::Receiver rx(&platform, "RxA", 9, radar::OperationMode::CW_MODE);

	rx.prepareStreamingData(3);
	REQUIRE(rx.getStreamingData().size() == 3);

	rx.setStreamingSample(1, ComplexType(1.0, 2.0));
	rx.setStreamingSample(1, ComplexType(0.5, -1.0));

	REQUIRE_THAT(rx.getStreamingData()[1].real(), WithinAbs(1.5, 1e-12));
	REQUIRE_THAT(rx.getStreamingData()[1].imag(), WithinAbs(1.0, 1e-12));

	rx.setStreamingSample(10, ComplexType(3.0, 3.0));
	REQUIRE_THAT(rx.getStreamingData()[1].real(), WithinAbs(1.5, 1e-12));
}

TEST_CASE("Receiver IF resampling preserves absolute output positions across schedule gaps", "[radar][receiver][fmcw]")
{
	ParamGuard guard;
	params::setTime(0.0, 1.0);
	params::setRate(10.0);
	params::setOversampleRatio(1);

	radar::Platform platform("RxPlatform");
	radar::Receiver rx(&platform, "IfRx", 9, radar::OperationMode::FMCW_MODE);
	rx.setDechirpMode(radar::Receiver::DechirpMode::Physical);

	const fers_signal::FmcwIfResamplerRequest request{.input_sample_rate_hz = 10.0,
													  .output_sample_rate_hz = 10.0,
													  .filter_bandwidth_hz = 2.0,
													  .filter_transition_width_hz = 1.0};
	rx.initializeFmcwIfResampling(fers_signal::planFmcwIfResampler(request));

	const std::array first_segment{ComplexType{1.0, 0.0}, ComplexType{2.0, 0.0}};
	const std::array second_segment{ComplexType{3.0, 0.0}};
	rx.beginFmcwIfResamplingSegment(0.2);
	rx.consumeFmcwIfBlock(first_segment, 0.2);
	rx.endFmcwIfResamplingSegment();
	rx.beginFmcwIfResamplingSegment(0.6);
	rx.consumeFmcwIfBlock(second_segment, 0.6);
	rx.endFmcwIfResamplingSegment();
	rx.flushFmcwIfResampling();

	const auto& samples = rx.getStreamingData();
	REQUIRE(samples.size() == 10);
	REQUIRE_THAT(samples[0].real(), WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(samples[1].real(), WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(samples[2].real(), WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(samples[3].real(), WithinAbs(2.0, 1e-12));
	REQUIRE_THAT(samples[4].real(), WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(samples[5].real(), WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(samples[6].real(), WithinAbs(3.0, 1e-12));
	REQUIRE_THAT(samples[7].real(), WithinAbs(0.0, 1e-12));
}

TEST_CASE("Receiver IF resampling aligns off-grid segments to the global IF grid", "[radar][receiver][fmcw]")
{
	ParamGuard guard;
	params::setTime(0.0, 1.0);
	params::setRate(100.0);
	params::setOversampleRatio(1);

	const fers_signal::FmcwIfResamplerRequest request{.input_sample_rate_hz = 100.0,
													  .output_sample_rate_hz = 20.0,
													  .filter_bandwidth_hz = 4.0,
													  .filter_transition_width_hz = 5.0,
													  .stopband_attenuation_db = 30.0};
	const auto plan = fers_signal::planFmcwIfResampler(request);

	std::vector<ComplexType> full_timeline(100, ComplexType{0.0, 0.0});
	for (std::size_t i = 13; i < 55; ++i)
	{
		full_timeline[i] = ComplexType{1.0 + static_cast<RealType>(i) * 0.01, -0.5 + static_cast<RealType>(i) * 0.005};
	}

	radar::Platform platform("RxPlatform");
	radar::Receiver reference(&platform, "ReferenceIfRx", 10, radar::OperationMode::FMCW_MODE);
	reference.setDechirpMode(radar::Receiver::DechirpMode::Physical);
	reference.initializeFmcwIfResampling(plan);
	reference.consumeFmcwIfBlock(full_timeline, 0.0);
	reference.flushFmcwIfResampling();

	radar::Receiver scheduled(&platform, "ScheduledIfRx", 11, radar::OperationMode::FMCW_MODE);
	scheduled.setDechirpMode(radar::Receiver::DechirpMode::Physical);
	scheduled.initializeFmcwIfResampling(plan);
	scheduled.beginFmcwIfResamplingSegment(0.13);
	scheduled.consumeFmcwIfBlock(std::span<const ComplexType>(full_timeline.data() + 13, 42), 0.13);
	scheduled.endFmcwIfResamplingSegment();
	scheduled.flushFmcwIfResampling();

	requireComplexVectorsNear(scheduled.getStreamingData(), reference.getStreamingData(), 1e-10);
}

TEST_CASE("Receiver IF resampling flushes filtered tail samples after scheduled segment ends",
		  "[radar][receiver][fmcw]")
{
	ParamGuard guard;
	params::setTime(0.0, 1.0);
	params::setRate(100.0);
	params::setOversampleRatio(1);

	const fers_signal::FmcwIfResamplerRequest request{.input_sample_rate_hz = 100.0,
													  .output_sample_rate_hz = 20.0,
													  .filter_bandwidth_hz = 4.0,
													  .filter_transition_width_hz = 5.0,
													  .stopband_attenuation_db = 30.0};
	const auto plan = fers_signal::planFmcwIfResampler(request);

	std::vector<ComplexType> full_timeline(100, ComplexType{0.0, 0.0});
	for (std::size_t i = 0; i < 50; ++i)
	{
		full_timeline[i] =
			ComplexType{std::sin(static_cast<RealType>(i) * 0.13), std::cos(static_cast<RealType>(i) * 0.17)};
	}

	radar::Platform platform("RxPlatform");
	radar::Receiver reference(&platform, "ReferenceTailIfRx", 12, radar::OperationMode::FMCW_MODE);
	reference.setDechirpMode(radar::Receiver::DechirpMode::Physical);
	reference.initializeFmcwIfResampling(plan);
	reference.consumeFmcwIfBlock(full_timeline, 0.0);
	reference.flushFmcwIfResampling();

	radar::Receiver scheduled(&platform, "ScheduledTailIfRx", 13, radar::OperationMode::FMCW_MODE);
	scheduled.setDechirpMode(radar::Receiver::DechirpMode::Physical);
	scheduled.initializeFmcwIfResampling(plan);
	scheduled.beginFmcwIfResamplingSegment(0.0);
	scheduled.consumeFmcwIfBlock(std::span<const ComplexType>(full_timeline.data(), 50), 0.0);
	scheduled.endFmcwIfResamplingSegment();
	scheduled.flushFmcwIfResampling();

	requireComplexVectorsNear(scheduled.getStreamingData(), reference.getStreamingData(), 1e-10);
}

TEST_CASE("Receiver IF resampling skips long inactive gaps without zero-sample work", "[radar][receiver][fmcw]")
{
	ParamGuard guard;
	params::setTime(0.0, 1000.0);
	params::setRate(100000.0);
	params::setOversampleRatio(1);

	const fers_signal::FmcwIfResamplerRequest request{.input_sample_rate_hz = 100000.0,
													  .output_sample_rate_hz = 100.0,
													  .filter_bandwidth_hz = 20.0,
													  .filter_transition_width_hz = 25.0,
													  .stopband_attenuation_db = 30.0};
	const auto plan = fers_signal::planFmcwIfResampler(request);

	radar::Platform platform("RxPlatform");
	radar::Receiver rx(&platform, "SparseIfRx", 14, radar::OperationMode::FMCW_MODE);
	rx.setDechirpMode(radar::Receiver::DechirpMode::Physical);
	rx.initializeFmcwIfResampling(plan);

	std::vector<ComplexType> active_block(10000, ComplexType{1.0, 0.0});
	rx.beginFmcwIfResamplingSegment(500.0);
	rx.consumeFmcwIfBlock(active_block, 500.0);
	rx.endFmcwIfResamplingSegment();
	rx.flushFmcwIfResampling();

	const auto& samples = rx.getStreamingData();
	REQUIRE(samples.size() == 100000);
	REQUIRE_THAT(std::abs(samples[49900]), WithinAbs(0.0, 1e-12));

	RealType local_peak = 0.0;
	for (std::size_t i = 49990; i < 50030; ++i)
	{
		local_peak = std::max(local_peak, std::abs(samples[i]));
	}
	REQUIRE(local_peak > 0.01);
}

TEST_CASE("Receiver schedule determines next window time", "[radar][receiver]")
{
	radar::Platform platform("RxPlatform");
	radar::Receiver rx(&platform, "RxA", 11, radar::OperationMode::PULSED_MODE);

	SECTION("No schedule means always on")
	{
		const auto next = rx.getNextWindowTime(2.0);
		REQUIRE(next.has_value());
		REQUIRE_THAT(*next, WithinAbs(2.0, 1e-12));
	}

	SECTION("Schedule enforces active windows")
	{
		std::vector<radar::SchedulePeriod> schedule = {{1.0, 2.0}, {4.0, 5.0}};
		rx.setSchedule(schedule);
		REQUIRE(rx.getSchedule().size() == 2);

		const auto inside = rx.getNextWindowTime(1.25);
		REQUIRE(inside.has_value());
		REQUIRE_THAT(*inside, WithinAbs(1.25, 1e-12));

		const auto before = rx.getNextWindowTime(3.5);
		REQUIRE(before.has_value());
		REQUIRE_THAT(*before, WithinAbs(4.0, 1e-12));

		const auto after = rx.getNextWindowTime(6.0);
		REQUIRE_FALSE(after.has_value());
	}
}

TEST_CASE("Receiver finalizer queue processes jobs", "[radar][receiver]")
{
	radar::Platform platform("RxPlatform");
	radar::Receiver rx(&platform, "RxA", 13, radar::OperationMode::PULSED_MODE);

	core::RenderingJob job_in{};
	job_in.ideal_start_time = 1.0;
	job_in.duration = 0.5;

	rx.enqueueFinalizerJob(std::move(job_in));

	core::RenderingJob job_out{};
	REQUIRE(rx.waitAndDequeueFinalizerJob(job_out));
	REQUIRE_THAT(job_out.ideal_start_time, WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(job_out.duration, WithinAbs(0.5, 1e-12));

	core::RenderingJob shutdown_job{};
	shutdown_job.duration = -1.0;
	rx.enqueueFinalizerJob(std::move(shutdown_job));
	REQUIRE_FALSE(rx.waitAndDequeueFinalizerJob(job_out));
}
