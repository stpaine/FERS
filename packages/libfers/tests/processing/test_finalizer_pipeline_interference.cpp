#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/parameters.h"
#include "math/coord.h"
#include "processing/finalizer_pipeline.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/response.h"
#include "signal/radar_signal.h"
#include "simulation/channel_model.h"
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

	struct FixedSignal final : public fers_signal::Signal
	{
		std::vector<ComplexType> data;

		std::vector<ComplexType> render(const std::vector<interp::InterpPoint>&, unsigned& size,
										RealType) const override
		{
			size = static_cast<unsigned>(data.size());
			return data;
		}
	};

	void setupPlatform(radar::Platform& platform, const math::Vec3& position)
	{
		platform.getMotionPath()->addCoord(math::Coord{position, 0.0});
		platform.getMotionPath()->finalize();
		platform.getRotationPath()->addCoord(math::RotationCoord{0.0, 0.0, 0.0});
		platform.getRotationPath()->finalize();
	}

	std::shared_ptr<timing::Timing> makeQuietTiming(const std::string& name, unsigned seed)
	{
		timing::PrototypeTiming prototype(name);
		prototype.setFrequency(1.0);

		auto timing_model = std::make_shared<timing::Timing>(name, seed);
		timing_model->initializeModel(&prototype);
		return timing_model;
	}

	std::unique_ptr<serial::Response>
	makeFixedResponse(const radar::Transmitter* transmitter,
					  std::vector<std::unique_ptr<fers_signal::RadarSignal>>& wave_store,
					  const std::vector<ComplexType>& samples, RealType sample_rate, RealType start_time)
	{
		auto signal = std::make_unique<FixedSignal>();
		signal->data = samples;
		signal->load(samples, static_cast<unsigned>(samples.size()), sample_rate);

		auto wave = std::make_unique<fers_signal::RadarSignal>(
			"wave", 1.0, 1.0e9, static_cast<RealType>(samples.size()) / sample_rate, std::move(signal));
		const auto* wave_ptr = wave.get();
		wave_store.push_back(std::move(wave));

		auto response = std::make_unique<serial::Response>(wave_ptr, transmitter);
		response->addInterpPoint({1.0, start_time, 0.0, 0.0});
		response->addInterpPoint({1.0, start_time + static_cast<RealType>(samples.size() - 1) / sample_rate, 0.0, 0.0});
		return response;
	}
}

TEST_CASE("applyStreamingInterference adds direct-path streaming energy sample by sample",
		  "[processing][finalizer][interference]")
{
	ParamGuard guard;
	params::params.reset();

	radar::Platform tx_platform("TxPlatform");
	setupPlatform(tx_platform, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_platform("RxPlatform");
	setupPlatform(rx_platform, math::Vec3{1000.0, 0.0, 0.0});

	antenna::Isotropic antenna("iso");
	auto timing_model = makeQuietTiming("clk", 11);

	radar::Transmitter transmitter(&tx_platform, "TxA", radar::OperationMode::CW_MODE, 101);
	transmitter.setAntenna(&antenna);
	transmitter.setTiming(timing_model);
	auto signal = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("cw", 25.0, 1.0e9, 1.0, std::move(signal), 301);
	transmitter.setSignal(&wave);

	radar::Receiver receiver(&rx_platform, "RxA", 99, radar::OperationMode::CW_MODE, 202);
	receiver.setAntenna(&antenna);
	receiver.setTiming(timing_model);

	std::vector<ComplexType> window(3, ComplexType{0.25, -0.5});
	const std::vector<ComplexType> baseline = window;
	const std::vector<core::ActiveStreamingSource> streaming_sources = {
		core::makeActiveSource(&transmitter, params::startTime(), std::numeric_limits<RealType>::max())};
	const std::vector<std::unique_ptr<radar::Target>> targets;
	const RealType start = 0.0;
	const RealType dt = 0.25;
	core::ReceiverTrackerCache tracker_cache;

	processing::pipeline::applyStreamingInterference(window, start, dt, &receiver, streaming_sources, &targets,
													 tracker_cache);

	for (size_t i = 0; i < window.size(); ++i)
	{
		const ComplexType expected = simulation::calculateStreamingDirectPathContribution(
			streaming_sources.front(), &receiver, start + static_cast<RealType>(i) * dt);
		const ComplexType actual = window[i] - baseline[i];
		REQUIRE_THAT(actual.real(), WithinAbs(expected.real(), 1e-12));
		REQUIRE_THAT(actual.imag(), WithinAbs(expected.imag(), 1e-12));
	}
}

TEST_CASE("applyStreamingInterference respects FLAG_NODIRECT and keeps only physically reflected energy",
		  "[processing][finalizer][interference]")
{
	ParamGuard guard;
	params::params.reset();

	radar::Platform tx_platform("TxPlatform");
	setupPlatform(tx_platform, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_platform("RxPlatform");
	setupPlatform(rx_platform, math::Vec3{1000.0, 1000.0, 0.0});
	radar::Platform target_platform("TargetPlatform");
	setupPlatform(target_platform, math::Vec3{1000.0, 0.0, 0.0});

	antenna::Isotropic antenna("iso");
	auto timing_model = makeQuietTiming("clk", 12);

	radar::Transmitter transmitter(&tx_platform, "TxA", radar::OperationMode::CW_MODE, 102);
	transmitter.setAntenna(&antenna);
	transmitter.setTiming(timing_model);
	auto signal = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("cw", 9.0, 1.0e9, 1.0, std::move(signal), 302);
	transmitter.setSignal(&wave);

	radar::Receiver receiver(&rx_platform, "RxA", 100, radar::OperationMode::CW_MODE, 203);
	receiver.setAntenna(&antenna);
	receiver.setTiming(timing_model);
	receiver.setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);

	std::vector<std::unique_ptr<radar::Target>> targets;
	targets.push_back(radar::createIsoTarget(&target_platform, "TargetA", 4.0, 7, 501));

	std::vector<ComplexType> window(3, ComplexType{});
	const std::vector<core::ActiveStreamingSource> streaming_sources = {
		core::makeActiveSource(&transmitter, params::startTime(), std::numeric_limits<RealType>::max())};
	const RealType start = 0.0;
	const RealType dt = 0.2;
	core::ReceiverTrackerCache tracker_cache;

	processing::pipeline::applyStreamingInterference(window, start, dt, &receiver, streaming_sources, &targets,
													 tracker_cache);

	for (size_t i = 0; i < window.size(); ++i)
	{
		const ComplexType expected = simulation::calculateStreamingReflectedPathContribution(
			streaming_sources.front(), &receiver, targets.front().get(), start + static_cast<RealType>(i) * dt);
		REQUIRE_THAT(window[i].real(), WithinAbs(expected.real(), 1e-12));
		REQUIRE_THAT(window[i].imag(), WithinAbs(expected.imag(), 1e-12));
	}
}

TEST_CASE("applyPulsedInterference clips rendered pulses to LO-active sample spans",
		  "[processing][finalizer][interference][dechirp]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(10.0);
	params::setTime(0.0, 1.0);

	radar::Platform tx_platform("TxPlatform");
	setupPlatform(tx_platform, math::Vec3{0.0, 0.0, 0.0});
	antenna::Isotropic antenna("iso");
	auto timing_model = makeQuietTiming("clk", 17);

	radar::Transmitter transmitter(&tx_platform, "PulseTx", radar::OperationMode::PULSED_MODE, 101);
	transmitter.setAntenna(&antenna);
	transmitter.setTiming(timing_model);

	std::vector<ComplexType> iq_buffer(16, ComplexType{0.0, 0.0});
	std::vector<std::unique_ptr<fers_signal::RadarSignal>> wave_store;
	std::vector<std::unique_ptr<serial::Response>> interference_log;
	interference_log.push_back(makeFixedResponse(&transmitter, wave_store,
												 {ComplexType{1.0, 0.0}, ComplexType{2.0, 0.0}, ComplexType{3.0, 0.0},
												  ComplexType{4.0, 0.0}, ComplexType{5.0, 0.0}},
												 params::rate(), 0.2));

	const std::vector<processing::pipeline::SampleSpan> active_spans = {
		processing::pipeline::SampleSpan{.start = 2, .end_exclusive = 3},
		processing::pipeline::SampleSpan{.start = 6, .end_exclusive = 7}};

	processing::pipeline::applyPulsedInterference(iq_buffer, interference_log, active_spans, params::rate());

	for (std::size_t i = 0; i < iq_buffer.size(); ++i)
	{
		const RealType expected = i == 2 ? 1.0 : (i == 6 ? 5.0 : 0.0);
		REQUIRE_THAT(iq_buffer[i].real(), WithinAbs(expected, 1e-12));
	}
}

TEST_CASE("applyPulsedInterference rejects pulse resampling across mismatched rates",
		  "[processing][finalizer][interference][dechirp]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(10.0);
	params::setTime(0.0, 1.0);

	radar::Platform tx_platform("TxPlatform");
	setupPlatform(tx_platform, math::Vec3{0.0, 0.0, 0.0});
	antenna::Isotropic antenna("iso");
	auto timing_model = makeQuietTiming("clk", 17);

	radar::Transmitter transmitter(&tx_platform, "PulseTx", radar::OperationMode::PULSED_MODE, 101);
	transmitter.setAntenna(&antenna);
	transmitter.setTiming(timing_model);

	std::vector<ComplexType> iq_buffer(16, ComplexType{0.0, 0.0});
	std::vector<std::unique_ptr<fers_signal::RadarSignal>> wave_store;
	std::vector<std::unique_ptr<serial::Response>> interference_log;
	interference_log.push_back(makeFixedResponse(&transmitter, wave_store,
												 {ComplexType{1.0, 0.0}, ComplexType{2.0, 0.0}}, params::rate(), 0.2));

	const std::vector<processing::pipeline::SampleSpan> active_spans = {
		processing::pipeline::SampleSpan{.start = 4, .end_exclusive = 5}};

	REQUIRE_THROWS_AS(
		processing::pipeline::applyPulsedInterference(iq_buffer, interference_log, active_spans, 2.0 * params::rate()),
		std::runtime_error);
}

TEST_CASE("applyStreamingInterference adds FMCW energy to pulsed receiver windows",
		  "[processing][finalizer][interference][fmcw]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1.0e6);
	params::setOversampleRatio(1);

	radar::Platform tx_platform("TxPlatform");
	setupPlatform(tx_platform, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_platform("RxPlatform");
	setupPlatform(rx_platform, math::Vec3{300.0, 0.0, 0.0});

	antenna::Isotropic antenna("iso");
	auto timing_model = makeQuietTiming("clk", 13);

	radar::Transmitter transmitter(&tx_platform, "FmcwTx", radar::OperationMode::FMCW_MODE, 103);
	transmitter.setAntenna(&antenna);
	transmitter.setTiming(timing_model);
	auto signal = std::make_unique<fers_signal::FmcwChirpSignal>(1.0e6, 50.0e-6, 100.0e-6);
	fers_signal::RadarSignal wave("fmcw", 16.0, 1.0e9, 50.0e-6, std::move(signal), 303);
	transmitter.setSignal(&wave);

	radar::Receiver receiver(&rx_platform, "PulsedRx", 101, radar::OperationMode::PULSED_MODE, 204);
	receiver.setAntenna(&antenna);
	receiver.setTiming(timing_model);
	receiver.setWindowProperties(150.0e-6, 1.0e3, 0.0);

	std::vector<ComplexType> window(3, ComplexType{});
	const std::vector<core::ActiveStreamingSource> streaming_sources = {
		core::makeActiveSource(&transmitter, 0.0, 300.0e-6)};
	const std::vector<std::unique_ptr<radar::Target>> targets;
	core::ReceiverTrackerCache tracker_cache;

	processing::pipeline::applyStreamingInterference(window, 10.0e-6, 50.0e-6, &receiver, streaming_sources, &targets,
													 tracker_cache);

	REQUIRE(std::abs(window[0]) > 0.0);
	REQUIRE_THAT(std::abs(window[1]), WithinAbs(0.0, 1.0e-18));
	REQUIRE(std::abs(window[2]) > 0.0);
}

TEST_CASE("applyStreamingInterference supports FMCW transmitter with CW streaming receiver",
		  "[processing][finalizer][interference][fmcw]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1.0e6);
	params::setOversampleRatio(1);

	radar::Platform tx_platform("TxPlatform");
	setupPlatform(tx_platform, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_platform("RxPlatform");
	setupPlatform(rx_platform, math::Vec3{300.0, 0.0, 0.0});

	antenna::Isotropic antenna("iso");
	auto timing_model = makeQuietTiming("clk", 14);

	radar::Transmitter transmitter(&tx_platform, "FmcwTx", radar::OperationMode::FMCW_MODE, 104);
	transmitter.setAntenna(&antenna);
	transmitter.setTiming(timing_model);
	auto signal = std::make_unique<fers_signal::FmcwChirpSignal>(1.0e6, 50.0e-6, 100.0e-6);
	fers_signal::RadarSignal wave("fmcw", 16.0, 1.0e9, 50.0e-6, std::move(signal), 304);
	transmitter.setSignal(&wave);

	radar::Receiver receiver(&rx_platform, "CwRx", 102, radar::OperationMode::CW_MODE, 205);
	receiver.setAntenna(&antenna);
	receiver.setTiming(timing_model);

	std::vector<ComplexType> window(3, ComplexType{0.1, -0.2});
	const std::vector<ComplexType> baseline = window;
	const std::vector<core::ActiveStreamingSource> streaming_sources = {
		core::makeActiveSource(&transmitter, 0.0, 300.0e-6)};
	const std::vector<std::unique_ptr<radar::Target>> targets;
	core::ReceiverTrackerCache tracker_cache;

	processing::pipeline::applyStreamingInterference(window, 10.0e-6, 50.0e-6, &receiver, streaming_sources, &targets,
													 tracker_cache);

	for (std::size_t i = 0; i < window.size(); ++i)
	{
		const ComplexType expected = simulation::calculateStreamingDirectPathContribution(
			streaming_sources.front(), &receiver, 10.0e-6 + static_cast<RealType>(i) * 50.0e-6);
		const ComplexType actual = window[i] - baseline[i];
		REQUIRE_THAT(actual.real(), WithinAbs(expected.real(), 1.0e-12));
		REQUIRE_THAT(actual.imag(), WithinAbs(expected.imag(), 1.0e-12));
	}
}

TEST_CASE("applyStreamingInterference superposes up- and down-chirp FMCW transmitters",
		  "[processing][finalizer][interference][fmcw]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1.0e6);
	params::setOversampleRatio(1);

	radar::Platform tx_up_platform("TxUpPlatform");
	setupPlatform(tx_up_platform, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform tx_down_platform("TxDownPlatform");
	setupPlatform(tx_down_platform, math::Vec3{150.0, 0.0, 0.0});
	radar::Platform rx_platform("RxPlatform");
	setupPlatform(rx_platform, math::Vec3{300.0, 0.0, 0.0});

	antenna::Isotropic antenna("iso");
	auto timing_model = makeQuietTiming("clk", 15);

	radar::Transmitter up_tx(&tx_up_platform, "UpTx", radar::OperationMode::FMCW_MODE, 105);
	up_tx.setAntenna(&antenna);
	up_tx.setTiming(timing_model);
	auto up_signal = std::make_unique<fers_signal::FmcwChirpSignal>(1.0e6, 50.0e-6, 100.0e-6);
	fers_signal::RadarSignal up_wave("up_fmcw", 16.0, 1.0e9, 50.0e-6, std::move(up_signal), 305);
	up_tx.setSignal(&up_wave);

	radar::Transmitter down_tx(&tx_down_platform, "DownTx", radar::OperationMode::FMCW_MODE, 106);
	down_tx.setAntenna(&antenna);
	down_tx.setTiming(timing_model);
	auto down_signal = std::make_unique<fers_signal::FmcwChirpSignal>(1.0e6, 50.0e-6, 100.0e-6, 0.0, std::nullopt,
																	  fers_signal::FmcwChirpDirection::Down);
	fers_signal::RadarSignal down_wave("down_fmcw", 16.0, 1.0e9, 50.0e-6, std::move(down_signal), 306);
	down_tx.setSignal(&down_wave);

	radar::Receiver receiver(&rx_platform, "CwRx", 103, radar::OperationMode::CW_MODE, 206);
	receiver.setAntenna(&antenna);
	receiver.setTiming(timing_model);

	std::vector<ComplexType> window(3, ComplexType{0.1, -0.2});
	const std::vector<ComplexType> baseline = window;
	const std::vector<core::ActiveStreamingSource> streaming_sources = {
		core::makeActiveSource(&up_tx, 0.0, 300.0e-6), core::makeActiveSource(&down_tx, 0.0, 300.0e-6)};
	const std::vector<std::unique_ptr<radar::Target>> targets;
	core::ReceiverTrackerCache tracker_cache;

	processing::pipeline::applyStreamingInterference(window, 10.0e-6, 50.0e-6, &receiver, streaming_sources, &targets,
													 tracker_cache);

	for (std::size_t i = 0; i < window.size(); ++i)
	{
		const RealType sample_time = 10.0e-6 + static_cast<RealType>(i) * 50.0e-6;
		const ComplexType expected =
			simulation::calculateStreamingDirectPathContribution(streaming_sources[0], &receiver, sample_time) +
			simulation::calculateStreamingDirectPathContribution(streaming_sources[1], &receiver, sample_time);
		const ComplexType actual = window[i] - baseline[i];
		REQUIRE_THAT(actual.real(), WithinAbs(expected.real(), 1.0e-12));
		REQUIRE_THAT(actual.imag(), WithinAbs(expected.imag(), 1.0e-12));
	}
}

TEST_CASE("applyStreamingInterference reuses tracker cache without carrying window state",
		  "[processing][finalizer][interference][fmcw]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1.0e6);
	params::setOversampleRatio(1);

	radar::Platform tx_platform("TxPlatform");
	setupPlatform(tx_platform, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_platform("RxPlatform");
	setupPlatform(rx_platform, math::Vec3{300.0, 0.0, 0.0});
	radar::Platform target_platform("TargetPlatform");
	setupPlatform(target_platform, math::Vec3{150.0, 100.0, 0.0});

	antenna::Isotropic antenna("iso");
	auto timing_model = makeQuietTiming("clk", 16);

	radar::Transmitter transmitter(&tx_platform, "FmcwTx", radar::OperationMode::FMCW_MODE, 107);
	transmitter.setAntenna(&antenna);
	transmitter.setTiming(timing_model);
	auto signal = std::make_unique<fers_signal::FmcwChirpSignal>(1.0e6, 50.0e-6, 100.0e-6);
	fers_signal::RadarSignal wave("fmcw", 16.0, 1.0e9, 50.0e-6, std::move(signal), 307);
	transmitter.setSignal(&wave);

	radar::Receiver receiver(&rx_platform, "PulsedRx", 104, radar::OperationMode::PULSED_MODE, 207);
	receiver.setAntenna(&antenna);
	receiver.setTiming(timing_model);
	receiver.setWindowProperties(150.0e-6, 1.0e3, 0.0);

	std::vector<std::unique_ptr<radar::Target>> targets;
	targets.push_back(radar::createIsoTarget(&target_platform, "TargetA", 2.0, 8, 502));

	const std::vector<core::ActiveStreamingSource> streaming_sources = {
		core::makeActiveSource(&transmitter, 0.0, 300.0e-6)};
	core::ReceiverTrackerCache tracker_cache;
	std::vector<ComplexType> first_window(3, ComplexType{});
	std::vector<ComplexType> second_window(3, ComplexType{});

	processing::pipeline::applyStreamingInterference(first_window, 10.0e-6, 50.0e-6, &receiver, streaming_sources,
													 &targets, tracker_cache);

	REQUIRE(tracker_cache.direct.size() == 1);
	REQUIRE(tracker_cache.reflected.size() == 1);
	REQUIRE(tracker_cache.reflected.front().size() == 1);
	const std::size_t direct_capacity = tracker_cache.direct.capacity();
	const std::size_t reflected_capacity = tracker_cache.reflected.capacity();
	const std::size_t reflected_row_capacity = tracker_cache.reflected.front().capacity();

	processing::pipeline::applyStreamingInterference(second_window, 10.0e-6, 50.0e-6, &receiver, streaming_sources,
													 &targets, tracker_cache);

	REQUIRE(tracker_cache.direct.capacity() == direct_capacity);
	REQUIRE(tracker_cache.reflected.capacity() == reflected_capacity);
	REQUIRE(tracker_cache.reflected.front().capacity() == reflected_row_capacity);
	REQUIRE(std::abs(first_window.front()) > 0.0);
	for (std::size_t i = 0; i < first_window.size(); ++i)
	{
		REQUIRE_THAT(second_window[i].real(), WithinAbs(first_window[i].real(), 1.0e-12));
		REQUIRE_THAT(second_window[i].imag(), WithinAbs(first_window[i].imag(), 1.0e-12));
	}
}

TEST_CASE("applyPulsedInterference maps pulse start times to simulation sample indices and clips overflow",
		  "[processing][finalizer][interference]")
{
	ParamGuard guard;
	params::setTime(10.0, 11.0);
	params::setRate(4.0);
	params::setOversampleRatio(1);

	radar::Platform tx_platform("TxPlatform");
	radar::Transmitter transmitter(&tx_platform, "TxA", radar::OperationMode::PULSED_MODE, 401);

	std::vector<std::unique_ptr<fers_signal::RadarSignal>> wave_store;
	std::vector<std::unique_ptr<serial::Response>> interference_log;
	interference_log.push_back(makeFixedResponse(
		&transmitter, wave_store, {ComplexType{1.0, 0.5}, ComplexType{2.0, -0.5}, ComplexType{3.0, 1.0}}, 4.0, 10.25));
	interference_log.push_back(
		makeFixedResponse(&transmitter, wave_store,
						  {ComplexType{10.0, 0.0}, ComplexType{20.0, 1.0}, ComplexType{30.0, 2.0}}, 4.0, 10.75));

	std::vector<ComplexType> iq_buffer(5, ComplexType{});

	processing::pipeline::applyPulsedInterference(iq_buffer, interference_log);

	const std::vector<ComplexType> expected = {
		ComplexType{0.0, 0.0},	ComplexType{1.0, 0.5},	ComplexType{2.0, -0.5},
		ComplexType{13.0, 1.0}, ComplexType{20.0, 1.0},
	};

	REQUIRE(iq_buffer.size() == expected.size());
	for (size_t i = 0; i < expected.size(); ++i)
	{
		REQUIRE_THAT(iq_buffer[i].real(), WithinAbs(expected[i].real(), 1e-12));
		REQUIRE_THAT(iq_buffer[i].imag(), WithinAbs(expected[i].imag(), 1e-12));
	}
}

TEST_CASE("applyPulsedInterference uses RF simulation rate for oversampled full-buffer path",
		  "[processing][finalizer][interference]")
{
	ParamGuard guard;
	params::setTime(10.0, 11.0);
	params::setRate(4.0);
	params::setOversampleRatio(2);

	radar::Platform tx_platform("TxPlatform");
	radar::Transmitter transmitter(&tx_platform, "TxA", radar::OperationMode::PULSED_MODE, 402);

	std::vector<std::unique_ptr<fers_signal::RadarSignal>> wave_store;
	std::vector<std::unique_ptr<serial::Response>> interference_log;
	interference_log.push_back(makeFixedResponse(&transmitter, wave_store,
												 {ComplexType{1.0, 0.5}, ComplexType{2.0, -0.5}, ComplexType{3.0, 1.0}},
												 params::rate(), 10.25));

	std::vector<ComplexType> iq_buffer(8, ComplexType{});

	processing::pipeline::applyPulsedInterference(iq_buffer, interference_log);

	const std::vector<ComplexType> expected = {
		ComplexType{0.0, 0.0}, ComplexType{0.0, 0.0}, ComplexType{1.0, 0.5}, ComplexType{2.0, -0.5},
		ComplexType{3.0, 1.0}, ComplexType{0.0, 0.0}, ComplexType{0.0, 0.0}, ComplexType{0.0, 0.0},
	};

	REQUIRE(iq_buffer.size() == expected.size());
	for (size_t i = 0; i < expected.size(); ++i)
	{
		REQUIRE_THAT(iq_buffer[i].real(), WithinAbs(expected[i].real(), 1e-12));
		REQUIRE_THAT(iq_buffer[i].imag(), WithinAbs(expected[i].imag(), 1e-12));
	}
}
