// Tests for channel_model reflected-path calculations: solveRe and
// calculateReflectedPathContribution. Verifies bistatic radar range equation,
// propagation delay, phase shift, and RCS scaling against hand-calculated values.

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <cmath>
#include <memory>
#include <optional>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/parameters.h"
#include "math/coord.h"
#include "math/geometry_ops.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "signal/radar_signal.h"
#include "simulation/channel_model.h"
#include "timing/prototype_timing.h"
#include "timing/timing.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
	};

	void setupPlatform(radar::Platform& plat, const math::Vec3& pos)
	{
		plat.getMotionPath()->addCoord(math::Coord{pos, 0.0});
		plat.getMotionPath()->finalize();
		plat.getRotationPath()->addCoord(math::RotationCoord{0.0, 0.0, 0.0});
		plat.getRotationPath()->finalize();
	}

	simulation::CwPhaseNoiseLookup makeLookup(const std::shared_ptr<timing::Timing>& timing)
	{
		const std::vector<std::shared_ptr<timing::Timing>> timings = {timing};
		return simulation::CwPhaseNoiseLookup::build(timings, params::startTime(), params::endTime());
	}

	RealType unwrapDelta(RealType delta)
	{
		while (delta > PI)
		{
			delta -= 2.0 * PI;
		}
		while (delta < -PI)
		{
			delta += 2.0 * PI;
		}
		return delta;
	}
}

// =============================================================================
// solveRe: Bistatic radar range equation
// =============================================================================

TEST_CASE("solveRe computes correct bistatic power for isotropic antennas", "[simulation][channel_model][reflected]")
{
	ParamGuard guard;
	params::params.reset();

	// Setup: Tx at (0,0,0), Target at (500,0,0), Rx at (500,500,0)
	// Carrier: 1 GHz, RCS: 10 m^2
	// Bistatic: Pr/Pt = Gt * Gr * sigma * lambda^2 / (64 * pi^3 * R1^2 * R2^2)

	const RealType c = params::c();
	const RealType carrier = 1.0e9;
	const RealType lambda = c / carrier;
	const RealType rcs = 10.0;

	const math::Vec3 tx_pos{0.0, 0.0, 0.0};
	const math::Vec3 tgt_pos{500.0, 0.0, 0.0};
	const math::Vec3 rx_pos{500.0, 500.0, 0.0};

	const RealType r1 = (tgt_pos - tx_pos).length(); // 500 m
	const RealType r2 = (rx_pos - tgt_pos).length(); // 500 m

	const RealType expected_power = (rcs * lambda * lambda) / (64.0 * PI * PI * PI * r1 * r1 * r2 * r2);
	const RealType expected_delay = (r1 + r2) / c;
	const RealType expected_phase = -expected_delay * 2.0 * PI * carrier;

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, tx_pos);

	radar::Platform tgt_plat("tgt_plat");
	setupPlatform(tgt_plat, tgt_pos);

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, rx_pos);

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("sig", 1.0, carrier, 1e-6, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	radar::IsoTarget tgt(&tgt_plat, "tgt", rcs, 42);

	simulation::ReResults results{};
	const auto time = std::chrono::duration<RealType>(0.0);

	simulation::solveRe(&tx, &rx, &tgt, time, &wave, results);

	REQUIRE_THAT(results.power, WithinRel(expected_power, 1e-6));
	REQUIRE_THAT(results.delay, WithinRel(expected_delay, 1e-9));
	REQUIRE_THAT(results.phase, WithinRel(expected_phase, 1e-9));
}

TEST_CASE("CW streaming reflected path gates schedules by retarded transmit time",
		  "[simulation][channel_model][reflected]")
{
	ParamGuard guard;
	params::params.reset();
	params::setC(100.0);

	const RealType tx_target_dist = 5.0;
	const RealType target_rx_dist = 5.0;
	const RealType tau = (tx_target_dist + target_rx_dist) / params::c();
	const RealType segment_start = 0.2;
	const RealType segment_end = 0.5;
	const RealType eps = 1.0e-6;

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform target_plat("target_plat");
	setupPlatform(target_plat, math::Vec3{tx_target_dist, 0.0, 0.0});
	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{tx_target_dist + target_rx_dist, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::CW_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);
	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("cw", 1.0, 1.0e6, 1.0, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::CW_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	radar::IsoTarget target(&target_plat, "target", 1.0, 7);
	const core::ActiveStreamingSource source = core::makeActiveSource(&tx, segment_start, segment_end);

	REQUIRE(std::abs(simulation::calculateStreamingReflectedPathContribution(source, &rx, &target,
																			 segment_start + tau - eps)) == 0.0);
	REQUIRE(std::abs(simulation::calculateStreamingReflectedPathContribution(source, &rx, &target,
																			 segment_start + tau + eps)) > 0.0);
	REQUIRE(std::abs(simulation::calculateStreamingReflectedPathContribution(source, &rx, &target,
																			 segment_end + 0.5 * tau)) > 0.0);
	REQUIRE(std::abs(simulation::calculateStreamingReflectedPathContribution(source, &rx, &target,
																			 segment_end + tau + eps)) == 0.0);
}

TEST_CASE("FMCW monostatic reflected path dechirps to expected stationary-target beat frequency",
		  "[simulation][channel_model][reflected][fmcw]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1.0e6);
	params::setSimSamplingRate(1.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 1.0e-3);

	const RealType target_range = 150.0;
	const RealType chirp_bandwidth = 1.0e6;
	const RealType chirp_duration = 1.0e-3;
	const RealType chirp_rate = chirp_bandwidth / chirp_duration;
	const RealType tau = (2.0 * target_range) / params::c();

	radar::Platform radar_platform("radar");
	setupPlatform(radar_platform, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform target_platform("target_platform");
	setupPlatform(target_platform, math::Vec3{target_range, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&radar_platform, "tx", radar::OperationMode::FMCW_MODE, 101);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);
	auto fmcw_signal = std::make_unique<fers_signal::FmcwChirpSignal>(chirp_bandwidth, chirp_duration, chirp_duration);
	fers_signal::RadarSignal wave("fmcw", 1.0, 10.0e6, chirp_duration, std::move(fmcw_signal), 301);
	tx.setSignal(&wave);

	radar::Receiver rx(&radar_platform, "rx", 43, radar::OperationMode::FMCW_MODE, 201);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);
	rx.setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
	tx.setAttached(&rx);
	rx.setAttached(&tx);

	auto target = radar::createIsoTarget(&target_platform, "target", 1.0, 7, 401);
	const core::ActiveStreamingSource source = core::makeActiveSource(&tx, 0.0, chirp_duration);
	core::FmcwChirpBoundaryTracker tracker;

	const RealType dt = 1.0 / params::simSamplingRate();
	const RealType first_time = 20.0e-6;
	const std::size_t sample_count = 500;
	std::vector<RealType> dechirped_phase;
	dechirped_phase.reserve(sample_count);

	for (std::size_t i = 0; i < sample_count; ++i)
	{
		const RealType t = first_time + static_cast<RealType>(i) * dt;
		const ComplexType sample =
			simulation::calculateStreamingReflectedPathContribution(source, &rx, target.get(), t, nullptr, &tracker);
		const ComplexType dechirped =
			sample * std::polar(1.0, -wave.getFmcwChirpSignal()->basebandPhaseForChirpTime(t));
		dechirped_phase.push_back(std::arg(dechirped));
	}

	RealType unwrapped_span = 0.0;
	for (std::size_t i = 1; i < dechirped_phase.size(); ++i)
	{
		unwrapped_span += unwrapDelta(dechirped_phase[i] - dechirped_phase[i - 1]);
	}

	const RealType measured_beat_hz =
		unwrapped_span / (2.0 * PI * dt * static_cast<RealType>(dechirped_phase.size() - 1));
	const RealType expected_beat_hz = -chirp_rate * tau;
	REQUIRE_THAT(measured_beat_hz, WithinRel(expected_beat_hz, 1.0e-3));
}

TEST_CASE("FMCW native dechirp convention produces positive up-chirp beat frequency",
		  "[simulation][channel_model][reflected][fmcw][dechirp]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1.0e6);
	params::setSimSamplingRate(1.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 1.0e-3);

	const RealType target_range = 150.0;
	const RealType chirp_bandwidth = 1.0e6;
	const RealType chirp_duration = 1.0e-3;
	const RealType chirp_rate = chirp_bandwidth / chirp_duration;
	const RealType tau = (2.0 * target_range) / params::c();

	radar::Platform radar_platform("radar");
	setupPlatform(radar_platform, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform target_platform("target_platform");
	setupPlatform(target_platform, math::Vec3{target_range, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&radar_platform, "tx", radar::OperationMode::FMCW_MODE, 101);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);
	auto fmcw_signal = std::make_unique<fers_signal::FmcwChirpSignal>(chirp_bandwidth, chirp_duration, chirp_duration);
	fers_signal::RadarSignal wave("fmcw", 1.0, 10.0e6, chirp_duration, std::move(fmcw_signal), 301);
	tx.setSignal(&wave);

	radar::Receiver rx(&radar_platform, "rx", 43, radar::OperationMode::FMCW_MODE, 201);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);
	rx.setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
	tx.setAttached(&rx);
	rx.setAttached(&tx);

	auto target = radar::createIsoTarget(&target_platform, "target", 1.0, 7, 401);
	const core::ActiveStreamingSource source = core::makeActiveSource(&tx, 0.0, chirp_duration);
	core::FmcwChirpBoundaryTracker channel_tracker;
	core::FmcwChirpBoundaryTracker reference_tracker;

	const RealType dt = 1.0 / params::simSamplingRate();
	const RealType first_time = 20.0e-6;
	const std::size_t sample_count = 500;
	std::vector<RealType> dechirped_phase;
	dechirped_phase.reserve(sample_count);

	for (std::size_t i = 0; i < sample_count; ++i)
	{
		const RealType t = first_time + static_cast<RealType>(i) * dt;
		const ComplexType sample = simulation::calculateStreamingReflectedPathContribution(
			source, &rx, target.get(), t, nullptr, &channel_tracker, simulation::StreamingTimingPhaseMode::None);
		RealType reference_phase = 0.0;
		REQUIRE(simulation::calculateStreamingReferencePhase(source, t, &reference_tracker, reference_phase));
		const ComplexType dechirped = std::polar(1.0, reference_phase) * std::conj(sample);
		dechirped_phase.push_back(std::arg(dechirped));
	}

	RealType unwrapped_span = 0.0;
	for (std::size_t i = 1; i < dechirped_phase.size(); ++i)
	{
		unwrapped_span += unwrapDelta(dechirped_phase[i] - dechirped_phase[i - 1]);
	}

	const RealType measured_beat_hz =
		unwrapped_span / (2.0 * PI * dt * static_cast<RealType>(dechirped_phase.size() - 1));
	const RealType expected_beat_hz = chirp_rate * tau;
	REQUIRE_THAT(measured_beat_hz, WithinRel(expected_beat_hz, 1.0e-3));
}

TEST_CASE("FMCW physical dechirp preserves timing decorrelation while ideal mode removes it",
		  "[simulation][channel_model][reflected][fmcw][dechirp]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(2.0e6);
	params::setSimSamplingRate(2.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 1.0e-3);

	const RealType target_range = 150.0;
	const RealType chirp_bandwidth = 1.0e6;
	const RealType chirp_duration = 1.0e-3;
	const RealType tau = (2.0 * target_range) / params::c();
	const RealType freq_offset = 250.0e3;
	const RealType sample_time = 40.0e-6;

	radar::Platform radar_platform("radar");
	setupPlatform(radar_platform, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform target_platform("target_platform");
	setupPlatform(target_platform, math::Vec3{target_range, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	timing::PrototypeTiming prototype("clk");
	prototype.setFrequency(10.0e6);
	prototype.setFreqOffset(freq_offset);
	auto timing = std::make_shared<timing::Timing>("clk", 42);
	timing->initializeModel(&prototype);

	radar::Transmitter tx(&radar_platform, "tx", radar::OperationMode::FMCW_MODE, 101);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);
	auto fmcw_signal = std::make_unique<fers_signal::FmcwChirpSignal>(chirp_bandwidth, chirp_duration, chirp_duration);
	fers_signal::RadarSignal wave("fmcw", 1.0, 10.0e6, chirp_duration, std::move(fmcw_signal), 301);
	tx.setSignal(&wave);

	radar::Receiver rx(&radar_platform, "rx", 43, radar::OperationMode::FMCW_MODE, 201);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);
	rx.setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
	tx.setAttached(&rx);
	rx.setAttached(&tx);

	auto target = radar::createIsoTarget(&target_platform, "target", 1.0, 7, 401);
	const core::ActiveStreamingSource source = core::makeActiveSource(&tx, 0.0, chirp_duration);
	core::FmcwChirpBoundaryTracker physical_tracker;
	core::FmcwChirpBoundaryTracker ideal_tracker;
	core::FmcwChirpBoundaryTracker reference_tracker;

	const ComplexType physical_channel = simulation::calculateStreamingReflectedPathContribution(
		source, &rx, target.get(), sample_time, nullptr, &physical_tracker,
		simulation::StreamingTimingPhaseMode::TransmitterOnly);
	const ComplexType ideal_channel = simulation::calculateStreamingReflectedPathContribution(
		source, &rx, target.get(), sample_time, nullptr, &ideal_tracker, simulation::StreamingTimingPhaseMode::None);

	RealType reference_phase = 0.0;
	REQUIRE(simulation::calculateStreamingReferencePhase(source, sample_time, &reference_tracker, reference_phase));

	const RealType rx_phase = 2.0 * PI * freq_offset * sample_time;
	const ComplexType physical_dechirped = std::polar(1.0, reference_phase + rx_phase) * std::conj(physical_channel);
	const ComplexType ideal_dechirped = std::polar(1.0, reference_phase) * std::conj(ideal_channel);
	const RealType measured_delta = std::arg(physical_dechirped * std::conj(ideal_dechirped));
	const RealType expected_delta = 2.0 * PI * freq_offset * tau;

	REQUIRE_THAT(measured_delta, WithinAbs(expected_delta, 1.0e-6));
}

TEST_CASE("FMCW down-chirp monostatic reflected path reverses stationary-target beat sign",
		  "[simulation][channel_model][reflected][fmcw]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1.0e6);
	params::setSimSamplingRate(1.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 1.0e-3);

	const RealType target_range = 150.0;
	const RealType chirp_bandwidth = 1.0e6;
	const RealType chirp_duration = 1.0e-3;
	const RealType chirp_rate = chirp_bandwidth / chirp_duration;
	const RealType tau = (2.0 * target_range) / params::c();

	radar::Platform radar_platform("radar");
	setupPlatform(radar_platform, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform target_platform("target_platform");
	setupPlatform(target_platform, math::Vec3{target_range, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&radar_platform, "tx", radar::OperationMode::FMCW_MODE, 101);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);
	auto fmcw_signal = std::make_unique<fers_signal::FmcwChirpSignal>(
		chirp_bandwidth, chirp_duration, chirp_duration, 0.0, std::nullopt, fers_signal::FmcwChirpDirection::Down);
	fers_signal::RadarSignal wave("fmcw", 1.0, 10.0e6, chirp_duration, std::move(fmcw_signal), 301);
	tx.setSignal(&wave);

	radar::Receiver rx(&radar_platform, "rx", 43, radar::OperationMode::FMCW_MODE, 201);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);
	rx.setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
	tx.setAttached(&rx);
	rx.setAttached(&tx);

	auto target = radar::createIsoTarget(&target_platform, "target", 1.0, 7, 401);
	const core::ActiveStreamingSource source = core::makeActiveSource(&tx, 0.0, chirp_duration);
	core::FmcwChirpBoundaryTracker tracker;

	const RealType dt = 1.0 / params::simSamplingRate();
	const RealType first_time = 20.0e-6;
	const std::size_t sample_count = 500;
	std::vector<RealType> dechirped_phase;
	dechirped_phase.reserve(sample_count);

	for (std::size_t i = 0; i < sample_count; ++i)
	{
		const RealType t = first_time + static_cast<RealType>(i) * dt;
		const ComplexType sample =
			simulation::calculateStreamingReflectedPathContribution(source, &rx, target.get(), t, nullptr, &tracker);
		const ComplexType dechirped =
			sample * std::polar(1.0, -wave.getFmcwChirpSignal()->basebandPhaseForChirpTime(t));
		dechirped_phase.push_back(std::arg(dechirped));
	}

	RealType unwrapped_span = 0.0;
	for (std::size_t i = 1; i < dechirped_phase.size(); ++i)
	{
		unwrapped_span += unwrapDelta(dechirped_phase[i] - dechirped_phase[i - 1]);
	}

	const RealType measured_beat_hz =
		unwrapped_span / (2.0 * PI * dt * static_cast<RealType>(dechirped_phase.size() - 1));
	const RealType expected_beat_hz = chirp_rate * tau;
	REQUIRE_THAT(measured_beat_hz, WithinRel(expected_beat_hz, 1.0e-3));
}

TEST_CASE("FMCW chirp-boundary tracker matches cold-path direct contribution across chirps",
		  "[simulation][channel_model][direct][fmcw]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1.0e6);
	params::setSimSamplingRate(1.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 3.0e-4);

	radar::Platform tx_platform("tx_platform");
	setupPlatform(tx_platform, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_platform("rx_platform");
	setupPlatform(rx_platform, math::Vec3{300.0, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 44);

	radar::Transmitter tx(&tx_platform, "tx", radar::OperationMode::FMCW_MODE, 102);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);
	auto fmcw_signal = std::make_unique<fers_signal::FmcwChirpSignal>(2.0e6, 2.0e-5, 5.0e-5);
	fers_signal::RadarSignal wave("fmcw", 5.0, 20.0e6, 2.0e-5, std::move(fmcw_signal), 302);
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_platform, "rx", 45, radar::OperationMode::CW_MODE, 202);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	const core::ActiveStreamingSource source = core::makeActiveSource(&tx, 0.0, 3.0e-4);
	core::FmcwChirpBoundaryTracker tracker;
	const RealType dt = 1.0 / params::simSamplingRate();

	for (std::size_t i = 0; i < 260; ++i)
	{
		const RealType t = static_cast<RealType>(i) * dt;
		const ComplexType tracked =
			simulation::calculateStreamingDirectPathContribution(source, &rx, t, nullptr, &tracker);
		const ComplexType cold = simulation::calculateStreamingDirectPathContribution(source, &rx, t);
		REQUIRE_THAT(tracked.real(), WithinAbs(cold.real(), 1.0e-12));
		REQUIRE_THAT(tracked.imag(), WithinAbs(cold.imag(), 1.0e-12));
	}
	REQUIRE(tracker.initialized);
	REQUIRE(tracker.n_current > 3);

	BENCHMARK("FMCW tracker hot path direct contribution")
	{
		core::FmcwChirpBoundaryTracker bench_tracker;
		ComplexType acc{0.0, 0.0};
		for (std::size_t i = 0; i < 260; ++i)
		{
			acc += simulation::calculateStreamingDirectPathContribution(source, &rx, static_cast<RealType>(i) * dt,
																		nullptr, &bench_tracker);
		}
		return acc;
	};
}

TEST_CASE("solveRe power scales linearly with RCS", "[simulation][channel_model][reflected]")
{
	ParamGuard guard;
	params::params.reset();

	// RCS is a linear multiplier in the bistatic equation
	// Doubling RCS should double the power

	const RealType carrier = 1.0e9;

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	// RCS = 5 m^2
	radar::Platform tx_plat1("tx1");
	setupPlatform(tx_plat1, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform tgt_plat1("tgt1");
	setupPlatform(tgt_plat1, math::Vec3{1000.0, 0.0, 0.0});
	radar::Platform rx_plat1("rx1");
	setupPlatform(rx_plat1, math::Vec3{2000.0, 0.0, 0.0});

	radar::Transmitter tx1(&tx_plat1, "tx1", radar::OperationMode::PULSED_MODE);
	tx1.setAntenna(&iso_ant);
	tx1.setTiming(timing);
	auto sig1 = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave1("sig1", 1.0, carrier, 1e-6, std::move(sig1));
	tx1.setSignal(&wave1);

	radar::Receiver rx1(&rx_plat1, "rx1", 42, radar::OperationMode::PULSED_MODE);
	rx1.setAntenna(&iso_ant);
	rx1.setTiming(timing);

	radar::IsoTarget tgt1(&tgt_plat1, "tgt1", 5.0, 42);

	simulation::ReResults r1{};
	simulation::solveRe(&tx1, &rx1, &tgt1, std::chrono::duration<RealType>(0.0), &wave1, r1);

	// RCS = 10 m^2 (doubled)
	radar::Platform tx_plat2("tx2");
	setupPlatform(tx_plat2, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform tgt_plat2("tgt2");
	setupPlatform(tgt_plat2, math::Vec3{1000.0, 0.0, 0.0});
	radar::Platform rx_plat2("rx2");
	setupPlatform(rx_plat2, math::Vec3{2000.0, 0.0, 0.0});

	radar::Transmitter tx2(&tx_plat2, "tx2", radar::OperationMode::PULSED_MODE);
	tx2.setAntenna(&iso_ant);
	tx2.setTiming(timing);
	auto sig2 = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave2("sig2", 1.0, carrier, 1e-6, std::move(sig2));
	tx2.setSignal(&wave2);

	radar::Receiver rx2(&rx_plat2, "rx2", 42, radar::OperationMode::PULSED_MODE);
	rx2.setAntenna(&iso_ant);
	rx2.setTiming(timing);

	radar::IsoTarget tgt2(&tgt_plat2, "tgt2", 10.0, 42);

	simulation::ReResults r2{};
	simulation::solveRe(&tx2, &rx2, &tgt2, std::chrono::duration<RealType>(0.0), &wave2, r2);

	REQUIRE_THAT(r2.power / r1.power, WithinRel(2.0, 1e-9));
}

TEST_CASE("solveRe delay is sum of Tx-Tgt and Tgt-Rx distances divided by c", "[simulation][channel_model][reflected]")
{
	ParamGuard guard;
	params::params.reset();

	const RealType c = params::c();
	const RealType carrier = 1.0e9;

	// Collinear case: Tx at origin, Tgt at 1000 m, Rx at 3000 m
	// r1 = 1000, r2 = 2000
	const RealType r1 = 1000.0;
	const RealType r2 = 2000.0;
	const RealType expected_delay = (r1 + r2) / c;

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform tgt_plat("tgt_plat");
	setupPlatform(tgt_plat, math::Vec3{r1, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{r1 + r2, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("sig", 1.0, carrier, 1e-6, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	radar::IsoTarget tgt(&tgt_plat, "tgt", 1.0, 42);

	simulation::ReResults results{};
	simulation::solveRe(&tx, &rx, &tgt, std::chrono::duration<RealType>(0.0), &wave, results);

	REQUIRE_THAT(results.delay, WithinRel(expected_delay, 1e-9));
}

TEST_CASE("solveRe with noproploss flag ignores distance in power", "[simulation][channel_model][reflected]")
{
	ParamGuard guard;
	params::params.reset();

	const RealType c = params::c();
	const RealType carrier = 1.0e9;
	const RealType lambda = c / carrier;
	const RealType rcs = 10.0;

	// With no_prop_loss: Pr/Pt = Gt * Gr * sigma * lambda^2 / (64 * pi^3) (no R terms)
	const RealType expected_power = (rcs * lambda * lambda) / (64.0 * PI * PI * PI);

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform tgt_plat("tgt_plat");
	setupPlatform(tgt_plat, math::Vec3{2000.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{2000.0, 3000.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("sig", 1.0, carrier, 1e-6, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);
	rx.setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);

	radar::IsoTarget tgt(&tgt_plat, "tgt", rcs, 42);

	simulation::ReResults results{};
	simulation::solveRe(&tx, &rx, &tgt, std::chrono::duration<RealType>(0.0), &wave, results);

	REQUIRE_THAT(results.power, WithinRel(expected_power, 1e-6));
}

TEST_CASE("solveRe power follows R^-4 for monostatic geometry", "[simulation][channel_model][reflected]")
{
	ParamGuard guard;
	params::params.reset();

	// For monostatic-like geometry (Tx and Rx co-located, target along x-axis):
	// Pr/Pt = sigma * lambda^2 / (64 * pi^3 * R^4) for r1 = r2 = R
	// Doubling R should reduce power by factor of 16

	const RealType carrier = 1.0e9;

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	// R = 500 m
	radar::Platform tx_plat1("tx1");
	setupPlatform(tx_plat1, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_plat1("rx1");
	setupPlatform(rx_plat1, math::Vec3{0.0, 0.1, 0.0}); // Slightly offset to avoid RangeError with target
	radar::Platform tgt_plat1("tgt1");
	setupPlatform(tgt_plat1, math::Vec3{500.0, 0.0, 0.0});

	radar::Transmitter tx1(&tx_plat1, "tx1", radar::OperationMode::PULSED_MODE);
	tx1.setAntenna(&iso_ant);
	tx1.setTiming(timing);
	auto sig1 = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave1("sig1", 1.0, carrier, 1e-6, std::move(sig1));
	tx1.setSignal(&wave1);

	radar::Receiver rx1(&rx_plat1, "rx1", 42, radar::OperationMode::PULSED_MODE);
	rx1.setAntenna(&iso_ant);
	rx1.setTiming(timing);

	radar::IsoTarget tgt1(&tgt_plat1, "tgt1", 1.0, 42);

	simulation::ReResults r1{};
	simulation::solveRe(&tx1, &rx1, &tgt1, std::chrono::duration<RealType>(0.0), &wave1, r1);

	// R = 1000 m (doubled)
	radar::Platform tx_plat2("tx2");
	setupPlatform(tx_plat2, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_plat2("rx2");
	setupPlatform(rx_plat2, math::Vec3{0.0, 0.1, 0.0});
	radar::Platform tgt_plat2("tgt2");
	setupPlatform(tgt_plat2, math::Vec3{1000.0, 0.0, 0.0});

	radar::Transmitter tx2(&tx_plat2, "tx2", radar::OperationMode::PULSED_MODE);
	tx2.setAntenna(&iso_ant);
	tx2.setTiming(timing);
	auto sig2 = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave2("sig2", 1.0, carrier, 1e-6, std::move(sig2));
	tx2.setSignal(&wave2);

	radar::Receiver rx2(&rx_plat2, "rx2", 42, radar::OperationMode::PULSED_MODE);
	rx2.setAntenna(&iso_ant);
	rx2.setTiming(timing);

	radar::IsoTarget tgt2(&tgt_plat2, "tgt2", 1.0, 42);

	simulation::ReResults r2{};
	simulation::solveRe(&tx2, &rx2, &tgt2, std::chrono::duration<RealType>(0.0), &wave2, r2);

	// r1_near ~ 500, r2_near ~ 500 => power1 ~ sigma*lambda^2 / (64*pi^3 * 500^2 * 500^2)
	// r1_far ~ 1000, r2_far ~ 1000 => power2 ~ sigma*lambda^2 / (64*pi^3 * 1000^2 * 1000^2)
	// Ratio ≈ (1000^4) / (500^4) = 16
	// Note: The Rx is offset by 0.1m, so this is approximately R^-4 but not exact
	// We use a looser tolerance to account for the 0.1m offset
	REQUIRE_THAT(r1.power / r2.power, WithinRel(16.0, 0.01));
}

// =============================================================================
// calculateReflectedPathContribution: CW reflected path complex sample
// =============================================================================

TEST_CASE("calculateReflectedPathContribution amplitude matches bistatic equation with signal power",
		  "[simulation][channel_model][reflected]")
{
	ParamGuard guard;
	params::params.reset();

	const RealType c = params::c();
	const RealType carrier = 1.0e9;
	const RealType lambda = c / carrier;
	const RealType rcs = 10.0;
	const RealType power = 100.0;

	const math::Vec3 tx_pos{0.0, 0.0, 0.0};
	const math::Vec3 tgt_pos{1000.0, 0.0, 0.0};
	const math::Vec3 rx_pos{1000.0, 1000.0, 0.0};

	const RealType r1 = (tgt_pos - tx_pos).length();
	const RealType r2 = (rx_pos - tgt_pos).length();

	const RealType bistatic = (rcs * lambda * lambda) / (64.0 * PI * PI * PI * r1 * r1 * r2 * r2);
	const RealType expected_amplitude = std::sqrt(power * bistatic);

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, tx_pos);

	radar::Platform tgt_plat("tgt_plat");
	setupPlatform(tgt_plat, tgt_pos);

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, rx_pos);

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::CW_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("sig", power, carrier, 1e-3, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::CW_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	radar::IsoTarget tgt(&tgt_plat, "tgt", rcs, 42);

	const ComplexType result = simulation::calculateReflectedPathContribution(&tx, &rx, &tgt, 0.0);

	REQUIRE_THAT(std::abs(result), WithinRel(expected_amplitude, 1e-6));
}

TEST_CASE("calculateReflectedPathContribution phase matches bistatic propagation delay",
		  "[simulation][channel_model][reflected]")
{
	ParamGuard guard;
	params::params.reset();

	const RealType c = params::c();
	const RealType carrier = 1.0e9;

	const math::Vec3 tx_pos{0.0, 0.0, 0.0};
	const math::Vec3 tgt_pos{1000.0, 0.0, 0.0};
	const math::Vec3 rx_pos{2000.0, 0.0, 0.0};

	const RealType r1 = (tgt_pos - tx_pos).length();
	const RealType r2 = (rx_pos - tgt_pos).length();
	const RealType tau = (r1 + r2) / c;

	// Expected phase = -2*pi*f*tau + timing_phase (0 for default timing)
	const RealType expected_phase = -2.0 * PI * carrier * tau;

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, tx_pos);

	radar::Platform tgt_plat("tgt_plat");
	setupPlatform(tgt_plat, tgt_pos);

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, rx_pos);

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::CW_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("sig", 1.0, carrier, 1e-3, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::CW_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	radar::IsoTarget tgt(&tgt_plat, "tgt", 1.0, 42);

	const ComplexType result = simulation::calculateReflectedPathContribution(&tx, &rx, &tgt, 0.0);
	const RealType result_phase = std::arg(result);

	// Wrap expected phase to [-pi, pi] for comparison
	RealType wrapped = std::fmod(expected_phase, 2.0 * PI);
	if (wrapped > PI)
		wrapped -= 2.0 * PI;
	if (wrapped < -PI)
		wrapped += 2.0 * PI;

	REQUIRE_THAT(result_phase, WithinAbs(wrapped, 1e-6));
}

TEST_CASE("calculateReflectedPathContribution preserves stronger phase-noise cancellation for shorter delays",
		  "[simulation][channel_model][reflected]")
{
	ParamGuard guard;
	params::params.reset();
	params::setTime(0.0, 1.0);
	params::setRate(10.0);
	params::setOversampleRatio(1);
	params::setC(10.0);

	const RealType carrier = 1.0;
	const RealType sample_time = 0.5;

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{-0.5, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{0.5, 0.0, 0.0});

	radar::Platform near_tgt_plat("near_tgt");
	setupPlatform(near_tgt_plat, math::Vec3{0.0, 0.0, 0.0}); // tau = 0.1

	radar::Platform far_tgt_plat("far_tgt");
	setupPlatform(far_tgt_plat, math::Vec3{1.0, 0.0, 0.0}); // tau = 0.2

	antenna::Isotropic iso_ant("iso");
	timing::PrototypeTiming prototype("clk");
	prototype.setFrequency(carrier);
	prototype.setFreqOffset(1.0);
	auto timing = std::make_shared<timing::Timing>("clk", 42);
	timing->initializeModel(&prototype);
	const auto lookup = makeLookup(timing);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::CW_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("sig", 1.0, carrier, 1.0, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::CW_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	radar::IsoTarget near_tgt(&near_tgt_plat, "near", 1.0, 42);
	radar::IsoTarget far_tgt(&far_tgt_plat, "far", 1.0, 43);

	const ComplexType near_ideal = simulation::calculateReflectedPathContribution(&tx, &rx, &near_tgt, sample_time);
	const ComplexType near_delayed =
		simulation::calculateReflectedPathContribution(&tx, &rx, &near_tgt, sample_time, &lookup);
	const ComplexType far_ideal = simulation::calculateReflectedPathContribution(&tx, &rx, &far_tgt, sample_time);
	const ComplexType far_delayed =
		simulation::calculateReflectedPathContribution(&tx, &rx, &far_tgt, sample_time, &lookup);

	const RealType near_extra_phase = std::abs(std::arg(near_delayed * std::conj(near_ideal)));
	const RealType far_extra_phase = std::abs(std::arg(far_delayed * std::conj(far_ideal)));

	REQUIRE_THAT(near_extra_phase, WithinAbs(0.2 * PI, 1e-6));
	REQUIRE_THAT(far_extra_phase, WithinAbs(0.4 * PI, 1e-6));
	REQUIRE(near_extra_phase < far_extra_phase);
}

TEST_CASE("solveRe with 3D geometry computes correct bistatic range", "[simulation][channel_model][reflected]")
{
	ParamGuard guard;
	params::params.reset();

	const RealType c = params::c();
	const RealType carrier = 2.0e9;
	const RealType lambda = c / carrier;
	const RealType rcs = 5.0;

	// 3D positions
	const math::Vec3 tx_pos{0.0, 0.0, 0.0};
	const math::Vec3 tgt_pos{300.0, 400.0, 0.0}; // r1 = 500
	const math::Vec3 rx_pos{300.0, 400.0, 600.0}; // r2 = 600

	const RealType r1 = (tgt_pos - tx_pos).length();
	const RealType r2 = (rx_pos - tgt_pos).length();

	REQUIRE_THAT(r1, WithinAbs(500.0, 1e-9));
	REQUIRE_THAT(r2, WithinAbs(600.0, 1e-9));

	const RealType expected_power = (rcs * lambda * lambda) / (64.0 * PI * PI * PI * r1 * r1 * r2 * r2);
	const RealType expected_delay = (r1 + r2) / c;

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, tx_pos);

	radar::Platform tgt_plat("tgt_plat");
	setupPlatform(tgt_plat, tgt_pos);

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, rx_pos);

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("sig", 1.0, carrier, 1e-6, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	radar::IsoTarget tgt(&tgt_plat, "tgt", rcs, 42);

	simulation::ReResults results{};
	simulation::solveRe(&tx, &rx, &tgt, std::chrono::duration<RealType>(0.0), &wave, results);

	REQUIRE_THAT(results.power, WithinRel(expected_power, 1e-6));
	REQUIRE_THAT(results.delay, WithinRel(expected_delay, 1e-9));
}
