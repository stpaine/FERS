// Tests for channel_model direct-path calculations: solveReDirect and
// calculateDirectPathContribution. Verifies Friis transmission equation,
// propagation delay, phase shift, and timing offsets against hand-calculated values.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <cmath>
#include <memory>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/parameters.h"
#include "math/coord.h"
#include "math/geometry_ops.h"
#include "radar/platform.h"
#include "radar/receiver.h"
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

	RealType rms(const std::vector<ComplexType>& samples)
	{
		RealType sum_square = 0.0;
		for (const auto& sample : samples)
		{
			sum_square += std::norm(sample);
		}
		return std::sqrt(sum_square / static_cast<RealType>(samples.size()));
	}
}

// =============================================================================
// solveReDirect: Friis transmission equation verification
// =============================================================================

TEST_CASE("solveReDirect computes correct Friis power for isotropic antennas", "[simulation][channel_model][direct]")
{
	ParamGuard guard;
	params::params.reset();

	// Setup: Tx at origin, Rx at (1000, 0, 0)
	// Carrier frequency: 1 GHz => lambda = c / 1e9
	// Distance: 1000 m
	// Isotropic antennas: Gt = Gr = 1
	// Expected Friis: Pr/Pt = Gt * Gr * lambda^2 / (16 * pi^2 * R^2)

	const RealType c = params::c();
	const RealType carrier = 1.0e9;
	const RealType lambda = c / carrier;
	const RealType dist = 1000.0;

	// Hand-calculated:
	// Pr/Pt = 1 * 1 * lambda^2 / (16 * pi^2 * 1000^2)
	const RealType expected_power = (lambda * lambda) / (16.0 * PI * PI * dist * dist);
	const RealType expected_delay = dist / c;
	const RealType expected_phase = -expected_delay * 2.0 * PI * carrier;

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{dist, 0.0, 0.0});

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

	simulation::ReResults results{};
	const auto time = std::chrono::duration<RealType>(0.0);

	simulation::solveReDirect(&tx, &rx, time, &wave, results);

	REQUIRE_THAT(results.power, WithinRel(expected_power, 1e-9));
	REQUIRE_THAT(results.delay, WithinRel(expected_delay, 1e-9));
	REQUIRE_THAT(results.phase, WithinRel(expected_phase, 1e-9));
}

TEST_CASE("solveReDirect Friis equation with different distances", "[simulation][channel_model][direct]")
{
	ParamGuard guard;
	params::params.reset();

	// Verify inverse-square law: doubling distance should quarter the power
	const RealType carrier = 3.0e9; // 3 GHz

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	auto sig1 = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave1("sig1", 1.0, carrier, 1e-6, std::move(sig1));

	auto sig2 = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave2("sig2", 1.0, carrier, 1e-6, std::move(sig2));

	// Distance 500 m
	radar::Platform tx_plat1("tx1");
	setupPlatform(tx_plat1, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_plat1("rx1");
	setupPlatform(rx_plat1, math::Vec3{500.0, 0.0, 0.0});

	radar::Transmitter tx1(&tx_plat1, "tx1", radar::OperationMode::PULSED_MODE);
	tx1.setAntenna(&iso_ant);
	tx1.setTiming(timing);
	tx1.setSignal(&wave1);

	radar::Receiver rx1(&rx_plat1, "rx1", 42, radar::OperationMode::PULSED_MODE);
	rx1.setAntenna(&iso_ant);
	rx1.setTiming(timing);

	simulation::ReResults results1{};
	simulation::solveReDirect(&tx1, &rx1, std::chrono::duration<RealType>(0.0), &wave1, results1);

	// Distance 1000 m (doubled)
	radar::Platform tx_plat2("tx2");
	setupPlatform(tx_plat2, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_plat2("rx2");
	setupPlatform(rx_plat2, math::Vec3{1000.0, 0.0, 0.0});

	radar::Transmitter tx2(&tx_plat2, "tx2", radar::OperationMode::PULSED_MODE);
	tx2.setAntenna(&iso_ant);
	tx2.setTiming(timing);
	tx2.setSignal(&wave2);

	radar::Receiver rx2(&rx_plat2, "rx2", 42, radar::OperationMode::PULSED_MODE);
	rx2.setAntenna(&iso_ant);
	rx2.setTiming(timing);

	simulation::ReResults results2{};
	simulation::solveReDirect(&tx2, &rx2, std::chrono::duration<RealType>(0.0), &wave2, results2);

	// Power ratio should be 4:1 (inverse square)
	REQUIRE_THAT(results1.power / results2.power, WithinRel(4.0, 1e-9));

	// Delay ratio should be 1:2
	REQUIRE_THAT(results2.delay / results1.delay, WithinRel(2.0, 1e-9));
}

TEST_CASE("solveReDirect with noproploss flag ignores distance in power", "[simulation][channel_model][direct]")
{
	ParamGuard guard;
	params::params.reset();

	const RealType c = params::c();
	const RealType carrier = 1.0e9;
	const RealType lambda = c / carrier;
	const RealType dist = 5000.0;

	// With no_prop_loss: Pr/Pt = Gt * Gr * lambda^2 / (16 * pi^2) (no R^2 in denominator)
	const RealType expected_power = (lambda * lambda) / (16.0 * PI * PI);

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{dist, 0.0, 0.0});

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

	simulation::ReResults results{};
	const auto time = std::chrono::duration<RealType>(0.0);

	simulation::solveReDirect(&tx, &rx, time, &wave, results);

	REQUIRE_THAT(results.power, WithinRel(expected_power, 1e-9));

	// Delay is still computed based on actual distance
	const RealType expected_delay = dist / c;
	REQUIRE_THAT(results.delay, WithinRel(expected_delay, 1e-9));
}

TEST_CASE("solveReDirect phase is proportional to carrier frequency", "[simulation][channel_model][direct]")
{
	ParamGuard guard;
	params::params.reset();

	// Phase = -delay * 2 * pi * carrier
	// For same distance, doubling carrier should double the phase magnitude
	const RealType dist = 300.0;

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	// Carrier at 1 GHz
	auto sig1 = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave1("sig1", 1.0, 1.0e9, 1e-6, std::move(sig1));

	// Carrier at 2 GHz
	auto sig2 = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave2("sig2", 1.0, 2.0e9, 1e-6, std::move(sig2));

	radar::Platform tx_plat1("tx1");
	setupPlatform(tx_plat1, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_plat1("rx1");
	setupPlatform(rx_plat1, math::Vec3{dist, 0.0, 0.0});

	radar::Transmitter tx1(&tx_plat1, "tx1", radar::OperationMode::PULSED_MODE);
	tx1.setAntenna(&iso_ant);
	tx1.setTiming(timing);
	tx1.setSignal(&wave1);

	radar::Receiver rx1(&rx_plat1, "rx1", 42, radar::OperationMode::PULSED_MODE);
	rx1.setAntenna(&iso_ant);
	rx1.setTiming(timing);

	simulation::ReResults r1{};
	simulation::solveReDirect(&tx1, &rx1, std::chrono::duration<RealType>(0.0), &wave1, r1);

	radar::Platform tx_plat2("tx2");
	setupPlatform(tx_plat2, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_plat2("rx2");
	setupPlatform(rx_plat2, math::Vec3{dist, 0.0, 0.0});

	radar::Transmitter tx2(&tx_plat2, "tx2", radar::OperationMode::PULSED_MODE);
	tx2.setAntenna(&iso_ant);
	tx2.setTiming(timing);
	tx2.setSignal(&wave2);

	radar::Receiver rx2(&rx_plat2, "rx2", 42, radar::OperationMode::PULSED_MODE);
	rx2.setAntenna(&iso_ant);
	rx2.setTiming(timing);

	simulation::ReResults r2{};
	simulation::solveReDirect(&tx2, &rx2, std::chrono::duration<RealType>(0.0), &wave2, r2);

	// Same delay (same distance)
	REQUIRE_THAT(r1.delay, WithinRel(r2.delay, 1e-12));

	// Phase ratio should be 1:2
	REQUIRE_THAT(r2.phase / r1.phase, WithinRel(2.0, 1e-9));
}

// =============================================================================
// calculateDirectPathContribution: CW direct path complex sample
// =============================================================================

TEST_CASE("calculateDirectPathContribution amplitude matches Friis equation with signal power",
		  "[simulation][channel_model][direct]")
{
	ParamGuard guard;
	params::params.reset();

	// Setup: Tx at (0,0,0), Rx at (1000,0,0)
	// Carrier: 1 GHz, Signal Power: 100 W
	// Expected amplitude = sqrt(Pt * Friis_factor)
	// Friis: Gt * Gr * lambda^2 / (16*pi^2*R^2)
	const RealType c = params::c();
	const RealType carrier = 1.0e9;
	const RealType lambda = c / carrier;
	const RealType dist = 1000.0;
	const RealType power = 100.0;

	const RealType friis = (lambda * lambda) / (16.0 * PI * PI * dist * dist);
	const RealType expected_amplitude = std::sqrt(power * friis);

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{dist, 0.0, 0.0});

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

	const ComplexType result = simulation::calculateDirectPathContribution(&tx, &rx, 0.0);

	REQUIRE_THAT(std::abs(result), WithinRel(expected_amplitude, 1e-6));
}

TEST_CASE("calculateDirectPathContribution phase matches propagation delay", "[simulation][channel_model][direct]")
{
	ParamGuard guard;
	params::params.reset();

	const RealType c = params::c();
	const RealType carrier = 1.0e9;
	const RealType dist = 1000.0;

	const RealType tau = dist / c;
	// Expected phase = -2*pi*f*tau (from carrier) + timing_phase (0 for default timing)
	const RealType expected_phase = -2.0 * PI * carrier * tau;

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{dist, 0.0, 0.0});

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

	const ComplexType result = simulation::calculateDirectPathContribution(&tx, &rx, 0.0);
	const RealType result_phase = std::arg(result);

	// Compare via complex unit vectors to avoid branch-cut wrapping issues
	// If the phases match, cos(result_phase) == cos(expected_phase) and
	// sin(result_phase) == sin(expected_phase)
	REQUIRE_THAT(std::cos(result_phase), WithinAbs(std::cos(expected_phase), 1e-6));
	REQUIRE_THAT(std::sin(result_phase), WithinAbs(std::sin(expected_phase), 1e-6));
}

TEST_CASE("FMCW streaming direct path preserves in-flight segment-end tail",
		  "[simulation][channel_model][direct][fmcw]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(80.0e6);
	params::setSimSamplingRate(80.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 0.00502);

	const RealType dist = 1500.0;
	const RealType tau = dist / params::c();
	const RealType segment_end = 0.005;
	const RealType sample_rate = 80.0e6;
	const RealType dt = 1.0 / sample_rate;
	const RealType chirp_bandwidth = 20.0e6;
	const RealType chirp_duration = 250.0e-6;
	const RealType chirp_period = chirp_duration;
	const RealType chirp_rate = chirp_bandwidth / chirp_duration;
	const std::size_t sample_count = static_cast<std::size_t>(std::ceil(tau / dt));

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{dist, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::FMCW_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::FmcwChirpSignal>(chirp_bandwidth, chirp_duration, chirp_period, 0.0, 20);
	fers_signal::RadarSignal wave("fmcw", 1000.0, 10.0e9, chirp_duration, std::move(sig));
	tx.setSignal(&wave);
	const auto* fmcw = wave.getFmcwChirpSignal();

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::CW_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	const core::ActiveStreamingSource source{.transmitter = &tx, .segment_start = 0.0, .segment_end = segment_end};
	core::FmcwChirpBoundaryTracker reference_tracker;
	core::FmcwChirpBoundaryTracker tail_tracker;
	std::vector<ComplexType> reference_samples;
	std::vector<ComplexType> tail_samples;
	reference_samples.reserve(sample_count);
	tail_samples.reserve(sample_count);

	for (std::size_t i = 0; i < sample_count; ++i)
	{
		const RealType offset = static_cast<RealType>(i) * dt;
		reference_samples.push_back(simulation::calculateStreamingDirectPathContribution(
			source, &rx, segment_end - tau + offset, nullptr, &reference_tracker));
		tail_samples.push_back(simulation::calculateStreamingDirectPathContribution(source, &rx, segment_end + offset,
																					nullptr, &tail_tracker));
	}

	const RealType reference_rms = rms(reference_samples);
	const RealType tail_rms = rms(tail_samples);
	REQUIRE(reference_samples.size() >= 400);
	REQUIRE(reference_rms > 0.0);
	REQUIRE(tail_rms > 0.0);
	REQUIRE_THAT(tail_rms, WithinRel(reference_rms, 1.0e-12));

	const RealType last_chirp_start = segment_end - chirp_period;
	RealType unwrapped_span = 0.0;
	RealType previous_phase = 0.0;
	for (std::size_t i = 0; i < tail_samples.size(); ++i)
	{
		const RealType t = segment_end + static_cast<RealType>(i) * dt;
		const RealType local_time = t - last_chirp_start;
		const ComplexType dechirped = tail_samples[i] * std::polar(1.0, -fmcw->basebandPhaseForChirpTime(local_time));
		const RealType phase = std::arg(dechirped);
		if (i > 0)
		{
			unwrapped_span += unwrapDelta(phase - previous_phase);
		}
		previous_phase = phase;
	}

	const RealType measured_beat_hz = unwrapped_span / (2.0 * PI * dt * static_cast<RealType>(tail_samples.size() - 1));
	const RealType expected_beat_hz = -chirp_rate * tau;
	REQUIRE_THAT(measured_beat_hz, WithinRel(expected_beat_hz, 1.0e-6));
}

TEST_CASE("CW streaming direct path gates schedules by retarded transmit time", "[simulation][channel_model][direct]")
{
	ParamGuard guard;
	params::params.reset();
	params::setC(100.0);

	const RealType dist = 10.0;
	const RealType tau = dist / params::c();
	const RealType segment_start = 0.2;
	const RealType segment_end = 0.5;
	const RealType eps = 1.0e-6;

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{dist, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::CW_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("cw", 1.0, 1.0e9, 1.0, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::CW_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	const core::ActiveStreamingSource source{
		.transmitter = &tx, .segment_start = segment_start, .segment_end = segment_end};

	REQUIRE(std::abs(simulation::calculateStreamingDirectPathContribution(source, &rx, segment_start + tau - eps)) ==
			0.0);
	REQUIRE(std::abs(simulation::calculateStreamingDirectPathContribution(source, &rx, segment_start + tau + eps)) >
			0.0);
	REQUIRE(std::abs(simulation::calculateStreamingDirectPathContribution(source, &rx, segment_end + 0.5 * tau)) > 0.0);
	REQUIRE(std::abs(simulation::calculateStreamingDirectPathContribution(source, &rx, segment_end + tau + eps)) ==
			0.0);
}

TEST_CASE("calculateDirectPathContribution with noproploss gives distance-independent amplitude",
		  "[simulation][channel_model][direct]")
{
	ParamGuard guard;
	params::params.reset();

	const RealType c = params::c();
	const RealType carrier = 1.0e9;
	const RealType lambda = c / carrier;
	const RealType power = 50.0;

	// With noproploss: Friis = lambda^2 / (16*pi^2) (no R^2)
	const RealType friis_noloss = (lambda * lambda) / (16.0 * PI * PI);
	const RealType expected_amplitude = std::sqrt(power * friis_noloss);

	// Test with two different distances - both should give same amplitude
	for (const RealType dist : {500.0, 2000.0})
	{
		radar::Platform tx_plat("tx_plat");
		setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

		radar::Platform rx_plat("rx_plat");
		setupPlatform(rx_plat, math::Vec3{dist, 0.0, 0.0});

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
		rx.setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);

		const ComplexType result = simulation::calculateDirectPathContribution(&tx, &rx, 0.0);

		REQUIRE_THAT(std::abs(result), WithinRel(expected_amplitude, 1e-6));
	}
}

TEST_CASE("calculateDirectPathContribution applies buffered delayed timing phase with interpolation",
		  "[simulation][channel_model][direct]")
{
	ParamGuard guard;
	params::params.reset();
	params::setTime(0.0, 1.0);
	params::setRate(10.0);
	params::setOversampleRatio(1);
	params::setC(10.0);

	const RealType carrier = 1.0;
	const RealType dist = 1.5; // tau = 0.15 s
	const RealType time = 0.5;
	const RealType tau = dist / params::c();
	const RealType expected_extra_phase = -2.0 * PI * tau;

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{dist, 0.0, 0.0});

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

	const ComplexType ideal = simulation::calculateDirectPathContribution(&tx, &rx, time);
	const ComplexType delayed = simulation::calculateDirectPathContribution(&tx, &rx, time, &lookup);
	const RealType extra_phase = std::arg(delayed * std::conj(ideal));

	REQUIRE_THAT(std::cos(extra_phase), WithinAbs(std::cos(expected_extra_phase), 1e-6));
	REQUIRE_THAT(std::sin(extra_phase), WithinAbs(std::sin(expected_extra_phase), 1e-6));
}

TEST_CASE("solveReDirect at 3D separation produces correct distance and delay", "[simulation][channel_model][direct]")
{
	ParamGuard guard;
	params::params.reset();

	// Tx at (100, 200, 300), Rx at (400, 600, 800)
	// Distance = sqrt(300^2 + 400^2 + 500^2) = sqrt(90000 + 160000 + 250000) = sqrt(500000)
	const RealType c = params::c();
	const RealType carrier = 2.0e9;
	const RealType lambda = c / carrier;

	const math::Vec3 tx_pos{100.0, 200.0, 300.0};
	const math::Vec3 rx_pos{400.0, 600.0, 800.0};
	const RealType dist = std::sqrt(300.0 * 300.0 + 400.0 * 400.0 + 500.0 * 500.0);

	const RealType expected_delay = dist / c;
	const RealType expected_power = (lambda * lambda) / (16.0 * PI * PI * dist * dist);
	const RealType expected_phase = -expected_delay * 2.0 * PI * carrier;

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, tx_pos);

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

	simulation::ReResults results{};
	simulation::solveReDirect(&tx, &rx, std::chrono::duration<RealType>(0.0), &wave, results);

	REQUIRE_THAT(results.delay, WithinRel(expected_delay, 1e-9));
	REQUIRE_THAT(results.power, WithinRel(expected_power, 1e-9));
	REQUIRE_THAT(results.phase, WithinRel(expected_phase, 1e-9));
}

TEST_CASE("solveReDirect power scales with lambda squared", "[simulation][channel_model][direct]")
{
	ParamGuard guard;
	params::params.reset();

	// Friis: Pr/Pt ∝ lambda^2
	// At frequency f1, lambda1 = c/f1
	// At frequency f2 = 2*f1, lambda2 = c/(2*f1) = lambda1/2
	// Power ratio: (lambda1/lambda2)^2 = 4
	const RealType dist = 1000.0;

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	// Test at 1 GHz
	auto sig1 = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave1("sig1", 1.0, 1.0e9, 1e-6, std::move(sig1));

	radar::Platform tx_plat1("tx1");
	setupPlatform(tx_plat1, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_plat1("rx1");
	setupPlatform(rx_plat1, math::Vec3{dist, 0.0, 0.0});

	radar::Transmitter tx1(&tx_plat1, "tx1", radar::OperationMode::PULSED_MODE);
	tx1.setAntenna(&iso_ant);
	tx1.setTiming(timing);
	tx1.setSignal(&wave1);

	radar::Receiver rx1(&rx_plat1, "rx1", 42, radar::OperationMode::PULSED_MODE);
	rx1.setAntenna(&iso_ant);
	rx1.setTiming(timing);

	simulation::ReResults r1{};
	simulation::solveReDirect(&tx1, &rx1, std::chrono::duration<RealType>(0.0), &wave1, r1);

	// Test at 2 GHz
	auto sig2 = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave2("sig2", 1.0, 2.0e9, 1e-6, std::move(sig2));

	radar::Platform tx_plat2("tx2");
	setupPlatform(tx_plat2, math::Vec3{0.0, 0.0, 0.0});
	radar::Platform rx_plat2("rx2");
	setupPlatform(rx_plat2, math::Vec3{dist, 0.0, 0.0});

	radar::Transmitter tx2(&tx_plat2, "tx2", radar::OperationMode::PULSED_MODE);
	tx2.setAntenna(&iso_ant);
	tx2.setTiming(timing);
	tx2.setSignal(&wave2);

	radar::Receiver rx2(&rx_plat2, "rx2", 42, radar::OperationMode::PULSED_MODE);
	rx2.setAntenna(&iso_ant);
	rx2.setTiming(timing);

	simulation::ReResults r2{};
	simulation::solveReDirect(&tx2, &rx2, std::chrono::duration<RealType>(0.0), &wave2, r2);

	// Power at 1 GHz should be 4x power at 2 GHz
	REQUIRE_THAT(r1.power / r2.power, WithinRel(4.0, 1e-9));
}
