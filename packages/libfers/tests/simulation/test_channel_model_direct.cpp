// Tests for channel_model direct-path calculations: solveReDirect and
// calculateDirectPathContribution. Verifies Friis transmission equation,
// propagation delay, phase shift, and timing offsets against hand-calculated values.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <cmath>
#include <memory>

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
	const RealType c = params::c();
	const RealType carrier = 3.0e9; // 3 GHz
	const RealType lambda = c / carrier;

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
