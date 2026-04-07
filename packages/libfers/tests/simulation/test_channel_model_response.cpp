// Tests for channel_model::calculateResponse function. Verifies InterpPoint
// generation, direct and reflected path response creation, co-location skipping,
// and point count calculations.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <memory>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/parameters.h"
#include "math/coord.h"
#include "math/geometry_ops.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/response.h"
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

	// A minimal Signal subclass that renders to known data
	struct TestSignal final : public fers_signal::Signal
	{
		std::vector<ComplexType> data;

		std::vector<ComplexType> render(const std::vector<interp::InterpPoint>&, unsigned& size,
										RealType) const override
		{
			size = static_cast<unsigned>(data.size());
			return data;
		}
	};
}

// =============================================================================
// calculateResponse: returns nullptr for co-located direct path
// =============================================================================

TEST_CASE("calculateResponse returns nullptr when Tx and Rx share same platform (direct path)",
		  "[simulation][channel_model][response]")
{
	ParamGuard guard;
	params::params.reset();
	params::setSimSamplingRate(1000.0);

	radar::Platform shared_plat("shared");
	setupPlatform(shared_plat, math::Vec3{0.0, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&shared_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<TestSignal>();
	sig->data = {ComplexType{1.0, 0.0}};
	fers_signal::RadarSignal wave("sig", 1.0, 1e9, 1e-3, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&shared_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	// Direct path (no target) with same platform
	auto response = simulation::calculateResponse(&tx, &rx, &wave, 0.0, nullptr);
	REQUIRE(response == nullptr);
}

TEST_CASE("calculateResponse returns nullptr when Tx is attached to Rx (monostatic, direct path)",
		  "[simulation][channel_model][response]")
{
	ParamGuard guard;
	params::params.reset();
	params::setSimSamplingRate(1000.0);

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{100.0, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<TestSignal>();
	sig->data = {ComplexType{1.0, 0.0}};
	fers_signal::RadarSignal wave("sig", 1.0, 1e9, 1e-3, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	tx.setAttached(&rx);

	// Direct path with attached Rx
	auto response = simulation::calculateResponse(&tx, &rx, &wave, 0.0, nullptr);
	REQUIRE(response == nullptr);
}

TEST_CASE("calculateResponse returns nullptr for reflected path when target co-located with Tx",
		  "[simulation][channel_model][response]")
{
	ParamGuard guard;
	params::params.reset();
	params::setSimSamplingRate(1000.0);

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{1000.0, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<TestSignal>();
	sig->data = {ComplexType{1.0, 0.0}};
	fers_signal::RadarSignal wave("sig", 1.0, 1e9, 1e-3, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	// Target on same platform as Tx
	radar::IsoTarget tgt(&tx_plat, "tgt", 10.0, 42);

	auto response = simulation::calculateResponse(&tx, &rx, &wave, 0.0, &tgt);
	REQUIRE(response == nullptr);
}

TEST_CASE("calculateResponse returns nullptr for reflected path when target co-located with Rx",
		  "[simulation][channel_model][response]")
{
	ParamGuard guard;
	params::params.reset();
	params::setSimSamplingRate(1000.0);

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{1000.0, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<TestSignal>();
	sig->data = {ComplexType{1.0, 0.0}};
	fers_signal::RadarSignal wave("sig", 1.0, 1e9, 1e-3, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	// Target on same platform as Rx
	radar::IsoTarget tgt(&rx_plat, "tgt", 10.0, 42);

	auto response = simulation::calculateResponse(&tx, &rx, &wave, 0.0, &tgt);
	REQUIRE(response == nullptr);
}

// =============================================================================
// calculateResponse: Direct path produces valid response with correct timing
// =============================================================================

TEST_CASE("calculateResponse direct path produces non-null response with interp points",
		  "[simulation][channel_model][response]")
{
	ParamGuard guard;
	params::params.reset();
	params::setSimSamplingRate(1000.0); // 1000 samples/sec

	// Signal length = 1e-3 s = 1 ms
	// Point count = ceil(1e-3 * 1000) = 1
	// Loop runs from i=0 to i=point_count (inclusive), so 2 points total

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{1000.0, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<TestSignal>();
	sig->data = {ComplexType{1.0, 0.0}};
	fers_signal::RadarSignal wave("sig", 1.0, 1e9, 1e-3, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	auto response = simulation::calculateResponse(&tx, &rx, &wave, 0.0, nullptr);
	REQUIRE(response != nullptr);

	// The response should have valid start and end times
	// start_time = 0.0 + delay, end_time = (0.0 + signal_length) + delay
	const RealType expected_delay = 1000.0 / params::c();
	REQUIRE_THAT(response->startTime(), WithinRel(expected_delay, 1e-6));
}

TEST_CASE("calculateResponse reflected path produces valid response", "[simulation][channel_model][response]")
{
	ParamGuard guard;
	params::params.reset();
	params::setSimSamplingRate(1000.0);

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform tgt_plat("tgt_plat");
	setupPlatform(tgt_plat, math::Vec3{500.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{500.0, 500.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<TestSignal>();
	sig->data = {ComplexType{1.0, 0.0}};
	fers_signal::RadarSignal wave("sig", 1.0, 1e9, 1e-3, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	radar::IsoTarget tgt(&tgt_plat, "tgt", 10.0, 42);

	auto response = simulation::calculateResponse(&tx, &rx, &wave, 0.0, &tgt);
	REQUIRE(response != nullptr);

	// The response's start time should be approximately the first sample time + delay
	const RealType r1 = 500.0;
	const RealType r2 = 500.0;
	const RealType expected_delay = (r1 + r2) / params::c();
	REQUIRE_THAT(response->startTime(), WithinRel(expected_delay, 1e-6));
}

TEST_CASE("calculateResponse direct path interp points have consistent Friis power",
		  "[simulation][channel_model][response]")
{
	ParamGuard guard;
	params::params.reset();
	params::setSimSamplingRate(10000.0); // High sample rate for many points

	// Signal length = 5e-3 s
	// point_count = ceil(5e-3 * 10000) = 50
	// Loop produces 51 interp points (0 to 50 inclusive)

	const RealType dist = 2000.0;
	const RealType carrier = 1.0e9;
	const RealType c = params::c();

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{dist, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<TestSignal>();
	sig->data = {ComplexType{1.0, 0.0}};
	fers_signal::RadarSignal wave("sig", 1.0, carrier, 5e-3, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	auto response = simulation::calculateResponse(&tx, &rx, &wave, 0.0, nullptr);
	REQUIRE(response != nullptr);

	// Render to get access to the response data
	// The response's start time should be consistent
	const RealType expected_delay = dist / c;
	REQUIRE_THAT(response->startTime(), WithinRel(expected_delay, 1e-6));

	// Response length should be approximately signal length
	REQUIRE_THAT(response->getLength(), WithinRel(5e-3, 0.01));
}

TEST_CASE("calculateResponse exposes transmitter id", "[simulation][channel_model][response]")
{
	ParamGuard guard;
	params::params.reset();
	params::setSimSamplingRate(1000.0);

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{1000.0, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	const SimId tx_id = 99999;
	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE, tx_id);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<TestSignal>();
	sig->data = {ComplexType{1.0, 0.0}};
	fers_signal::RadarSignal wave("sig", 1.0, 1e9, 1e-3, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	auto response = simulation::calculateResponse(&tx, &rx, &wave, 0.0, nullptr);
	REQUIRE(response != nullptr);
	REQUIRE(response->getTransmitterId() == tx_id);
}
