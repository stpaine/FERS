// Tests for channel_model.h/cpp: RangeError, ReResults, and edge cases that exercise
// internal helper functions (computeLink, wattsToDbm, wattsToDb, isComponentActive, etc.)
// indirectly through the public API.

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
#include "radar/target.h"
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

	// Helper to set up a platform at a static position
	void setupPlatform(radar::Platform& plat, const math::Vec3& pos)
	{
		plat.getMotionPath()->addCoord(math::Coord{pos, 0.0});
		plat.getMotionPath()->finalize();
		plat.getRotationPath()->addCoord(math::RotationCoord{0.0, 0.0, 0.0});
		plat.getRotationPath()->finalize();
	}
}

// =============================================================================
// RangeError exception
// =============================================================================

TEST_CASE("RangeError provides descriptive message", "[simulation][channel_model][helpers]")
{
	const simulation::RangeError error;
	REQUIRE(std::string(error.what()).find("Range error") != std::string::npos);
}

TEST_CASE("RangeError inherits from std::exception", "[simulation][channel_model][helpers]")
{
	const simulation::RangeError error;
	const std::exception& base = error;
	REQUIRE(base.what() != nullptr);
}

// =============================================================================
// ReResults struct
// =============================================================================

TEST_CASE("ReResults default-initializes to zero", "[simulation][channel_model][helpers]")
{
	simulation::ReResults results{};
	REQUIRE_THAT(results.power, WithinAbs(0.0, 0.0));
	REQUIRE_THAT(results.delay, WithinAbs(0.0, 0.0));
	REQUIRE_THAT(results.phase, WithinAbs(0.0, 0.0));
}

TEST_CASE("ReResults holds assigned values", "[simulation][channel_model][helpers]")
{
	simulation::ReResults results{.power = 1.5e-10, .delay = 3.3e-6, .phase = -1.2};
	REQUIRE_THAT(results.power, WithinAbs(1.5e-10, 1e-20));
	REQUIRE_THAT(results.delay, WithinAbs(3.3e-6, 1e-16));
	REQUIRE_THAT(results.phase, WithinAbs(-1.2, 1e-12));
}

// =============================================================================
// RangeError thrown for co-located objects (tests computeLink indirectly)
// =============================================================================

TEST_CASE("solveReDirect throws RangeError when Tx and Rx are at the same position",
		  "[simulation][channel_model][helpers]")
{
	ParamGuard guard;
	params::params.reset();

	// Two platforms at the exact same position
	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{100.0, 200.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{100.0, 200.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("sig", 1.0, 1e9, 1e-6, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	simulation::ReResults results{};
	const auto time = std::chrono::duration<RealType>(0.0);

	REQUIRE_THROWS_AS(simulation::solveReDirect(&tx, &rx, time, &wave, results), simulation::RangeError);
}

TEST_CASE("solveRe throws RangeError when target is at transmitter position", "[simulation][channel_model][helpers]")
{
	ParamGuard guard;
	params::params.reset();

	// Tx and target at same position; Rx elsewhere
	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform tgt_plat("tgt_plat");
	setupPlatform(tgt_plat, math::Vec3{0.0, 0.0, 0.0}); // same as Tx

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{1000.0, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("sig", 1.0, 1e9, 1e-6, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	radar::IsoTarget tgt(&tgt_plat, "tgt", 10.0, 42);

	simulation::ReResults results{};
	const auto time = std::chrono::duration<RealType>(0.0);

	REQUIRE_THROWS_AS(simulation::solveRe(&tx, &rx, &tgt, time, &wave, results), simulation::RangeError);
}

TEST_CASE("solveRe throws RangeError when target is at receiver position", "[simulation][channel_model][helpers]")
{
	ParamGuard guard;
	params::params.reset();

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{1000.0, 0.0, 0.0});

	radar::Platform tgt_plat("tgt_plat");
	setupPlatform(tgt_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{0.0, 0.0, 0.0}); // same as target

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::PULSED_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("sig", 1.0, 1e9, 1e-6, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::PULSED_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	radar::IsoTarget tgt(&tgt_plat, "tgt", 10.0, 42);

	simulation::ReResults results{};
	const auto time = std::chrono::duration<RealType>(0.0);

	REQUIRE_THROWS_AS(simulation::solveRe(&tx, &rx, &tgt, time, &wave, results), simulation::RangeError);
}

// =============================================================================
// Enum value tests
// =============================================================================

TEST_CASE("LinkType enum values are distinct", "[simulation][channel_model][helpers]")
{
	REQUIRE(simulation::LinkType::Monostatic != simulation::LinkType::BistaticTxTgt);
	REQUIRE(simulation::LinkType::BistaticTxTgt != simulation::LinkType::BistaticTgtRx);
	REQUIRE(simulation::LinkType::BistaticTgtRx != simulation::LinkType::DirectTxRx);
	REQUIRE(simulation::LinkType::Monostatic != simulation::LinkType::DirectTxRx);
}

TEST_CASE("LinkQuality enum values are distinct", "[simulation][channel_model][helpers]")
{
	REQUIRE(simulation::LinkQuality::Strong != simulation::LinkQuality::Weak);
}

// =============================================================================
// PreviewLink struct
// =============================================================================

TEST_CASE("PreviewLink holds assigned values", "[simulation][channel_model][helpers]")
{
	simulation::PreviewLink link{.type = simulation::LinkType::Monostatic,
								 .quality = simulation::LinkQuality::Strong,
								 .label = "Test Label",
								 .display_value = -42.0,
								 .source_id = 100,
								 .dest_id = 200,
								 .origin_id = 100};

	REQUIRE(link.type == simulation::LinkType::Monostatic);
	REQUIRE(link.quality == simulation::LinkQuality::Strong);
	REQUIRE(link.label == "Test Label");
	REQUIRE(link.display_value == -42.0);
	REQUIRE(link.source_id == 100);
	REQUIRE(link.dest_id == 200);
	REQUIRE(link.origin_id == 100);
}

// =============================================================================
// CW path returns zero for co-located Tx/Rx (same platform)
// =============================================================================

TEST_CASE("calculateDirectPathContribution returns zero for co-located Tx and Rx on same platform",
		  "[simulation][channel_model][helpers]")
{
	ParamGuard guard;
	params::params.reset();

	radar::Platform shared_plat("shared");
	setupPlatform(shared_plat, math::Vec3{500.0, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&shared_plat, "tx", radar::OperationMode::CW_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("sig", 1.0, 1e9, 1e-6, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&shared_plat, "rx", 42, radar::OperationMode::CW_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	const ComplexType result = simulation::calculateDirectPathContribution(&tx, &rx, 0.0);
	REQUIRE_THAT(std::abs(result), WithinAbs(0.0, 1e-15));
}

TEST_CASE("calculateReflectedPathContribution returns zero when target shares platform with Tx",
		  "[simulation][channel_model][helpers]")
{
	ParamGuard guard;
	params::params.reset();

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{1000.0, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::CW_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("sig", 1.0, 1e9, 1e-6, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::CW_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	// Target on same platform as Tx
	radar::IsoTarget tgt(&tx_plat, "tgt", 10.0, 42);

	const ComplexType result = simulation::calculateReflectedPathContribution(&tx, &rx, &tgt, 0.0);
	REQUIRE_THAT(std::abs(result), WithinAbs(0.0, 1e-15));
}

TEST_CASE("calculateReflectedPathContribution returns zero when target shares platform with Rx",
		  "[simulation][channel_model][helpers]")
{
	ParamGuard guard;
	params::params.reset();

	radar::Platform tx_plat("tx_plat");
	setupPlatform(tx_plat, math::Vec3{0.0, 0.0, 0.0});

	radar::Platform rx_plat("rx_plat");
	setupPlatform(rx_plat, math::Vec3{1000.0, 0.0, 0.0});

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	radar::Transmitter tx(&tx_plat, "tx", radar::OperationMode::CW_MODE);
	tx.setAntenna(&iso_ant);
	tx.setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal wave("sig", 1.0, 1e9, 1e-6, std::move(sig));
	tx.setSignal(&wave);

	radar::Receiver rx(&rx_plat, "rx", 42, radar::OperationMode::CW_MODE);
	rx.setAntenna(&iso_ant);
	rx.setTiming(timing);

	// Target on same platform as Rx
	radar::IsoTarget tgt(&rx_plat, "tgt", 10.0, 42);

	const ComplexType result = simulation::calculateReflectedPathContribution(&tx, &rx, &tgt, 0.0);
	REQUIRE_THAT(std::abs(result), WithinAbs(0.0, 1e-15));
}
