// Tests for channel_model::calculatePreviewLinks. Verifies link generation for
// monostatic, bistatic, and direct-path scenarios, including schedule filtering,
// link type categorization, and signal strength classification.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/parameters.h"
#include "core/world.h"
#include "math/coord.h"
#include "math/geometry_ops.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/schedule_period.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "signal/radar_signal.h"
#include "simulation/channel_model.h"
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

	void setupPlatform(radar::Platform& plat, const math::Vec3& pos)
	{
		plat.getMotionPath()->addCoord(math::Coord{pos, 0.0});
		plat.getMotionPath()->finalize();
		plat.getRotationPath()->addCoord(math::RotationCoord{0.0, 0.0, 0.0});
		plat.getRotationPath()->finalize();
	}

	// Count links of a specific type in a vector
	std::size_t countLinkType(const std::vector<simulation::PreviewLink>& links, simulation::LinkType type)
	{
		std::size_t count = 0;
		for (const auto& link : links)
		{
			if (link.type == type)
			{
				++count;
			}
		}
		return count;
	}

	void requireSingleValuePreviewLabel(const std::string& label)
	{
		REQUIRE(label.find("avg ") == std::string::npos);
		REQUIRE(label.find("peak ") == std::string::npos);
		REQUIRE(label.find('(') == std::string::npos);
		REQUIRE(label.find(')') == std::string::npos);
	}
}

// =============================================================================
// Empty world produces no links
// =============================================================================

TEST_CASE("calculatePreviewLinks returns empty for empty world", "[simulation][channel_model][preview]")
{
	ParamGuard guard;
	params::params.reset();

	core::World world;
	const auto links = simulation::calculatePreviewLinks(world, 0.0);
	REQUIRE(links.empty());
}

// =============================================================================
// Monostatic scenario: Tx with attached Rx and one target
// =============================================================================

TEST_CASE("calculatePreviewLinks monostatic produces Monostatic and BistaticTxTgt links",
		  "[simulation][channel_model][preview]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1e6);

	// Setup: Tx/Rx co-located at origin, target at 1000 m
	auto tx_plat = std::make_unique<radar::Platform>("tx_plat");
	setupPlatform(*tx_plat, math::Vec3{0.0, 0.0, 0.0});
	auto* tx_plat_ptr = tx_plat.get();

	auto tgt_plat = std::make_unique<radar::Platform>("tgt_plat");
	setupPlatform(*tgt_plat, math::Vec3{1000.0, 0.0, 0.0});
	auto* tgt_plat_ptr = tgt_plat.get();

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	auto tx = std::make_unique<radar::Transmitter>(tx_plat_ptr, "tx", radar::OperationMode::PULSED_MODE, 100);
	tx->setAntenna(&iso_ant);
	tx->setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	auto wave = std::make_unique<fers_signal::RadarSignal>("sig", 1000.0, 1e9, 1e-3, std::move(sig), 300);
	tx->setSignal(wave.get());

	auto rx = std::make_unique<radar::Receiver>(tx_plat_ptr, "rx", 42, radar::OperationMode::PULSED_MODE, 200);
	rx->setAntenna(&iso_ant);
	rx->setTiming(timing);
	rx->setNoiseTemperature(290.0);

	tx->setAttached(rx.get());

	auto tgt = radar::createIsoTarget(tgt_plat_ptr, "tgt", 10.0, 42, 400);

	core::World world;
	world.add(std::move(tx_plat));
	world.add(std::move(tgt_plat));
	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt));
	world.add(std::move(wave));

	const auto links = simulation::calculatePreviewLinks(world, 0.0);

	// Should have exactly 1 BistaticTxTgt (Tx -> Target illuminator link)
	// and 1 Monostatic (Tx/Rx -> Target round trip)
	REQUIRE(countLinkType(links, simulation::LinkType::BistaticTxTgt) == 1);
	REQUIRE(countLinkType(links, simulation::LinkType::Monostatic) == 1);

	// No direct Tx->Rx or BistaticTgtRx links in monostatic
	REQUIRE(countLinkType(links, simulation::LinkType::DirectTxRx) == 0);
	REQUIRE(countLinkType(links, simulation::LinkType::BistaticTgtRx) == 0);
}

TEST_CASE("calculatePreviewLinks monostatic link references correct IDs", "[simulation][channel_model][preview]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1e6);

	const SimId tx_id = 100;
	const SimId rx_id = 200;
	const SimId tgt_id = 400;

	auto tx_plat = std::make_unique<radar::Platform>("tx_plat");
	setupPlatform(*tx_plat, math::Vec3{0.0, 0.0, 0.0});
	auto* tx_plat_ptr = tx_plat.get();

	auto tgt_plat = std::make_unique<radar::Platform>("tgt_plat");
	setupPlatform(*tgt_plat, math::Vec3{1000.0, 0.0, 0.0});
	auto* tgt_plat_ptr = tgt_plat.get();

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	auto tx = std::make_unique<radar::Transmitter>(tx_plat_ptr, "tx", radar::OperationMode::PULSED_MODE, tx_id);
	tx->setAntenna(&iso_ant);
	tx->setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	auto wave = std::make_unique<fers_signal::RadarSignal>("sig", 1000.0, 1e9, 1e-3, std::move(sig), 300);
	tx->setSignal(wave.get());

	auto rx = std::make_unique<radar::Receiver>(tx_plat_ptr, "rx", 42, radar::OperationMode::PULSED_MODE, rx_id);
	rx->setAntenna(&iso_ant);
	rx->setTiming(timing);
	rx->setNoiseTemperature(290.0);

	tx->setAttached(rx.get());

	auto tgt = radar::createIsoTarget(tgt_plat_ptr, "tgt", 10.0, 42, tgt_id);

	core::World world;
	world.add(std::move(tx_plat));
	world.add(std::move(tgt_plat));
	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt));
	world.add(std::move(wave));

	const auto links = simulation::calculatePreviewLinks(world, 0.0);

	// Find the Monostatic link
	for (const auto& link : links)
	{
		if (link.type == simulation::LinkType::Monostatic)
		{
			REQUIRE(link.source_id == tx_id);
			REQUIRE(link.dest_id == tgt_id);
			REQUIRE(link.origin_id == tx_id);
		}
		if (link.type == simulation::LinkType::BistaticTxTgt)
		{
			REQUIRE(link.source_id == tx_id);
			REQUIRE(link.dest_id == tgt_id);
			REQUIRE(link.origin_id == tx_id);
		}
	}
}

// =============================================================================
// Bistatic scenario: Separate Tx and Rx, one target
// =============================================================================

TEST_CASE("calculatePreviewLinks bistatic produces correct link types", "[simulation][channel_model][preview]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1e6);

	auto tx_plat = std::make_unique<radar::Platform>("tx_plat");
	setupPlatform(*tx_plat, math::Vec3{0.0, 0.0, 0.0});
	auto* tx_plat_ptr = tx_plat.get();

	auto rx_plat = std::make_unique<radar::Platform>("rx_plat");
	setupPlatform(*rx_plat, math::Vec3{2000.0, 0.0, 0.0});
	auto* rx_plat_ptr = rx_plat.get();

	auto tgt_plat = std::make_unique<radar::Platform>("tgt_plat");
	setupPlatform(*tgt_plat, math::Vec3{1000.0, 1000.0, 0.0});
	auto* tgt_plat_ptr = tgt_plat.get();

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	auto tx = std::make_unique<radar::Transmitter>(tx_plat_ptr, "tx", radar::OperationMode::PULSED_MODE, 100);
	tx->setAntenna(&iso_ant);
	tx->setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	auto wave = std::make_unique<fers_signal::RadarSignal>("sig", 1000.0, 1e9, 1e-3, std::move(sig), 300);
	tx->setSignal(wave.get());

	auto rx = std::make_unique<radar::Receiver>(rx_plat_ptr, "rx", 42, radar::OperationMode::PULSED_MODE, 200);
	rx->setAntenna(&iso_ant);
	rx->setTiming(timing);
	rx->setNoiseTemperature(290.0);
	// NOT attached - this is bistatic

	auto tgt = radar::createIsoTarget(tgt_plat_ptr, "tgt", 10.0, 42, 400);

	core::World world;
	world.add(std::move(tx_plat));
	world.add(std::move(rx_plat));
	world.add(std::move(tgt_plat));
	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt));
	world.add(std::move(wave));

	const auto links = simulation::calculatePreviewLinks(world, 0.0);

	// Bistatic scenario should have:
	// 1 BistaticTxTgt (Tx -> Target, in outer loop)
	// 1 DirectTxRx (Tx -> Rx direct path)
	// 1 BistaticTgtRx (Target -> Rx scattered)
	// 0 Monostatic (not attached)
	REQUIRE(countLinkType(links, simulation::LinkType::BistaticTxTgt) == 1);
	REQUIRE(countLinkType(links, simulation::LinkType::DirectTxRx) == 1);
	REQUIRE(countLinkType(links, simulation::LinkType::BistaticTgtRx) == 1);
	REQUIRE(countLinkType(links, simulation::LinkType::Monostatic) == 0);
}

// =============================================================================
// NODIRECT flag suppresses direct Tx->Rx link
// =============================================================================

TEST_CASE("calculatePreviewLinks bistatic with NODIRECT flag omits direct link", "[simulation][channel_model][preview]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1e6);

	auto tx_plat = std::make_unique<radar::Platform>("tx_plat");
	setupPlatform(*tx_plat, math::Vec3{0.0, 0.0, 0.0});
	auto* tx_plat_ptr = tx_plat.get();

	auto rx_plat = std::make_unique<radar::Platform>("rx_plat");
	setupPlatform(*rx_plat, math::Vec3{2000.0, 0.0, 0.0});
	auto* rx_plat_ptr = rx_plat.get();

	auto tgt_plat = std::make_unique<radar::Platform>("tgt_plat");
	setupPlatform(*tgt_plat, math::Vec3{1000.0, 1000.0, 0.0});
	auto* tgt_plat_ptr = tgt_plat.get();

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	auto tx = std::make_unique<radar::Transmitter>(tx_plat_ptr, "tx", radar::OperationMode::PULSED_MODE, 100);
	tx->setAntenna(&iso_ant);
	tx->setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	auto wave = std::make_unique<fers_signal::RadarSignal>("sig", 1000.0, 1e9, 1e-3, std::move(sig), 300);
	tx->setSignal(wave.get());

	auto rx = std::make_unique<radar::Receiver>(rx_plat_ptr, "rx", 42, radar::OperationMode::PULSED_MODE, 200);
	rx->setAntenna(&iso_ant);
	rx->setTiming(timing);
	rx->setNoiseTemperature(290.0);
	rx->setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);

	auto tgt = radar::createIsoTarget(tgt_plat_ptr, "tgt", 10.0, 42, 400);

	core::World world;
	world.add(std::move(tx_plat));
	world.add(std::move(rx_plat));
	world.add(std::move(tgt_plat));
	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt));
	world.add(std::move(wave));

	const auto links = simulation::calculatePreviewLinks(world, 0.0);

	// DirectTxRx should be suppressed
	REQUIRE(countLinkType(links, simulation::LinkType::DirectTxRx) == 0);
	// Other links should still be present
	REQUIRE(countLinkType(links, simulation::LinkType::BistaticTxTgt) == 1);
	REQUIRE(countLinkType(links, simulation::LinkType::BistaticTgtRx) == 1);
}

// =============================================================================
// Schedule filtering
// =============================================================================

TEST_CASE("calculatePreviewLinks respects transmitter schedule", "[simulation][channel_model][preview]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1e6);

	auto tx_plat = std::make_unique<radar::Platform>("tx_plat");
	setupPlatform(*tx_plat, math::Vec3{0.0, 0.0, 0.0});
	auto* tx_plat_ptr = tx_plat.get();

	auto rx_plat = std::make_unique<radar::Platform>("rx_plat");
	setupPlatform(*rx_plat, math::Vec3{1000.0, 0.0, 0.0});
	auto* rx_plat_ptr = rx_plat.get();

	auto tgt_plat = std::make_unique<radar::Platform>("tgt_plat");
	setupPlatform(*tgt_plat, math::Vec3{500.0, 0.0, 0.0});
	auto* tgt_plat_ptr = tgt_plat.get();

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	auto tx = std::make_unique<radar::Transmitter>(tx_plat_ptr, "tx", radar::OperationMode::PULSED_MODE, 100);
	tx->setAntenna(&iso_ant);
	tx->setTiming(timing);
	// Schedule: active from 1.0 to 5.0 seconds only
	tx->setSchedule({radar::SchedulePeriod{1.0, 5.0}});

	auto sig = std::make_unique<fers_signal::CwSignal>();
	auto wave = std::make_unique<fers_signal::RadarSignal>("sig", 1000.0, 1e9, 1e-3, std::move(sig), 300);
	tx->setSignal(wave.get());

	auto rx = std::make_unique<radar::Receiver>(rx_plat_ptr, "rx", 42, radar::OperationMode::PULSED_MODE, 200);
	rx->setAntenna(&iso_ant);
	rx->setTiming(timing);
	rx->setNoiseTemperature(290.0);

	auto tgt = radar::createIsoTarget(tgt_plat_ptr, "tgt", 10.0, 42, 400);

	core::World world;
	world.add(std::move(tx_plat));
	world.add(std::move(rx_plat));
	world.add(std::move(tgt_plat));
	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt));
	world.add(std::move(wave));

	// At time 0.0, Tx is not active
	auto links_inactive = simulation::calculatePreviewLinks(world, 0.0);
	REQUIRE(links_inactive.empty());

	// At time 2.0, Tx is active
	auto links_active = simulation::calculatePreviewLinks(world, 2.0);
	REQUIRE_FALSE(links_active.empty());
}

TEST_CASE("calculatePreviewLinks respects receiver schedule", "[simulation][channel_model][preview]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1e6);

	auto tx_plat = std::make_unique<radar::Platform>("tx_plat");
	setupPlatform(*tx_plat, math::Vec3{0.0, 0.0, 0.0});
	auto* tx_plat_ptr = tx_plat.get();

	auto rx_plat = std::make_unique<radar::Platform>("rx_plat");
	setupPlatform(*rx_plat, math::Vec3{1000.0, 0.0, 0.0});
	auto* rx_plat_ptr = rx_plat.get();

	auto tgt_plat = std::make_unique<radar::Platform>("tgt_plat");
	setupPlatform(*tgt_plat, math::Vec3{500.0, 0.0, 0.0});
	auto* tgt_plat_ptr = tgt_plat.get();

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	auto tx = std::make_unique<radar::Transmitter>(tx_plat_ptr, "tx", radar::OperationMode::PULSED_MODE, 100);
	tx->setAntenna(&iso_ant);
	tx->setTiming(timing);
	// Tx always active (no schedule)

	auto sig = std::make_unique<fers_signal::CwSignal>();
	auto wave = std::make_unique<fers_signal::RadarSignal>("sig", 1000.0, 1e9, 1e-3, std::move(sig), 300);
	tx->setSignal(wave.get());

	auto rx = std::make_unique<radar::Receiver>(rx_plat_ptr, "rx", 42, radar::OperationMode::PULSED_MODE, 200);
	rx->setAntenna(&iso_ant);
	rx->setTiming(timing);
	rx->setNoiseTemperature(290.0);
	// Receiver active from 3.0 to 7.0 only
	rx->setSchedule({radar::SchedulePeriod{3.0, 7.0}});

	auto tgt = radar::createIsoTarget(tgt_plat_ptr, "tgt", 10.0, 42, 400);

	core::World world;
	world.add(std::move(tx_plat));
	world.add(std::move(rx_plat));
	world.add(std::move(tgt_plat));
	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt));
	world.add(std::move(wave));

	// At time 2.0: Tx active, Rx not active
	// The illuminator path (Tx->Tgt) is in the outer Tx loop, so it should still appear
	auto links_rx_inactive = simulation::calculatePreviewLinks(world, 2.0);
	// We should have the illuminator link (BistaticTxTgt), but no Rx-dependent links
	REQUIRE(countLinkType(links_rx_inactive, simulation::LinkType::BistaticTxTgt) == 1);
	REQUIRE(countLinkType(links_rx_inactive, simulation::LinkType::DirectTxRx) == 0);
	REQUIRE(countLinkType(links_rx_inactive, simulation::LinkType::BistaticTgtRx) == 0);

	// At time 5.0: both active
	auto links_active = simulation::calculatePreviewLinks(world, 5.0);
	REQUIRE(countLinkType(links_active, simulation::LinkType::BistaticTxTgt) == 1);
	REQUIRE(countLinkType(links_active, simulation::LinkType::DirectTxRx) == 1);
	REQUIRE(countLinkType(links_active, simulation::LinkType::BistaticTgtRx) == 1);
}

// =============================================================================
// Multiple targets produce multiple links
// =============================================================================

TEST_CASE("calculatePreviewLinks scales link count with number of targets", "[simulation][channel_model][preview]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1e6);

	auto tx_plat = std::make_unique<radar::Platform>("tx_plat");
	setupPlatform(*tx_plat, math::Vec3{0.0, 0.0, 0.0});
	auto* tx_plat_ptr = tx_plat.get();

	auto tgt_plat1 = std::make_unique<radar::Platform>("tgt_plat1");
	setupPlatform(*tgt_plat1, math::Vec3{1000.0, 0.0, 0.0});
	auto* tgt_plat1_ptr = tgt_plat1.get();

	auto tgt_plat2 = std::make_unique<radar::Platform>("tgt_plat2");
	setupPlatform(*tgt_plat2, math::Vec3{0.0, 1000.0, 0.0});
	auto* tgt_plat2_ptr = tgt_plat2.get();

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	auto tx = std::make_unique<radar::Transmitter>(tx_plat_ptr, "tx", radar::OperationMode::PULSED_MODE, 100);
	tx->setAntenna(&iso_ant);
	tx->setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	auto wave = std::make_unique<fers_signal::RadarSignal>("sig", 1000.0, 1e9, 1e-3, std::move(sig), 300);
	tx->setSignal(wave.get());

	auto rx = std::make_unique<radar::Receiver>(tx_plat_ptr, "rx", 42, radar::OperationMode::PULSED_MODE, 200);
	rx->setAntenna(&iso_ant);
	rx->setTiming(timing);
	rx->setNoiseTemperature(290.0);
	tx->setAttached(rx.get());

	auto tgt1 = radar::createIsoTarget(tgt_plat1_ptr, "tgt1", 10.0, 42, 401);
	auto tgt2 = radar::createIsoTarget(tgt_plat2_ptr, "tgt2", 5.0, 42, 402);

	core::World world;
	world.add(std::move(tx_plat));
	world.add(std::move(tgt_plat1));
	world.add(std::move(tgt_plat2));
	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt1));
	world.add(std::move(tgt2));
	world.add(std::move(wave));

	const auto links = simulation::calculatePreviewLinks(world, 0.0);

	// 2 targets: 2 BistaticTxTgt + 2 Monostatic
	REQUIRE(countLinkType(links, simulation::LinkType::BistaticTxTgt) == 2);
	REQUIRE(countLinkType(links, simulation::LinkType::Monostatic) == 2);
}

// =============================================================================
// Signal strength classification (Strong vs Weak)
// =============================================================================

TEST_CASE("calculatePreviewLinks classifies very weak signals as Weak", "[simulation][channel_model][preview]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1e6); // 1 MHz bandwidth for noise floor calculation

	// Place target very far away with low power and small RCS
	// to ensure received power is below noise floor

	auto tx_plat = std::make_unique<radar::Platform>("tx_plat");
	setupPlatform(*tx_plat, math::Vec3{0.0, 0.0, 0.0});
	auto* tx_plat_ptr = tx_plat.get();

	auto tgt_plat = std::make_unique<radar::Platform>("tgt_plat");
	setupPlatform(*tgt_plat, math::Vec3{1e7, 0.0, 0.0}); // 10000 km!
	auto* tgt_plat_ptr = tgt_plat.get();

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	auto tx = std::make_unique<radar::Transmitter>(tx_plat_ptr, "tx", radar::OperationMode::PULSED_MODE, 100);
	tx->setAntenna(&iso_ant);
	tx->setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	// Very low transmit power: 1e-6 W (1 uW)
	auto wave = std::make_unique<fers_signal::RadarSignal>("sig", 1e-6, 1e9, 1e-3, std::move(sig), 300);
	tx->setSignal(wave.get());

	auto rx = std::make_unique<radar::Receiver>(tx_plat_ptr, "rx", 42, radar::OperationMode::PULSED_MODE, 200);
	rx->setAntenna(&iso_ant);
	rx->setTiming(timing);
	rx->setNoiseTemperature(290.0); // Standard noise temperature
	tx->setAttached(rx.get());

	// Very small RCS: 0.001 m^2
	auto tgt = radar::createIsoTarget(tgt_plat_ptr, "tgt", 0.001, 42, 400);

	core::World world;
	world.add(std::move(tx_plat));
	world.add(std::move(tgt_plat));
	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt));
	world.add(std::move(wave));

	const auto links = simulation::calculatePreviewLinks(world, 0.0);

	// Find the Monostatic link and check it's Weak
	bool found_mono = false;
	for (const auto& link : links)
	{
		if (link.type == simulation::LinkType::Monostatic)
		{
			REQUIRE(link.quality == simulation::LinkQuality::Weak);
			found_mono = true;
		}
	}
	REQUIRE(found_mono);
}

TEST_CASE("calculatePreviewLinks classifies strong signals as Strong", "[simulation][channel_model][preview]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1e6);

	// Close target, high power, large RCS => strong signal

	auto tx_plat = std::make_unique<radar::Platform>("tx_plat");
	setupPlatform(*tx_plat, math::Vec3{0.0, 0.0, 0.0});
	auto* tx_plat_ptr = tx_plat.get();

	auto tgt_plat = std::make_unique<radar::Platform>("tgt_plat");
	setupPlatform(*tgt_plat, math::Vec3{100.0, 0.0, 0.0}); // 100 m
	auto* tgt_plat_ptr = tgt_plat.get();

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	auto tx = std::make_unique<radar::Transmitter>(tx_plat_ptr, "tx", radar::OperationMode::PULSED_MODE, 100);
	tx->setAntenna(&iso_ant);
	tx->setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	auto wave = std::make_unique<fers_signal::RadarSignal>("sig", 1e6, 1e9, 1e-3, std::move(sig), 300);
	tx->setSignal(wave.get());

	auto rx = std::make_unique<radar::Receiver>(tx_plat_ptr, "rx", 42, radar::OperationMode::PULSED_MODE, 200);
	rx->setAntenna(&iso_ant);
	rx->setTiming(timing);
	rx->setNoiseTemperature(290.0);
	tx->setAttached(rx.get());

	auto tgt = radar::createIsoTarget(tgt_plat_ptr, "tgt", 100.0, 42, 400);

	core::World world;
	world.add(std::move(tx_plat));
	world.add(std::move(tgt_plat));
	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt));
	world.add(std::move(wave));

	const auto links = simulation::calculatePreviewLinks(world, 0.0);

	// Find the Monostatic link and check it's Strong
	bool found_mono = false;
	for (const auto& link : links)
	{
		if (link.type == simulation::LinkType::Monostatic)
		{
			REQUIRE(link.quality == simulation::LinkQuality::Strong);
			found_mono = true;
		}
	}
	REQUIRE(found_mono);
}

// =============================================================================
// No waveform attached: uses default lambda and power=0
// =============================================================================

TEST_CASE("calculatePreviewLinks handles transmitter with no waveform", "[simulation][channel_model][preview]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1e6);

	auto tx_plat = std::make_unique<radar::Platform>("tx_plat");
	setupPlatform(*tx_plat, math::Vec3{0.0, 0.0, 0.0});
	auto* tx_plat_ptr = tx_plat.get();

	auto tgt_plat = std::make_unique<radar::Platform>("tgt_plat");
	setupPlatform(*tgt_plat, math::Vec3{1000.0, 0.0, 0.0});
	auto* tgt_plat_ptr = tgt_plat.get();

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	auto tx = std::make_unique<radar::Transmitter>(tx_plat_ptr, "tx", radar::OperationMode::PULSED_MODE, 100);
	tx->setAntenna(&iso_ant);
	tx->setTiming(timing);
	// NOT setting a signal: tx->setSignal(...)

	auto rx = std::make_unique<radar::Receiver>(tx_plat_ptr, "rx", 42, radar::OperationMode::PULSED_MODE, 200);
	rx->setAntenna(&iso_ant);
	rx->setTiming(timing);
	rx->setNoiseTemperature(290.0);
	tx->setAttached(rx.get());

	auto tgt = radar::createIsoTarget(tgt_plat_ptr, "tgt", 10.0, 42, 400);

	core::World world;
	world.add(std::move(tx_plat));
	world.add(std::move(tgt_plat));
	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt));
	// No waveform added — testing the no-signal case

	// Should not crash; uses default lambda = 0.3m, power = 0
	const auto links = simulation::calculatePreviewLinks(world, 0.0);

	// Links should still be generated (geometric visualization)
	REQUIRE(countLinkType(links, simulation::LinkType::BistaticTxTgt) == 1);
}

// =============================================================================
// Label format verification (exercises wattsToDbm/wattsToDb indirectly)
// =============================================================================

TEST_CASE("calculatePreviewLinks monostatic label contains dBm", "[simulation][channel_model][preview]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1e6);

	auto tx_plat = std::make_unique<radar::Platform>("tx_plat");
	setupPlatform(*tx_plat, math::Vec3{0.0, 0.0, 0.0});
	auto* tx_plat_ptr = tx_plat.get();

	auto tgt_plat = std::make_unique<radar::Platform>("tgt_plat");
	setupPlatform(*tgt_plat, math::Vec3{1000.0, 0.0, 0.0});
	auto* tgt_plat_ptr = tgt_plat.get();

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	auto tx = std::make_unique<radar::Transmitter>(tx_plat_ptr, "tx", radar::OperationMode::PULSED_MODE, 100);
	tx->setAntenna(&iso_ant);
	tx->setTiming(timing);

	auto sig = std::make_unique<fers_signal::CwSignal>();
	auto wave = std::make_unique<fers_signal::RadarSignal>("sig", 1000.0, 1e9, 1e-3, std::move(sig), 300);
	tx->setSignal(wave.get());

	auto rx = std::make_unique<radar::Receiver>(tx_plat_ptr, "rx", 42, radar::OperationMode::PULSED_MODE, 200);
	rx->setAntenna(&iso_ant);
	rx->setTiming(timing);
	rx->setNoiseTemperature(290.0);
	tx->setAttached(rx.get());

	auto tgt = radar::createIsoTarget(tgt_plat_ptr, "tgt", 10.0, 42, 400);

	core::World world;
	world.add(std::move(tx_plat));
	world.add(std::move(tgt_plat));
	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt));
	world.add(std::move(wave));

	const auto links = simulation::calculatePreviewLinks(world, 0.0);
	bool saw_monostatic = false;
	bool saw_tx_tgt = false;

	for (const auto& link : links)
	{
		if (link.type == simulation::LinkType::Monostatic)
		{
			saw_monostatic = true;
			REQUIRE(link.label.find("dBm") != std::string::npos);
			REQUIRE_THAT(link.rcs, WithinAbs(10.0, 1e-12));
			REQUIRE(link.actual_power_dbm > -999.0);

			const double label_power_dbm = std::stod(link.label);
			REQUIRE_THAT(link.actual_power_dbm, WithinAbs(label_power_dbm + 10.0 * std::log10(link.rcs), 0.1));
		}
		if (link.type == simulation::LinkType::BistaticTxTgt)
		{
			saw_tx_tgt = true;
			REQUIRE(link.label.find("dBW") != std::string::npos);
		}
	}

	REQUIRE(saw_monostatic);
	REQUIRE(saw_tx_tgt);
}

TEST_CASE("calculatePreviewLinks FMCW labels use single-value preview style", "[simulation][channel_model][preview]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1e6);

	auto tx_plat = std::make_unique<radar::Platform>("tx_plat");
	setupPlatform(*tx_plat, math::Vec3{0.0, 0.0, 0.0});
	auto* tx_plat_ptr = tx_plat.get();

	auto tgt_plat = std::make_unique<radar::Platform>("tgt_plat");
	setupPlatform(*tgt_plat, math::Vec3{1000.0, 0.0, 0.0});
	auto* tgt_plat_ptr = tgt_plat.get();

	antenna::Isotropic iso_ant("iso");
	auto timing = std::make_shared<timing::Timing>("clk", 42);

	auto tx = std::make_unique<radar::Transmitter>(tx_plat_ptr, "tx", radar::OperationMode::FMCW_MODE, 100);
	tx->setAntenna(&iso_ant);
	tx->setTiming(timing);

	auto sig = std::make_unique<fers_signal::FmcwChirpSignal>(1.0e6, 1.0e-4, 4.0e-4);
	auto wave = std::make_unique<fers_signal::RadarSignal>("sig", 1000.0, 1e9, 1.0e-4, std::move(sig), 300);
	tx->setSignal(wave.get());

	auto rx = std::make_unique<radar::Receiver>(tx_plat_ptr, "rx", 42, radar::OperationMode::FMCW_MODE, 200);
	rx->setAntenna(&iso_ant);
	rx->setTiming(timing);
	rx->setNoiseTemperature(290.0);
	tx->setAttached(rx.get());

	auto tgt = radar::createIsoTarget(tgt_plat_ptr, "tgt", 1.0, 42, 400);

	core::World world;
	world.add(std::move(tx_plat));
	world.add(std::move(tgt_plat));
	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt));
	world.add(std::move(wave));

	const auto links = simulation::calculatePreviewLinks(world, 0.0);
	bool saw_fmcw_power_label = false;

	for (const auto& link : links)
	{
		if (link.type == simulation::LinkType::Monostatic)
		{
			saw_fmcw_power_label = true;
			requireSingleValuePreviewLabel(link.label);
			REQUIRE(link.label.find("dBm") != std::string::npos);
			(void)std::stod(link.label);
		}
		if (link.type == simulation::LinkType::BistaticTxTgt)
		{
			requireSingleValuePreviewLabel(link.label);
			REQUIRE(link.label.find("dBW/m") != std::string::npos);
			(void)std::stod(link.label);
		}
	}

	REQUIRE(saw_fmcw_power_label);
}
