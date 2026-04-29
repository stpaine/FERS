#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <memory>
#include <string>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/parameters.h"
#include "core/world.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "signal/radar_signal.h"
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

	struct EventSnapshot
	{
		RealType timestamp{};
		core::EventType type{};
		std::string source_name;
	};

	std::vector<EventSnapshot> drainQueue(core::World& world)
	{
		auto& queue = world.getEventQueue();
		std::vector<EventSnapshot> events;
		while (!queue.empty())
		{
			const auto [timestamp, type, source] = queue.top();
			queue.pop();
			EventSnapshot snapshot;
			snapshot.timestamp = timestamp;
			snapshot.type = type;
			snapshot.source_name = source ? source->getName() : std::string{};
			events.push_back(std::move(snapshot));
		}
		return events;
	}

	std::vector<EventSnapshot> eventsFor(const std::vector<EventSnapshot>& events, const std::string& name)
	{
		std::vector<EventSnapshot> filtered;
		for (const auto& event : events)
		{
			if (event.source_name == name)
			{
				filtered.push_back(event);
			}
		}
		return filtered;
	}
}

TEST_CASE("World stores and retrieves added objects", "[core][world]")
{
	core::World world;

	auto platform = std::make_unique<radar::Platform>("Platform-A", 11);
	auto target = std::make_unique<radar::IsoTarget>(platform.get(), "Target-A", 1.5, 9, 22);
	world.add(std::move(platform));
	world.add(std::move(target));

	auto tx_platform = std::make_unique<radar::Platform>("Platform-Tx", 12);
	auto rx_platform = std::make_unique<radar::Platform>("Platform-Rx", 13);
	auto tx = std::make_unique<radar::Transmitter>(tx_platform.get(), "Tx-A", radar::OperationMode::PULSED_MODE, 101);
	auto rx = std::make_unique<radar::Receiver>(rx_platform.get(), "Rx-A", 123, radar::OperationMode::CW_MODE, 202);
	world.add(std::move(tx_platform));
	world.add(std::move(rx_platform));
	world.add(std::move(tx));
	world.add(std::move(rx));

	auto signal = std::make_unique<fers_signal::RadarSignal>("Wave-A", 2.0, 1.0e9, 0.2,
															 std::make_unique<fers_signal::CwSignal>(), 301);
	world.add(std::move(signal));

	auto antenna = std::make_unique<antenna::Isotropic>("Ant-A", 401);
	world.add(std::move(antenna));

	auto timing = std::make_unique<timing::PrototypeTiming>("Timing-A", 501);
	world.add(std::move(timing));

	REQUIRE(world.getPlatforms().size() == 3);
	REQUIRE(world.getTargets().size() == 1);
	REQUIRE(world.getTransmitters().size() == 1);
	REQUIRE(world.getReceivers().size() == 1);
	REQUIRE(world.getWaveforms().size() == 1);
	REQUIRE(world.getAntennas().size() == 1);
	REQUIRE(world.getTimings().size() == 1);

	REQUIRE(world.findWaveform(301) != nullptr);
	REQUIRE(world.findAntenna(401) != nullptr);
	REQUIRE(world.findTiming(501) != nullptr);
	REQUIRE(world.findTransmitterByName("Tx-A") == world.getTransmitters().front().get());
	REQUIRE(world.findWaveformByName("Wave-A") == world.findWaveform(301));

	REQUIRE(world.findWaveform(9999) == nullptr);
	REQUIRE(world.findAntenna(9999) == nullptr);
	REQUIRE(world.findTiming(9999) == nullptr);
	REQUIRE(world.findTransmitterByName("missing") == nullptr);
	REQUIRE(world.findWaveformByName("missing") == nullptr);
}

TEST_CASE("World finds platform by ID", "[core][world]")
{
	core::World world;
	auto platform = std::make_unique<radar::Platform>("FindMe", 777);
	auto* plat_ptr = platform.get();
	world.add(std::move(platform));

	REQUIRE(world.findPlatform(777) == plat_ptr);
	REQUIRE(world.findPlatform(999) == nullptr);
}

TEST_CASE("World finds radar components by ID", "[core][world]")
{
	core::World world;
	auto platform = std::make_unique<radar::Platform>("Plat", 1);
	auto* plat_ptr = platform.get();
	world.add(std::move(platform));

	auto tx = std::make_unique<radar::Transmitter>(plat_ptr, "Tx", radar::OperationMode::PULSED_MODE, 101);
	auto rx = std::make_unique<radar::Receiver>(plat_ptr, "Rx", 42, radar::OperationMode::PULSED_MODE, 202);
	auto tgt = std::make_unique<radar::IsoTarget>(plat_ptr, "Tgt", 1.0, 42, 303);

	auto* tx_ptr = tx.get();
	auto* rx_ptr = rx.get();
	auto* tgt_ptr = tgt.get();

	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt));

	REQUIRE(world.findTransmitter(101) == tx_ptr);
	REQUIRE(world.findReceiver(202) == rx_ptr);
	REQUIRE(world.findTarget(303) == tgt_ptr);

	REQUIRE(world.findTransmitter(999) == nullptr);
	REQUIRE(world.findReceiver(999) == nullptr);
	REQUIRE(world.findTarget(999) == nullptr);
}

TEST_CASE("World replaces target", "[core][world]")
{
	core::World world;
	auto platform = std::make_unique<radar::Platform>("Plat", 1);
	auto* plat_ptr = platform.get();
	world.add(std::move(platform));

	// 1. Add initial target
	auto old_tgt = std::make_unique<radar::IsoTarget>(plat_ptr, "OldTgt", 1.0, 42, 100);
	world.add(std::move(old_tgt));

	REQUIRE(world.getTargets().size() == 1);
	REQUIRE(world.findTarget(100)->getName() == "OldTgt");

	// 2. Replace existing target
	auto new_tgt = std::make_unique<radar::IsoTarget>(plat_ptr, "NewTgt", 5.0, 42, 100);
	auto* new_tgt_ptr = new_tgt.get();
	world.replace(std::move(new_tgt));

	REQUIRE(world.getTargets().size() == 1);
	REQUIRE(world.findTarget(100) == new_tgt_ptr);
	REQUIRE(world.findTarget(100)->getName() == "NewTgt");

	// 3. Replace (add) non-existent target
	auto another_tgt = std::make_unique<radar::IsoTarget>(plat_ptr, "AnotherTgt", 10.0, 42, 200);
	auto* another_tgt_ptr = another_tgt.get();
	world.replace(std::move(another_tgt));

	REQUIRE(world.getTargets().size() == 2);
	REQUIRE(world.findTarget(200) == another_tgt_ptr);
}

TEST_CASE("World replaces antenna and updates dependent components", "[core][world]")
{
	core::World world;

	// 1. Setup initial antenna
	auto old_ant = std::make_unique<antenna::Isotropic>("OldAnt", 100);
	auto* old_ant_ptr = old_ant.get();
	world.add(std::move(old_ant));

	// 2. Setup platform and components using the antenna
	auto plat = std::make_unique<radar::Platform>("Plat", 1);
	auto tx = std::make_unique<radar::Transmitter>(plat.get(), "Tx", radar::OperationMode::CW_MODE, 2);
	auto rx = std::make_unique<radar::Receiver>(plat.get(), "Rx", 42, radar::OperationMode::CW_MODE, 3);

	tx->setAntenna(old_ant_ptr);
	rx->setAntenna(old_ant_ptr);

	auto* tx_ptr = tx.get();
	auto* rx_ptr = rx.get();

	world.add(std::move(plat));
	world.add(std::move(tx));
	world.add(std::move(rx));

	REQUIRE(tx_ptr->getAntenna() == old_ant_ptr);
	REQUIRE(rx_ptr->getAntenna() == old_ant_ptr);

	// 3. Replace the antenna
	auto new_ant = std::make_unique<antenna::Isotropic>("NewAnt", 100); // Same ID
	auto* new_ant_ptr = new_ant.get();
	world.replace(std::move(new_ant));

	// 4. Verify replacement and pointer updates
	REQUIRE(world.findAntenna(100) == new_ant_ptr);
	REQUIRE(world.findAntenna(100)->getName() == "NewAnt");
	REQUIRE(tx_ptr->getAntenna() == new_ant_ptr);
	REQUIRE(rx_ptr->getAntenna() == new_ant_ptr);
}

TEST_CASE("World replaces waveform and updates dependent components", "[core][world]")
{
	core::World world;

	// 1. Setup initial waveform
	auto old_wf = std::make_unique<fers_signal::RadarSignal>("OldWf", 1.0, 1e9, 1.0,
															 std::make_unique<fers_signal::CwSignal>(), 200);
	auto* old_wf_ptr = old_wf.get();
	world.add(std::move(old_wf));

	// 2. Setup platform and transmitter using the waveform
	auto plat = std::make_unique<radar::Platform>("Plat", 1);
	auto tx = std::make_unique<radar::Transmitter>(plat.get(), "Tx", radar::OperationMode::CW_MODE, 2);

	tx->setSignal(old_wf_ptr);
	auto* tx_ptr = tx.get();

	world.add(std::move(plat));
	world.add(std::move(tx));

	REQUIRE(tx_ptr->getSignal() == old_wf_ptr);

	// 3. Replace the waveform
	auto new_wf = std::make_unique<fers_signal::RadarSignal>("NewWf", 2.0, 2e9, 1.0,
															 std::make_unique<fers_signal::CwSignal>(), 200); // Same ID
	auto* new_wf_ptr = new_wf.get();
	world.replace(std::move(new_wf));

	// 4. Verify replacement and pointer updates
	REQUIRE(world.findWaveform(200) == new_wf_ptr);
	REQUIRE(world.findWaveform(200)->getName() == "NewWf");
	REQUIRE(world.findWaveformByName("OldWf") == nullptr);
	REQUIRE(world.findWaveformByName("NewWf") == new_wf_ptr);
	REQUIRE(tx_ptr->getSignal() == new_wf_ptr);
}

TEST_CASE("World replaces timing and refreshes dependent radar timing models", "[core][world]")
{
	core::World world;

	auto old_timing = std::make_unique<timing::PrototypeTiming>("OldTiming", 300);
	old_timing->setFrequency(1.0e6);
	old_timing->setFreqOffset(2.0);
	old_timing->setPhaseOffset(0.1);
	world.add(std::move(old_timing));

	auto plat = std::make_unique<radar::Platform>("Plat", 1);
	auto tx = std::make_unique<radar::Transmitter>(plat.get(), "Tx", radar::OperationMode::CW_MODE, 2);
	auto rx = std::make_unique<radar::Receiver>(plat.get(), "Rx", 42, radar::OperationMode::CW_MODE, 3);

	auto tx_timing = std::make_shared<timing::Timing>("OldTiming", 12345, 300);
	tx_timing->initializeModel(world.findTiming(300));
	auto rx_timing = std::make_shared<timing::Timing>("OldTiming", 54321, 300);
	rx_timing->initializeModel(world.findTiming(300));

	tx->setTiming(tx_timing);
	rx->setTiming(rx_timing);

	auto* tx_ptr = tx.get();
	auto* rx_ptr = rx.get();

	world.add(std::move(plat));
	world.add(std::move(tx));
	world.add(std::move(rx));

	auto new_timing = std::make_unique<timing::PrototypeTiming>("NewTiming", 300);
	new_timing->setFrequency(2.5e6);
	new_timing->setSyncOnPulse();
	new_timing->setFreqOffset(7.5);
	new_timing->setPhaseOffset(0.75);
	new_timing->setAlpha(1.0, 0.5);
	world.replace(std::move(new_timing));

	auto* replaced = world.findTiming(300);
	REQUIRE(replaced != nullptr);
	REQUIRE(replaced->getName() == "NewTiming");
	REQUIRE_THAT(replaced->getFrequency(), WithinAbs(2.5e6, 1e-9));

	REQUIRE(tx_ptr->getTiming().get() != tx_timing.get());
	REQUIRE(rx_ptr->getTiming().get() != rx_timing.get());
	REQUIRE(tx_ptr->getTiming().get() != rx_ptr->getTiming().get());
	REQUIRE(tx_ptr->getTiming()->getSeed() == 12345);
	REQUIRE(rx_ptr->getTiming()->getSeed() == 54321);
	REQUIRE(tx_ptr->getTiming()->getName() == "NewTiming");
	REQUIRE(rx_ptr->getTiming()->getName() == "NewTiming");
	REQUIRE_THAT(tx_ptr->getTiming()->getFrequency(), WithinAbs(2.5e6, 1e-9));
	REQUIRE_THAT(rx_ptr->getTiming()->getFrequency(), WithinAbs(2.5e6, 1e-9));
	REQUIRE(tx_ptr->getTiming()->getSyncOnPulse());
	REQUIRE(rx_ptr->getTiming()->getSyncOnPulse());
	REQUIRE_THAT(tx_ptr->getTiming()->getFreqOffset(), WithinAbs(7.5, 1e-9));
	REQUIRE_THAT(rx_ptr->getTiming()->getFreqOffset(), WithinAbs(7.5, 1e-9));
	REQUIRE_THAT(tx_ptr->getTiming()->getPhaseOffset(), WithinAbs(0.75, 1e-9));
	REQUIRE_THAT(rx_ptr->getTiming()->getPhaseOffset(), WithinAbs(0.75, 1e-9));
}

TEST_CASE("World replacement preserves shared timing instances", "[core][world]")
{
	core::World world;

	auto old_timing = std::make_unique<timing::PrototypeTiming>("OldTiming", 301);
	old_timing->setFrequency(1.0e6);
	world.add(std::move(old_timing));

	auto plat = std::make_unique<radar::Platform>("Plat", 11);
	auto tx = std::make_unique<radar::Transmitter>(plat.get(), "Tx", radar::OperationMode::CW_MODE, 12);
	auto rx = std::make_unique<radar::Receiver>(plat.get(), "Rx", 42, radar::OperationMode::CW_MODE, 13);

	auto shared_timing = std::make_shared<timing::Timing>("OldTiming", 777, 301);
	shared_timing->initializeModel(world.findTiming(301));
	tx->setTiming(shared_timing);
	rx->setTiming(shared_timing);

	auto* tx_ptr = tx.get();
	auto* rx_ptr = rx.get();

	world.add(std::move(plat));
	world.add(std::move(tx));
	world.add(std::move(rx));

	auto new_timing = std::make_unique<timing::PrototypeTiming>("NewTiming", 301);
	new_timing->setFrequency(2.0e6);
	world.replace(std::move(new_timing));

	REQUIRE(tx_ptr->getTiming().get() == rx_ptr->getTiming().get());
	REQUIRE(tx_ptr->getTiming().get() != shared_timing.get());
	REQUIRE(tx_ptr->getTiming()->getSeed() == 777);
	REQUIRE_THAT(tx_ptr->getTiming()->getFrequency(), WithinAbs(2.0e6, 1e-9));
}

TEST_CASE("World enforces unique ids for assets", "[core][world]")
{
	core::World world;

	SECTION("Waveform ids are unique")
	{
		world.add(std::make_unique<fers_signal::RadarSignal>("Wave-1", 1.0, 1.0e9, 0.1,
															 std::make_unique<fers_signal::CwSignal>(), 77));
		REQUIRE_THROWS_AS(world.add(std::make_unique<fers_signal::RadarSignal>(
							  "Wave-2", 2.0, 2.0e9, 0.2, std::make_unique<fers_signal::CwSignal>(), 77)),
						  std::runtime_error);
	}

	SECTION("Antenna ids are unique")
	{
		world.add(std::make_unique<antenna::Isotropic>("Ant-1", 88));
		REQUIRE_THROWS_AS(world.add(std::make_unique<antenna::Isotropic>("Ant-2", 88)), std::runtime_error);
	}

	SECTION("Timing ids are unique")
	{
		world.add(std::make_unique<timing::PrototypeTiming>("Timing-1", 99));
		REQUIRE_THROWS_AS(world.add(std::make_unique<timing::PrototypeTiming>("Timing-2", 99)), std::runtime_error);
	}
}

TEST_CASE("World clear resets storage and state", "[core][world]")
{
	core::World world;
	world.add(std::make_unique<radar::Platform>("Platform-A", 11));
	world.add(std::make_unique<radar::Platform>("Platform-B", 12));
	world.add(std::make_unique<radar::Transmitter>(world.getPlatforms().front().get(), "Tx-A",
												   radar::OperationMode::CW_MODE, 101));
	world.add(std::make_unique<radar::Receiver>(world.getPlatforms().front().get(), "Rx-A", 1,
												radar::OperationMode::CW_MODE, 202));
	world.add(std::make_unique<radar::IsoTarget>(world.getPlatforms().front().get(), "Target-A", 1.0, 5, 303));
	world.add(std::make_unique<fers_signal::RadarSignal>("Wave-A", 1.0, 1.0e9, 0.1,
														 std::make_unique<fers_signal::CwSignal>(), 404));
	world.add(std::make_unique<antenna::Isotropic>("Ant-A", 505));
	world.add(std::make_unique<timing::PrototypeTiming>("Timing-A", 606));

	world.getEventQueue().push({1.0, core::EventType::TX_STREAMING_START, world.getTransmitters().front().get()});
	world.getSimulationState().t_current = 42.0;
	core::ActiveStreamingSource source{};
	source.transmitter = world.getTransmitters().front().get();
	world.getSimulationState().active_streaming_transmitters.push_back(source);

	world.clear();

	REQUIRE(world.getPlatforms().empty());
	REQUIRE(world.getTargets().empty());
	REQUIRE(world.getTransmitters().empty());
	REQUIRE(world.getReceivers().empty());
	REQUIRE(world.getWaveforms().empty());
	REQUIRE(world.getAntennas().empty());
	REQUIRE(world.getTimings().empty());
	REQUIRE(world.getEventQueue().empty());
	REQUIRE(world.findTransmitterByName("Tx-A") == nullptr);
	REQUIRE(world.findWaveformByName("Wave-A") == nullptr);
	REQUIRE(world.getSimulationState().t_current == 0.0);
	REQUIRE(world.getSimulationState().active_streaming_transmitters.empty());
}

TEST_CASE("World schedules pulsed transmitter events", "[core][world]")
{
	ParamGuard guard;
	params::setTime(0.0, 10.0);

	core::World world;
	auto platform = std::make_unique<radar::Platform>("Platform-Tx", 1);
	auto tx = std::make_unique<radar::Transmitter>(platform.get(), "Tx-Pulsed", radar::OperationMode::PULSED_MODE, 2);
	std::vector<radar::SchedulePeriod> schedule = {{2.0, 3.0}, {6.0, 7.0}};
	tx->setSchedule(schedule);
	world.add(std::move(platform));
	world.add(std::move(tx));

	world.scheduleInitialEvents();

	const auto events = drainQueue(world);
	REQUIRE(events.size() == 1);
	REQUIRE(events[0].type == core::EventType::TX_PULSED_START);
	REQUIRE(events[0].source_name == "Tx-Pulsed");
	REQUIRE_THAT(events[0].timestamp, WithinAbs(2.0, 1e-12));
}

TEST_CASE("World schedules CW transmitter events with bounds", "[core][world]")
{
	ParamGuard guard;
	params::setTime(0.0, 10.0);

	core::World world;
	auto platform = std::make_unique<radar::Platform>("Platform-Tx", 1);
	auto tx = std::make_unique<radar::Transmitter>(platform.get(), "Tx-CW", radar::OperationMode::CW_MODE, 2);
	std::vector<radar::SchedulePeriod> schedule = {{-5.0, 1.0}, {3.0, 5.0}, {9.0, 12.0}};
	tx->setSchedule(schedule);
	world.add(std::move(platform));
	world.add(std::move(tx));

	world.scheduleInitialEvents();

	const auto events = drainQueue(world);
	REQUIRE(events.size() == 6);

	const auto tx_events = eventsFor(events, "Tx-CW");
	REQUIRE(tx_events.size() == 6);
	REQUIRE(tx_events[0].type == core::EventType::TX_STREAMING_START);
	REQUIRE(tx_events[1].type == core::EventType::TX_STREAMING_END);
	REQUIRE(tx_events[2].type == core::EventType::TX_STREAMING_START);
	REQUIRE(tx_events[3].type == core::EventType::TX_STREAMING_END);
	REQUIRE(tx_events[4].type == core::EventType::TX_STREAMING_START);
	REQUIRE(tx_events[5].type == core::EventType::TX_STREAMING_END);

	REQUIRE_THAT(tx_events[0].timestamp, WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(tx_events[1].timestamp, WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(tx_events[2].timestamp, WithinAbs(3.0, 1e-12));
	REQUIRE_THAT(tx_events[3].timestamp, WithinAbs(5.0, 1e-12));
	REQUIRE_THAT(tx_events[4].timestamp, WithinAbs(9.0, 1e-12));
	REQUIRE_THAT(tx_events[5].timestamp, WithinAbs(10.0, 1e-12));
}

TEST_CASE("World schedules CW transmitter always-on when schedule empty", "[core][world]")
{
	ParamGuard guard;
	params::setTime(1.5, 4.5);

	core::World world;
	auto platform = std::make_unique<radar::Platform>("Platform-Tx", 1);
	auto tx = std::make_unique<radar::Transmitter>(platform.get(), "Tx-CW", radar::OperationMode::CW_MODE, 2);
	world.add(std::move(platform));
	world.add(std::move(tx));

	world.scheduleInitialEvents();

	const auto events = drainQueue(world);
	REQUIRE(events.size() == 2);
	REQUIRE(events[0].type == core::EventType::TX_STREAMING_START);
	REQUIRE(events[1].type == core::EventType::TX_STREAMING_END);
	REQUIRE_THAT(events[0].timestamp, WithinAbs(1.5, 1e-12));
	REQUIRE_THAT(events[1].timestamp, WithinAbs(4.5, 1e-12));
}

TEST_CASE("World schedules pulsed receiver events", "[core][world]")
{
	ParamGuard guard;
	params::setTime(0.0, 10.0);
	params::setRate(1000.0);
	params::setOversampleRatio(1);

	core::World world;
	auto platform = std::make_unique<radar::Platform>("Platform-Rx", 1);
	auto rx = std::make_unique<radar::Receiver>(platform.get(), "Rx-Pulsed", 3, radar::OperationMode::PULSED_MODE, 4);
	rx->setWindowProperties(1.0, 2.0, 0.5);
	rx->setTiming(std::make_shared<timing::Timing>("RxClock", 11));
	std::vector<radar::SchedulePeriod> schedule = {{2.0, 3.0}, {6.0, 7.0}};
	rx->setSchedule(schedule);
	world.add(std::move(platform));
	world.add(std::move(rx));

	world.scheduleInitialEvents();

	const auto events = drainQueue(world);
	REQUIRE(events.size() == 1);
	REQUIRE(events[0].type == core::EventType::RX_PULSED_WINDOW_START);
	REQUIRE(events[0].source_name == "Rx-Pulsed");
	REQUIRE_THAT(events[0].timestamp, WithinAbs(2.0, 1e-12));
}

TEST_CASE("World schedules CW receiver events with bounds", "[core][world]")
{
	ParamGuard guard;
	params::setTime(0.0, 10.0);

	core::World world;
	auto platform = std::make_unique<radar::Platform>("Platform-Rx", 1);
	auto rx = std::make_unique<radar::Receiver>(platform.get(), "Rx-CW", 3, radar::OperationMode::CW_MODE, 4);
	std::vector<radar::SchedulePeriod> schedule = {{-1.0, 2.0}, {4.0, 6.0}, {9.0, 11.0}};
	rx->setSchedule(schedule);
	world.add(std::move(platform));
	world.add(std::move(rx));

	world.scheduleInitialEvents();

	const auto events = drainQueue(world);
	REQUIRE(events.size() == 6);

	const auto rx_events = eventsFor(events, "Rx-CW");
	REQUIRE(rx_events.size() == 6);
	REQUIRE(rx_events[0].type == core::EventType::RX_STREAMING_START);
	REQUIRE(rx_events[1].type == core::EventType::RX_STREAMING_END);
	REQUIRE(rx_events[2].type == core::EventType::RX_STREAMING_START);
	REQUIRE(rx_events[3].type == core::EventType::RX_STREAMING_END);
	REQUIRE(rx_events[4].type == core::EventType::RX_STREAMING_START);
	REQUIRE(rx_events[5].type == core::EventType::RX_STREAMING_END);

	REQUIRE_THAT(rx_events[0].timestamp, WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(rx_events[1].timestamp, WithinAbs(2.0, 1e-12));
	REQUIRE_THAT(rx_events[2].timestamp, WithinAbs(4.0, 1e-12));
	REQUIRE_THAT(rx_events[3].timestamp, WithinAbs(6.0, 1e-12));
	REQUIRE_THAT(rx_events[4].timestamp, WithinAbs(9.0, 1e-12));
	REQUIRE_THAT(rx_events[5].timestamp, WithinAbs(10.0, 1e-12));
}

TEST_CASE("World schedules CW receiver always-on when schedule empty", "[core][world]")
{
	ParamGuard guard;
	params::setTime(2.0, 5.0);

	core::World world;
	auto platform = std::make_unique<radar::Platform>("Platform-Rx", 1);
	auto rx = std::make_unique<radar::Receiver>(platform.get(), "Rx-CW", 3, radar::OperationMode::CW_MODE, 4);
	world.add(std::move(platform));
	world.add(std::move(rx));

	world.scheduleInitialEvents();

	const auto events = drainQueue(world);
	REQUIRE(events.size() == 2);
	REQUIRE(events[0].type == core::EventType::RX_STREAMING_START);
	REQUIRE(events[1].type == core::EventType::RX_STREAMING_END);
	REQUIRE_THAT(events[0].timestamp, WithinAbs(2.0, 1e-12));
	REQUIRE_THAT(events[1].timestamp, WithinAbs(5.0, 1e-12));
}

TEST_CASE("World dumpEventQueue formats content", "[core][world]")
{
	core::World world;
	REQUIRE(world.dumpEventQueue() == "Event Queue is empty.\n");

	auto platform = std::make_unique<radar::Platform>("Platform-Tx", 1);
	auto tx = std::make_unique<radar::Transmitter>(platform.get(), "Tx-Format", radar::OperationMode::CW_MODE, 2);
	world.add(std::move(platform));
	world.add(std::move(tx));

	world.getEventQueue().push({1.25, core::EventType::TX_STREAMING_START, world.getTransmitters().front().get()});

	const auto output = world.dumpEventQueue();
	REQUIRE(output.find("Event Queue Contents") != std::string::npos);
	REQUIRE(output.find("TxStreamingStart") != std::string::npos);
	REQUIRE(output.find("Tx-Format") != std::string::npos);
	REQUIRE(output.find("1.250000") != std::string::npos);
}
