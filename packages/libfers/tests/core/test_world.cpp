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

	REQUIRE(world.findWaveform(9999) == nullptr);
	REQUIRE(world.findAntenna(9999) == nullptr);
	REQUIRE(world.findTiming(9999) == nullptr);
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

	world.getEventQueue().push({1.0, core::EventType::TX_CW_START, world.getTransmitters().front().get()});
	world.getSimulationState().t_current = 42.0;
	world.getSimulationState().active_cw_transmitters.push_back(world.getTransmitters().front().get());

	world.clear();

	REQUIRE(world.getPlatforms().empty());
	REQUIRE(world.getTargets().empty());
	REQUIRE(world.getTransmitters().empty());
	REQUIRE(world.getReceivers().empty());
	REQUIRE(world.getWaveforms().empty());
	REQUIRE(world.getAntennas().empty());
	REQUIRE(world.getTimings().empty());
	REQUIRE(world.getEventQueue().empty());
	REQUIRE(world.getSimulationState().t_current == 0.0);
	REQUIRE(world.getSimulationState().active_cw_transmitters.empty());
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
	REQUIRE(tx_events[0].type == core::EventType::TX_CW_START);
	REQUIRE(tx_events[1].type == core::EventType::TX_CW_END);
	REQUIRE(tx_events[2].type == core::EventType::TX_CW_START);
	REQUIRE(tx_events[3].type == core::EventType::TX_CW_END);
	REQUIRE(tx_events[4].type == core::EventType::TX_CW_START);
	REQUIRE(tx_events[5].type == core::EventType::TX_CW_END);

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
	REQUIRE(events[0].type == core::EventType::TX_CW_START);
	REQUIRE(events[1].type == core::EventType::TX_CW_END);
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
	REQUIRE(rx_events[0].type == core::EventType::RX_CW_START);
	REQUIRE(rx_events[1].type == core::EventType::RX_CW_END);
	REQUIRE(rx_events[2].type == core::EventType::RX_CW_START);
	REQUIRE(rx_events[3].type == core::EventType::RX_CW_END);
	REQUIRE(rx_events[4].type == core::EventType::RX_CW_START);
	REQUIRE(rx_events[5].type == core::EventType::RX_CW_END);

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
	REQUIRE(events[0].type == core::EventType::RX_CW_START);
	REQUIRE(events[1].type == core::EventType::RX_CW_END);
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

	world.getEventQueue().push({1.25, core::EventType::TX_CW_START, world.getTransmitters().front().get()});

	const auto output = world.dumpEventQueue();
	REQUIRE(output.find("Event Queue Contents") != std::string::npos);
	REQUIRE(output.find("TxCwStart") != std::string::npos);
	REQUIRE(output.find("Tx-Format") != std::string::npos);
	REQUIRE(output.find("1.250000") != std::string::npos);
}
