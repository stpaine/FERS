#include <catch2/catch_test_macros.hpp>
#include <queue>

#include "core/sim_events.h"

TEST_CASE("EventComparator orders earliest timestamps first", "[core][events]")
{
	core::EventComparator comparator;

	const core::Event early{1.0, core::EventType::TX_PULSED_START, nullptr};
	const core::Event late{5.0, core::EventType::RX_CW_END, nullptr};

	REQUIRE(comparator(late, early));
	REQUIRE_FALSE(comparator(early, late));

	std::priority_queue<core::Event, std::vector<core::Event>, core::EventComparator> queue;
	queue.push(late);
	queue.push(early);

	REQUIRE(queue.top().timestamp == early.timestamp);
}

TEST_CASE("EventType toString covers all values", "[core][events]")
{
	REQUIRE(core::toString(core::EventType::TX_PULSED_START) == "TxPulsedStart");
	REQUIRE(core::toString(core::EventType::RX_PULSED_WINDOW_START) == "RxPulsedWindowStart");
	REQUIRE(core::toString(core::EventType::RX_PULSED_WINDOW_END) == "RxPulsedWindowEnd");
	REQUIRE(core::toString(core::EventType::TX_CW_START) == "TxCwStart");
	REQUIRE(core::toString(core::EventType::TX_CW_END) == "TxCwEnd");
	REQUIRE(core::toString(core::EventType::RX_CW_START) == "RxCwStart");
	REQUIRE(core::toString(core::EventType::RX_CW_END) == "RxCwEnd");

	const auto unknown = static_cast<core::EventType>(-1);
	REQUIRE(core::toString(unknown) == "UnknownEvent");
}
