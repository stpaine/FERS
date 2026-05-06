#include <catch2/catch_test_macros.hpp>
#include <queue>

#include "core/sim_events.h"

TEST_CASE("EventComparator orders earliest timestamps first", "[core][events]")
{
	core::EventComparator comparator;

	const core::Event early{1.0, core::EventType::TX_PULSED_START, nullptr};
	const core::Event late{5.0, core::EventType::RX_STREAMING_END, nullptr};

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
	REQUIRE(core::toString(core::EventType::TX_STREAMING_START) == "TxStreamingStart");
	REQUIRE(core::toString(core::EventType::TX_STREAMING_END) == "TxStreamingEnd");
	REQUIRE(core::toString(core::EventType::RX_STREAMING_START) == "RxStreamingStart");
	REQUIRE(core::toString(core::EventType::RX_STREAMING_END) == "RxStreamingEnd");

	const auto unknown = static_cast<core::EventType>(-1);
	REQUIRE(core::toString(unknown) == "UnknownEvent");
}
