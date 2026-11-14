#include <catch2/catch_test_macros.hpp>
#include <queue>

#include "core/sim_events.h"

TEST_CASE("EventComparator orders earliest timestamps first", "[core][events]")
{
	core::EventComparator const comparator;

	const core::Event early{1.0, core::EventType::TX_PULSED_START, nullptr};
	const core::Event late{5.0, core::EventType::RX_STREAMING_END, nullptr};

	REQUIRE(comparator(late, early));
	REQUIRE_FALSE(comparator(early, late));

	std::priority_queue<core::Event, std::vector<core::Event>, core::EventComparator> queue;
	queue.push(late);
	queue.push(early);

	REQUIRE(queue.top().timestamp == early.timestamp);
}

TEST_CASE("EventComparator has no event-type tie-break for equal timestamps", "[core][events][equal-time]")
{
	core::EventComparator const comparator;

	const core::Event transmit{2.0, core::EventType::TX_STREAMING_START, nullptr};
	const core::Event receive{2.0, core::EventType::RX_STREAMING_END, nullptr};

	// The current comparator is timestamp-only. Distinct event types at the same
	// timestamp are therefore equivalent under the priority-queue comparator; the
	// EventType declaration order is not a secondary priority rule.
	REQUIRE_FALSE(comparator(transmit, receive));
	REQUIRE_FALSE(comparator(receive, transmit));
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

	const auto unknown = static_cast<core::EventType>(255);
	REQUIRE(core::toString(unknown) == "UnknownEvent");
}
