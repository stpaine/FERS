// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file sim_events.h
 * @brief Defines the core structures for the event-driven simulation engine.
 */

#pragma once

#include "config.h"

namespace radar
{
	class Radar;
}

namespace core
{
	/**
	 * @enum EventType
	 * @brief Enumerates the types of events that can occur in the simulation.
	 */
	enum class EventType
	{
		TX_PULSED_START, ///< A pulsed transmitter begins emitting a pulse.
		RX_PULSED_WINDOW_START, ///< A pulsed receiver opens its listening window.
		RX_PULSED_WINDOW_END, ///< A pulsed receiver closes its listening window.
		TX_STREAMING_START, ///< A streaming transmitter starts transmitting.
		TX_STREAMING_END, ///< A streaming transmitter stops transmitting.
		RX_STREAMING_START, ///< A streaming receiver starts listening.
		RX_STREAMING_END, ///< A streaming receiver stops listening.
	};

	/**
	 * @struct Event
	 * @brief Represents a single event in the simulation's time-ordered queue.
	 */
	struct Event
	{
		RealType timestamp; ///< The simulation time at which the event occurs.
		EventType type; ///< The type of the event.
		radar::Radar* source_object; ///< Pointer to the object that generated the event.
	};

	/**
	 * @struct EventComparator
	 * @brief A custom comparator for the event priority queue.
	 *
	 * This comparator creates a min-heap, ensuring that the event with the
	 * smallest timestamp is always at the top of the queue.
	 */
	struct EventComparator
	{
		/**
		 * @brief Compares two events based on their timestamps.
		 * @param a The first event.
		 * @param b The second event.
		 * @return True if event 'a' should occur after event 'b'.
		 */
		bool operator()(const Event& a, const Event& b) const noexcept { return a.timestamp > b.timestamp; }
	};

	/**
	 * @brief Converts an EventType enum to its string representation.
	 * @param type The event type.
	 * @return A string representing the event type.
	 */
	inline std::string toString(const EventType type)
	{
		switch (type)
		{
		case EventType::TX_PULSED_START:
			return "TxPulsedStart";
		case EventType::RX_PULSED_WINDOW_START:
			return "RxPulsedWindowStart";
		case EventType::RX_PULSED_WINDOW_END:
			return "RxPulsedWindowEnd";
		case EventType::TX_STREAMING_START:
			return "TxStreamingStart";
		case EventType::TX_STREAMING_END:
			return "TxStreamingEnd";
		case EventType::RX_STREAMING_START:
			return "RxStreamingStart";
		case EventType::RX_STREAMING_END:
			return "RxStreamingEnd";
		default:
			return "UnknownEvent";
		}
	}
}
