// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#pragma once

#include <atomic>
#include <cassert>
#include <cstdint>
#include <stdexcept>

/**
 * @brief 64-bit Unique Simulation ID.
 * Structure: [16-bit Type][48-bit Counter]
 */
using SimId = uint64_t;

/**
 * @enum ObjectType
 * @brief Categorizes objects for ID generation.
 */
enum class ObjectType : uint16_t
{
	Unknown = 0,
	Platform = 1,
	Transmitter = 2,
	Receiver = 3,
	Target = 4,
	Antenna = 5,
	Waveform = 6,
	Timing = 7,
	Debug = 0xFFFF
};

/**
 * @class SimIdGenerator
 * @brief Thread-safe Meyers singleton for generating unique object IDs.
 */
class SimIdGenerator
{
public:
	/**
	 * @brief Get the singleton instance of SimIdGenerator.
	 * @return Reference to the SimIdGenerator instance.
	 */
	static SimIdGenerator& instance()
	{
		static SimIdGenerator instance;
		return instance;
	}

	/**
	 * @brief Generate a unique SimId for a given object type.
	 * @param type The ObjectType for which to generate the ID.
	 * @return A unique SimId.
	 * @throws std::overflow_error if the 48-bit counter overflows.
	 */
	SimId generateId(ObjectType type)
	{
		assert(type != ObjectType::Debug && "generateId called with reserved ObjectType::Debug");

		// Increment 48-bit counter
		uint64_t count = _counter.fetch_add(1, std::memory_order_relaxed);
		constexpr uint64_t MAX_COUNTER_VALUE = 0x0000FFFFFFFFFFFF;
		if (count >= MAX_COUNTER_VALUE)
		{
			throw std::overflow_error("FERS object ID counter has overflowed the 48-bit space.");
		}

		// The mask is technically redundant due to the check above
		count &= MAX_COUNTER_VALUE;

		// Shift type into upper 16 bits
		const uint64_t type_bits = static_cast<uint64_t>(type) << 48;
		return type_bits | count;
	}

	/**
	 * @brief Generate a debug SimId.
	 * @return A unique SimId with ObjectType::Debug.
	 */
	SimId generateDebugId()
	{
		// Increment 48-bit counter
		uint64_t count = _counter.fetch_add(1, std::memory_order_relaxed);
		constexpr uint64_t MAX_COUNTER_VALUE = 0x0000FFFFFFFFFFFF;
		if (count >= MAX_COUNTER_VALUE)
		{
			throw std::overflow_error("FERS object ID counter has overflowed the 48-bit space.");
		}

		// The mask is technically redundant due to the check above
		count &= MAX_COUNTER_VALUE;

		// Shift debug type into upper 16 bits
		constexpr uint64_t type_bits = static_cast<uint64_t>(ObjectType::Debug) << 48;
		return type_bits | count;
	}

	/**
	 * @brief Extract object type from SimId
	 */
	static ObjectType getType(const SimId id) { return static_cast<ObjectType>(id >> 48); }

	/**
	 * @brief Extract counter from SimId
	 */
	static uint64_t getCounter(const SimId id) { return id & 0x0000FFFFFFFFFFFF; }

private:
	SimIdGenerator() = default; ///< Private constructor for singleton
	std::atomic<uint64_t> _counter{1}; ///< 48-bit counter for unique IDs
};
