// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#pragma once

#include <string>
#include <vector>

#include "core/config.h"

namespace radar
{
	/**
	 * @struct SchedulePeriod
	 * @brief Represents a time period during which the transmitter is active.
	 */
	struct SchedulePeriod
	{
		RealType start{};
		RealType end{};
	};

	/**
	 * @brief Processes a raw list of schedule periods.
	 *
	 * This function performs the following operations:
	 * 1. Filters invalid periods (start >= end).
	 * 2. Filters periods completely outside simulation bounds.
	 * 3. Sorts periods by start time.
	 * 4. Merges overlapping or adjacent periods.
	 * 5. Checks against PRI constraints (if pulsed).
	 *
	 * @param periods The raw vector of periods.
	 * @param ownerName The name of the object owning this schedule (for logging).
	 * @param isPulsed Whether the object is operating in pulsed mode.
	 * @param pri The Pulse Repetition Interval (only used if isPulsed is true).
	 * @return A sorted, merged, and validated vector of periods.
	 */
	std::vector<SchedulePeriod> processRawSchedule(std::vector<SchedulePeriod> periods, const std::string& ownerName,
												   bool isPulsed, RealType pri);
}
