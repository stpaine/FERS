// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "schedule_period.h"

#include <algorithm>

#include "core/logging.h"
#include "core/parameters.h"

namespace radar
{

	std::vector<SchedulePeriod> processRawSchedule(std::vector<SchedulePeriod> periods, const std::string& ownerName,
												   const bool isPulsed, const RealType pri)
	{
		if (periods.empty())
		{
			return {};
		}

		std::vector<SchedulePeriod> valid_periods;
		valid_periods.reserve(periods.size());

		const RealType sim_start = params::startTime();
		const RealType sim_end = params::endTime();

		// 1. Filter invalid and out-of-bounds periods
		for (const auto& p : periods)
		{
			if (p.start >= p.end)
			{
				LOG(logging::Level::WARNING,
					"Object '{}' has a schedule period with start ({}) >= end ({}). Ignoring period.", ownerName,
					p.start, p.end);
				continue;
			}

			if (p.end <= sim_start || p.start >= sim_end)
			{
				LOG(logging::Level::WARNING,
					"Object '{}' has a schedule period [{}, {}] completely outside simulation time. Ignoring.",
					ownerName, p.start, p.end);
				continue;
			}
			valid_periods.push_back(p);
		}

		if (valid_periods.empty())
		{
			return {};
		}

		// 2. Sort by start time
		std::sort(valid_periods.begin(), valid_periods.end(),
				  [](const auto& a, const auto& b) { return a.start < b.start; });

		// 3. Merge overlapping intervals
		std::vector<SchedulePeriod> merged;
		merged.reserve(valid_periods.size());
		merged.push_back(valid_periods.front());

		for (size_t i = 1; i < valid_periods.size(); ++i)
		{
			if (valid_periods[i].start <= merged.back().end)
			{
				// Overlap or adjacency: extend the previous period
				merged.back().end = std::max(merged.back().end, valid_periods[i].end);
			}
			else
			{
				merged.push_back(valid_periods[i]);
			}
		}

		// 4. Check PRI constraints
		if (isPulsed)
		{
			for (const auto& period : merged)
			{
				if (period.end - period.start < pri)
				{
					LOG(logging::Level::WARNING, "Object '{}' has a schedule period [{}, {}] shorter than PRI ({}s).",
						ownerName, period.start, period.end, pri);
				}
			}
		}

		return merged;
	}

} // namespace radar
