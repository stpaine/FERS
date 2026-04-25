// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#pragma once

#include <functional>
#include <string>
#include <vector>

#include "core/config.h"
#include "radar/radar_obj.h"
#include "radar/schedule_period.h"

namespace fers_signal
{
	class FmcwChirpSignal;
	class RadarSignal;
}

namespace serial::fmcw_validation
{
	using Thrower = std::function<void(const std::string&)>;

	void validateWaveform(const fers_signal::RadarSignal& wave, const std::string& owner, const Thrower& throw_error);

	void validateWaveformModeMatch(const fers_signal::RadarSignal& wave, radar::OperationMode mode,
								   const std::string& owner, const Thrower& throw_error);

	void validateSchedule(const std::vector<radar::SchedulePeriod>& schedule, const fers_signal::FmcwChirpSignal& fmcw,
						  const std::string& owner, const Thrower& throw_error);
}
