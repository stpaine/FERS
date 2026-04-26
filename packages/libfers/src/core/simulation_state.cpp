// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "simulation_state.h"

#include <algorithm>
#include <cmath>

#include "signal/radar_signal.h"

namespace core
{
	ActiveStreamingSource makeActiveSource(const radar::Transmitter* const tx, const RealType segment_start,
										   const RealType segment_end)
	{
		ActiveStreamingSource source{};
		source.transmitter = tx;
		source.segment_start = segment_start;
		source.segment_end = segment_end;
		if (tx == nullptr)
		{
			return source;
		}

		const auto* const signal = tx->getSignal();
		if (signal == nullptr)
		{
			return source;
		}

		source.carrier_freq = signal->getCarrier();
		source.amplitude = std::sqrt(signal->getPower());
		source.is_fmcw = signal->isFmcwUpChirp();
		if (!source.is_fmcw)
		{
			return source;
		}

		source.fmcw = signal->getFmcwChirpSignal();
		if (source.fmcw == nullptr)
		{
			source.is_fmcw = false;
			return source;
		}

		source.chirp_duration = source.fmcw->getChirpDuration();
		source.chirp_period = source.fmcw->getChirpPeriod();
		source.chirp_rate = source.fmcw->getChirpRate();
		source.start_freq_off = source.fmcw->getStartFrequencyOffset();
		source.two_pi_f0 = 2.0 * PI * source.start_freq_off;
		source.pi_alpha = PI * source.chirp_rate;
		source.chirp_count = source.fmcw->getChirpCount();
		if (source.chirp_count.has_value())
		{
			source.segment_end = std::min(
				source.segment_end, segment_start + static_cast<RealType>(*source.chirp_count) * source.chirp_period);
		}

		return source;
	}
}
