// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "simulation_state.h"

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "signal/radar_signal.h"

namespace core
{
	namespace
	{
		/// Rounds a non-negative floating-point value up to an unsigned integer.
		std::uint64_t ceilToUint(const RealType value)
		{
			if (value <= 0.0)
			{
				return 0;
			}

			const RealType nearest = std::round(value);
			const RealType tolerance = 1.0e-12 * std::max<RealType>(1.0, std::abs(nearest));
			if (std::abs(value - nearest) <= tolerance)
			{
				return static_cast<std::uint64_t>(nearest);
			}
			return static_cast<std::uint64_t>(std::ceil(value));
		}

		/// Returns the first FMCW chirp index that can contribute inside an interval.
		std::optional<std::uint64_t> firstFmcwChirpIndex(const ActiveStreamingSource& source,
														 const RealType active_start, const RealType active_end)
		{
			if (!source.is_fmcw || source.chirp_period <= 0.0)
			{
				return std::nullopt;
			}

			const RealType clipped_start = std::max(active_start, source.segment_start);
			const RealType clipped_end = std::min(active_end, source.segment_end);
			if (clipped_end <= clipped_start)
			{
				return std::nullopt;
			}

			const auto first_index = clipped_start <= source.segment_start
				? std::uint64_t{0}
				: ceilToUint((clipped_start - source.segment_start) / source.chirp_period);
			if (source.chirp_count.has_value() && first_index >= *source.chirp_count)
			{
				return std::nullopt;
			}

			const RealType first_start =
				source.segment_start + static_cast<RealType>(first_index) * source.chirp_period;
			if (first_start >= clipped_end)
			{
				return std::nullopt;
			}
			return first_index;
		}
	}

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
		source.is_fmcw = signal->isFmcwChirp();
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
		source.signed_chirp_rate = source.fmcw->getSignedChirpRate();
		source.start_freq_off = source.fmcw->getStartFrequencyOffset();
		source.two_pi_f0 = 2.0 * PI * source.start_freq_off;
		source.s_pi_alpha = PI * source.signed_chirp_rate;
		source.chirp_count = source.fmcw->getChirpCount();
		if (source.chirp_count.has_value())
		{
			source.segment_end = std::min(
				source.segment_end, segment_start + static_cast<RealType>(*source.chirp_count) * source.chirp_period);
		}

		return source;
	}

	std::optional<RealType> firstFmcwChirpStart(const ActiveStreamingSource& source, const RealType active_start,
												const RealType active_end)
	{
		const auto first_index = firstFmcwChirpIndex(source, active_start, active_end);
		if (!first_index.has_value())
		{
			return std::nullopt;
		}
		return source.segment_start + static_cast<RealType>(*first_index) * source.chirp_period;
	}

	std::uint64_t countFmcwChirpStarts(const ActiveStreamingSource& source, const RealType active_start,
									   const RealType active_end)
	{
		const auto first_index = firstFmcwChirpIndex(source, active_start, active_end);
		if (!first_index.has_value())
		{
			return 0;
		}

		const RealType clipped_end = std::min(active_end, source.segment_end);
		const RealType first_start = source.segment_start + static_cast<RealType>(*first_index) * source.chirp_period;
		const auto starts_in_interval = ceilToUint((clipped_end - first_start) / source.chirp_period);
		if (!source.chirp_count.has_value())
		{
			return starts_in_interval;
		}

		const auto configured = static_cast<std::uint64_t>(*source.chirp_count);
		return std::min(starts_in_interval, configured - *first_index);
	}
}
