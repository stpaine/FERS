// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "fmcw_validation.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <format>

#include "core/logging.h"
#include "core/parameters.h"
#include "signal/radar_signal.h"

namespace serial::fmcw_validation
{
	namespace
	{
		/// Emits a warning when a streaming scenario may allocate a large I/Q buffer.
		void warnStreamingMemory(const std::string& owner)
		{
			const RealType duration = std::max<RealType>(0.0, params::endTime() - params::startTime());
			const RealType effective_rate = params::rate() * params::oversampleRatio();
			const auto estimated_samples = static_cast<std::uint64_t>(std::ceil(duration * effective_rate));
			const std::uint64_t estimated_bytes = estimated_samples * sizeof(ComplexType);
			constexpr std::uint64_t warning_threshold_bytes = 2ULL * 1024ULL * 1024ULL * 1024ULL;
			if (estimated_bytes < warning_threshold_bytes)
			{
				return;
			}

			const double estimated_gib = static_cast<double>(estimated_bytes) / (1024.0 * 1024.0 * 1024.0);
			LOG(logging::Level::WARNING,
				"{} will allocate about {:.2f} GiB of FMCW streaming IQ data for simulation_duration={} s at "
				"effective_rate={} Hz. Large simulation_time * bandwidth/rate products need chunking.",
				owner, estimated_gib, duration, effective_rate);
		}
	}

	void validateWaveform(const fers_signal::RadarSignal& wave, const std::string& owner, const Thrower& throw_error)
	{
		const auto* fmcw = wave.getFmcwChirpSignal();
		if (fmcw == nullptr)
		{
			return;
		}

		if (fmcw->getChirpDuration() <= 0.0)
		{
			throw_error(owner + " has FMCW chirp_duration <= 0.");
		}
		if (fmcw->getChirpBandwidth() <= 0.0)
		{
			throw_error(owner + " has FMCW chirp_bandwidth <= 0.");
		}
		if (fmcw->getChirpPeriod() < fmcw->getChirpDuration())
		{
			throw_error(owner + " has chirp_period shorter than chirp_duration; FMCW requires T_rep >= T_c.");
		}

		const RealType effective_rate = params::rate() * params::oversampleRatio();
		const RealType max_baseband = std::max(std::abs(fmcw->getStartFrequencyOffset()),
											   std::abs(fmcw->getStartFrequencyOffset() + fmcw->getChirpBandwidth()));
		if (effective_rate <= max_baseband)
		{
			throw_error(owner + " violates FMCW baseband aliasing constraint.");
		}
		if (max_baseband > 0.0 && effective_rate < 1.1 * max_baseband)
		{
			LOG(logging::Level::WARNING,
				"{} is within 10% of the FMCW aliasing limit: effective_rate={} Hz, max_baseband={} Hz.", owner,
				effective_rate, max_baseband);
		}

		const RealType min_rf = wave.getCarrier() +
			std::min(fmcw->getStartFrequencyOffset(), fmcw->getStartFrequencyOffset() + fmcw->getChirpBandwidth());
		if (min_rf <= 0.0)
		{
			throw_error(owner + " yields a non-positive RF-equivalent frequency.");
		}

		warnStreamingMemory(owner);
	}

	void validateWaveformModeMatch(const fers_signal::RadarSignal& wave, const radar::OperationMode mode,
								   const std::string& owner, const Thrower& throw_error)
	{
		const bool is_pulsed = !wave.isCw() && !wave.isFmcwUpChirp();
		const bool matches = (mode == radar::OperationMode::PULSED_MODE && is_pulsed) ||
			(mode == radar::OperationMode::CW_MODE && wave.isCw()) ||
			(mode == radar::OperationMode::FMCW_MODE && wave.isFmcwUpChirp());
		if (!matches)
		{
			throw_error(owner + " mode does not match waveform '" + wave.getName() + "'.");
		}
	}

	void validateSchedule(const std::vector<radar::SchedulePeriod>& schedule, const fers_signal::FmcwChirpSignal& fmcw,
						  const std::string& owner, const Thrower& throw_error)
	{
		for (const auto& period : schedule)
		{
			const RealType duration = period.end - period.start;
			if (duration < fmcw.getChirpDuration())
			{
				throw_error(std::format(
					"{} has schedule period [{}, {}] duration {} s shorter than FMCW chirp_duration T_c={} s.", owner,
					period.start, period.end, duration, fmcw.getChirpDuration()));
			}
			if (duration < fmcw.getChirpPeriod())
			{
				LOG(logging::Level::WARNING, "{} has a schedule period [{}, {}] shorter than FMCW chirp_period ({}s).",
					owner, period.start, period.end, fmcw.getChirpPeriod());
			}
		}
	}
}
