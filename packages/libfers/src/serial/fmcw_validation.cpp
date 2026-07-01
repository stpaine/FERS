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
		void warnStreamingMemory(const std::string& owner, const std::string_view mode_name)
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
				"{} will buffer about {:.2f} GiB of {} streaming output for simulation_duration={} s at "
				"effective_rate={} Hz. Large simulation_time * bandwidth/rate products need chunking.",
				owner, estimated_gib, mode_name, duration, effective_rate);
		}

		/// Returns schedule periods, or the full simulation as one implicit active period.
		std::vector<radar::SchedulePeriod> effectiveSchedule(const std::vector<radar::SchedulePeriod>& schedule)
		{
			if (!schedule.empty())
			{
				return schedule;
			}
			return {radar::SchedulePeriod{.start = params::startTime(), .end = params::endTime()}};
		}

		/// Validates common FMCW sweep constraints.
		void validateCommonSweep(const fers_signal::RadarSignal& wave, const std::string& owner,
								 const RealType chirp_duration, const RealType chirp_bandwidth,
								 const RealType sweep_start, const RealType sweep_end, const Thrower& throw_error)
		{
			if (chirp_duration <= 0.0)
			{
				throw_error(owner + " has FMCW chirp_duration <= 0.");
			}
			if (chirp_bandwidth <= 0.0)
			{
				throw_error(owner + " has FMCW chirp_bandwidth <= 0.");
			}

			const RealType f_low = std::min(sweep_start, sweep_end);
			const RealType f_high = std::max(sweep_start, sweep_end);
			const RealType effective_rate = params::rate() * params::oversampleRatio();
			const RealType max_baseband = std::max(std::abs(f_low), std::abs(f_high));
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

			const RealType min_rf = wave.getCarrier() + f_low;
			if (min_rf <= 0.0)
			{
				throw_error(owner + " yields a non-positive RF-equivalent frequency.");
			}
		}

		void validateSfcwWaveform(const fers_signal::RadarSignal& wave, const std::string& owner,
								  const fers_signal::SteppedFrequencySignal& sfcw, const Thrower& throw_error)
		{
			if (sfcw.getStepCount() < 1U)
			{
				throw_error(owner + " has SFCW step_count < 1.");
			}
			if (sfcw.getStepSize() == 0.0 || !std::isfinite(sfcw.getStepSize()))
			{
				throw_error(owner + " has invalid SFCW step_size.");
			}
			if (sfcw.getDwellTime() <= 0.0 || !std::isfinite(sfcw.getDwellTime()))
			{
				throw_error(owner + " has SFCW dwell_time <= 0.");
			}
			if (sfcw.getStepPeriod() <= 0.0 || !std::isfinite(sfcw.getStepPeriod()))
			{
				throw_error(owner + " has SFCW step_period <= 0.");
			}
			if (sfcw.getDwellTime() > sfcw.getStepPeriod())
			{
				throw_error(owner + " has SFCW dwell_time longer than step_period.");
			}

			const RealType f_low =
				std::min(sfcw.firstFrequency(wave.getCarrier()), sfcw.lastFrequency(wave.getCarrier()));
			if (f_low <= 0.0)
			{
				throw_error(owner + " yields a non-positive SFCW RF step frequency.");
			}

			const RealType output_samples_per_dwell = params::rate() * sfcw.getDwellTime();
			if (output_samples_per_dwell < 1.0)
			{
				LOG(logging::Level::WARNING, "{} has fewer than one output sample per SFCW dwell.", owner);
			}
			const RealType internal_samples_per_dwell =
				params::rate() * static_cast<RealType>(params::oversampleRatio()) * sfcw.getDwellTime();
			if (internal_samples_per_dwell < 1.0)
			{
				LOG(logging::Level::WARNING, "{} has fewer than one internal sample per SFCW dwell.", owner);
			}
			warnStreamingMemory(owner, "SFCW");
		}
	}

	void validateWaveform(const fers_signal::RadarSignal& wave, const std::string& owner, const Thrower& throw_error)
	{
		if (const auto* fmcw = wave.getFmcwChirpSignal(); fmcw != nullptr)
		{
			if (fmcw->getChirpPeriod() < fmcw->getChirpDuration())
			{
				throw_error(owner + " has chirp_period shorter than chirp_duration; FMCW requires T_rep >= T_c.");
			}

			const RealType sweep_start = fmcw->getStartFrequencyOffset();
			const RealType sweep_end =
				sweep_start + (fmcw->isDownChirp() ? -fmcw->getChirpBandwidth() : fmcw->getChirpBandwidth());
			validateCommonSweep(wave, owner, fmcw->getChirpDuration(), fmcw->getChirpBandwidth(), sweep_start,
								sweep_end, throw_error);
			warnStreamingMemory(owner, "FMCW");
			return;
		}

		if (const auto* triangle = wave.getFmcwTriangleSignal(); triangle != nullptr)
		{
			const RealType sweep_start = triangle->getStartFrequencyOffset();
			const RealType sweep_end = sweep_start + triangle->getChirpBandwidth();
			validateCommonSweep(wave, owner, triangle->getChirpDuration(), triangle->getChirpBandwidth(), sweep_start,
								sweep_end, throw_error);
			warnStreamingMemory(owner, "FMCW");
			return;
		}

		if (const auto* sfcw = wave.getSteppedFrequencySignal(); sfcw != nullptr)
		{
			validateSfcwWaveform(wave, owner, *sfcw, throw_error);
		}
	}

	void validateWaveformModeMatch(const fers_signal::RadarSignal& wave, const radar::OperationMode mode,
								   const std::string& owner, const Thrower& throw_error)
	{
		const bool is_pulsed = !wave.isCw() && !wave.isFmcwFamily() && !wave.isSteppedFrequency();
		const bool matches = (mode == radar::OperationMode::PULSED_MODE && is_pulsed) ||
			(mode == radar::OperationMode::CW_MODE && wave.isCw()) ||
			(mode == radar::OperationMode::FMCW_MODE && wave.isFmcwFamily()) ||
			(mode == radar::OperationMode::SFCW_MODE && wave.isSteppedFrequency());
		if (!matches)
		{
			throw_error(owner + " mode does not match waveform '" + wave.getName() + "'.");
		}
	}

	void validateSchedule(const std::vector<radar::SchedulePeriod>& schedule, const fers_signal::FmcwChirpSignal& fmcw,
						  const std::string& owner, const Thrower& throw_error)
	{
		for (const auto& period : effectiveSchedule(schedule))
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

	void validateSchedule(const std::vector<radar::SchedulePeriod>& schedule, const fers_signal::RadarSignal& wave,
						  const std::string& owner, const Thrower& throw_error)
	{
		if (const auto* fmcw = wave.getFmcwChirpSignal(); fmcw != nullptr)
		{
			validateSchedule(schedule, *fmcw, owner, throw_error);
			return;
		}

		if (const auto* sfcw = wave.getSteppedFrequencySignal(); sfcw != nullptr)
		{
			for (const auto& period : effectiveSchedule(schedule))
			{
				const RealType duration = period.end - period.start;
				if (duration < sfcw->getDwellTime())
				{
					LOG(logging::Level::WARNING, "{} has schedule period [{}, {}] shorter than one SFCW dwell ({}s).",
						owner, period.start, period.end, sfcw->getDwellTime());
				}
				if (duration < sfcw->getSweepPeriod())
				{
					LOG(logging::Level::WARNING, "{} has schedule period [{}, {}] shorter than one SFCW sweep ({}s).",
						owner, period.start, period.end, sfcw->getSweepPeriod());
				}
			}
			return;
		}

		const auto* triangle = wave.getFmcwTriangleSignal();
		if (triangle == nullptr)
		{
			return;
		}

		const RealType T_tri = triangle->getTrianglePeriod();
		for (const auto& period : effectiveSchedule(schedule))
		{
			const RealType duration = period.end - period.start;
			if (duration < T_tri)
			{
				throw_error(std::format(
					"{} has schedule period [{}, {}] duration {} s shorter than FMCW triangle_period T_tri={} s.",
					owner, period.start, period.end, duration, T_tri));
			}
			const RealType full_triangles = std::floor(duration / T_tri);
			const RealType used = full_triangles * T_tri;
			const RealType leftover = duration - used;
			if (leftover > 1.0e-12)
			{
				LOG(logging::Level::WARNING,
					"{} has schedule period [{}, {}] that is not an integer multiple of FMCW triangle_period ({}s); "
					"{}s will be silent after the last complete triangle.",
					owner, period.start, period.end, T_tri, leftover);
			}
		}
	}
}
