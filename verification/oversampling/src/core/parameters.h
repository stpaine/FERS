// SPDX-License-Identifier: GPL-2.0-only
//
// Local harness support code derived from:
//   packages/libfers/src/core/parameters.h

#pragma once

#include <optional>
#include <random>
#include <stdexcept>

#include "config.h"

namespace params
{
	struct Parameters
	{
		constexpr static RealType DEFAULT_C = 299792458.0;
		constexpr static RealType DEFAULT_BOLTZMANN_K = 1.3806503e-23;

		RealType c = DEFAULT_C;
		RealType boltzmann_k = DEFAULT_BOLTZMANN_K;
		RealType start = 0.0;
		RealType end = 0.0;
		RealType sim_sampling_rate = 1000.0;
		RealType rate = 0.0;
		std::optional<unsigned> random_seed;
		unsigned adc_bits = 0;
		unsigned filter_length = 33;
		unsigned render_threads = 1;
		unsigned oversample_ratio = 1;

		void reset() noexcept { *this = Parameters{}; }
	};

	inline Parameters params;

	inline RealType c() noexcept { return params.c; }
	inline RealType boltzmannK() noexcept { return params.boltzmann_k; }
	inline RealType startTime() noexcept { return params.start; }
	inline RealType endTime() noexcept { return params.end; }
	inline RealType simSamplingRate() noexcept { return params.sim_sampling_rate; }
	inline RealType rate() noexcept { return params.rate; }
	inline unsigned randomSeed() noexcept { return params.random_seed.value_or(0); }
	inline unsigned adcBits() noexcept { return params.adc_bits; }
	inline unsigned renderFilterLength() noexcept { return params.filter_length; }
	inline unsigned renderThreads() noexcept { return params.render_threads; }
	inline unsigned oversampleRatio() noexcept { return params.oversample_ratio; }

	inline void setTime(const RealType start_time, const RealType end_time) noexcept
	{
		params.start = start_time;
		params.end = end_time;
	}

	inline void setRate(const RealType rate_value)
	{
		if (rate_value <= 0.0)
		{
			throw std::runtime_error("Sampling rate must be > 0");
		}
		params.rate = rate_value;
	}

	inline void setRandomSeed(const unsigned seed) noexcept { params.random_seed = seed; }
	inline void setAdcBits(const unsigned bits) noexcept { params.adc_bits = bits; }

	inline void setOversampleRatio(const unsigned ratio)
	{
		if (ratio == 0)
		{
			throw std::runtime_error("Oversample ratio must be >= 1");
		}
		params.oversample_ratio = ratio;
	}
}

