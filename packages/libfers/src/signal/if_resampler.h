// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file if_resampler.h
 * @brief Internal FMCW IF rational resampler planning and streaming sink.
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <vector>

#include "core/config.h"

namespace fers_signal
{
	enum class FmcwIfResamplerStageKind
	{
		HalfBandDecimateBy2,
		RationalPolyphase
	};

	struct FmcwIfResamplerLimits
	{
		std::size_t max_taps_per_stage = 4096;
		RealType max_macs_per_output_sample = 2048.0;
		std::size_t max_phase_refinement = 64;
		std::uint64_t max_ratio_denominator = 1'000'000;
		RealType ratio_relative_tolerance = 1.0e-12;
	};

	struct FmcwIfRateRatio
	{
		std::uint64_t numerator = 1;
		std::uint64_t denominator = 1;
		RealType requested_ratio = 1.0;
		RealType actual_ratio = 1.0;
		RealType relative_error = 0.0;
	};

	struct FmcwIfResamplerRequest
	{
		RealType input_sample_rate_hz = 0.0;
		RealType output_sample_rate_hz = 0.0;
		RealType filter_bandwidth_hz = 0.0;
		std::optional<RealType> filter_transition_width_hz;
		RealType stopband_attenuation_db = 80.0;
		FmcwIfResamplerLimits limits{};
	};

	struct FmcwIfResamplerStagePlan
	{
		FmcwIfResamplerStageKind kind = FmcwIfResamplerStageKind::RationalPolyphase;
		RealType input_sample_rate_hz = 0.0;
		RealType output_sample_rate_hz = 0.0;
		std::uint64_t up_factor = 1;
		std::uint64_t down_factor = 1;
		std::size_t tap_count = 0;
		std::size_t phase_refinement = 1;
		std::size_t phase_count = 1;
		std::int64_t initial_input_advance = 0;
		std::size_t initial_phase_accumulator = 0;
		RealType initial_branch_interpolation_fraction = 0.0;
		bool applies_fractional_delay = false;
		RealType filter_bandwidth_hz = 0.0;
		RealType transition_width_hz = 0.0;
		RealType cutoff_hz = 0.0;
		RealType stopband_attenuation_db = 0.0;
		RealType group_delay_seconds = 0.0;
		RealType estimated_macs_per_stage_output = 0.0;
	};

	struct FmcwIfResamplerPlan
	{
		RealType input_sample_rate_hz = 0.0;
		RealType requested_output_sample_rate_hz = 0.0;
		RealType actual_output_sample_rate_hz = 0.0;
		RealType filter_bandwidth_hz = 0.0;
		RealType filter_transition_width_hz = 0.0;
		RealType stopband_attenuation_db = 0.0;
		FmcwIfRateRatio overall_ratio{};
		std::vector<FmcwIfResamplerStagePlan> stages;
		RealType estimated_macs_per_output_sample = 0.0;
		RealType group_delay_seconds = 0.0;
		RealType group_delay_output_samples = 0.0;
		std::uint64_t warmup_discard_samples = 0;
		RealType fractional_output_delay_samples = 0.0;
		std::size_t phase_refinement = 1;
		RealType fractional_phase_offset = 0.0;
		RealType branch_interpolation_fraction = 0.0;
		RealType estimated_timing_error_seconds = 0.0;
		RealType estimated_phase_error_radians = 0.0;
		bool group_delay_compensated = true;
		FmcwIfResamplerLimits limits{};
	};

	[[nodiscard]] FmcwIfRateRatio approximateFmcwIfRateRatio(RealType output_sample_rate_hz,
															 RealType input_sample_rate_hz,
															 const FmcwIfResamplerLimits& limits = {});

	[[nodiscard]] FmcwIfResamplerPlan planFmcwIfResampler(const FmcwIfResamplerRequest& request);

	struct FmcwIfZeroInputResult
	{
		std::vector<ComplexType> emitted;
		std::size_t skipped_output_samples = 0;
	};

	class FmcwIfResamplingSink
	{
	public:
		explicit FmcwIfResamplingSink(FmcwIfResamplerPlan plan);
		~FmcwIfResamplingSink();

		FmcwIfResamplingSink(const FmcwIfResamplingSink&) = delete;
		FmcwIfResamplingSink& operator=(const FmcwIfResamplingSink&) = delete;
		FmcwIfResamplingSink(FmcwIfResamplingSink&&) noexcept;
		FmcwIfResamplingSink& operator=(FmcwIfResamplingSink&&) noexcept;

		void consume(std::span<const ComplexType> block);
		[[nodiscard]] FmcwIfZeroInputResult consumeZeroInput(std::size_t input_count);
		[[nodiscard]] std::vector<ComplexType> takeOutput();
		[[nodiscard]] std::vector<ComplexType> finish();
		void reset();

		[[nodiscard]] const FmcwIfResamplerPlan& plan() const noexcept { return _plan; }

	private:
		class Stage;

		FmcwIfResamplerPlan _plan;
		std::vector<std::unique_ptr<Stage>> _stages;
		std::vector<ComplexType> _output;
		bool _finished = false;
	};
}
