// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file if_resampler.cpp
 * @brief Internal FMCW IF rational resampler planning and streaming sink.
 */

#include "if_resampler.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{
	constexpr RealType kHalfBandTransitionFraction = 0.05;
	constexpr RealType kHalfBandPassbandFraction = 0.20;

	[[nodiscard]] RealType sinc(const RealType x) noexcept { return x == 0.0 ? 1.0 : std::sin(PI * x) / (PI * x); }

	[[nodiscard]] RealType besselI0(const RealType x)
	{
		if (x < 0.0)
		{
			throw std::invalid_argument("Kaiser Bessel approximation requires x >= 0");
		}

		RealType t = x / 3.75;
		if (t <= 1.0)
		{
			t *= t;
			return 1.0 +
				t * (3.5156229 + t * (3.0899424 + t * (1.2067492 + t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))));
		}

		const RealType inv = 3.75 / x;
		const RealType poly = 0.39894228 +
			inv *
				(0.01328592 +
				 inv *
					 (0.00225319 +
					  inv *
						  (-0.00157565 +
						   inv *
							   (0.00916281 +
								inv * (-0.02057706 + inv * (0.02635537 + inv * (-0.01647633 + inv * 0.00392377)))))));
		return poly * std::exp(x) / std::sqrt(x);
	}

	[[nodiscard]] RealType kaiserBeta(const RealType attenuation_db) noexcept
	{
		if (attenuation_db > 50.0)
		{
			return 0.1102 * (attenuation_db - 8.7);
		}
		if (attenuation_db >= 21.0)
		{
			return 0.5842 * std::pow(attenuation_db - 21.0, 0.4) + 0.07886 * (attenuation_db - 21.0);
		}
		return 0.0;
	}

	[[nodiscard]] RealType kaiserWindow(const std::size_t index, const std::size_t taps, const RealType beta,
										const RealType bessel_beta)
	{
		if (taps <= 1)
		{
			return 1.0;
		}
		const RealType ratio = (2.0 * static_cast<RealType>(index)) / static_cast<RealType>(taps - 1) - 1.0;
		const RealType arg = beta * std::sqrt(std::max<RealType>(0.0, 1.0 - ratio * ratio));
		return besselI0(arg) / bessel_beta;
	}

	[[nodiscard]] std::size_t oddTapCount(std::size_t taps) noexcept
	{
		taps = std::max<std::size_t>(3, taps);
		return (taps % 2 == 0) ? taps + 1 : taps;
	}

	[[nodiscard]] std::size_t estimateKaiserTapCount(const RealType transition_width_hz,
													 const RealType design_sample_rate_hz,
													 const RealType attenuation_db)
	{
		if (!std::isfinite(transition_width_hz) || transition_width_hz <= 0.0 ||
			!std::isfinite(design_sample_rate_hz) || design_sample_rate_hz <= 0.0)
		{
			throw std::invalid_argument("IF resampler transition width and design sample rate must be positive");
		}

		const RealType normalized_transition = transition_width_hz / design_sample_rate_hz;
		if (normalized_transition <= 0.0 || normalized_transition >= 0.5)
		{
			throw std::invalid_argument("IF resampler normalized transition width is outside (0, 0.5)");
		}

		const RealType width_rad = 2.0 * PI * normalized_transition;
		const RealType estimate = std::max<RealType>(3.0, (attenuation_db - 8.0) / (2.285 * width_rad));
		return oddTapCount(static_cast<std::size_t>(std::ceil(estimate)));
	}

	void validateLimit(const bool condition, const std::string& message)
	{
		if (!condition)
		{
			throw std::invalid_argument(message);
		}
	}

	[[nodiscard]] RealType derivedTransitionWidth(const fers_signal::FmcwIfResamplerRequest& request)
	{
		const RealType nyquist = request.output_sample_rate_hz * 0.5;
		const RealType available = nyquist - request.filter_bandwidth_hz;
		if (request.filter_transition_width_hz.has_value())
		{
			validateLimit(*request.filter_transition_width_hz > 0.0,
						  "IF resampler filter transition width must be positive");
			validateLimit(*request.filter_transition_width_hz <= available,
						  "IF resampler filter transition width must fit below output Nyquist");
			return *request.filter_transition_width_hz;
		}
		return available;
	}

	[[nodiscard]] std::size_t ceilBranchMacs(const std::size_t taps, const std::uint64_t up_factor,
											 const bool interpolates_adjacent_phase) noexcept
	{
		const auto up = std::max<std::uint64_t>(1, up_factor);
		const auto branch = static_cast<std::size_t>((taps + up - 1) / up);
		return interpolates_adjacent_phase ? branch * 2 : branch;
	}

	[[nodiscard]] fers_signal::FmcwIfResamplerStagePlan
	makeStage(const fers_signal::FmcwIfResamplerStageKind kind, const RealType input_rate_hz,
			  const std::uint64_t up_factor, const std::uint64_t down_factor, const RealType passband_hz,
			  const RealType requested_transition_hz, const RealType stopband_attenuation_db,
			  const fers_signal::FmcwIfResamplerLimits& limits)
	{
		const RealType output_rate_hz =
			input_rate_hz * static_cast<RealType>(up_factor) / static_cast<RealType>(down_factor);
		const RealType limiting_nyquist = 0.5 * std::min(input_rate_hz, output_rate_hz);
		const RealType available_transition = limiting_nyquist - passband_hz;
		if (available_transition <= 0.0)
		{
			throw std::runtime_error("IF resampler passband does not fit stage Nyquist limit");
		}

		const RealType transition_hz = std::min(requested_transition_hz, available_transition);
		if (transition_hz <= 0.0)
		{
			throw std::runtime_error("IF resampler transition width is too narrow for stage");
		}

		const RealType design_rate_hz = kind == fers_signal::FmcwIfResamplerStageKind::RationalPolyphase
			? input_rate_hz * static_cast<RealType>(up_factor)
			: input_rate_hz;
		const auto taps = estimateKaiserTapCount(transition_hz, design_rate_hz, stopband_attenuation_db);
		if (taps > limits.max_taps_per_stage)
		{
			throw std::runtime_error(
				"IF resampler FIR tap count " + std::to_string(taps) + " exceeds maximum " +
				std::to_string(limits.max_taps_per_stage) +
				"; increase if_sample_rate, reduce if_filter_bandwidth, or increase transition width");
		}

		fers_signal::FmcwIfResamplerStagePlan stage;
		stage.kind = kind;
		stage.input_sample_rate_hz = input_rate_hz;
		stage.output_sample_rate_hz = output_rate_hz;
		stage.up_factor = up_factor;
		stage.down_factor = down_factor;
		stage.tap_count = taps;
		stage.phase_refinement =
			kind == fers_signal::FmcwIfResamplerStageKind::RationalPolyphase ? limits.max_phase_refinement : 1;
		stage.phase_count = static_cast<std::size_t>(up_factor) * stage.phase_refinement;
		stage.filter_bandwidth_hz = passband_hz;
		stage.transition_width_hz = transition_hz;
		stage.cutoff_hz = passband_hz + transition_hz * 0.5;
		stage.stopband_attenuation_db = stopband_attenuation_db;
		stage.group_delay_seconds = static_cast<RealType>(taps - 1) / (2.0 * design_rate_hz);
		if (kind == fers_signal::FmcwIfResamplerStageKind::HalfBandDecimateBy2)
		{
			stage.estimated_macs_per_stage_output = static_cast<RealType>((taps + 3) / 4);
		}
		else
		{
			stage.estimated_macs_per_stage_output =
				static_cast<RealType>(ceilBranchMacs(taps, up_factor, stage.phase_refinement > 1));
		}
		return stage;
	}

	[[nodiscard]] RealType checkedRatio(const RealType output_sample_rate_hz, const RealType input_sample_rate_hz)
	{
		validateLimit(std::isfinite(input_sample_rate_hz) && input_sample_rate_hz > 0.0,
					  "IF resampler input sample rate must be positive");
		validateLimit(std::isfinite(output_sample_rate_hz) && output_sample_rate_hz > 0.0,
					  "IF resampler output sample rate must be positive");
		validateLimit(output_sample_rate_hz <= input_sample_rate_hz,
					  "IF resampler output sample rate must be <= input sample rate");
		return output_sample_rate_hz / input_sample_rate_hz;
	}

	[[nodiscard]] std::size_t ceilOutputCount(const std::size_t input_count, const std::uint64_t up_factor,
											  const std::uint64_t down_factor) noexcept
	{
		const long double scaled = static_cast<long double>(input_count) * static_cast<long double>(up_factor) /
			static_cast<long double>(down_factor);
		return static_cast<std::size_t>(std::ceil(scaled - 1.0e-12L));
	}

	[[nodiscard]] RealType fractionalStageMacs(const fers_signal::FmcwIfResamplerStagePlan& stage,
											   const std::size_t refinement) noexcept
	{
		return static_cast<RealType>(ceilBranchMacs(stage.tap_count, stage.up_factor, refinement > 0));
	}

	void recomputePlanCost(fers_signal::FmcwIfResamplerPlan& plan)
	{
		plan.estimated_macs_per_output_sample = 0.0;
		plan.phase_refinement = 1;
		for (const auto& stage : plan.stages)
		{
			plan.estimated_macs_per_output_sample += stage.estimated_macs_per_stage_output;
			plan.phase_refinement = std::max(plan.phase_refinement, stage.phase_refinement);
		}
	}

	void configureFractionalDelay(fers_signal::FmcwIfResamplerPlan& plan)
	{
		constexpr RealType kFractionalEpsilon = 1.0e-12;
		if (plan.stages.empty() || plan.fractional_output_delay_samples <= kFractionalEpsilon)
		{
			return;
		}

		auto& stage = plan.stages.back();
		const RealType base_cost = plan.estimated_macs_per_output_sample - stage.estimated_macs_per_stage_output;
		std::size_t selected_refinement = 1;
		for (std::size_t refinement = plan.limits.max_phase_refinement; refinement >= 1; --refinement)
		{
			const RealType candidate_cost = base_cost + fractionalStageMacs(stage, refinement);
			if (candidate_cost <= plan.limits.max_macs_per_output_sample)
			{
				selected_refinement = refinement;
				break;
			}
			if (refinement == 1)
			{
				break;
			}
		}

		stage.phase_refinement = selected_refinement;
		stage.phase_count = static_cast<std::size_t>(stage.up_factor) * selected_refinement;
		stage.estimated_macs_per_stage_output = fractionalStageMacs(stage, selected_refinement);

		const RealType u0 = plan.fractional_output_delay_samples * static_cast<RealType>(stage.down_factor) *
			static_cast<RealType>(selected_refinement);
		const auto modulus = static_cast<RealType>(stage.phase_count);
		stage.initial_input_advance = static_cast<std::int64_t>(std::floor(u0 / modulus));
		const RealType residual = u0 - static_cast<RealType>(stage.initial_input_advance) * modulus;
		stage.initial_phase_accumulator = static_cast<std::size_t>(std::floor(residual));
		stage.initial_branch_interpolation_fraction = residual - static_cast<RealType>(stage.initial_phase_accumulator);
		stage.applies_fractional_delay = true;

		plan.fractional_phase_offset = plan.fractional_output_delay_samples * static_cast<RealType>(stage.down_factor);
		plan.branch_interpolation_fraction = stage.initial_branch_interpolation_fraction;
		plan.phase_refinement = selected_refinement;
		plan.estimated_timing_error_seconds =
			0.5 / (plan.actual_output_sample_rate_hz * static_cast<RealType>(selected_refinement));
		plan.estimated_phase_error_radians = 2.0 * PI * plan.filter_bandwidth_hz * plan.estimated_timing_error_seconds;

		recomputePlanCost(plan);
		if (plan.estimated_macs_per_output_sample > plan.limits.max_macs_per_output_sample)
		{
			throw std::runtime_error("IF resampler fractional-delay MAC cost " +
									 std::to_string(plan.estimated_macs_per_output_sample) + " exceeds maximum " +
									 std::to_string(plan.limits.max_macs_per_output_sample));
		}
	}
}

namespace fers_signal
{
	FmcwIfRateRatio approximateFmcwIfRateRatio(const RealType output_sample_rate_hz,
											   const RealType input_sample_rate_hz, const FmcwIfResamplerLimits& limits)
	{
		const RealType target = checkedRatio(output_sample_rate_hz, input_sample_rate_hz);
		if (std::abs(target - 1.0) <= limits.ratio_relative_tolerance)
		{
			return {.numerator = 1,
					.denominator = 1,
					.requested_ratio = target,
					.actual_ratio = 1.0,
					.relative_error = std::abs(target - 1.0)};
		}

		std::uint64_t prev_num = 0;
		std::uint64_t num = 1;
		std::uint64_t prev_den = 1;
		std::uint64_t den = 0;
		RealType x = target;
		FmcwIfRateRatio best{.requested_ratio = target};

		for (std::size_t iter = 0; iter < 128; ++iter)
		{
			const auto a = static_cast<std::uint64_t>(std::floor(x));
			const auto next_num = a * num + prev_num;
			const auto next_den = a * den + prev_den;
			if (next_den == 0 || next_den > limits.max_ratio_denominator)
			{
				break;
			}

			const RealType actual = static_cast<RealType>(next_num) / static_cast<RealType>(next_den);
			const RealType relative_error = std::abs(actual - target) / target;
			best = {.numerator = next_num,
					.denominator = next_den,
					.requested_ratio = target,
					.actual_ratio = actual,
					.relative_error = relative_error};
			if (relative_error <= limits.ratio_relative_tolerance)
			{
				return best;
			}

			const RealType remainder = x - static_cast<RealType>(a);
			if (std::abs(remainder) <= std::numeric_limits<RealType>::epsilon())
			{
				break;
			}
			x = 1.0 / remainder;
			prev_num = num;
			num = next_num;
			prev_den = den;
			den = next_den;
		}

		throw std::runtime_error("IF resampler cannot approximate output/input sample-rate ratio within tolerance; "
								 "increase max denominator or choose a representable IF sample rate");
	}

	FmcwIfResamplerPlan planFmcwIfResampler(const FmcwIfResamplerRequest& request)
	{
		const auto ratio =
			approximateFmcwIfRateRatio(request.output_sample_rate_hz, request.input_sample_rate_hz, request.limits);
		validateLimit(std::isfinite(request.filter_bandwidth_hz) && request.filter_bandwidth_hz > 0.0,
					  "IF resampler filter bandwidth must be positive");
		validateLimit(request.filter_bandwidth_hz < request.output_sample_rate_hz * 0.5,
					  "IF resampler filter bandwidth must be below output Nyquist");
		validateLimit(request.limits.max_taps_per_stage > 0, "IF resampler max taps per stage must be positive");
		validateLimit(request.limits.max_macs_per_output_sample > 0.0, "IF resampler max MAC budget must be positive");
		validateLimit(request.limits.max_phase_refinement > 0, "IF resampler phase refinement limit must be positive");

		const RealType transition_hz = derivedTransitionWidth(request);

		FmcwIfResamplerPlan plan;
		plan.input_sample_rate_hz = request.input_sample_rate_hz;
		plan.requested_output_sample_rate_hz = request.output_sample_rate_hz;
		plan.actual_output_sample_rate_hz = request.input_sample_rate_hz * ratio.actual_ratio;
		plan.filter_bandwidth_hz = request.filter_bandwidth_hz;
		plan.filter_transition_width_hz = transition_hz;
		plan.stopband_attenuation_db = request.stopband_attenuation_db;
		plan.overall_ratio = ratio;
		plan.limits = request.limits;

		auto remaining_up = ratio.numerator;
		auto remaining_down = ratio.denominator;
		RealType current_rate_hz = request.input_sample_rate_hz;

		while (remaining_down % 2 == 0 && remaining_down > 1)
		{
			if (request.filter_bandwidth_hz > current_rate_hz * kHalfBandPassbandFraction)
			{
				break;
			}

			const RealType available_transition = current_rate_hz * 0.25 - request.filter_bandwidth_hz;
			if (available_transition <= 0.0)
			{
				break;
			}

			const RealType halfband_transition =
				std::min(current_rate_hz * kHalfBandTransitionFraction, available_transition);
			plan.stages.push_back(makeStage(FmcwIfResamplerStageKind::HalfBandDecimateBy2, current_rate_hz, 1, 2,
											request.filter_bandwidth_hz, halfband_transition,
											request.stopband_attenuation_db, request.limits));
			current_rate_hz *= 0.5;
			remaining_down /= 2;
		}

		if (remaining_up != remaining_down)
		{
			plan.stages.push_back(makeStage(FmcwIfResamplerStageKind::RationalPolyphase, current_rate_hz, remaining_up,
											remaining_down, request.filter_bandwidth_hz, transition_hz,
											request.stopband_attenuation_db, request.limits));
		}

		for (const auto& stage : plan.stages)
		{
			plan.group_delay_seconds += stage.group_delay_seconds;
		}
		recomputePlanCost(plan);

		if (plan.estimated_macs_per_output_sample > request.limits.max_macs_per_output_sample)
		{
			throw std::runtime_error(
				"IF resampler estimated MAC cost " + std::to_string(plan.estimated_macs_per_output_sample) +
				" exceeds maximum " + std::to_string(request.limits.max_macs_per_output_sample) +
				"; increase if_sample_rate, reduce if_filter_bandwidth, or increase transition width");
		}

		plan.group_delay_output_samples = plan.group_delay_seconds * plan.actual_output_sample_rate_hz;
		plan.warmup_discard_samples = static_cast<std::uint64_t>(std::floor(plan.group_delay_output_samples));
		plan.fractional_output_delay_samples =
			plan.group_delay_output_samples - static_cast<RealType>(plan.warmup_discard_samples);
		plan.fractional_phase_offset = plan.fractional_output_delay_samples * static_cast<RealType>(ratio.denominator);
		plan.branch_interpolation_fraction =
			plan.fractional_phase_offset * static_cast<RealType>(plan.phase_refinement) -
			std::floor(plan.fractional_phase_offset * static_cast<RealType>(plan.phase_refinement));
		plan.estimated_timing_error_seconds =
			0.5 / (plan.actual_output_sample_rate_hz * static_cast<RealType>(plan.phase_refinement));
		plan.estimated_phase_error_radians =
			2.0 * PI * request.filter_bandwidth_hz * plan.estimated_timing_error_seconds;
		configureFractionalDelay(plan);

		return plan;
	}

	class FmcwIfResamplingSink::Stage
	{
	public:
		explicit Stage(FmcwIfResamplerStagePlan plan) : _plan(plan) { buildPhaseTable(); }

		struct ZeroInputResult
		{
			std::vector<ComplexType> emitted;
			std::size_t skipped_output_samples = 0;
		};

		void consume(std::span<const ComplexType> block)
		{
			if (_finished)
			{
				throw std::logic_error("Cannot consume IF resampler input after finish");
			}
			_input.insert(_input.end(), block.begin(), block.end());
			_input_total += block.size();
			process(false);
		}

		[[nodiscard]] ZeroInputResult consumeZeroInput(std::size_t input_count)
		{
			if (_finished)
			{
				throw std::logic_error("Cannot consume IF resampler input after finish");
			}

			ZeroInputResult result;
			if (input_count == 0)
			{
				return result;
			}

			const auto drained = std::min(input_count, zeroDrainInputCount());
			if (drained > 0)
			{
				_input.resize(_input.size() + drained, ComplexType{0.0, 0.0});
				_input_total += drained;
				input_count -= drained;
				process(false);
				result.emitted = take();
			}

			if (input_count == 0)
			{
				return result;
			}

			_input_total += input_count;
			const auto available_outputs = availableOutputCount(_input_total);
			if (available_outputs > _next_output_index)
			{
				result.skipped_output_samples = available_outputs - _next_output_index;
				_next_output_index = available_outputs;
			}
			keepTrailingZeroHistory();
			return result;
		}

		[[nodiscard]] std::vector<ComplexType> take()
		{
			std::vector<ComplexType> out;
			out.swap(_pending);
			return out;
		}

		[[nodiscard]] std::vector<ComplexType> finish()
		{
			if (!_finished)
			{
				_finished = true;
				process(true);
				_input.clear();
			}
			return take();
		}

		void reset()
		{
			_input.clear();
			_pending.clear();
			_input_total = 0;
			_input_start_index = 0;
			_next_output_index = 0;
			_finished = false;
		}

	private:
		[[nodiscard]] std::size_t zeroDrainInputCount() const noexcept { return _plan.tap_count + 2; }

		[[nodiscard]] std::size_t zeroHistoryKeepCount() const noexcept { return _plan.tap_count + 2; }

		void buildPhaseTable()
		{
			const std::size_t taps = _plan.tap_count;
			const std::size_t phase_count = std::max<std::size_t>(1, _plan.phase_count);
			_phase_table.assign(phase_count * taps, 0.0);

			const RealType beta = kaiserBeta(_plan.stopband_attenuation_db);
			const RealType bessel_beta = besselI0(beta);
			const RealType limiting_nyquist = 0.5 * std::min(_plan.input_sample_rate_hz, _plan.output_sample_rate_hz);
			const RealType cutoff_hz = std::min(_plan.cutoff_hz, limiting_nyquist * 0.999);
			const RealType cutoff = cutoff_hz / _plan.input_sample_rate_hz;
			const auto half = static_cast<RealType>((taps - 1) / 2);

			for (std::size_t phase = 0; phase < phase_count; ++phase)
			{
				const RealType frac = static_cast<RealType>(phase) / static_cast<RealType>(phase_count);
				RealType sum = 0.0;
				for (std::size_t i = 0; i < taps; ++i)
				{
					const RealType offset = frac - (static_cast<RealType>(i) - half);
					const RealType window = kaiserWindow(i, taps, beta, bessel_beta);
					const RealType coeff = 2.0 * cutoff * sinc(2.0 * cutoff * offset) * window;
					_phase_table[phase * taps + i] = coeff;
					sum += coeff;
				}
				if (std::abs(sum) > std::numeric_limits<RealType>::epsilon())
				{
					for (std::size_t i = 0; i < taps; ++i)
					{
						_phase_table[phase * taps + i] /= sum;
					}
				}
			}
		}

		[[nodiscard]] bool hasSamplesFor(const std::int64_t base_index, const std::int64_t input_offset,
										 const bool final) const noexcept
		{
			const auto half = static_cast<std::int64_t>((_plan.tap_count - 1) / 2);
			const std::int64_t last_needed = base_index + input_offset + half;
			return final || last_needed < static_cast<std::int64_t>(_input_total);
		}

		[[nodiscard]] ComplexType sampleAt(const std::int64_t global_index) const
		{
			if (global_index < 0 || global_index >= static_cast<std::int64_t>(_input_total))
			{
				return {0.0, 0.0};
			}
			const auto local = global_index - _input_start_index;
			if (local < 0 || local >= static_cast<std::int64_t>(_input.size()))
			{
				throw std::logic_error("IF resampler discarded input still needed by FIR branch");
			}
			return _input[static_cast<std::size_t>(local)];
		}

		[[nodiscard]] ComplexType evaluateBranch(const std::size_t phase, const std::int64_t base_index) const
		{
			const auto half = static_cast<std::int64_t>((_plan.tap_count - 1) / 2);
			const auto* coeffs = _phase_table.data() + phase * _plan.tap_count;
			ComplexType result{0.0, 0.0};
			for (std::size_t i = 0; i < _plan.tap_count; ++i)
			{
				const auto sample_index = base_index + static_cast<std::int64_t>(i) - half;
				result += sampleAt(sample_index) * coeffs[i];
			}
			return result;
		}

		[[nodiscard]] long double stagePosition(const std::size_t output_index) const noexcept
		{
			const std::size_t phase_count = std::max<std::size_t>(1, _plan.phase_count);
			const long double base = static_cast<long double>(output_index) *
				static_cast<long double>(_plan.down_factor) / static_cast<long double>(_plan.up_factor);
			if (!_plan.applies_fractional_delay)
			{
				return base;
			}

			const long double initial_phase = (static_cast<long double>(_plan.initial_phase_accumulator) +
											   static_cast<long double>(_plan.initial_branch_interpolation_fraction)) /
				static_cast<long double>(phase_count);
			return base + static_cast<long double>(_plan.initial_input_advance) + initial_phase;
		}

		struct BranchPosition
		{
			std::int64_t base_index = 0;
			std::size_t lower_phase = 0;
			std::size_t upper_phase = 0;
			std::int64_t upper_offset = 0;
			RealType mu = 0.0;
		};

		[[nodiscard]] BranchPosition branchPosition(const std::size_t output_index) const noexcept
		{
			const std::size_t phase_count = std::max<std::size_t>(1, _plan.phase_count);
			const long double position = stagePosition(output_index);
			const auto base_index = static_cast<std::int64_t>(std::floor(position + 1.0e-12L));
			long double frac = position - static_cast<long double>(base_index);
			if (frac < 0.0L && frac > -1.0e-12L)
			{
				frac = 0.0L;
			}

			const long double phase_position = frac * static_cast<long double>(phase_count);
			auto lower_phase = static_cast<std::size_t>(std::floor(phase_position + 1.0e-12L));
			auto mu = static_cast<RealType>(phase_position - static_cast<long double>(lower_phase));
			if (lower_phase >= phase_count)
			{
				lower_phase = 0;
				mu = 0.0;
			}

			std::size_t upper_phase = lower_phase + 1;
			std::int64_t upper_offset = 0;
			if (upper_phase == phase_count)
			{
				upper_phase = 0;
				upper_offset = 1;
			}

			return {.base_index = base_index,
					.lower_phase = lower_phase,
					.upper_phase = upper_phase,
					.upper_offset = upper_offset,
					.mu = mu};
		}

		[[nodiscard]] bool outputHasSamplesForTotal(const std::size_t output_index,
													const std::size_t input_total) const noexcept
		{
			const auto branch = branchPosition(output_index);
			const auto half = static_cast<std::int64_t>((_plan.tap_count - 1) / 2);
			const std::int64_t last_needed = branch.base_index + branch.upper_offset + half;
			return last_needed < static_cast<std::int64_t>(input_total);
		}

		[[nodiscard]] std::size_t availableOutputCount(const std::size_t input_total) const noexcept
		{
			std::size_t low = 0;
			std::size_t high = ceilOutputCount(input_total, _plan.up_factor, _plan.down_factor);
			while (low < high)
			{
				const auto mid = low + (high - low) / 2;
				if (outputHasSamplesForTotal(mid, input_total))
				{
					low = mid + 1;
				}
				else
				{
					high = mid;
				}
			}
			return low;
		}

		void process(const bool final)
		{
			const std::size_t final_output_count =
				final ? ceilOutputCount(_input_total, _plan.up_factor, _plan.down_factor) : 0;

			while (!final || _next_output_index < final_output_count)
			{
				const auto branch = branchPosition(_next_output_index);
				if (!hasSamplesFor(branch.base_index, 0, final) ||
					!hasSamplesFor(branch.base_index, branch.upper_offset, final))
				{
					break;
				}

				const ComplexType lower = evaluateBranch(branch.lower_phase, branch.base_index);
				const ComplexType upper = evaluateBranch(branch.upper_phase, branch.base_index + branch.upper_offset);
				_pending.push_back((1.0 - branch.mu) * lower + branch.mu * upper);
				++_next_output_index;
			}

			discardUnneededInput();
		}

		void discardUnneededInput()
		{
			if (_input.empty())
			{
				return;
			}

			const long double next_position = stagePosition(_next_output_index);
			const auto half = static_cast<std::int64_t>((_plan.tap_count - 1) / 2);
			const auto first_needed = static_cast<std::int64_t>(std::floor(next_position + 1.0e-12L)) - half;
			const auto new_start = std::max<std::int64_t>(_input_start_index, first_needed);
			const auto discard = new_start - _input_start_index;
			if (discard <= 0)
			{
				return;
			}

			const auto count = std::min<std::size_t>(static_cast<std::size_t>(discard), _input.size());
			_input.erase(_input.begin(), _input.begin() + static_cast<std::ptrdiff_t>(count));
			_input_start_index += static_cast<std::int64_t>(count);
		}

		void keepTrailingZeroHistory()
		{
			const auto keep = std::min<std::size_t>(zeroHistoryKeepCount(), _input_total);
			_input.assign(keep, ComplexType{0.0, 0.0});
			_input_start_index = static_cast<std::int64_t>(_input_total - keep);
		}

		FmcwIfResamplerStagePlan _plan;
		std::vector<RealType> _phase_table;
		std::vector<ComplexType> _input;
		std::vector<ComplexType> _pending;
		std::size_t _input_total = 0;
		std::int64_t _input_start_index = 0;
		std::size_t _next_output_index = 0;
		bool _finished = false;
	};

	FmcwIfResamplingSink::FmcwIfResamplingSink(FmcwIfResamplerPlan plan) : _plan(std::move(plan))
	{
		for (const auto& stage : _plan.stages)
		{
			_stages.push_back(std::make_unique<Stage>(stage));
		}
	}

	FmcwIfResamplingSink::~FmcwIfResamplingSink() = default;
	FmcwIfResamplingSink::FmcwIfResamplingSink(FmcwIfResamplingSink&&) noexcept = default;
	FmcwIfResamplingSink& FmcwIfResamplingSink::operator=(FmcwIfResamplingSink&&) noexcept = default;

	void FmcwIfResamplingSink::consume(std::span<const ComplexType> block)
	{
		if (_finished)
		{
			throw std::logic_error("Cannot consume IF resampler input after finish");
		}

		if (_stages.empty())
		{
			_output.insert(_output.end(), block.begin(), block.end());
			return;
		}

		std::vector<ComplexType> current(block.begin(), block.end());
		for (auto& stage : _stages)
		{
			stage->consume(current);
			current = stage->take();
		}
		_output.insert(_output.end(), current.begin(), current.end());
	}

	FmcwIfZeroInputResult FmcwIfResamplingSink::consumeZeroInput(const std::size_t input_count)
	{
		if (_finished)
		{
			throw std::logic_error("Cannot consume IF resampler input after finish");
		}

		FmcwIfZeroInputResult result;
		if (input_count == 0)
		{
			return result;
		}

		if (_stages.empty())
		{
			result.emitted = takeOutput();
			result.skipped_output_samples = input_count;
			return result;
		}

		std::vector<ComplexType> current;
		std::size_t zero_count = input_count;
		for (auto& stage : _stages)
		{
			if (!current.empty())
			{
				stage->consume(current);
				current = stage->take();
			}
			if (zero_count == 0)
			{
				continue;
			}

			auto zero_result = stage->consumeZeroInput(zero_count);
			current.insert(current.end(), zero_result.emitted.begin(), zero_result.emitted.end());
			zero_count = zero_result.skipped_output_samples;
		}

		_output.insert(_output.end(), current.begin(), current.end());
		result.emitted = takeOutput();
		result.skipped_output_samples = zero_count;
		return result;
	}

	std::vector<ComplexType> FmcwIfResamplingSink::takeOutput()
	{
		std::vector<ComplexType> out;
		out.swap(_output);
		return out;
	}

	std::vector<ComplexType> FmcwIfResamplingSink::finish()
	{
		if (_finished)
		{
			return takeOutput();
		}
		_finished = true;

		if (_stages.empty())
		{
			return takeOutput();
		}

		std::vector<ComplexType> current;
		for (std::size_t i = 0; i < _stages.size(); ++i)
		{
			if (i > 0)
			{
				_stages[i]->consume(current);
			}
			current = _stages[i]->finish();
		}
		_output.insert(_output.end(), current.begin(), current.end());
		return takeOutput();
	}

	void FmcwIfResamplingSink::reset()
	{
		for (auto& stage : _stages)
		{
			stage->reset();
		}
		_output.clear();
		_finished = false;
	}
}
