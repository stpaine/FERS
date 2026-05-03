#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <optional>
#include <vector>

#include "signal/if_resampler.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace
{
	[[nodiscard]] std::vector<ComplexType> complexTone(const RealType sample_rate_hz, const RealType tone_hz,
													   const std::size_t count)
	{
		std::vector<ComplexType> data(count);
		for (std::size_t i = 0; i < data.size(); ++i)
		{
			const RealType phase = 2.0 * PI * tone_hz * static_cast<RealType>(i) / sample_rate_hz;
			data[i] = {std::cos(phase), std::sin(phase)};
		}
		return data;
	}
}

TEST_CASE("FMCW IF resampler plans rational 5/512 conversion", "[signal][if_resampler][plan]")
{
	const fers_signal::FmcwIfResamplerRequest request{.input_sample_rate_hz = 1.024e9,
													  .output_sample_rate_hz = 10.0e6,
													  .filter_bandwidth_hz = 4.0e6,
													  .filter_transition_width_hz = std::nullopt};

	const auto plan = fers_signal::planFmcwIfResampler(request);

	REQUIRE(plan.overall_ratio.numerator == 5);
	REQUIRE(plan.overall_ratio.denominator == 512);
	REQUIRE_THAT(plan.actual_output_sample_rate_hz, WithinAbs(10.0e6, 1e-6));
	REQUIRE(plan.stages.size() == 7);
	for (std::size_t i = 0; i < 6; ++i)
	{
		REQUIRE(plan.stages[i].kind == fers_signal::FmcwIfResamplerStageKind::HalfBandDecimateBy2);
		REQUIRE(plan.stages[i].down_factor == 2);
	}
	const auto& final_stage = plan.stages.back();
	REQUIRE(final_stage.kind == fers_signal::FmcwIfResamplerStageKind::RationalPolyphase);
	REQUIRE(final_stage.up_factor == 5);
	REQUIRE(final_stage.down_factor == 8);
	REQUIRE(final_stage.phase_refinement == request.limits.max_phase_refinement);
	REQUIRE(final_stage.phase_count == 5 * request.limits.max_phase_refinement);
	REQUIRE(plan.group_delay_seconds > 0.0);
	REQUIRE(plan.group_delay_output_samples > 0.0);
	REQUIRE(plan.warmup_discard_samples > 0);
	REQUIRE(plan.fractional_output_delay_samples >= 0.0);
	REQUIRE(plan.fractional_output_delay_samples < 1.0);
	REQUIRE(plan.estimated_macs_per_output_sample <= request.limits.max_macs_per_output_sample);
}

TEST_CASE("FMCW IF resampler rejects excessive near-Nyquist filter cost", "[signal][if_resampler][plan]")
{
	const fers_signal::FmcwIfResamplerRequest request{.input_sample_rate_hz = 1.024e9,
													  .output_sample_rate_hz = 10.0e6,
													  .filter_bandwidth_hz = 4.999e6,
													  .filter_transition_width_hz = std::nullopt};

	REQUIRE_THROWS_WITH(fers_signal::planFmcwIfResampler(request), ContainsSubstring("FIR tap count"));
}

TEST_CASE("FMCW IF streaming sink resamples a complex tone by blocks", "[signal][if_resampler][stream]")
{
	constexpr RealType input_rate_hz = 1024.0;
	constexpr RealType output_rate_hz = 1000.0;
	constexpr RealType tone_hz = 40.0;
	constexpr std::size_t input_count = 4096;

	const fers_signal::FmcwIfResamplerRequest request{.input_sample_rate_hz = input_rate_hz,
													  .output_sample_rate_hz = output_rate_hz,
													  .filter_bandwidth_hz = 120.0,
													  .filter_transition_width_hz = 100.0};
	const auto plan = fers_signal::planFmcwIfResampler(request);
	fers_signal::FmcwIfResamplingSink sink(plan);

	const auto input = complexTone(input_rate_hz, tone_hz, input_count);
	std::vector<ComplexType> output;
	for (std::size_t offset = 0; offset < input.size(); offset += 257)
	{
		const auto count = std::min<std::size_t>(257, input.size() - offset);
		sink.consume(std::span<const ComplexType>{input.data() + offset, count});
		auto chunk = sink.takeOutput();
		output.insert(output.end(), chunk.begin(), chunk.end());
	}
	auto tail = sink.finish();
	output.insert(output.end(), tail.begin(), tail.end());

	REQUIRE(output.size() == 4000);

	ComplexType correlation{0.0, 0.0};
	RealType output_energy = 0.0;
	RealType reference_energy = 0.0;
	for (std::size_t i = 600; i < output.size() - 600; ++i)
	{
		const RealType phase = 2.0 * PI * tone_hz * static_cast<RealType>(i) / output_rate_hz;
		const ComplexType expected{std::cos(phase), std::sin(phase)};
		correlation += output[i] * std::conj(expected);
		output_energy += std::norm(output[i]);
		reference_energy += std::norm(expected);
	}

	const RealType normalized_correlation = std::abs(correlation) / std::sqrt(output_energy * reference_energy);
	const RealType gain = std::sqrt(output_energy / reference_energy);
	REQUIRE(normalized_correlation > 0.995);
	REQUIRE(gain > 0.95);
	REQUIRE(gain < 1.05);
}
