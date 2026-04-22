#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <vector>

#include "core/parameters.h"
#include "signal/dsp_filters.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
	};
}

TEST_CASE("IirFilter acts as identity with unity coefficients", "[signal][dsp][iir]")
{
	const RealType a[] = {1.0, 0.0};
	const RealType b[] = {1.0, 0.0};

	fers_signal::IirFilter filter(a, b, 2);
	std::vector<RealType> samples = {1.0, -2.5, 3.25, 0.0};

	filter.filter(samples);

	REQUIRE(samples.size() == 4);
	REQUIRE_THAT(samples[0], WithinAbs(1.0, 0.0));
	REQUIRE_THAT(samples[1], WithinAbs(-2.5, 0.0));
	REQUIRE_THAT(samples[2], WithinAbs(3.25, 0.0));
	REQUIRE_THAT(samples[3], WithinAbs(0.0, 0.0));
}

TEST_CASE("IirFilter feedback produces expected decay", "[signal][dsp][iir]")
{
	const RealType a[] = {1.0, -0.5};
	const RealType b[] = {1.0, 0.0};

	fers_signal::IirFilter filter(a, b, 2);
	std::vector<RealType> input = {1.0, 0.0, 0.0, 0.0};
	std::vector<RealType> expected = {1.0, 0.5, 0.25, 0.125};

	for (size_t i = 0; i < input.size(); ++i)
	{
		const auto out = filter.filter(input[i]);
		REQUIRE_THAT(out, WithinAbs(expected[i], 1e-12));
	}
}

TEST_CASE("FirFilter applies coefficients to complex samples", "[signal][dsp][fir]")
{
	const std::vector<RealType> coeffs = {1.0, 2.0};
	std::vector<ComplexType> samples = {ComplexType{1.0, 1.0}, ComplexType{0.0, 0.0}, ComplexType{0.0, 0.0}};

	fers_signal::FirFilter filter(coeffs);
	filter.filter(samples);

	REQUIRE(samples.size() == 3);
	REQUIRE_THAT(samples[0].real(), WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(samples[0].imag(), WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(samples[1].real(), WithinAbs(2.0, 1e-12));
	REQUIRE_THAT(samples[1].imag(), WithinAbs(2.0, 1e-12));
	REQUIRE_THAT(samples[2].real(), WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(samples[2].imag(), WithinAbs(0.0, 1e-12));
}

TEST_CASE("FirFilter scalar overload returns zero", "[signal][dsp][fir]")
{
	const std::vector<RealType> coeffs = {1.0, 2.0, 3.0};
	fers_signal::FirFilter filter(coeffs);

	REQUIRE_THAT(filter.filter(0.0), WithinAbs(0.0, 0.0));
	REQUIRE_THAT(filter.filter(1.0), WithinAbs(0.0, 0.0));
	REQUIRE_THAT(filter.filter(-7.5), WithinAbs(0.0, 0.0));
}

TEST_CASE("FirFilter real-span overload is a no-op", "[signal][dsp][fir]")
{
	const std::vector<RealType> coeffs = {1.0, 2.0, 3.0};
	fers_signal::FirFilter filter(coeffs);

	std::vector<RealType> samples = {1.0, -2.0, 3.5, 0.0};
	const auto expected = samples;
	filter.filter(samples);
	REQUIRE(samples == expected);

	std::vector<RealType> empty;
	filter.filter(empty);
	REQUIRE(empty.empty());
}

TEST_CASE("DecadeUpsampler rejects incorrect output size", "[signal][dsp][upsample]")
{
	fers_signal::DecadeUpsampler upsampler;
	std::vector<RealType> out(9, 0.0);
	REQUIRE_THROWS_AS(upsampler.upsample(1.0, out), std::invalid_argument);
}

TEST_CASE("DecadeUpsampler produces zeros for zero input", "[signal][dsp][upsample]")
{
	fers_signal::DecadeUpsampler upsampler;
	std::vector<RealType> out(10, 1.0);

	upsampler.upsample(0.0, out);

	for (auto value : out)
	{
		REQUIRE_THAT(value, WithinAbs(0.0, 1e-12));
	}
}

TEST_CASE("Upsample and downsample preserve low-frequency content", "[signal][dsp][resample]")
{
	ParamGuard guard;
	params::params.filter_length = 8;
	params::setOversampleRatio(2);

	constexpr size_t sample_count = 64;
	constexpr RealType frequency = 0.05;
	std::vector<ComplexType> input(sample_count);
	for (size_t i = 0; i < input.size(); ++i)
	{
		const RealType phase = 2.0 * PI * frequency * static_cast<RealType>(i);
		input[i] = ComplexType{std::cos(phase), std::sin(phase)};
	}

	std::vector<ComplexType> upsampled(sample_count * params::oversampleRatio());
	fers_signal::upsample(input, static_cast<unsigned>(input.size()), upsampled);

	auto downsampled = fers_signal::downsample(upsampled);
	REQUIRE(downsampled.size() == input.size());

	constexpr RealType tolerance = 1e-2;
	for (size_t i = 10; i < input.size() - 10; ++i)
	{
		const RealType error = std::abs(downsampled[i] - input[i]);
		REQUIRE(error < tolerance);
	}
}

TEST_CASE("Upsample and downsample preserve boundary-ratio low-frequency content", "[signal][dsp][resample]")
{
	ParamGuard guard;
	params::params.filter_length = 33;
	params::setOversampleRatio(8);

	constexpr size_t sample_count = 128;
	constexpr RealType frequency = 0.01;
	std::vector<ComplexType> input(sample_count);
	for (size_t i = 0; i < input.size(); ++i)
	{
		const RealType phase = 2.0 * PI * frequency * static_cast<RealType>(i);
		input[i] = ComplexType{std::cos(phase), std::sin(phase)};
	}

	std::vector<ComplexType> upsampled(sample_count * params::oversampleRatio());
	fers_signal::upsample(input, static_cast<unsigned>(input.size()), upsampled);

	auto downsampled = fers_signal::downsample(upsampled);
	REQUIRE(downsampled.size() == input.size());

	ComplexType correlation{0.0, 0.0};
	RealType output_energy = 0.0;
	RealType reference_energy = 0.0;
	for (size_t i = 20; i < input.size() - 20; ++i)
	{
		correlation += downsampled[i] * std::conj(input[i]);
		output_energy += std::norm(downsampled[i]);
		reference_energy += std::norm(input[i]);
	}

	const RealType normalized_correlation = std::abs(correlation) / std::sqrt(output_energy * reference_energy);
	const RealType gain = std::sqrt(output_energy / reference_energy);

	REQUIRE(normalized_correlation > 0.995);
	REQUIRE(gain > 0.99);
	REQUIRE(gain < 1.01);
}

TEST_CASE("Downsample truncates non-divisible oversampled buffers", "[signal][dsp][downsample]")
{
	ParamGuard guard;
	params::setOversampleRatio(8);

	const std::vector<ComplexType> oversampled(17, ComplexType{0.0, 0.0});
	const auto downsampled = fers_signal::downsample(oversampled);

	REQUIRE(downsampled.size() == oversampled.size() / params::oversampleRatio());
}

TEST_CASE("Resampler fails fast when oversample ratio exceeds fixed-filter limit", "[signal][dsp][resample]")
{
	ParamGuard guard;
	params::params.oversample_ratio = 9;

	const std::vector<ComplexType> input = {ComplexType{1.0, 0.0}, ComplexType{0.0, 1.0}};
	std::vector<ComplexType> upsampled(input.size() * params::params.oversample_ratio);

	REQUIRE_THROWS_WITH(fers_signal::upsample(input, static_cast<unsigned>(input.size()), upsampled),
						ContainsSubstring("Oversampling ratios > 8 are not supported"));
	REQUIRE_THROWS_WITH(fers_signal::downsample(upsampled),
						ContainsSubstring("Oversampling ratios > 8 are not supported"));
}

TEST_CASE("Downsample rejects empty input", "[signal][dsp][downsample]")
{
	const std::vector<ComplexType> empty;
	REQUIRE_THROWS_AS(fers_signal::downsample(empty), std::invalid_argument);
}
