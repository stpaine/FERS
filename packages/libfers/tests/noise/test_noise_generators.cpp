#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

#include "core/config.h"
#include "core/parameters.h"
#include "noise/falpha_branch.h"
#include "noise/noise_generators.h"

using Catch::Matchers::WithinAbs;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
	};

	RealType mean(const std::vector<RealType>& values)
	{
		RealType sum = 0.0;
		for (const auto value : values)
		{
			sum += value;
		}
		return sum / static_cast<RealType>(values.size());
	}

	RealType variance(const std::vector<RealType>& values, const RealType mean_value)
	{
		RealType sum = 0.0;
		for (const auto value : values)
		{
			const auto diff = value - mean_value;
			sum += diff * diff;
		}
		return sum / static_cast<RealType>(values.size() - 1);
	}
}

TEST_CASE("WgnGenerator produces expected mean and variance", "[noise][wgn]")
{
	std::mt19937 rng(12345);
	constexpr RealType stddev = 2.0;
	noise::WgnGenerator generator(rng, stddev);

	constexpr size_t sample_count = 20000;
	std::vector<RealType> samples(sample_count);
	for (auto& value : samples)
	{
		value = generator.getSample();
	}

	const auto sample_mean = mean(samples);
	const auto sample_var = variance(samples, sample_mean);

	REQUIRE_THAT(sample_mean, WithinAbs(0.0, 0.1));
	REQUIRE_THAT(sample_var, WithinAbs(stddev * stddev, 0.25));
}

TEST_CASE("GammaGenerator matches theoretical mean and variance", "[noise][gamma]")
{
	std::mt19937 rng(54321);
	constexpr RealType shape = 2.0;
	noise::GammaGenerator generator(rng, shape);

	constexpr size_t sample_count = 30000;
	std::vector<RealType> samples(sample_count);
	RealType min_sample = std::numeric_limits<RealType>::max();
	for (auto& value : samples)
	{
		value = generator.getSample();
		min_sample = std::min(min_sample, value);
	}

	const auto sample_mean = mean(samples);
	const auto sample_var = variance(samples, sample_mean);

	REQUIRE(min_sample >= 0.0);
	REQUIRE_THAT(sample_mean, WithinAbs(shape, 0.1));
	REQUIRE_THAT(sample_var, WithinAbs(shape, 0.15));
}

TEST_CASE("MultirateGenerator applies scaling based on alpha", "[noise][multirate]")
{
	std::mt19937 rng_actual(2021);
	std::mt19937 rng_expected(2021);

	constexpr RealType alpha = 1.0;
	constexpr unsigned branches = 1;

	noise::MultirateGenerator generator(rng_actual, alpha, branches);
	noise::FAlphaBranch reference_branch(rng_expected, 0.5, 0, nullptr, true);

	const RealType expected_scale = 1.0 / std::pow(10.0, (-alpha + 2.0) * 2.0);
	const RealType expected_sample = reference_branch.getSample() * expected_scale;

	REQUIRE_THAT(generator.getSample(), WithinAbs(expected_sample, 1e-12));
}

TEST_CASE("MultirateGenerator skipSamples handles short and long skips", "[noise][multirate]")
{
	std::mt19937 rng(7);
	noise::MultirateGenerator generator(rng, 0.0, 2);

	const auto first = generator.getSample();
	REQUIRE(std::isfinite(first));

	generator.skipSamples(10);
	const auto after_short_skip = generator.getSample();
	REQUIRE(std::isfinite(after_short_skip));

	generator.skipSamples(1000);
	const auto after_long_skip = generator.getSample();
	REQUIRE(std::isfinite(after_long_skip));

	generator.reset();
	const auto after_reset = generator.getSample();
	REQUIRE(std::isfinite(after_reset));
}

TEST_CASE("ClockModelGenerator applies offsets and respects sample count", "[noise][clock]")
{
	ParamGuard guard;
	params::params.rate = 100.0;

	std::mt19937 rng(11);
	const std::vector<RealType> alpha;
	const std::vector<RealType> weights;

	noise::ClockModelGenerator generator(rng, alpha, weights, 1.0, 1.5, 2.0, 1);
	REQUIRE(generator.enabled());

	const RealType sample0 = generator.getSample();
	const RealType sample1 = generator.getSample();

	REQUIRE_THAT(sample0, WithinAbs(1.5, 1e-12));
	REQUIRE_THAT(sample1, WithinAbs(1.5 + 2.0 * PI * 2.0 / params::rate(), 1e-12));

	generator.skipSamples(8);
	const RealType sample_after_skip = generator.getSample();
	REQUIRE_THAT(sample_after_skip, WithinAbs(1.5 + 2.0 * PI * 2.0 * 10.0 / params::rate(), 1e-12));

	generator.reset();
	const RealType sample_after_reset = generator.getSample();
	REQUIRE_THAT(sample_after_reset, WithinAbs(1.5, 1e-12));
}

TEST_CASE("ClockModelGenerator handles weight scaling branches", "[noise][clock]")
{
	ParamGuard guard;
	params::params.rate = 100.0;

	std::mt19937 rng(29);
	const std::vector<RealType> alpha = {2.0, 1.0, 0.0, -1.0, -2.0};
	const std::vector<RealType> weights = {0.0, 0.0, 0.0, 0.0, 0.0};

	noise::ClockModelGenerator generator(rng, alpha, weights, 1.0, 0.0, 0.0, 1);
	REQUIRE(generator.enabled());

	const RealType sample = generator.getSample();
	REQUIRE_THAT(sample, WithinAbs(0.0, 1e-12));
}
