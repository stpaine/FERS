#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <random>
#include <vector>

#include "core/config.h"
#include "noise/falpha_branch.h"
#include "signal/dsp_filters.h"

using Catch::Matchers::WithinAbs;

namespace
{
	constexpr std::array<RealType, 12> kHighpassNum = {
		3.817871081981451e-01, -4.093384095523618e+00, 2.005300512623078e+01, -5.924672881811163e+01,
		1.172948159891025e+02, -1.633810410083022e+02, 1.633810410083034e+02, -1.172948159891052e+02,
		5.924672881811390e+01, -2.005300512623186e+01, 4.093384095523903e+00, -3.817871081981776e-01,
	};

	constexpr std::array<RealType, 12> kHighpassDen = {
		1.000000000000000e+00, -8.829695665523831e+00, 3.583068809011030e+01, -8.811479652970442e+01,
		1.457874067329429e+02, -1.702715637111961e+02, 1.431504350055831e+02, -8.656925883534657e+01,
		3.687395592491803e+01, -1.052413841411803e+01, 1.808292123637038e+00, -1.412932578340511e-01,
	};

	constexpr std::array<RealType, 16> kShapeNum = {
		5.210373977738306e-03,	-7.694671394585578e-03, 1.635979377907092e-03,	9.852449140857658e-05,
		-2.080553126780113e-03, 4.088764157029523e-03,	-1.549082440084623e-03, 9.054734252370680e-04,
		-3.467369912368729e-04, 4.516383087838856e-04,	-1.063356106118517e-03, 1.330008998057684e-04,
		6.556909567323943e-04,	-4.839476350293955e-04, 6.664936170526832e-05,	1.528520559763056e-05,
	};

	constexpr std::array<RealType, 16> kShapeDen = {
		1.000000000000000e+00,	-2.065565041154101e+00, 1.130909190864681e+00,	-1.671244644503288e-01,
		-3.331474931013877e-01, 9.952625337612708e-01,	-7.123036343635182e-01, 3.297062696290504e-01,
		-1.925691520710595e-01, 1.301247006176314e-01,	-2.702016290409912e-01, 1.455380885858886e-01,
		1.091921868353888e-01,	-1.524953111510459e-01, 5.667716332023935e-02,	-2.890314873767405e-03,
	};

	constexpr RealType kShapeGain = 5.210373977738306e-03;
}

TEST_CASE("FAlphaBranch rejects unsupported fractional integrator values", "[noise][falpha]")
{
	std::mt19937 rng(1);
	REQUIRE_THROWS_AS(noise::FAlphaBranch(rng, 0.25, 0, nullptr, true), std::runtime_error);
}

TEST_CASE("FAlphaBranch rejects unsupported integer integrator values", "[noise][falpha]")
{
	std::mt19937 rng(2);
	REQUIRE_THROWS_AS(noise::FAlphaBranch(rng, 0.0, 3, nullptr, true), std::runtime_error);
}

TEST_CASE("FAlphaBranch non-last branch matches upsampled pre-offset buffer", "[noise][falpha]")
{
	constexpr unsigned fint = 0;
	constexpr RealType ffrac = 0.0;
	const RealType upsample_scale = std::pow(10.0, ffrac + fint + 0.5);

	std::mt19937 rng_branch(77);
	std::mt19937 rng_manual(77);

	auto pre = std::make_unique<noise::FAlphaBranch>(rng_branch, 0.0, 0, nullptr, true);
	noise::FAlphaBranch branch(rng_branch, ffrac, fint, std::move(pre), false);

	fers_signal::IirFilter highpass(kHighpassDen.data(), kHighpassNum.data(), kHighpassDen.size());
	fers_signal::DecadeUpsampler upsampler;
	std::normal_distribution<> dist_main(0.0, 1.0);
	std::normal_distribution<> dist_pre(0.0, 1.0);

	std::vector<RealType> expected(20);
	RealType offset_sample = 0.0;
	for (size_t refill = 0; refill < 2; ++refill)
	{
		const RealType main_sample = dist_main(rng_manual);
		const RealType pre_sample = dist_pre(rng_manual);
		if (refill == 0)
		{
			offset_sample = pre_sample;
		}
		const RealType filtered = highpass.filter(main_sample);
		const RealType combined = (refill == 0) ? filtered : filtered + pre_sample - offset_sample;

		std::vector<RealType> buffer(10, 0.0);
		upsampler.upsample(combined, buffer);
		const size_t base_index = refill * buffer.size();
		for (size_t i = 0; i < buffer.size(); ++i)
		{
			buffer[i] = buffer[i] * upsample_scale + offset_sample;
			expected[base_index + i] = buffer[i];
		}
	}

	for (size_t i = 0; i < expected.size(); ++i)
	{
		REQUIRE_THAT(branch.getSample(), WithinAbs(expected[i], 1e-12));
	}
}

TEST_CASE("FAlphaBranch fractional shaping matches reference filter", "[noise][falpha]")
{
	std::mt19937 rng_branch(123);
	std::mt19937 rng_manual(123);

	noise::FAlphaBranch branch(rng_branch, 0.5, 0, nullptr, true);
	std::normal_distribution<> dist(0.0, 1.0);
	fers_signal::IirFilter shape_filter(kShapeDen.data(), kShapeNum.data(), kShapeDen.size());

	for (int i = 0; i < 12; ++i)
	{
		const RealType raw = dist(rng_manual);
		const RealType expected = shape_filter.filter(raw) / kShapeGain;
		REQUIRE_THAT(branch.getSample(), WithinAbs(expected, 1e-12));
	}
}

TEST_CASE("FAlphaBranch integer integration matches reference filter", "[noise][falpha]")
{
	struct Case
	{
		unsigned fint;
		std::vector<RealType> den;
		std::vector<RealType> num;
	};

	const std::vector<Case> cases = {
		{1, {1.0, -1.0}, {1.0, 0.0}},
		{2, {1.0, -2.0, 1.0}, {1.0, 0.0, 0.0}},
	};

	for (const auto& test_case : cases)
	{
		std::mt19937 rng_branch(456 + test_case.fint);
		std::mt19937 rng_manual(456 + test_case.fint);

		noise::FAlphaBranch branch(rng_branch, 0.0, test_case.fint, nullptr, true);
		std::normal_distribution<> dist(0.0, 1.0);
		fers_signal::IirFilter integ_filter(test_case.den.data(), test_case.num.data(),
											static_cast<unsigned>(test_case.den.size()));

		for (int i = 0; i < 8; ++i)
		{
			const RealType raw = dist(rng_manual);
			const RealType expected = integ_filter.filter(raw);
			REQUIRE_THAT(branch.getSample(), WithinAbs(expected, 1e-12));
		}
	}
}
