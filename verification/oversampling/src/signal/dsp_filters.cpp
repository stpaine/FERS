// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2007-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// Upstream reference:
//   packages/libfers/src/signal/dsp_filters.cpp

#include "dsp_filters.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <numeric>
#include <ranges>
#include <set>
#include <stdexcept>

#include "core/logging.h"
#include "core/parameters.h"

constexpr RealType BLACKMAN_A0 = 0.42;
constexpr RealType BLACKMAN_A1 = 0.5;
constexpr RealType BLACKMAN_A2 = 0.08;

namespace
{
	constexpr RealType sinc(const RealType x) noexcept { return x == 0 ? 1.0 : std::sin(x * PI) / (x * PI); }

	struct BlackmanFirDesign
	{
		std::vector<RealType> coeffs;
		RealType coeff_sum = 0.0;
	};

	void warnIfUnderspecified(const unsigned ratio, const unsigned filter_length, const RealType coeff_sum) noexcept
	{
		if (ratio <= 1)
		{
			return;
		}

		const RealType relative_error = std::abs(coeff_sum - static_cast<RealType>(ratio)) / static_cast<RealType>(ratio);
		if (relative_error <= 0.01)
		{
			return;
		}

		static std::set<std::pair<unsigned, unsigned>> warned_pairs;
		const auto [_, inserted] = warned_pairs.emplace(ratio, filter_length);
		if (!inserted)
		{
			return;
		}

		const RealType dc_gain = coeff_sum / static_cast<RealType>(ratio);
		const RealType roundtrip_gain = dc_gain * dc_gain;
		LOG(logging::Level::WARNING,
			"Oversampling FIR under-spec'd for ratio=", ratio,
			", filter_length=", filter_length,
			", fir_sum=", coeff_sum,
			", estimated_roundtrip_gain=", roundtrip_gain,
			". Practical envelope roughly filter_length >= 4 * ratio.");
	}

	BlackmanFirDesign blackmanFir(const RealType cutoff, const unsigned ratio, const unsigned filter_length) noexcept
	{
		const unsigned filt_length = filter_length * 2;
		std::vector<RealType> coeffs(filt_length);
		const RealType n = filt_length / 2.0;
		const RealType pi_n = PI / n;

		std::ranges::for_each(coeffs,
							  [cutoff, n, pi_n, i = 0u](RealType& coeff) mutable
							  {
								  const RealType sinc_val = sinc(cutoff * (i - n));
								  const RealType window = BLACKMAN_A0 - BLACKMAN_A1 * std::cos(pi_n * i) +
									  BLACKMAN_A2 * std::cos(2 * pi_n * i);
									  coeff = sinc_val * window;
									  ++i;
								  });

		BlackmanFirDesign design{std::move(coeffs)};
		design.coeff_sum = std::accumulate(design.coeffs.begin(), design.coeffs.end(), 0.0);
		warnIfUnderspecified(ratio, filter_length, design.coeff_sum);
		return design;
	}
}

namespace fers_signal
{
	void upsample(const std::span<const ComplexType> in, const unsigned size, std::span<ComplexType> out)
	{
		const unsigned ratio = params::oversampleRatio();
		const unsigned filter_length = params::renderFilterLength();
		const auto design = blackmanFir(1 / static_cast<RealType>(ratio), ratio, filter_length);
		const unsigned filt_length = static_cast<unsigned>(design.coeffs.size());

		std::vector tmp(static_cast<size_t>(size * ratio + filt_length), ComplexType{0.0, 0.0});

		for (unsigned i = 0; i < size; ++i)
		{
			tmp[static_cast<size_t>(i * ratio)] = in[i];
		}

		const FirFilter filt(design.coeffs);
		filt.filter(tmp);

		const auto delay = filt_length / 2;
		std::ranges::copy_n(tmp.begin() + delay, size * ratio, out.begin());
	}

	std::vector<ComplexType> downsample(std::span<const ComplexType> in)
	{
		if (in.empty())
		{
			throw std::invalid_argument("Input span is empty in Downsample");
		}

		const unsigned ratio = params::oversampleRatio();
		const unsigned filter_length = params::renderFilterLength();
		const auto design = blackmanFir(1 / static_cast<RealType>(ratio), ratio, filter_length);
		const unsigned filt_length = static_cast<unsigned>(design.coeffs.size());

		std::vector tmp(in.size() + filt_length, ComplexType{0, 0});
		std::ranges::copy(in, tmp.begin());

		const FirFilter filt(design.coeffs);
		filt.filter(tmp);

		const auto downsampled_size = in.size() / ratio;
		std::vector<ComplexType> out(downsampled_size);
		for (unsigned i = 0; i < downsampled_size; ++i)
		{
			out[i] = tmp[static_cast<size_t>(i * ratio + filt_length / 2)] / static_cast<RealType>(ratio);
		}

		return out;
	}

	IirFilter::IirFilter(const RealType* denCoeffs, const RealType* numCoeffs, const unsigned order) noexcept :
		_a(denCoeffs, denCoeffs + order), _b(numCoeffs, numCoeffs + order), _w(order, 0.0), _order(order)
	{
	}

	RealType IirFilter::filter(const RealType sample) noexcept
	{
		std::ranges::rotate(_w, _w.end() - 1);
		_w[0] = sample;

		for (unsigned j = 1; j < _order; ++j)
		{
			_w[0] -= _a[j] * _w[j];
		}

		return std::inner_product(_b.begin(), _b.end(), _w.begin(), 0.0);
	}

	void IirFilter::filter(std::span<RealType> samples) noexcept
	{
		for (auto& sample : samples)
		{
			std::ranges::rotate(_w, _w.end() - 1);
			_w[0] = sample;

			for (unsigned j = 1; j < _order; ++j)
			{
				_w[0] -= _a[j] * _w[j];
			}

			sample = std::inner_product(_b.begin(), _b.end(), _w.begin(), 0.0);
		}
	}

	void FirFilter::filter(std::vector<ComplexType>& samples) const
	{
		std::vector line(_order, ComplexType{0.0, 0.0});

		for (auto& sample : samples)
		{
			line[0] = sample;
			ComplexType result{0.0, 0.0};

			result = std::transform_reduce(line.begin(), line.end(), _filter.begin(), ComplexType{0.0, 0.0},
										   std::plus<ComplexType>{},
										   [](const ComplexType& x, const RealType coeff) { return x * coeff; });

			sample = result;
			std::ranges::rotate(std::views::reverse(line), line.rbegin() + 1);
		}
	}

	DecadeUpsampler::DecadeUpsampler()
	{
		constexpr std::array den_coeffs = {1.0,
										   -10.301102119865,
										   48.5214567642597,
										   -137.934509572412,
										   262.914952985445,
										   -352.788381841481,
										   340.027874008585,
										   -235.39260470286,
										   114.698499845697,
										   -37.4634653062448,
										   7.38208765922137,
										   -0.664807695826097};

		constexpr std::array num_coeffs = {2.7301694322809e-06,	 -1.8508123430239e-05, 5.75739466753894e-05,
										   -0.000104348734423658, 0.000111949190289715,  -4.9384188225528e-05,
										   -4.9384188225522e-05,  0.00011194919028971,   -0.000104348734423656,
										   5.75739466753884e-05,  -1.85081234302388e-05, 2.73016943228086e-06};

		_filter = std::make_unique<IirFilter>(den_coeffs.data(), num_coeffs.data(), den_coeffs.size());
	}

	void DecadeUpsampler::upsample(const RealType sample, std::span<RealType> out) const
	{
		if (out.size() != 10)
		{
			throw std::invalid_argument("Output span must have a size of 10.");
		}
		out[0] = sample;
		std::fill(out.begin() + 1, out.end(), 0);
		_filter->filter(out);
	}
}
