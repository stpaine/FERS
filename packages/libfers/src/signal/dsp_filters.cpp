// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2007-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file dsp_filters.cpp
 * @brief Implementation file for Digital Signal Processing (DSP) filters and upsampling/downsampling functionality.
 */

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
	/**
	 * @brief Sinc function for FIR filter design.
	 *
	 * @param x Input value.
	 * @return Sinc value at x.
	 */
	constexpr RealType sinc(const RealType x) noexcept { return x == 0 ? 1.0 : std::sin(x * PI) / (x * PI); }

	/**
	 * @brief Stores a generated Blackman-windowed FIR design and its coefficient sum.
	 *
	 * The coefficient sum approximates the interpolator DC gain and is used to detect
	 * under-specified oversampling configurations.
	 */
	struct BlackmanFirDesign
	{
		std::vector<RealType> coeffs;
		RealType coeff_sum = 0.0;
	};

	/**
	 * @brief Enforces production oversampling limits for the fixed-length FIR design.
	 *
	 * @param ratio Requested oversampling ratio.
	 * @throws std::runtime_error if the configured ratio is unsupported.
	 */
	void validateOversamplingConfig(const unsigned ratio) { params::validateOversampleRatio(ratio); }

	/**
	 * @brief Emits a one-time warning when the FIR tap budget is too small for the ratio.
	 *
	 * @param ratio Oversampling ratio used to design the FIR.
	 * @param filter_length Configured half-length of the production FIR kernel.
	 * @param coeff_sum Sum of the generated FIR coefficients.
	 */
	void warnIfUnderspecified(const unsigned ratio, const unsigned filter_length, const RealType coeff_sum)
	{
		if (ratio <= 1)
		{
			return;
		}

		const RealType relative_error =
			std::abs(coeff_sum - static_cast<RealType>(ratio)) / static_cast<RealType>(ratio);
		if (relative_error <= 0.01)
		{
			return;
		}

		// TODO: warning deduplication via static set has data race potential
		static std::set<std::pair<unsigned, unsigned>> warned_pairs;
		const auto [_, inserted] = warned_pairs.emplace(ratio, filter_length);
		if (!inserted)
		{
			return;
		}

		const RealType dc_gain = coeff_sum / static_cast<RealType>(ratio);
		const RealType roundtrip_gain = dc_gain * dc_gain;
		LOG(logging::Level::WARNING, "Oversampling FIR under-spec'd for ratio=", ratio,
			", filter_length=", filter_length, ", fir_sum=", coeff_sum, ", estimated_roundtrip_gain=", roundtrip_gain,
			". Practical envelope roughly filter_length >= 4 * ratio.");
	}

	/**
	 * @brief Generates a Blackman-windowed low-pass FIR for the oversampling stage.
	 *
	 * @param cutoff Normalized cutoff frequency for the anti-imaging / anti-alias filter.
	 * @param ratio Oversampling ratio associated with this design.
	 * @param filter_length Configured half-length of the FIR kernel.
	 * @return FIR coefficients plus their sum for DC-gain auditing.
	 */
	BlackmanFirDesign blackmanFir(const RealType cutoff, const unsigned ratio, const unsigned filter_length)
	{
		const unsigned filt_length = filter_length * 2;
		std::vector<RealType> coeffs(filt_length);
		const RealType n = filt_length / 2.0;
		const RealType pi_n = PI / n;

		// We use the Blackman window, for a suitable tradeoff between rolloff and stopband attenuation
		// Equivalent Kaiser beta = 7.04 (Oppenhiem and Schaffer, Hamming)
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
	/**
	 * @brief Upsamples a complex waveform with zero-stuffing followed by Blackman FIR filtering.
	 *
	 * The production path uses a fixed-length single-stage FIR, so unsupported ratios
	 * fail fast before filtering begins.
	 *
	 * @param in Input span of base-rate complex samples.
	 * @param size Number of input samples to process from @p in.
	 * @param out Output span receiving @p size * oversampleRatio() filtered samples.
	 * @throws std::runtime_error if the configured oversampling ratio is unsupported.
	 */
	void upsample(const std::span<const ComplexType> in, const unsigned size, std::span<ComplexType> out)
	{
		const unsigned ratio = params::oversampleRatio();
		// TODO: this would be better as a multirate upsampler
		// This implementation is functional but suboptimal.
		// Users requiring higher accuracy should oversample outside FERS until this is addressed.
		validateOversamplingConfig(ratio);
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

	/**
	 * @brief Low-pass filters and decimates an oversampled complex waveform back to base rate.
	 *
	 * The same fixed-length FIR design is reused for anti-alias filtering, and unsupported
	 * ratios fail fast before filtering begins.
	 *
	 * @param in Input span of oversampled complex samples.
	 * @return Base-rate complex samples truncated to floor(in.size() / oversampleRatio()).
	 * @throws std::invalid_argument if @p in is empty.
	 * @throws std::runtime_error if the configured oversampling ratio is unsupported.
	 */
	std::vector<ComplexType> downsample(std::span<const ComplexType> in)
	{
		if (in.empty())
		{
			throw std::invalid_argument("Input span is empty in Downsample");
		}

		const unsigned ratio = params::oversampleRatio();
		// TODO: Replace with a more efficient multirate downsampling implementation.
		validateOversamplingConfig(ratio);
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
		/// 11th order elliptic lowpass at 0.1fs
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

		constexpr std::array num_coeffs = {2.7301694322809e-06,	  -1.8508123430239e-05,	 5.75739466753894e-05,
										   -0.000104348734423658, 0.000111949190289715,	 -4.9384188225528e-05,
										   -4.9384188225522e-05,  0.00011194919028971,	 -0.000104348734423656,
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
