// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2007-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file dsp_filters.h
 * @brief Header file for Digital Signal Processing (DSP) filters and upsampling/downsampling functionality.
 */

#pragma once

#include <memory>
#include <span>
#include <vector>

#include "core/config.h"

namespace fers_signal
{
	/**
	 * @brief Upsamples a signal by a given ratio.
	 *
	 * @param in Input span of complex samples.
	 * @param size Size of the input signal.
	 * @param out Output span for upsampled complex samples.
	 * @throws std::invalid_argument if the input or output spans are empty or the ratio is zero.
	 */
	void upsample(std::span<const ComplexType> in, unsigned size, std::span<ComplexType> out);

	/**
	 * @brief Downsamples a signal by a given ratio.
	 *
	 * @param in Input span of complex samples.
	 * @throws std::invalid_argument if the input or output spans are empty or the ratio is zero.
	 */
	std::vector<ComplexType> downsample(std::span<const ComplexType> in);

	/**
	 * @class DspFilter
	 * @brief Abstract base class for digital filters.
	 */
	class DspFilter
	{
	public:
		DspFilter() = default;

		virtual ~DspFilter() = default;

		/**
		 * @brief Filters a single sample.
		 *
		 * @param sample A single real-valued sample to be filtered.
		 * @return The filtered sample.
		 */
		virtual RealType filter(RealType sample) = 0;

		/**
		 * @brief Filters a block of samples.
		 *
		 * @param samples Span of real-valued samples to be filtered.
		 */
		virtual void filter(std::span<RealType> samples) = 0;

		DspFilter(const DspFilter&) = delete;

		DspFilter& operator=(const DspFilter&) = delete;

		DspFilter(DspFilter&&) noexcept = default;

		DspFilter& operator=(DspFilter&&) noexcept = default;
	};

	/**
	 * @class IirFilter
	 * @brief Implements an Infinite Impulse Response (IIR) filter.
	 */
	class IirFilter final : public DspFilter
	{
	public:
		/**
		 * @brief Constructs an IIR filter with given numerator and denominator coefficients and order.
		 *
		 * @param denCoeffs Pointer to the denominator coefficients array.
		 * @param numCoeffs Pointer to the numerator coefficients array.
		 * @param order The order of the filter.
		 */
		IirFilter(const RealType* denCoeffs, const RealType* numCoeffs, unsigned order) noexcept;

		~IirFilter() override = default;

		/**
		 * @brief Filters a single sample.
		 *
		 * @param sample The sample to be filtered.
		 * @return The filtered sample.
		 */
		RealType filter(RealType sample) noexcept override;

		/**
		 * @brief Filters a block of samples.
		 *
		 * @param samples Span of samples to be filtered.
		 */
		void filter(std::span<RealType> samples) noexcept override;

	private:
		std::vector<RealType> _a; ///< Denominator coefficients
		std::vector<RealType> _b; ///< Numerator coefficients
		std::vector<RealType> _w; ///< Internal state
		unsigned _order{}; ///< Filter order
	};

	/**
	 * @class FirFilter
	 * @brief Implements a Finite Impulse Response (FIR) filter.
	 */
	class FirFilter final : public DspFilter
	{
	public:
		/**
		 * @brief Constructs an FIR filter with the given coefficients.
		 *
		 * @param coeffs Span of filter coefficients.
		 */
		explicit FirFilter(std::span<const RealType> coeffs) :
			_filter(coeffs.begin(), coeffs.end()), _w(coeffs.size()), _order(coeffs.size())
		{
		}

		~FirFilter() override = default;

		RealType filter(RealType) override { return 0; }

		void filter(std::span<RealType> /*samples*/) noexcept override {}

		/**
		 * @brief Filters a block of complex samples.
		 *
		 * @param samples Span of complex samples to be filtered.
		 */
		void filter(std::vector<ComplexType>& samples) const;

	private:
		std::vector<RealType> _filter; ///< Filter coefficients
		std::vector<RealType> _w; ///< Internal state
		unsigned _order{}; ///< Filter order
	};

	/**
	 * @class DecadeUpsampler
	 * @brief Implements a specialized upsampler with a fixed upsampling factor of 10.
	 */
	class DecadeUpsampler
	{
	public:
		DecadeUpsampler();

		~DecadeUpsampler() = default;

		/**
		 * @brief Upsamples a single sample.
		 *
		 * @param sample The sample to be upsampled.
		 * @param out Span of output samples.
		 */
		void upsample(RealType sample, std::span<RealType> out) const;

		DecadeUpsampler(const DecadeUpsampler&) = delete;

		DecadeUpsampler& operator=(const DecadeUpsampler&) = delete;

		DecadeUpsampler(DecadeUpsampler&&) noexcept = default;

		DecadeUpsampler& operator=(DecadeUpsampler&&) noexcept = default;

	private:
		std::unique_ptr<IirFilter> _filter; ///< IIR filter for upsampling
	};
}
