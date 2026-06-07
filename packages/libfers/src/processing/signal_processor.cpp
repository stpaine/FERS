// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file signal_processor.cpp
 * @brief Implementation for receiver-side signal processing and rendering.
 *
 * This file contains the implementation of functions that perform digital signal
 * processing on received radar signals, as declared in signal_processor.h.
 */

#include "signal_processor.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <future>
#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "core/parameters.h"
#include "core/thread_pool.h"
#include "noise/noise_generators.h"
#include "serial/response.h"

namespace
{
	[[nodiscard]] std::int16_t scaleComponentToInt16(const RealType value, const RealType fullscale,
													 bool& clipped) noexcept
	{
		if (value > fullscale)
		{
			clipped = true;
			return std::numeric_limits<std::int16_t>::max();
		}
		if (value < -fullscale)
		{
			clipped = true;
			return std::numeric_limits<std::int16_t>::min();
		}
		if (value == fullscale)
		{
			return std::numeric_limits<std::int16_t>::max();
		}
		if (value == -fullscale)
		{
			return std::numeric_limits<std::int16_t>::min();
		}

		const RealType scaled =
			std::round((value / fullscale) * static_cast<RealType>(std::numeric_limits<std::int16_t>::max()));
		return static_cast<std::int16_t>(
			std::clamp<RealType>(scaled, static_cast<RealType>(std::numeric_limits<std::int16_t>::min()),
								 static_cast<RealType>(std::numeric_limits<std::int16_t>::max())));
	}

	/**
	 * @brief Simulate an ADC quantization process on a window of complex samples.
	 * @param data A span of ComplexType objects representing the window to quantize.
	 * @param bits The number of bits used for quantization.
	 * @param fullscale The full-scale value used for quantization.
	 */
	void adcSimulate(std::span<ComplexType> data, const unsigned bits, const RealType fullscale) noexcept
	{
		const RealType levels = std::pow(2, bits - 1);

		for (auto& sample : data)
		{
			auto [i, q] = std::tuple{std::clamp(std::floor(levels * sample.real() / fullscale) / levels, -1.0, 1.0),
									 std::clamp(std::floor(levels * sample.imag() / fullscale) / levels, -1.0, 1.0)};
			sample = ComplexType(i, q);
		}
	}

	/**
	 * @brief Process a single response and add its rendered data to a local window buffer.
	 * @param resp A pointer to a Response object to process.
	 * @param localWindow A reference to a vector of ComplexType representing the local window.
	 * @param rate The sample rate in Hz.
	 * @param start The start time of the window in seconds.
	 * @param fracDelay The fractional delay of the window in seconds.
	 * @param localWindowSize The size of the local window in samples.
	 */
	void processResponse(const serial::Response* resp, std::vector<ComplexType>& localWindow, const RealType rate,
						 const RealType start, const RealType fracDelay, const std::size_t localWindowSize)
	{
		unsigned psize = 0;
		RealType prate = std::numeric_limits<RealType>::quiet_NaN();
		const auto array = resp->renderBinary(prate, psize, fracDelay);
		int const start_sample = static_cast<int>(std::round(rate * (resp->startTime() - start)));
		const auto roffset = start_sample < 0 ? static_cast<std::size_t>(-start_sample) : std::size_t{0};
		const auto start_offset = static_cast<std::size_t>(std::max(start_sample, 0));

		for (std::size_t i = roffset; i < psize && i + start_offset < localWindowSize; ++i)
		{
			localWindow[i + start_offset] += array[i];
		}
	}
}

namespace processing
{
	void renderWindow(std::vector<ComplexType>& window, const RealType length, const RealType start,
					  const RealType fracDelay, const std::span<const std::unique_ptr<serial::Response>> responses)
	{
		const RealType end = start + length;
		std::queue<serial::Response*> work_list;

		for (const auto& response : responses)
		{
			if (response->startTime() <= end && response->endTime() >= start)
			{
				work_list.push(response.get());
			}
		}

		const RealType rate = params::rate() * params::oversampleRatio();
		const auto local_window_size = static_cast<std::size_t>(std::ceil(length * rate));

		std::vector local_window(local_window_size, ComplexType{});

		while (!work_list.empty())
		{
			const auto* resp = work_list.front();
			work_list.pop();
			processResponse(resp, local_window, rate, start, fracDelay, local_window_size);
		}

		for (std::size_t i = 0; i < local_window_size; ++i)
		{
			window[i] += local_window[i];
		}
	}

	void applyThermalNoise(std::span<ComplexType> window, const RealType noiseTemperature, std::mt19937& rngEngine)
	{
		if (noiseTemperature == 0)
		{
			return;
		}

		const RealType b = params::rate() / (2.0 * params::oversampleRatio());
		const RealType total_power = params::boltzmannK() * noiseTemperature * b;

		// Split total power equally between the I and Q channels
		const RealType per_channel_power = total_power / 2.0;
		const RealType stddev = std::sqrt(per_channel_power);

		noise::WgnGenerator generator(rngEngine, stddev);
		for (auto& sample : window)
		{
			sample += ComplexType(generator.getSample(), generator.getSample());
		}
	}

	void applyThermalNoiseAtSampleRate(std::span<ComplexType> window, const RealType noiseTemperature,
									   std::mt19937& rngEngine, const RealType sampleRateHz)
	{
		if (noiseTemperature == 0 || sampleRateHz <= 0.0)
		{
			return;
		}

		const RealType complex_power = params::boltzmannK() * noiseTemperature * sampleRateHz;
		const RealType per_channel_power = complex_power / 2.0;
		const RealType stddev = std::sqrt(per_channel_power);

		noise::WgnGenerator generator(rngEngine, stddev);
		for (auto& sample : window)
		{
			sample += ComplexType(generator.getSample(), generator.getSample());
		}
	}

	RealType quantizeAndScaleWindow(std::span<ComplexType> window)
	{
		RealType max_value = 0;
		for (const auto& sample : window)
		{
			const RealType real_abs = std::fabs(sample.real());
			const RealType imag_abs = std::fabs(sample.imag());

			max_value = std::max({max_value, real_abs, imag_abs});
		}

		if (const auto adc_bits = params::adcBits(); adc_bits > 0)
		{
			adcSimulate(window, adc_bits, max_value);
		}
		else if (max_value != 0)
		{
			for (auto& sample : window)
			{
				sample /= max_value;
			}
		}
		return max_value;
	}

	FixedFullscaleScalingResult scaleToInt16FixedFullscale(const std::span<const ComplexType> samples,
														   const RealType fullscale)
	{
		if (fullscale <= 0.0 || !std::isfinite(fullscale))
		{
			throw std::invalid_argument("Fixed full-scale must be positive and finite.");
		}

		FixedFullscaleScalingResult result;
		result.samples.reserve(samples.size());
		for (const auto& sample : samples)
		{
			bool clipped = false;
			const auto i = scaleComponentToInt16(sample.real(), fullscale, clipped);
			const auto q = scaleComponentToInt16(sample.imag(), fullscale, clipped);
			if (clipped)
			{
				++result.clipped_sample_count;
			}
			result.samples.push_back(FixedFullscaleIqSample{.i = i, .q = q});
		}
		return result;
	}
}
