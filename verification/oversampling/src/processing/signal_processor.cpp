// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// Upstream reference:
//   packages/libfers/src/processing/signal_processor.cpp
//
// Local simplification:
//   thermal noise uses std::normal_distribution directly instead of noise::WgnGenerator.

#include "signal_processor.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <tuple>
#include <vector>

#include "core/parameters.h"
#include "serial/response.h"

namespace
{
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

	void processResponse(const serial::Response* resp, std::vector<ComplexType>& localWindow, const RealType rate,
						 const RealType start, const RealType fracDelay, const unsigned localWindowSize)
	{
		unsigned psize;
		RealType prate;
		const auto array = resp->renderBinary(prate, psize, fracDelay);
		int start_sample = static_cast<int>(std::round(rate * (resp->startTime() - start)));
		const unsigned roffset = start_sample < 0 ? static_cast<unsigned>(-start_sample) : 0u;
		start_sample = std::max(start_sample, 0);

		for (unsigned i = roffset; i < psize && i + static_cast<unsigned>(start_sample) < localWindowSize; ++i)
		{
			localWindow[i + static_cast<unsigned>(start_sample)] += array[i];
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
		const auto local_window_size = static_cast<unsigned>(std::ceil(length * rate));

		std::vector local_window(local_window_size, ComplexType{});

		while (!work_list.empty())
		{
			const auto* resp = work_list.front();
			work_list.pop();
			processResponse(resp, local_window, rate, start, fracDelay, local_window_size);
		}

		for (unsigned i = 0; i < local_window_size; ++i)
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
		const RealType per_channel_power = total_power / 2.0;
		const RealType stddev = std::sqrt(per_channel_power);
		std::normal_distribution<RealType> dist(0.0, stddev);

		for (auto& sample : window)
		{
			sample += ComplexType(dist(rngEngine), dist(rngEngine));
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
}

