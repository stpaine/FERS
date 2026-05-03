// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file radar_signal.cpp
 * @brief Classes for handling radar waveforms and signals.
 */

#include "radar_signal.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iterator>
#include <stdexcept>
#include <utility>

#include "core/parameters.h"
#include "dsp_filters.h"
#include "interpolation/interpolation_filter.h"
#include "interpolation/interpolation_point.h"

namespace fers_signal
{
	std::string_view fmcwChirpDirectionToken(const FmcwChirpDirection direction) noexcept
	{
		return direction == FmcwChirpDirection::Down ? "down" : "up";
	}

	FmcwChirpDirection parseFmcwChirpDirection(const std::string_view direction)
	{
		if (direction == "up")
		{
			return FmcwChirpDirection::Up;
		}
		if (direction == "down")
		{
			return FmcwChirpDirection::Down;
		}
		throw std::runtime_error("Unsupported FMCW chirp direction '" + std::string(direction) + "'.");
	}

	std::vector<ComplexType> CwSignal::render(const std::vector<interp::InterpPoint>& /*points*/, unsigned& size,
											  const RealType /*fracWinDelay*/) const
	{
		size = 0;
		return {};
	}

	FmcwChirpSignal::FmcwChirpSignal(const RealType chirp_bandwidth, const RealType chirp_duration,
									 const RealType chirp_period, const RealType start_frequency_offset,
									 std::optional<std::size_t> chirp_count, const FmcwChirpDirection direction) :
		_chirp_bandwidth(chirp_bandwidth), _chirp_duration(chirp_duration), _chirp_period(chirp_period),
		_start_frequency_offset(start_frequency_offset), _chirp_count(std::move(chirp_count)),
		_chirp_rate(chirp_bandwidth / chirp_duration), _direction(direction)
	{
	}

	std::optional<std::size_t>
	FmcwChirpSignal::activeChirpIndexAt(const RealType time_since_segment_start) const noexcept
	{
		if (time_since_segment_start < 0.0)
		{
			return std::nullopt;
		}

		const auto chirp_index = static_cast<std::size_t>(std::floor(time_since_segment_start / _chirp_period));
		if (_chirp_count.has_value() && chirp_index >= *_chirp_count)
		{
			return std::nullopt;
		}

		const RealType chirp_time = time_since_segment_start - static_cast<RealType>(chirp_index) * _chirp_period;
		// Exact arithmetic cannot make this negative, but floating-point boundary rounding can.
		if (chirp_time < 0.0 || chirp_time >= _chirp_duration)
		{
			return std::nullopt;
		}

		return chirp_index;
	}

	std::optional<RealType>
	FmcwChirpSignal::instantaneousBasebandPhase(const RealType time_since_segment_start) const noexcept
	{
		const auto chirp_index = activeChirpIndexAt(time_since_segment_start);
		if (!chirp_index.has_value())
		{
			return std::nullopt;
		}

		const RealType chirp_time = time_since_segment_start - static_cast<RealType>(*chirp_index) * _chirp_period;
		return basebandPhaseForChirpTime(chirp_time);
	}

	RealType FmcwChirpSignal::basebandPhaseForChirpTime(const RealType chirp_time) const noexcept
	{
		return 2.0 * PI * _start_frequency_offset * chirp_time + PI * getSignedChirpRate() * chirp_time * chirp_time;
	}

	std::vector<ComplexType> FmcwChirpSignal::render(const std::vector<interp::InterpPoint>& /*points*/, unsigned& size,
													 const RealType /*fracWinDelay*/) const
	{
		size = 0;
		return {};
	}

	FmcwTriangleSignal::FmcwTriangleSignal(const RealType chirp_bandwidth, const RealType chirp_duration,
										   const RealType start_frequency_offset,
										   std::optional<std::size_t> triangle_count) :
		_chirp_bandwidth(chirp_bandwidth), _chirp_duration(chirp_duration),
		_start_frequency_offset(start_frequency_offset), _triangle_count(std::move(triangle_count)),
		_chirp_rate(chirp_bandwidth / chirp_duration), _triangle_period(2.0 * chirp_duration),
		_delta_phi_up(2.0 * PI * start_frequency_offset * chirp_duration +
					  PI * _chirp_rate * chirp_duration * chirp_duration)
	{
	}

	RealType FmcwTriangleSignal::basebandPhaseForTriangleTime(const RealType triangle_time) const noexcept
	{
		if (triangle_time <= 0.0)
		{
			return 0.0;
		}

		const auto triangle_index = static_cast<std::size_t>(std::floor(triangle_time / _triangle_period));
		const RealType local_triangle_time = triangle_time - static_cast<RealType>(triangle_index) * _triangle_period;
		const bool down_leg = local_triangle_time >= _chirp_duration;
		const RealType u = down_leg ? local_triangle_time - _chirp_duration : local_triangle_time;
		const RealType phi_base =
			static_cast<RealType>(triangle_index) * 2.0 * _delta_phi_up + (down_leg ? _delta_phi_up : 0.0);
		if (!down_leg)
		{
			return phi_base + 2.0 * PI * _start_frequency_offset * u + PI * _chirp_rate * u * u;
		}
		return phi_base + 2.0 * PI * (_start_frequency_offset + _chirp_bandwidth) * u - PI * _chirp_rate * u * u;
	}

	std::optional<RealType>
	FmcwTriangleSignal::instantaneousBasebandPhase(const RealType time_since_segment_start) const noexcept
	{
		if (time_since_segment_start < 0.0)
		{
			return std::nullopt;
		}

		const auto triangle_index = static_cast<std::size_t>(std::floor(time_since_segment_start / _triangle_period));
		if (_triangle_count.has_value() && triangle_index >= *_triangle_count)
		{
			return std::nullopt;
		}

		const RealType local_triangle_time =
			time_since_segment_start - static_cast<RealType>(triangle_index) * _triangle_period;
		if (local_triangle_time < 0.0 || local_triangle_time >= _triangle_period)
		{
			return std::nullopt;
		}
		return basebandPhaseForTriangleTime(time_since_segment_start);
	}

	std::vector<ComplexType> FmcwTriangleSignal::render(const std::vector<interp::InterpPoint>& /*points*/,
														unsigned& size, const RealType /*fracWinDelay*/) const
	{
		size = 0;
		return {};
	}

	RadarSignal::RadarSignal(std::string name, const RealType power, const RealType carrierfreq, const RealType length,
							 std::unique_ptr<Signal> signal, const SimId id) :
		_name(std::move(name)), _id(id == 0 ? SimIdGenerator::instance().generateId(ObjectType::Waveform) : id),
		_power(power), _carrierfreq(carrierfreq), _length(length), _signal(std::move(signal))
	{
		if (!_signal)
		{
			throw std::runtime_error("Signal is empty");
		}
	}

	std::vector<ComplexType> RadarSignal::render(const std::vector<interp::InterpPoint>& points, unsigned& size,
												 const RealType fracWinDelay) const
	{
		auto data = _signal->render(points, size, fracWinDelay);
		const RealType scale = std::sqrt(_power);

		std::ranges::for_each(data, [scale](auto& value) { value *= scale; });

		return data;
	}

	std::vector<ComplexType> RadarSignal::renderSlice(const std::vector<interp::InterpPoint>& points,
													  const RealType outputStartTime, const RealType outputSampleRate,
													  const std::size_t sampleCount, const RealType fracWinDelay) const
	{
		auto data = _signal->renderSlice(points, outputStartTime, outputSampleRate, sampleCount, fracWinDelay);
		const RealType scale = std::sqrt(_power);

		std::ranges::for_each(data, [scale](auto& value) { value *= scale; });

		return data;
	}

	bool RadarSignal::isCw() const noexcept { return dynamic_cast<const CwSignal*>(_signal.get()) != nullptr; }

	bool RadarSignal::isFmcwChirp() const noexcept
	{
		return dynamic_cast<const FmcwChirpSignal*>(_signal.get()) != nullptr;
	}

	bool RadarSignal::isFmcwTriangle() const noexcept
	{
		return dynamic_cast<const FmcwTriangleSignal*>(_signal.get()) != nullptr;
	}

	bool RadarSignal::isFmcwFamily() const noexcept { return _signal->isFmcwFamily(); }

	const FmcwChirpSignal* RadarSignal::getFmcwChirpSignal() const noexcept
	{
		return dynamic_cast<const FmcwChirpSignal*>(_signal.get());
	}

	const FmcwTriangleSignal* RadarSignal::getFmcwTriangleSignal() const noexcept
	{
		return dynamic_cast<const FmcwTriangleSignal*>(_signal.get());
	}

	void Signal::clear() noexcept
	{
		_size = 0;
		_rate = 0;
	}

	void Signal::load(std::span<const ComplexType> inData, const unsigned samples, const RealType sampleRate)
	{
		clear();
		const unsigned ratio = params::oversampleRatio();
		_data.resize(samples * ratio);
		_size = samples * ratio;
		_rate = sampleRate * ratio;

		if (ratio == 1)
		{
			std::ranges::copy(inData, _data.begin());
		}
		else
		{
			upsample(inData, samples, _data);
		}
	}

	std::vector<ComplexType> Signal::render(const std::vector<interp::InterpPoint>& points, unsigned& size,
											const double fracWinDelay) const
	{
		size = _size;
		if (points.empty())
		{
			return std::vector<ComplexType>(_size);
		}
		return renderSlice(points, points.front().time, _rate, _size, fracWinDelay);
	}

	std::vector<ComplexType> Signal::renderSlice(const std::vector<interp::InterpPoint>& points,
												 const RealType outputStartTime, const RealType outputSampleRate,
												 const std::size_t sampleCount, const RealType fracWinDelay) const
	{
		auto out = std::vector<ComplexType>(sampleCount);
		if (_size == 0 || _rate <= 0.0 || outputSampleRate <= 0.0 || points.empty())
		{
			return out;
		}

		const RealType timestep = 1.0 / outputSampleRate;
		const int filt_length = static_cast<int>(params::renderFilterLength());
		const auto& interp = interp::InterpFilter::getInstance();

		auto iter = points.begin();
		auto next = points.size() > 1 ? std::next(iter) : iter;
		const RealType idelay = std::round(_rate * iter->delay);
		RealType sample_time = outputStartTime;

		for (std::size_t i = 0; i < sampleCount; ++i)
		{
			while (sample_time > next->time && next != iter)
			{
				iter = next;
				if (std::next(next) != points.end())
				{
					++next;
				}
				else
				{
					break;
				}
			}

			auto [amplitude, phase, fdelay, i_sample_unwrap] =
				calculateWeightsAndDelays(iter, next, sample_time, idelay, fracWinDelay);
			const auto& filt = interp.getFilter(fdelay);
			const auto source_index = static_cast<int>(std::llround((sample_time - points.front().time) * _rate));
			ComplexType accum = performConvolution(source_index, filt.data(), filt_length, amplitude, i_sample_unwrap);
			out[i] = std::exp(ComplexType(0.0, 1.0) * phase) * accum;

			sample_time += timestep;
		}

		return out;
	}

	std::tuple<RealType, RealType, RealType, int>
	Signal::calculateWeightsAndDelays(const std::vector<interp::InterpPoint>::const_iterator iter,
									  const std::vector<interp::InterpPoint>::const_iterator next,
									  const RealType sampleTime, const RealType idelay,
									  const RealType fracWinDelay) const noexcept
	{
		const RealType bw = iter < next ? (sampleTime - iter->time) / (next->time - iter->time) : 0.0;

		const RealType amplitude = std::lerp(std::sqrt(iter->power), std::sqrt(next->power), bw);
		const RealType phase = std::lerp(iter->phase, next->phase, bw);
		RealType fdelay = -(std::lerp(iter->delay, next->delay, bw) * _rate - idelay + fracWinDelay);

		const int i_sample_unwrap = static_cast<int>(std::floor(fdelay));
		fdelay -= i_sample_unwrap;

		return {amplitude, phase, fdelay, i_sample_unwrap};
	}

	ComplexType Signal::performConvolution(const int i, const RealType* filt, const int filtLength,
										   const RealType amplitude, const int iSampleUnwrap) const noexcept
	{
		const int start = std::max(-filtLength / 2, -i);
		const int end = std::min(filtLength / 2, static_cast<int>(_size) - i);

		ComplexType accum(0.0, 0.0);

		for (int j = start; j < end; ++j)
		{
			const int sample_idx = i + j + iSampleUnwrap;
			const int filt_idx = j + filtLength / 2;
			if (sample_idx >= 0 && sample_idx < static_cast<int>(_size) && filt_idx >= 0 && filt_idx < filtLength)
			{
				accum += amplitude * _data[static_cast<std::size_t>(sample_idx)] * filt[filt_idx];
			}
		}

		return accum;
	}
}
