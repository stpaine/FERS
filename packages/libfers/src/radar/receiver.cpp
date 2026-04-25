// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file receiver.cpp
 * @brief Implementation of the Receiver class.
 */

#include "receiver.h"

#include <algorithm>
#include <utility>

#include "core/parameters.h"
#include "serial/response.h"

namespace radar
{
	Receiver::Receiver(Platform* platform, std::string name, const unsigned seed, const OperationMode mode,
					   const SimId id) noexcept :
		Radar(platform, std::move(name), id == 0 ? SimIdGenerator::instance().generateId(ObjectType::Receiver) : id),
		_mode(mode), _rng(seed)
	{
	}

	void Receiver::addResponseToInbox(std::unique_ptr<serial::Response> response) noexcept
	{
		std::scoped_lock lock(_inbox_mutex);
		_inbox.push_back(std::move(response));
	}

	void Receiver::addInterferenceToLog(std::unique_ptr<serial::Response> response) noexcept
	{
		std::scoped_lock lock(_interference_log_mutex);
		_pulsed_interference_log.push_back(std::move(response));
	}

	std::vector<std::unique_ptr<serial::Response>> Receiver::drainInbox() noexcept
	{
		std::scoped_lock lock(_inbox_mutex);
		std::vector<std::unique_ptr<serial::Response>> drained_responses;
		drained_responses.swap(_inbox);
		return drained_responses;
	}

	void Receiver::enqueueFinalizerJob(core::RenderingJob&& job)
	{
		{
			std::scoped_lock lock(_finalizer_queue_mutex);
			_finalizer_queue.push(std::move(job));
		}
		_finalizer_queue_cv.notify_one();
	}

	bool Receiver::waitAndDequeueFinalizerJob(core::RenderingJob& job)
	{
		std::unique_lock lock(_finalizer_queue_mutex);
		_finalizer_queue_cv.wait(lock, [this] { return !_finalizer_queue.empty(); });

		job = std::move(_finalizer_queue.front());
		_finalizer_queue.pop();

		// Check for shutdown signal (negative duration)
		if (job.duration < 0.0)
		{
			return false; // Shutdown signal
		}
		return true;
	}

	RealType Receiver::getNoiseTemperature(const math::SVec3& angle) const noexcept
	{
		return _noise_temperature + Radar::getNoiseTemperature(angle);
	}

	void Receiver::setNoiseTemperature(const RealType temp)
	{
		if (temp < -EPSILON)
		{
			LOG(logging::Level::FATAL, "Noise temperature for receiver {} is negative", getName());
			throw std::runtime_error("Noise temperature must be positive");
		}
		_noise_temperature = temp;
	}

	void Receiver::setWindowProperties(const RealType length, const RealType prf, const RealType skip) noexcept
	{
		const auto rate = params::rate() * params::oversampleRatio();
		_window_length = length;
		_window_prf = 1 / (std::floor(rate / prf) / rate);
		_window_skip = std::floor(rate * skip) / rate;
	}

	unsigned Receiver::getWindowCount() const noexcept
	{
		const RealType time = params::endTime() - params::startTime();
		const RealType pulses = time * _window_prf;
		return static_cast<unsigned>(std::ceil(pulses));
	}

	RealType Receiver::getWindowStart(const unsigned window) const
	{
		const RealType stime = static_cast<RealType>(window) / _window_prf + _window_skip;
		if (!_timing)
		{
			LOG(logging::Level::FATAL, "Receiver must be associated with timing source");
			throw std::logic_error("Receiver must be associated with timing source");
		}
		return stime;
	}

	void Receiver::prepareStreamingData(const size_t numSamples)
	{
		std::scoped_lock lock(_cw_mutex);
		_streaming_iq_data.resize(numSamples);
	}

	void Receiver::setStreamingSample(const size_t index, const ComplexType sample)
	{
		if (index < _streaming_iq_data.size())
		{
			_streaming_iq_data[index] += sample;
		}
	}

	void Receiver::setSchedule(std::vector<SchedulePeriod> schedule) { _schedule = std::move(schedule); }

	std::optional<RealType> Receiver::getNextWindowTime(RealType time) const
	{
		// If no schedule is defined, assume always on.
		if (_schedule.empty())
		{
			return time;
		}
		for (const auto& period : _schedule)
		{
			// If time is within this period, it's valid.
			if (time >= period.start && time <= period.end)
			{
				return time;
			}
			// If time is before this period, skip to the start of this period.
			if (time < period.start)
			{
				return period.start;
			}
			// If time is after this period, continue to next period.
		}
		// Time is after the last scheduled period.
		return std::nullopt;
	}
}
