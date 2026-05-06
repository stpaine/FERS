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
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <utility>

#include "core/parameters.h"
#include "serial/response.h"

namespace radar
{
	namespace
	{
		[[nodiscard]] std::size_t ceilSampleIndexAtOrAfter(const RealType time, const RealType sample_rate)
		{
			if (sample_rate <= 0.0 || time <= params::startTime())
			{
				return 0;
			}
			const RealType raw_index = (time - params::startTime()) * sample_rate;
			const RealType nearest = std::round(raw_index);
			const RealType tolerance = 1.0e-12 * std::max<RealType>(1.0, std::abs(nearest));
			if (std::abs(raw_index - nearest) <= tolerance)
			{
				return static_cast<std::size_t>(nearest);
			}
			return static_cast<std::size_t>(std::ceil(raw_index));
		}
	}

	std::string_view dechirpModeToken(const Receiver::DechirpMode mode) noexcept
	{
		switch (mode)
		{
		case Receiver::DechirpMode::Physical:
			return "physical";
		case Receiver::DechirpMode::Ideal:
			return "ideal";
		case Receiver::DechirpMode::None:
			return "none";
		}
		return "none";
	}

	Receiver::DechirpMode parseDechirpModeToken(const std::string_view token)
	{
		if (token == "none")
		{
			return Receiver::DechirpMode::None;
		}
		if (token == "physical")
		{
			return Receiver::DechirpMode::Physical;
		}
		if (token == "ideal")
		{
			return Receiver::DechirpMode::Ideal;
		}
		throw std::runtime_error("Unsupported FMCW dechirp_mode '" + std::string(token) + "'.");
	}

	std::string_view dechirpReferenceSourceToken(const Receiver::DechirpReferenceSource source) noexcept
	{
		switch (source)
		{
		case Receiver::DechirpReferenceSource::Attached:
			return "attached";
		case Receiver::DechirpReferenceSource::Transmitter:
			return "transmitter";
		case Receiver::DechirpReferenceSource::Custom:
			return "custom";
		case Receiver::DechirpReferenceSource::None:
			return "none";
		}
		return "none";
	}

	Receiver::DechirpReferenceSource parseDechirpReferenceSourceToken(const std::string_view token)
	{
		if (token == "attached")
		{
			return Receiver::DechirpReferenceSource::Attached;
		}
		if (token == "transmitter")
		{
			return Receiver::DechirpReferenceSource::Transmitter;
		}
		if (token == "custom")
		{
			return Receiver::DechirpReferenceSource::Custom;
		}
		throw std::runtime_error("Unsupported dechirp_reference source '" + std::string(token) + "'.");
	}

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

	void Receiver::setMode(const OperationMode mode) noexcept
	{
		_mode = mode;
		if (_mode != OperationMode::FMCW_MODE)
		{
			setDechirpMode(DechirpMode::None);
			_fmcw_if_chain = {};
		}
	}

	void Receiver::setDechirpMode(const DechirpMode mode) noexcept
	{
		_dechirp_mode = mode;
		if (_dechirp_mode == DechirpMode::None)
		{
			_dechirp_reference = {};
			_dechirp_sources.clear();
			_fmcw_if_chain = {};
			_fmcw_if_plan.reset();
			_fmcw_if_sink.reset();
			_fmcw_if_samples_to_discard = 0;
			_fmcw_if_input_cursor = 0;
			_fmcw_if_output_cursor = 0;
			_fmcw_if_segment_active = false;
		}
	}

	void Receiver::setDechirpReference(DechirpReference reference)
	{
		_dechirp_reference = std::move(reference);
		_dechirp_sources.clear();
	}

	void Receiver::setFmcwIfChainRequest(FmcwIfChainRequest request) noexcept
	{
		_fmcw_if_chain = std::move(request);
		_fmcw_if_plan.reset();
		_fmcw_if_sink.reset();
		_fmcw_if_samples_to_discard = 0;
		_fmcw_if_input_cursor = 0;
		_fmcw_if_output_cursor = 0;
		_fmcw_if_segment_active = false;
	}

	void Receiver::initializeFmcwIfResampling(fers_signal::FmcwIfResamplerPlan plan)
	{
		_streaming_iq_data.clear();
		_fmcw_if_plan = std::move(plan);
		const auto expected_samples =
			static_cast<std::size_t>(std::ceil(std::max<RealType>(0.0, params::endTime() - params::startTime()) *
											   _fmcw_if_plan->actual_output_sample_rate_hz));
		_streaming_iq_data.reserve(expected_samples);
		_fmcw_if_samples_to_discard = _fmcw_if_plan->warmup_discard_samples;
		_fmcw_if_sink = std::make_unique<fers_signal::FmcwIfResamplingSink>(*_fmcw_if_plan);
		_fmcw_if_input_cursor = 0;
		_fmcw_if_output_cursor = 0;
		_fmcw_if_segment_active = false;
	}

	void Receiver::beginFmcwIfResamplingSegment(const RealType segment_start_time)
	{
		if (_fmcw_if_sink == nullptr || !_fmcw_if_plan.has_value())
		{
			return;
		}
		consumeFmcwIfZerosUntil(ceilSampleIndexAtOrAfter(segment_start_time, _fmcw_if_plan->input_sample_rate_hz));
		_fmcw_if_segment_active = true;
	}

	void Receiver::consumeFmcwIfBlock(const std::span<const ComplexType> block, const RealType block_start_time)
	{
		if (_fmcw_if_sink == nullptr || !_fmcw_if_plan.has_value())
		{
			throw std::logic_error("FMCW IF resampling sink has not been initialized.");
		}
		if (!_fmcw_if_segment_active)
		{
			beginFmcwIfResamplingSegment(block_start_time);
		}
		const auto block_start_index = ceilSampleIndexAtOrAfter(block_start_time, _fmcw_if_plan->input_sample_rate_hz);
		if (block_start_index < _fmcw_if_input_cursor)
		{
			throw std::logic_error("FMCW IF resampling input blocks must be supplied in chronological order.");
		}
		consumeFmcwIfZerosUntil(block_start_index);
		_fmcw_if_sink->consume(block);
		_fmcw_if_input_cursor += block.size();
		auto emitted = _fmcw_if_sink->takeOutput();
		appendFmcwIfOutput(std::move(emitted));
	}

	void Receiver::appendFmcwIfOutput(std::vector<ComplexType> emitted)
	{
		if (_fmcw_if_samples_to_discard > 0)
		{
			const auto discard = std::min<std::uint64_t>(_fmcw_if_samples_to_discard, emitted.size());
			emitted.erase(emitted.begin(), emitted.begin() + static_cast<std::ptrdiff_t>(discard));
			_fmcw_if_samples_to_discard -= discard;
		}
		if (emitted.empty())
		{
			return;
		}
		if (_streaming_iq_data.size() < _fmcw_if_output_cursor)
		{
			_streaming_iq_data.resize(_fmcw_if_output_cursor);
		}
		if (_streaming_iq_data.size() < _fmcw_if_output_cursor + emitted.size())
		{
			_streaming_iq_data.resize(_fmcw_if_output_cursor + emitted.size());
		}
		for (std::size_t i = 0; i < emitted.size(); ++i)
		{
			_streaming_iq_data[_fmcw_if_output_cursor + i] += emitted[i];
		}
		_fmcw_if_output_cursor += emitted.size();
	}

	void Receiver::advanceFmcwIfOutputZeros(std::size_t sample_count)
	{
		if (_fmcw_if_samples_to_discard > 0)
		{
			const auto discard = std::min<std::uint64_t>(_fmcw_if_samples_to_discard, sample_count);
			sample_count -= static_cast<std::size_t>(discard);
			_fmcw_if_samples_to_discard -= discard;
		}
		if (sample_count == 0)
		{
			return;
		}
		if (_streaming_iq_data.size() < _fmcw_if_output_cursor + sample_count)
		{
			_streaming_iq_data.resize(_fmcw_if_output_cursor + sample_count);
		}
		_fmcw_if_output_cursor += sample_count;
	}

	void Receiver::endFmcwIfResamplingSegment()
	{
		if (_fmcw_if_sink == nullptr)
		{
			return;
		}
		_fmcw_if_segment_active = false;
	}

	void Receiver::flushFmcwIfResampling()
	{
		if (_fmcw_if_sink == nullptr)
		{
			return;
		}
		endFmcwIfResamplingSegment();
		if (_fmcw_if_plan.has_value())
		{
			const RealType flush_until_time = params::endTime() + _fmcw_if_plan->group_delay_seconds +
				1.0 / _fmcw_if_plan->actual_output_sample_rate_hz;
			consumeFmcwIfZerosUntil(ceilSampleIndexAtOrAfter(flush_until_time, _fmcw_if_plan->input_sample_rate_hz));
			auto emitted = _fmcw_if_sink->finish();
			appendFmcwIfOutput(std::move(emitted));

			const auto expected_samples =
				static_cast<std::size_t>(std::ceil(std::max<RealType>(0.0, params::endTime() - params::startTime()) *
												   _fmcw_if_plan->actual_output_sample_rate_hz));
			if (_streaming_iq_data.size() > expected_samples)
			{
				_streaming_iq_data.resize(expected_samples);
			}
			else if (_streaming_iq_data.size() < expected_samples)
			{
				_streaming_iq_data.resize(expected_samples);
			}
		}
		_fmcw_if_sink.reset();
	}

	void Receiver::consumeFmcwIfZerosUntil(const std::size_t target_input_cursor)
	{
		if (_fmcw_if_sink == nullptr)
		{
			throw std::logic_error("FMCW IF resampling sink has not been initialized.");
		}
		if (target_input_cursor <= _fmcw_if_input_cursor)
		{
			return;
		}

		const auto count = target_input_cursor - _fmcw_if_input_cursor;
		auto result = _fmcw_if_sink->consumeZeroInput(count);
		_fmcw_if_input_cursor = target_input_cursor;
		appendFmcwIfOutput(std::move(result.emitted));
		advanceFmcwIfOutputZeros(result.skipped_output_samples);
	}

	void Receiver::setResolvedDechirpSources(std::vector<core::ActiveStreamingSource> sources)
	{
		_dechirp_sources = std::move(sources);
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
