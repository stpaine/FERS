// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "finalizer_pipeline.h"

#include <algorithm>
#include <cmath>
#include <highfive/highfive.hpp>
#include <limits>
#include <optional>
#include <ranges>
#include <unordered_map>

#include "core/logging.h"
#include "core/output_metadata.h"
#include "core/parameters.h"
#include "processing/signal_processor.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/hdf5_handler.h"
#include "serial/response.h"
#include "signal/dsp_filters.h"
#include "simulation/channel_model.h"
#include "timing/timing.h"

namespace processing::pipeline
{
	namespace
	{
		/// Prepares caller-owned tracker storage for one independent receive window without shrinking capacity.
		void prepareWindowTrackerCache(core::ReceiverTrackerCache& tracker_cache, const std::size_t source_count,
									   const std::size_t target_count)
		{
			if (tracker_cache.direct.size() < source_count)
			{
				tracker_cache.direct.resize(source_count);
			}
			if (tracker_cache.reflected.size() < source_count)
			{
				tracker_cache.reflected.resize(source_count);
			}

			for (std::size_t source_index = 0; source_index < source_count; ++source_index)
			{
				tracker_cache.direct[source_index] = {};

				auto& reflected_trackers = tracker_cache.reflected[source_index];
				if (reflected_trackers.size() < target_count)
				{
					reflected_trackers.resize(target_count);
				}
				for (std::size_t target_index = 0; target_index < target_count; ++target_index)
				{
					reflected_trackers[target_index] = {};
				}
			}
		}
	}

	void advanceTimingModel(timing::Timing* timing_model, const radar::Receiver* receiver, const RealType rate)
	{

		if ((timing_model == nullptr) || !timing_model->isEnabled())
		{
			return;
		}

		// Convert the duration-derived sample advance to a signed temporary first.
		// The physical calculation can legitimately produce a negative or zero result
		// (for example due to floor() or an inter-pulse gap that is not positive),
		// but skipSamples() models advancing by a non-negative sample count only.
		// After validating that the computed value is > 0, cast to std::size_t for
		// the skipSamples() call chain:
		//   Timing::skipSamples -> ClockModelGenerator::skipSamples ->
		//   MultirateGenerator::skipSamples.
		// This preserves the sign until validation and avoids passing invalid negative
		// counts into the timing/noise generators.
		if (timing_model->getSyncOnPulse())
		{
			timing_model->reset();
			const auto skip_samples = static_cast<long long>(std::floor(rate * receiver->getWindowSkip()));
			if (skip_samples > 0)
			{
				timing_model->skipSamples(static_cast<std::size_t>(skip_samples));
			}
		}
		else
		{
			const RealType inter_pulse_skip_duration = 1.0 / receiver->getWindowPrf() - receiver->getWindowLength();
			const auto samples_to_skip = static_cast<long long>(std::floor(rate * inter_pulse_skip_duration));
			if (samples_to_skip > 0)
			{
				timing_model->skipSamples(static_cast<std::size_t>(samples_to_skip));
			}
		}
	}

	std::tuple<RealType, RealType> calculateJitteredStart(const RealType ideal_start, const RealType first_phase_noise,
														  const RealType carrier_freq, const RealType rate)
	{
		const RealType actual_start = ideal_start + first_phase_noise / (2.0 * PI * carrier_freq);
		const RealType rounded_start = std::round(actual_start * rate) / rate;
		const RealType fractional_delay = actual_start * rate - std::round(actual_start * rate);
		return {rounded_start, fractional_delay};
	}

	void applyStreamingInterference(std::span<ComplexType> window, const RealType actual_start, const RealType dt,
									const radar::Receiver* receiver,
									const std::vector<core::ActiveStreamingSource>& streaming_sources,
									const std::vector<std::unique_ptr<radar::Target>>* targets,
									core::ReceiverTrackerCache& tracker_cache,
									const simulation::CwPhaseNoiseLookup* phase_noise_lookup)
	{
		const simulation::CwPhaseNoiseLookup* lookup = phase_noise_lookup;
		std::optional<simulation::CwPhaseNoiseLookup> owned_lookup;
		if (lookup == nullptr)
		{
			std::unordered_map<SimId, std::shared_ptr<timing::Timing>> unique_timings;
			unique_timings.try_emplace(receiver->getTiming()->getId(), receiver->getTiming());
			for (const auto& streaming_source : streaming_sources)
			{
				unique_timings.try_emplace(streaming_source.transmitter->getTiming()->getId(),
										   streaming_source.transmitter->getTiming());
			}

			std::vector<std::shared_ptr<timing::Timing>> timings;
			timings.reserve(unique_timings.size());
			for (const auto& entry : unique_timings)
			{
				timings.push_back(entry.second);
			}

			const RealType end_time =
				actual_start + dt * static_cast<RealType>(window.empty() ? 0 : (window.size() - 1));
			owned_lookup = simulation::CwPhaseNoiseLookup::build(timings, actual_start, end_time);
			lookup = &*owned_lookup;
		}

		prepareWindowTrackerCache(tracker_cache, streaming_sources.size(), targets->size());

		RealType t_sample = actual_start;
		for (auto& window_sample : window)
		{
			ComplexType streaming_interference_sample{0.0, 0.0};
			for (std::size_t source_index = 0; source_index < streaming_sources.size(); ++source_index)
			{
				const auto& streaming_source = streaming_sources[source_index];
				if (!receiver->checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT))
				{
					streaming_interference_sample += simulation::calculateStreamingDirectPathContribution(
						streaming_source, receiver, t_sample, lookup, &tracker_cache.direct[source_index]);
				}
				for (std::size_t target_index = 0; target_index < targets->size(); ++target_index)
				{
					const auto& target_ptr = (*targets)[target_index];
					streaming_interference_sample += simulation::calculateStreamingReflectedPathContribution(
						streaming_source, receiver, target_ptr.get(), t_sample, lookup,
						&tracker_cache.reflected[source_index][target_index]);
				}
			}
			window_sample += streaming_interference_sample;
			t_sample += dt;
		}
	}

	void applyPulsedInterference(std::vector<ComplexType>& iq_buffer,
								 const std::vector<std::unique_ptr<serial::Response>>& interference_log)
	{
		for (const auto& response : interference_log)
		{
			unsigned psize;
			RealType prate;
			const auto rendered_pulse = response->renderBinary(prate, psize, 0.0);

			const RealType dt_sim = 1.0 / prate;
			const auto start_index = static_cast<size_t>((response->startTime() - params::startTime()) / dt_sim);

			for (size_t i = 0; i < psize; ++i)
			{
				if (start_index + i < iq_buffer.size())
				{
					iq_buffer[start_index + i] += rendered_pulse[i];
				}
			}
		}
	}

	void addPhaseNoiseToWindow(std::span<const RealType> noise, std::span<ComplexType> window)
	{
		for (auto [n, w] : std::views::zip(noise, window))
		{
			w *= std::polar(1.0, n);
		}
	}

	RealType applyDownsamplingAndQuantization(std::vector<ComplexType>& buffer)
	{
		if (params::oversampleRatio() > 1)
		{
			buffer = std::move(fers_signal::downsample(buffer));
		}
		return quantizeAndScaleWindow(buffer);
	}

	void exportStreamingToHdf5(const std::string& filename, const std::vector<ComplexType>& iq_buffer,
							   const RealType fullscale, const RealType ref_freq,
							   const core::OutputFileMetadata* metadata)
	{
		std::scoped_lock lock(serial::hdf5_global_mutex);
		try
		{
			HighFive::File file(filename, HighFive::File::Truncate);

			std::vector<RealType> i_data(iq_buffer.size());
			std::vector<RealType> q_data(iq_buffer.size());
			std::ranges::transform(iq_buffer, i_data.begin(), [](const auto& c) { return c.real(); });
			std::ranges::transform(iq_buffer, q_data.begin(), [](const auto& c) { return c.imag(); });

			HighFive::DataSet i_dataset = file.createDataSet<RealType>("I_data", HighFive::DataSpace::From(i_data));
			i_dataset.write(i_data);
			HighFive::DataSet q_dataset = file.createDataSet<RealType>("Q_data", HighFive::DataSpace::From(q_data));
			q_dataset.write(q_data);

			file.createAttribute("sampling_rate", params::rate());
			file.createAttribute("start_time", params::startTime());
			file.createAttribute("fullscale", fullscale);
			file.createAttribute("reference_carrier_frequency", ref_freq);
			if (metadata != nullptr)
			{
				serial::writeOutputFileMetadataAttributes(file, *metadata);
			}

			LOG(logging::Level::INFO, "Successfully exported streaming data to '{}'", filename);
		}
		catch (const HighFive::Exception& err)
		{
			LOG(logging::Level::FATAL, "Error writing streaming data to HDF5 file '{}': {}", filename, err.what());
		}
	}

}
