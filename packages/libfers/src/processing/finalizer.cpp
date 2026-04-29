// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "finalizer.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <format>
#include <highfive/highfive.hpp>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "core/logging.h"
#include "core/output_metadata.h"
#include "core/parameters.h"
#include "core/rendering_job.h"
#include "core/sim_threading.h"
#include "processing/finalizer_pipeline.h"
#include "processing/signal_processor.h"
#include "radar/receiver.h"
#include "radar/transmitter.h"
#include "serial/hdf5_handler.h"
#include "signal/radar_signal.h"
#include "timing/timing.h"

namespace processing
{
	namespace
	{
		/// Finalizes aggregate pulsed metadata after all chunks have been collected.
		void finalizePulsedMetadata(core::OutputFileMetadata& metadata)
		{
			metadata.pulse_count = static_cast<std::uint64_t>(metadata.chunks.size());
			metadata.total_samples = 0;
			metadata.min_pulse_length_samples = metadata.chunks.empty() ? 0 : std::numeric_limits<std::uint64_t>::max();
			metadata.max_pulse_length_samples = 0;
			metadata.uniform_pulse_length = true;

			for (const auto& chunk : metadata.chunks)
			{
				metadata.total_samples += chunk.sample_count;
				metadata.min_pulse_length_samples = std::min(metadata.min_pulse_length_samples, chunk.sample_count);
				metadata.max_pulse_length_samples = std::max(metadata.max_pulse_length_samples, chunk.sample_count);
			}

			if (!metadata.chunks.empty())
			{
				const auto expected = metadata.chunks.front().sample_count;
				metadata.uniform_pulse_length = std::ranges::all_of(metadata.chunks, [expected](const auto& chunk)
																	{ return chunk.sample_count == expected; });
			}

			metadata.sample_start = 0;
			metadata.sample_end_exclusive = metadata.total_samples;
		}

		/// Converts a cached streaming source to reusable FMCW waveform metadata.
		core::FmcwMetadata buildFmcwMetadata(const core::ActiveStreamingSource& source)
		{
			if (source.kind == core::StreamingWaveformKind::FmcwTriangle)
			{
				return core::FmcwMetadata{
					.waveform_shape = "triangle",
					.chirp_bandwidth = source.triangle != nullptr ? source.triangle->getChirpBandwidth() : 0.0,
					.chirp_duration = source.chirp_duration,
					.chirp_rate = source.chirp_rate,
					.start_frequency_offset = source.start_freq_off,
					.triangle_period = source.triangle_period,
					.triangle_count = source.triangle_count.has_value()
						? std::optional<std::uint64_t>(static_cast<std::uint64_t>(*source.triangle_count))
						: std::nullopt};
			}

			return core::FmcwMetadata{
				.waveform_shape = "linear",
				.chirp_bandwidth = source.fmcw != nullptr ? source.fmcw->getChirpBandwidth() : 0.0,
				.chirp_duration = source.chirp_duration,
				.chirp_period = source.chirp_period,
				.chirp_rate = source.chirp_rate,
				.chirp_rate_signed = source.signed_chirp_rate,
				.chirp_direction = source.fmcw != nullptr
					? std::string(fers_signal::fmcwChirpDirectionToken(source.fmcw->getDirection()))
					: std::string("up"),
				.start_frequency_offset = source.start_freq_off,
				.chirp_count = source.chirp_count.has_value()
					? std::optional<std::uint64_t>(static_cast<std::uint64_t>(*source.chirp_count))
					: std::nullopt};
		}

		/// Builds one FMCW source schedule segment from an active source cache.
		core::FmcwSourceSegmentMetadata buildFmcwSourceSegment(const core::ActiveStreamingSource& source)
		{
			const RealType active_start = std::max(params::startTime(), source.segment_start);
			const RealType active_end = std::min(params::endTime(), source.segment_end);
			core::FmcwSourceSegmentMetadata segment{.start_time = source.segment_start, .end_time = source.segment_end};
			if (source.kind == core::StreamingWaveformKind::FmcwTriangle)
			{
				segment.first_triangle_start_time = core::firstFmcwTriangleStart(source, active_start, active_end);
				segment.emitted_triangle_count = core::countFmcwTriangleStarts(source, active_start, active_end);
			}
			else
			{
				segment.first_chirp_start_time = core::firstFmcwChirpStart(source, active_start, active_end);
				segment.emitted_chirp_count = core::countFmcwChirpStarts(source, active_start, active_end);
			}
			return segment;
		}

		/// Finds the first source metadata entry for a transmitter/waveform pair.
		std::vector<core::FmcwSourceMetadata>::iterator findFmcwSource(std::vector<core::FmcwSourceMetadata>& sources,
																	   const SimId transmitter_id,
																	   const SimId waveform_id)
		{
			return std::ranges::find_if(
				sources, [&](const core::FmcwSourceMetadata& source)
				{ return source.transmitter_id == transmitter_id && source.waveform_id == waveform_id; });
		}

		/// Builds explicit per-source FMCW metadata from active streaming transmitters.
		std::vector<core::FmcwSourceMetadata>
		buildFmcwSources(const std::vector<core::ActiveStreamingSource>& streaming_sources)
		{
			std::vector<core::FmcwSourceMetadata> fmcw_sources;
			for (const auto& streaming_source : streaming_sources)
			{
				if (!streaming_source.is_fmcw || streaming_source.transmitter == nullptr)
				{
					continue;
				}

				const auto* signal = streaming_source.transmitter->getSignal();
				if (signal == nullptr)
				{
					continue;
				}

				const auto transmitter_id = streaming_source.transmitter->getId();
				const auto waveform_id = signal->getId();
				auto existing = findFmcwSource(fmcw_sources, transmitter_id, waveform_id);
				if (existing == fmcw_sources.end())
				{
					core::FmcwSourceMetadata source{.transmitter_id = transmitter_id,
													.transmitter_name = streaming_source.transmitter->getName(),
													.waveform_id = waveform_id,
													.waveform_name = signal->getName(),
													.carrier_frequency = signal->getCarrier(),
													.waveform = buildFmcwMetadata(streaming_source)};
					source.segments.push_back(buildFmcwSourceSegment(streaming_source));
					fmcw_sources.push_back(std::move(source));
					continue;
				}

				existing->segments.push_back(buildFmcwSourceSegment(streaming_source));
			}
			return fmcw_sources;
		}

		/// Adds scalar compatibility chirp metadata to receiver streaming segments for one FMCW source.
		void annotateStreamingSegmentsForSingleFmcwSource(core::OutputFileMetadata& metadata,
														  const core::ActiveStreamingSource& source)
		{
			for (auto& segment : metadata.streaming_segments)
			{
				const RealType active_start = std::max(segment.start_time, source.segment_start);
				const RealType active_end = std::min(segment.end_time, source.segment_end);
				if (source.kind == core::StreamingWaveformKind::FmcwTriangle)
				{
					const auto first_triangle = core::firstFmcwTriangleStart(source, active_start, active_end);
					const auto emitted = core::countFmcwTriangleStarts(source, active_start, active_end);
					if (first_triangle.has_value() || emitted > 0)
					{
						segment.first_triangle_start_time = first_triangle;
						segment.emitted_triangle_count = emitted;
					}
				}
				else
				{
					const auto first_chirp = core::firstFmcwChirpStart(source, active_start, active_end);
					const auto emitted = core::countFmcwChirpStarts(source, active_start, active_end);
					if (first_chirp.has_value() || emitted > 0)
					{
						segment.first_chirp_start_time = first_chirp;
						segment.emitted_chirp_count = emitted;
					}
				}
			}
		}

		/// Half-open time interval in simulation seconds.
		using TimeSpan = std::pair<RealType, RealType>;

		/// Converts a time interval to a half-open sample span at a specific output rate.
		std::optional<pipeline::SampleSpan> timeSpanToSampleSpan(const TimeSpan& span, const std::size_t total_samples,
																 const RealType sample_rate)
		{
			const auto start_sample = static_cast<std::size_t>(std::min<RealType>(
				static_cast<RealType>(total_samples),
				std::max<RealType>(0.0, std::ceil((span.first - params::startTime()) * sample_rate))));
			const auto end_sample = static_cast<std::size_t>(std::min<RealType>(
				static_cast<RealType>(total_samples),
				std::max<RealType>(0.0, std::ceil((span.second - params::startTime()) * sample_rate))));
			if (start_sample >= end_sample)
			{
				return std::nullopt;
			}
			return pipeline::SampleSpan{.start = start_sample, .end_exclusive = end_sample};
		}

		/// Merges overlapping or adjacent time spans.
		void normalizeTimeSpans(std::vector<TimeSpan>& spans)
		{
			std::ranges::sort(spans, [](const TimeSpan& lhs, const TimeSpan& rhs) { return lhs.first < rhs.first; });
			std::vector<TimeSpan> merged;
			for (const auto& span : spans)
			{
				if (span.second <= span.first)
				{
					continue;
				}
				if (merged.empty() || span.first > merged.back().second)
				{
					merged.push_back(span);
					continue;
				}
				merged.back().second = std::max(merged.back().second, span.second);
			}
			spans = std::move(merged);
		}

		/// Returns receiver active intervals clipped to simulation time.
		std::vector<TimeSpan> receiverActiveTimeSpans(const radar::Receiver* receiver)
		{
			std::vector<TimeSpan> spans;
			if (receiver->getSchedule().empty())
			{
				spans.emplace_back(params::startTime(), params::endTime());
				return spans;
			}

			for (const auto& period : receiver->getSchedule())
			{
				const RealType start = std::max(params::startTime(), period.start);
				const RealType end = std::min(params::endTime(), period.end);
				if (start < end)
				{
					spans.emplace_back(start, end);
				}
			}
			return spans;
		}

		/// Adds LO-active intervals for one source intersected with a receiver-active interval.
		void appendDechirpSourceIntervals(const core::ActiveStreamingSource& source, const TimeSpan& receiver_span,
										  std::vector<TimeSpan>& output)
		{
			const RealType clipped_start = std::max({receiver_span.first, source.segment_start, params::startTime()});
			const RealType clipped_end = std::min({receiver_span.second, source.segment_end, params::endTime()});
			if (clipped_start >= clipped_end)
			{
				return;
			}

			if (source.kind == core::StreamingWaveformKind::FmcwLinear)
			{
				if (source.chirp_period <= 0.0 || source.chirp_duration <= 0.0)
				{
					return;
				}
				auto chirp_index = clipped_start <= source.segment_start
					? std::size_t{0}
					: static_cast<std::size_t>(
						  std::floor((clipped_start - source.segment_start) / source.chirp_period));
				while (true)
				{
					if (source.chirp_count.has_value() && chirp_index >= *source.chirp_count)
					{
						return;
					}
					const RealType chirp_start =
						source.segment_start + static_cast<RealType>(chirp_index) * source.chirp_period;
					if (chirp_start >= clipped_end)
					{
						return;
					}
					const RealType chirp_end = std::min(chirp_start + source.chirp_duration, source.segment_end);
					const RealType active_start = std::max(chirp_start, clipped_start);
					const RealType active_end = std::min(chirp_end, clipped_end);
					if (active_start < active_end)
					{
						output.emplace_back(active_start, active_end);
					}
					++chirp_index;
				}
			}

			output.emplace_back(clipped_start, clipped_end);
		}

		/// Returns exact LO-active time spans for a dechirped receiver.
		std::vector<TimeSpan> dechirpActiveTimeSpans(const radar::Receiver* receiver)
		{
			std::vector<TimeSpan> spans;
			const auto receiver_spans = receiverActiveTimeSpans(receiver);
			for (const auto& receiver_span : receiver_spans)
			{
				for (const auto& source : receiver->getDechirpSources())
				{
					appendDechirpSourceIntervals(source, receiver_span, spans);
				}
			}
			normalizeTimeSpans(spans);
			return spans;
		}

		/// Converts time spans to sample spans.
		std::vector<pipeline::SampleSpan> sampleSpansFromTimeSpans(const std::vector<TimeSpan>& spans,
																   const std::size_t total_samples,
																   const RealType sample_rate)
		{
			std::vector<pipeline::SampleSpan> sample_spans;
			for (const auto& span : spans)
			{
				if (const auto sample_span = timeSpanToSampleSpan(span, total_samples, sample_rate))
				{
					if (!sample_spans.empty() && sample_spans.back().end_exclusive >= sample_span->start)
					{
						sample_spans.back().end_exclusive =
							std::max(sample_spans.back().end_exclusive, sample_span->end_exclusive);
					}
					else
					{
						sample_spans.push_back(*sample_span);
					}
				}
			}
			return sample_spans;
		}

		/// Applies thermal noise only inside selected spans.
		void applyThermalNoiseToSpans(std::vector<ComplexType>& iq_buffer, const RealType noise_temperature,
									  std::mt19937& rng_engine, std::span<const pipeline::SampleSpan> spans)
		{
			for (const auto& span : spans)
			{
				if (span.start >= span.end_exclusive || span.start >= iq_buffer.size())
				{
					continue;
				}
				const auto end = std::min(span.end_exclusive, iq_buffer.size());
				applyThermalNoise(std::span(iq_buffer).subspan(span.start, end - span.start), noise_temperature,
								  rng_engine);
			}
		}

		/// Builds output metadata for a streaming receiver result file.
		core::OutputFileMetadata
		buildStreamingMetadata(const radar::Receiver* receiver, const std::string& hdf5_filename,
							   const std::size_t total_samples,
							   const std::vector<core::ActiveStreamingSource>& streaming_sources,
							   const RealType output_sample_rate, const std::vector<TimeSpan>& dechirp_time_spans = {})
		{
			core::OutputFileMetadata metadata{.receiver_id = receiver->getId(),
											  .receiver_name = receiver->getName(),
											  .mode = receiver->getMode() == radar::OperationMode::FMCW_MODE ? "fmcw"
																											 : "cw",
											  .path = hdf5_filename,
											  .sampling_rate = output_sample_rate,
											  .total_samples = static_cast<std::uint64_t>(total_samples),
											  .sample_start = 0,
											  .sample_end_exclusive = static_cast<std::uint64_t>(total_samples)};

			const auto append_segment = [&](const RealType start_time, const RealType end_time)
			{
				const auto start_sample = static_cast<std::uint64_t>(std::min<RealType>(
					static_cast<RealType>(total_samples),
					std::max<RealType>(0.0, std::ceil((start_time - params::startTime()) * output_sample_rate))));
				const auto end_sample = static_cast<std::uint64_t>(std::min<RealType>(
					static_cast<RealType>(total_samples),
					std::max<RealType>(0.0, std::ceil((end_time - params::startTime()) * output_sample_rate))));
				if (start_sample < end_sample)
				{
					core::StreamingSegmentMetadata segment{.start_time = start_time,
														   .end_time = end_time,
														   .sample_count = end_sample - start_sample,
														   .sample_start = start_sample,
														   .sample_end_exclusive = end_sample};
					metadata.streaming_segments.push_back(std::move(segment));
				}
			};

			if (receiver->isDechirpEnabled())
			{
				for (const auto& span : dechirp_time_spans)
				{
					append_segment(span.first, span.second);
				}
			}
			else
			{
				const auto& schedule = receiver->getSchedule();
				if (schedule.empty())
				{
					append_segment(params::startTime(), params::endTime());
				}
				else
				{
					for (const auto& period : schedule)
					{
						const RealType start = std::max(params::startTime(), period.start);
						const RealType end = std::min(params::endTime(), period.end);
						if (start < end)
						{
							append_segment(start, end);
						}
					}
				}
			}

			metadata.fmcw_sources = buildFmcwSources(streaming_sources);
			if (metadata.fmcw_sources.size() == 1)
			{
				metadata.fmcw = metadata.fmcw_sources.front().waveform;
				for (const auto& streaming_source : streaming_sources)
				{
					if (streaming_source.is_fmcw && streaming_source.transmitter != nullptr &&
						streaming_source.transmitter->getId() == metadata.fmcw_sources.front().transmitter_id)
					{
						annotateStreamingSegmentsForSingleFmcwSource(metadata, streaming_source);
					}
				}
			}

			metadata.fmcw_dechirp_mode = std::string(radar::dechirpModeToken(receiver->getDechirpMode()));
			if (receiver->isDechirpEnabled())
			{
				const auto& reference = receiver->getDechirpReference();
				metadata.fmcw_dechirp_reference_source =
					std::string(radar::dechirpReferenceSourceToken(reference.source));
				if (reference.source == radar::Receiver::DechirpReferenceSource::Attached ||
					reference.source == radar::Receiver::DechirpReferenceSource::Transmitter)
				{
					metadata.fmcw_dechirp_reference_transmitter_id = reference.transmitter_id;
					metadata.fmcw_dechirp_reference_transmitter_name = reference.transmitter_name;
				}
				else if (reference.source == radar::Receiver::DechirpReferenceSource::Custom)
				{
					metadata.fmcw_dechirp_reference_waveform_id = reference.waveform_id;
					metadata.fmcw_dechirp_reference_waveform_name = reference.waveform_name;
					if (!receiver->getDechirpSources().empty())
					{
						metadata.fmcw_dechirp_reference_waveform =
							buildFmcwMetadata(receiver->getDechirpSources().front());
					}
				}
			}

			return metadata;
		}

		/// Returns the UI-facing finalization label for a streaming receiver.
		std::string streamingFinalizerLabel(const radar::Receiver* receiver)
		{
			if (receiver->getMode() == radar::OperationMode::FMCW_MODE)
			{
				return "FMCW";
			}
			return "CW";
		}
	}

	void runPulsedFinalizer(radar::Receiver* receiver, const std::vector<std::unique_ptr<radar::Target>>* targets,
							std::shared_ptr<core::ProgressReporter> reporter, const std::string& output_dir,
							std::shared_ptr<core::OutputMetadataCollector> metadata_collector)
	{
		const auto timing_model = receiver->getTiming()->clone();
		if (!timing_model)
		{
			LOG(logging::Level::FATAL, "Failed to clone timing model for receiver '{}'", receiver->getName());
			return;
		}

		std::filesystem::path out_path(output_dir);
		if (!std::filesystem::exists(out_path))
		{
			std::filesystem::create_directories(out_path);
		}
		const auto hdf5_filename = (out_path / std::format("{}_results.h5", receiver->getName())).string();
		core::OutputFileMetadata file_metadata{.receiver_id = receiver->getId(),
											   .receiver_name = receiver->getName(),
											   .mode = "pulsed",
											   .path = hdf5_filename,
											   .sampling_rate = params::rate()};

		std::unique_ptr<HighFive::File> h5_file;
		{
			std::scoped_lock lock(serial::hdf5_global_mutex);
			h5_file = std::make_unique<HighFive::File>(hdf5_filename, HighFive::File::Truncate);
		}

		unsigned chunk_index = 0;

		LOG(logging::Level::INFO, "Finalizer thread started for receiver '{}'. Outputting to '{}'.",
			receiver->getName(), hdf5_filename);

		auto last_report_time = std::chrono::steady_clock::now();
		const auto report_interval = std::chrono::milliseconds(100);
		const RealType rate = params::rate() * params::oversampleRatio();
		const RealType dt = 1.0 / rate;
		core::ReceiverTrackerCache streaming_tracker_cache;

		while (true)
		{
			core::RenderingJob job;
			if (!receiver->waitAndDequeueFinalizerJob(job))
			{
				break; // Shutdown signal received
			}

			const auto window_samples = static_cast<unsigned>(std::ceil(job.duration * rate));
			std::vector pnoise(window_samples, 0.0);

			RealType actual_start = job.ideal_start_time;
			RealType frac_delay = 0.0;

			if (timing_model->isEnabled())
			{
				pipeline::advanceTimingModel(timing_model.get(), receiver, rate);
				std::ranges::generate(pnoise, [&] { return timing_model->getNextSample(); });
				std::tie(actual_start, frac_delay) = pipeline::calculateJitteredStart(
					job.ideal_start_time, pnoise[0], timing_model->getFrequency(), rate);
			}

			std::vector<ComplexType> window_buffer(window_samples);

			applyThermalNoise(window_buffer, receiver->getNoiseTemperature(receiver->getRotation(actual_start)),
							  receiver->getRngEngine());

			pipeline::applyStreamingInterference(window_buffer, actual_start, dt, receiver,
												 job.active_streaming_sources, targets, streaming_tracker_cache);

			renderWindow(window_buffer, job.duration, actual_start, frac_delay, job.responses);

			if (timing_model->isEnabled())
			{
				pipeline::addPhaseNoiseToWindow(pnoise, window_buffer);
			}

			const RealType fullscale = pipeline::applyDownsamplingAndQuantization(window_buffer);

			const auto current_chunk_index = chunk_index++;
			const auto sample_start = file_metadata.total_samples;
			core::PulseChunkMetadata chunk_metadata{.chunk_index = current_chunk_index,
													.i_dataset = std::format("chunk_{:06}_I", current_chunk_index),
													.q_dataset = std::format("chunk_{:06}_Q", current_chunk_index),
													.start_time = actual_start,
													.sample_count = static_cast<std::uint64_t>(window_buffer.size()),
													.sample_start = sample_start,
													.sample_end_exclusive = sample_start +
														static_cast<std::uint64_t>(window_buffer.size())};

			serial::addChunkToFile(*h5_file, window_buffer, actual_start, fullscale, current_chunk_index,
								   &chunk_metadata);
			file_metadata.chunks.push_back(std::move(chunk_metadata));
			file_metadata.total_samples = file_metadata.chunks.back().sample_end_exclusive;

			if (reporter)
			{
				const auto now = std::chrono::steady_clock::now();
				if ((now - last_report_time) >= report_interval)
				{
					reporter->report(std::format("Exporting {}: Chunk {}", receiver->getName(), chunk_index),
									 static_cast<int>(chunk_index), 0);
					last_report_time = now;
				}
			}
		}

		finalizePulsedMetadata(file_metadata);
		{
			std::scoped_lock lock(serial::hdf5_global_mutex);
			serial::writeOutputFileMetadataAttributes(*h5_file, file_metadata);
		}

		{
			// Safe destruction of the HDF5 object inside a lock
			std::scoped_lock lock(serial::hdf5_global_mutex);
			h5_file.reset();
		}

		if (metadata_collector)
		{
			metadata_collector->addFile(std::move(file_metadata));
		}

		if (reporter)
		{
			reporter->report(std::format("Finished Exporting {}", receiver->getName()), 100, 100);
		}
		LOG(logging::Level::INFO, "Finalizer thread for receiver '{}' finished.", receiver->getName());
	}

	void finalizeStreamingReceiver(radar::Receiver* receiver, pool::ThreadPool* /*pool*/,
								   std::shared_ptr<core::ProgressReporter> reporter, const std::string& output_dir,
								   std::shared_ptr<core::OutputMetadataCollector> metadata_collector,
								   std::vector<core::ActiveStreamingSource> streaming_sources)
	{
		LOG(logging::Level::INFO, "Finalization task started for streaming receiver '{}'.", receiver->getName());

		auto& iq_buffer = receiver->getMutableStreamingData();
		if (iq_buffer.empty())
		{
			LOG(logging::Level::INFO, "No streaming data to finalize for receiver '{}'.", receiver->getName());
			return;
		}

		if (reporter)
		{
			reporter->report(
				std::format("Finalizing {} Receiver {}", streamingFinalizerLabel(receiver), receiver->getName()), 0,
				100);
		}

		const bool dechirped = receiver->isDechirpEnabled();
		const RealType output_sample_rate =
			dechirped ? params::rate() * static_cast<RealType>(params::oversampleRatio()) : params::rate();
		const auto dechirp_time_spans = dechirped ? dechirpActiveTimeSpans(receiver) : std::vector<TimeSpan>{};
		const auto dechirp_sample_spans = dechirped
			? sampleSpansFromTimeSpans(dechirp_time_spans, iq_buffer.size(), output_sample_rate)
			: std::vector<pipeline::SampleSpan>{};

		if (reporter)
		{
			reporter->report(std::format("Rendering Interference for {}", receiver->getName()), 25, 100);
		}
		if (dechirped)
		{
			pipeline::applyPulsedInterference(iq_buffer, receiver->getPulsedInterferenceLog(), dechirp_sample_spans,
											  output_sample_rate);
		}
		else
		{
			pipeline::applyPulsedInterference(iq_buffer, receiver->getPulsedInterferenceLog());
		}

		if (reporter)
		{
			reporter->report(std::format("Applying Noise for {}", receiver->getName()), 50, 100);
		}
		if (dechirped)
		{
			applyThermalNoiseToSpans(iq_buffer, receiver->getNoiseTemperature(), receiver->getRngEngine(),
									 dechirp_sample_spans);
		}
		else
		{
			applyThermalNoise(iq_buffer, receiver->getNoiseTemperature(), receiver->getRngEngine());
		}

		const RealType fullscale =
			dechirped ? quantizeAndScaleWindow(iq_buffer) : pipeline::applyDownsamplingAndQuantization(iq_buffer);

		if (reporter)
		{
			reporter->report(std::format("Writing HDF5 for {}", receiver->getName()), 75, 100);
		}

		std::filesystem::path out_path(output_dir);
		if (!std::filesystem::exists(out_path))
		{
			std::filesystem::create_directories(out_path);
		}
		const auto hdf5_filename = (out_path / std::format("{}_results.h5", receiver->getName())).string();
		auto file_metadata = buildStreamingMetadata(receiver, hdf5_filename, iq_buffer.size(), streaming_sources,
													output_sample_rate, dechirp_time_spans);
		pipeline::exportStreamingToHdf5(hdf5_filename, iq_buffer, fullscale, receiver->getTiming()->getFrequency(),
										&file_metadata, output_sample_rate);
		if (metadata_collector)
		{
			metadata_collector->addFile(std::move(file_metadata));
		}

		if (reporter)
		{
			reporter->report(std::format("Finalized {}", receiver->getName()), 100, 100);
		}
	}

}
