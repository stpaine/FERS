// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "finalizer.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <format>
#include <optional>
#include <span>
#include <stdexcept>
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
#include "signal/radar_signal.h"
#include "timing/timing.h"

namespace processing
{
	namespace
	{
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

		core::SfcwMetadata buildSfcwMetadata(const core::ActiveStreamingSource& source)
		{
			if (source.sfcw == nullptr)
			{
				return {};
			}
			const RealType effective_bandwidth = source.sfcw->effectiveBandwidth();
			const RealType step_size = std::abs(source.sfcw->getStepSize());
			return core::SfcwMetadata{
				.carrier_frequency = source.carrier_freq,
				.start_frequency_offset = source.sfcw->getStartFrequencyOffset(),
				.step_size = source.sfcw->getStepSize(),
				.step_count = static_cast<std::uint64_t>(source.sfcw->getStepCount()),
				.dwell_time = source.sfcw->getDwellTime(),
				.step_period = source.sfcw->getStepPeriod(),
				.sweep_count = source.sfcw->getSweepCount().has_value()
					? std::optional<std::uint64_t>(static_cast<std::uint64_t>(*source.sfcw->getSweepCount()))
					: std::nullopt,
				.first_frequency = source.sfcw->firstFrequency(source.carrier_freq),
				.last_frequency = source.sfcw->lastFrequency(source.carrier_freq),
				.frequency_span = source.sfcw->frequencySpan(),
				.effective_bandwidth = effective_bandwidth,
				.range_resolution = effective_bandwidth > 0.0 ? params::c() / (2.0 * effective_bandwidth) : 0.0,
				.unambiguous_range = step_size > 0.0 ? params::c() / (2.0 * step_size) : 0.0};
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

		core::SfcwSourceSegmentMetadata buildSfcwSourceSegment(const core::ActiveStreamingSource& source)
		{
			const RealType active_start = std::max(params::startTime(), source.segment_start);
			const RealType active_end = std::min(params::endTime(), source.segment_end);
			return core::SfcwSourceSegmentMetadata{
				.start_time = source.segment_start,
				.end_time = source.segment_end,
				.first_step_start_time = core::firstSfcwStepStart(source, active_start, active_end),
				.emitted_step_count = core::countSfcwStepStarts(source, active_start, active_end)};
		}

		std::vector<core::SfcwSourceMetadata>::iterator findSfcwSource(std::vector<core::SfcwSourceMetadata>& sources,
																	   const SimId transmitter_id,
																	   const SimId waveform_id)
		{
			return std::ranges::find_if(
				sources, [&](const core::SfcwSourceMetadata& source)
				{ return source.transmitter_id == transmitter_id && source.waveform_id == waveform_id; });
		}

		std::vector<core::SfcwSourceMetadata>
		buildSfcwSources(const std::vector<core::ActiveStreamingSource>& streaming_sources)
		{
			std::vector<core::SfcwSourceMetadata> sfcw_sources;
			for (const auto& streaming_source : streaming_sources)
			{
				if (!streaming_source.is_sfcw || streaming_source.transmitter == nullptr)
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
				auto existing = findSfcwSource(sfcw_sources, transmitter_id, waveform_id);
				if (existing == sfcw_sources.end())
				{
					core::SfcwSourceMetadata source{.transmitter_id = transmitter_id,
													.transmitter_name = streaming_source.transmitter->getName(),
													.waveform_id = waveform_id,
													.waveform_name = signal->getName(),
													.waveform = buildSfcwMetadata(streaming_source)};
					source.segments.push_back(buildSfcwSourceSegment(streaming_source));
					sfcw_sources.push_back(std::move(source));
					continue;
				}

				existing->segments.push_back(buildSfcwSourceSegment(streaming_source));
			}
			return sfcw_sources;
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

		void annotateStreamingSegmentsForSingleSfcwSource(core::OutputFileMetadata& metadata,
														  const core::ActiveStreamingSource& source)
		{
			for (auto& segment : metadata.streaming_segments)
			{
				const RealType active_start = std::max(segment.start_time, source.segment_start);
				const RealType active_end = std::min(segment.end_time, source.segment_end);
				const auto first_step = core::firstSfcwStepStart(source, active_start, active_end);
				const auto emitted = core::countSfcwStepStarts(source, active_start, active_end);
				if (first_step.has_value() || emitted > 0)
				{
					segment.first_sfcw_step_start_time = first_step;
					segment.emitted_sfcw_step_count = emitted;
				}
			}
		}

		/// Half-open time interval in simulation seconds.
		using TimeSpan = std::pair<RealType, RealType>;

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

		void appendStreamingSegment(core::OutputFileMetadata& metadata, const std::size_t total_samples,
									const RealType output_sample_rate, const RealType start_time,
									const RealType end_time)
		{
			const auto start_sample = static_cast<std::uint64_t>(std::min<RealType>(
				static_cast<RealType>(total_samples),
				std::max<RealType>(0.0, std::ceil((start_time - params::startTime()) * output_sample_rate))));
			const auto end_sample = static_cast<std::uint64_t>(std::min<RealType>(
				static_cast<RealType>(total_samples),
				std::max<RealType>(0.0, std::ceil((end_time - params::startTime()) * output_sample_rate))));
			if (start_sample < end_sample)
			{
				const core::StreamingSegmentMetadata segment{.start_time = start_time,
															 .end_time = end_time,
															 .sample_count = end_sample - start_sample,
															 .sample_start = start_sample,
															 .sample_end_exclusive = end_sample};
				metadata.streaming_segments.push_back(segment);
			}
		}

		void appendScheduledStreamingSegments(core::OutputFileMetadata& metadata, const radar::Receiver* receiver,
											  const std::size_t total_samples, const RealType output_sample_rate)
		{
			const auto& schedule = receiver->getSchedule();
			if (schedule.empty())
			{
				appendStreamingSegment(metadata, total_samples, output_sample_rate, params::startTime(),
									   params::endTime());
				return;
			}

			for (const auto& period : schedule)
			{
				const RealType start = std::max(params::startTime(), period.start);
				const RealType end = std::min(params::endTime(), period.end);
				if (start < end)
				{
					appendStreamingSegment(metadata, total_samples, output_sample_rate, start, end);
				}
			}
		}

		void appendStreamingSegments(core::OutputFileMetadata& metadata, const radar::Receiver* receiver,
									 const std::size_t total_samples, const RealType output_sample_rate,
									 const std::vector<TimeSpan>& dechirp_time_spans)
		{
			if (!receiver->isDechirpEnabled())
			{
				appendScheduledStreamingSegments(metadata, receiver, total_samples, output_sample_rate);
				return;
			}

			for (const auto& span : dechirp_time_spans)
			{
				appendStreamingSegment(metadata, total_samples, output_sample_rate, span.first, span.second);
			}
		}

		void annotateSingleFmcwSource(core::OutputFileMetadata& metadata,
									  const std::vector<core::ActiveStreamingSource>& streaming_sources)
		{
			metadata.fmcw_sources = buildFmcwSources(streaming_sources);
			if (metadata.fmcw_sources.size() != 1)
			{
				return;
			}

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

		void annotateSingleSfcwSource(core::OutputFileMetadata& metadata,
									  const std::vector<core::ActiveStreamingSource>& streaming_sources)
		{
			metadata.sfcw_sources = buildSfcwSources(streaming_sources);
			if (metadata.sfcw_sources.size() != 1)
			{
				return;
			}

			metadata.sfcw = metadata.sfcw_sources.front().waveform;
			for (const auto& streaming_source : streaming_sources)
			{
				if (streaming_source.is_sfcw && streaming_source.transmitter != nullptr &&
					streaming_source.transmitter->getId() == metadata.sfcw_sources.front().transmitter_id)
				{
					annotateStreamingSegmentsForSingleSfcwSource(metadata, streaming_source);
				}
			}
		}

		void populateFmcwIfMetadata(core::OutputFileMetadata& metadata, const radar::Receiver* receiver)
		{
			const auto& if_request = receiver->getFmcwIfChainRequest();
			const auto& if_plan = receiver->getFmcwIfResamplerPlan();
			metadata.fmcw_if_legacy_full_rate = !if_request.sample_rate_hz.has_value();
			metadata.fmcw_if_decimation_enabled = if_plan.has_value();
			if (if_request.sample_rate_hz.has_value())
			{
				metadata.fmcw_if_requested_sample_rate = if_request.sample_rate_hz;
			}
			if (!if_plan.has_value())
			{
				return;
			}

			metadata.fmcw_if_sample_rate = if_plan->actual_output_sample_rate_hz;
			metadata.fmcw_if_input_sample_rate = if_plan->input_sample_rate_hz;
			metadata.fmcw_if_resample_numerator = static_cast<unsigned>(if_plan->overall_ratio.numerator);
			metadata.fmcw_if_resample_denominator = static_cast<unsigned>(if_plan->overall_ratio.denominator);
			metadata.fmcw_if_decimation_factor = if_plan->actual_output_sample_rate_hz > 0.0
				? if_plan->input_sample_rate_hz / if_plan->actual_output_sample_rate_hz
				: 0.0;
			metadata.fmcw_if_filter_bandwidth = if_plan->filter_bandwidth_hz;
			metadata.fmcw_if_filter_transition_width = if_plan->filter_transition_width_hz;
			metadata.fmcw_if_filter_stopband = if_plan->stopband_attenuation_db;
			metadata.fmcw_if_filter_group_delay_seconds = if_plan->group_delay_seconds;
			metadata.fmcw_if_compensated_integer_delay_samples = if_plan->warmup_discard_samples;
			metadata.fmcw_if_compensated_fractional_delay_samples = if_plan->fractional_output_delay_samples;
			metadata.fmcw_if_warmup_discard_samples = if_plan->warmup_discard_samples;
			metadata.fmcw_if_phase_refinement = static_cast<unsigned>(if_plan->phase_refinement);
			metadata.fmcw_if_timing_error_seconds = if_plan->estimated_timing_error_seconds;
			metadata.fmcw_if_phase_error_radians = if_plan->estimated_phase_error_radians;
			metadata.fmcw_if_noise_variance =
				params::boltzmannK() * receiver->getNoiseTemperature() * if_plan->actual_output_sample_rate_hz;
			metadata.fmcw_if_group_delay_compensated = if_plan->group_delay_compensated;
		}

		void populateDechirpReferenceMetadata(core::OutputFileMetadata& metadata, const radar::Receiver* receiver)
		{
			const auto& reference = receiver->getDechirpReference();
			metadata.fmcw_dechirp_reference_source = std::string(radar::dechirpReferenceSourceToken(reference.source));
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
					metadata.fmcw_dechirp_reference_waveform = buildFmcwMetadata(receiver->getDechirpSources().front());
				}
			}
		}

		/// Builds output metadata for a streaming receiver result file.
		core::OutputFileMetadata
		buildStreamingMetadata(const radar::Receiver* receiver, const std::string& hdf5_filename,
							   const std::size_t total_samples,
							   const std::vector<core::ActiveStreamingSource>& streaming_sources,
							   const RealType output_sample_rate, const std::vector<TimeSpan>& dechirp_time_spans = {})
		{
			core::OutputFileMetadata metadata{
				.receiver_id = receiver->getId(),
				.receiver_name = receiver->getName(),
				.mode = receiver->getMode() == radar::OperationMode::FMCW_MODE
					? "fmcw"
					: (receiver->getMode() == radar::OperationMode::SFCW_MODE ? "sfcw" : "cw"),
				.path = hdf5_filename,
				.sampling_rate = output_sample_rate,
				.total_samples = static_cast<std::uint64_t>(total_samples),
				.sample_start = 0,
				.sample_end_exclusive = static_cast<std::uint64_t>(total_samples)};

			appendStreamingSegments(metadata, receiver, total_samples, output_sample_rate, dechirp_time_spans);
			annotateSingleFmcwSource(metadata, streaming_sources);
			annotateSingleSfcwSource(metadata, streaming_sources);

			metadata.fmcw_dechirp_mode = std::string(radar::dechirpModeToken(receiver->getDechirpMode()));
			if (receiver->isDechirpEnabled())
			{
				populateFmcwIfMetadata(metadata, receiver);
				populateDechirpReferenceMetadata(metadata, receiver);
			}

			return metadata;
		}

		/// Converts a receiver mode to the stable sink descriptor token.
		[[nodiscard]] std::string receiverModeToken(const radar::Receiver* receiver)
		{
			switch (receiver->getMode())
			{
			case radar::OperationMode::PULSED_MODE:
				return "pulsed";
			case radar::OperationMode::FMCW_MODE:
				return "fmcw";
			case radar::OperationMode::SFCW_MODE:
				return "sfcw";
			case radar::OperationMode::CW_MODE:
				return "cw";
			}
			return "unknown";
		}

		[[nodiscard]] std::string coordinateFrameToken(const params::CoordinateFrame frame)
		{
			switch (frame)
			{
			case params::CoordinateFrame::ENU:
				return "ENU";
			case params::CoordinateFrame::UTM:
				return "UTM";
			case params::CoordinateFrame::ECEF:
				return "ECEF";
			}
			return "ENU";
		}

		[[nodiscard]] core::ReceiverStreamDescriptor::CoordinateContext buildCoordinateContext()
		{
			return core::ReceiverStreamDescriptor::CoordinateContext{
				.frame = coordinateFrameToken(params::coordinateFrame()),
				.origin_latitude = params::originLatitude(),
				.origin_longitude = params::originLongitude(),
				.origin_altitude = params::originAltitude(),
				.utm_zone = params::utmZone(),
				.utm_north_hemisphere = params::utmNorthHemisphere()};
		}

		[[nodiscard]] core::ReceiverStreamDescriptor::PlatformState
		buildInitialPlatformState(const radar::Receiver* receiver)
		{
			core::ReceiverStreamDescriptor::PlatformState state;
			const auto* platform = receiver->getPlatform();
			if (platform == nullptr)
			{
				return state;
			}

			const RealType t0 = params::startTime();
			state.platform_id = platform->getId();
			state.platform_name = platform->getName();
			try
			{
				const auto position = platform->getPosition(t0);
				state.position_x = position.x;
				state.position_y = position.y;
				state.position_z = position.z;
			}
			catch (...)
			{
			}
			try
			{
				const auto velocity = platform->getMotionPath()->getVelocity(t0);
				state.velocity_x = velocity.x;
				state.velocity_y = velocity.y;
				state.velocity_z = velocity.z;
			}
			catch (...)
			{
			}
			try
			{
				const auto rotation = platform->getRotation(t0);
				state.azimuth = rotation.azimuth;
				state.elevation = rotation.elevation;
			}
			catch (...)
			{
			}
			return state;
		}

		[[nodiscard]] const core::ActiveStreamingSource*
		findFmcwContextSource(const radar::Receiver* receiver,
							  const std::span<const core::ActiveStreamingSource> streaming_sources)
		{
			if (receiver->isDechirpEnabled() && !receiver->getDechirpSources().empty())
			{
				return &receiver->getDechirpSources().front();
			}
			const auto found = std::ranges::find_if(streaming_sources, [](const core::ActiveStreamingSource& source)
													{ return source.is_fmcw; });
			return found == streaming_sources.end() ? nullptr : &*found;
		}

		[[nodiscard]] const radar::Transmitter* attachedTransmitter(const radar::Receiver* receiver)
		{
			return receiver == nullptr ? nullptr : dynamic_cast<const radar::Transmitter*>(receiver->getAttached());
		}

		void populateWaveformIdentity(core::ReceiverStreamDescriptor::PulsedContext& context,
									  const fers_signal::RadarSignal* signal)
		{
			if (signal == nullptr)
			{
				return;
			}
			context.waveform_id = signal->getId();
			context.waveform_name = signal->getName();
			context.carrier_frequency = signal->getCarrier();
			context.power = signal->getPower();
			context.pulse_width = signal->getLength();
			context.native_sample_rate = signal->getRate();
			context.native_sample_count = signal->getSampleCount();
		}

		void populateWaveformIdentity(core::ReceiverStreamDescriptor::CwContext& context,
									  const fers_signal::RadarSignal* signal)
		{
			if (signal == nullptr)
			{
				return;
			}
			context.waveform_id = signal->getId();
			context.waveform_name = signal->getName();
			context.carrier_frequency = signal->getCarrier();
			context.power = signal->getPower();
		}

		void populateWaveformIdentity(core::ReceiverStreamDescriptor::SfcwContext& context,
									  const fers_signal::RadarSignal* signal)
		{
			if (signal == nullptr)
			{
				return;
			}
			context.waveform_id = signal->getId();
			context.waveform_name = signal->getName();
			context.carrier_frequency = signal->getCarrier();
			context.power = signal->getPower();
			const auto* sfcw = signal->getSteppedFrequencySignal();
			if (sfcw == nullptr)
			{
				return;
			}
			const RealType effective_bandwidth = sfcw->effectiveBandwidth();
			const RealType step_size = std::abs(sfcw->getStepSize());
			context.start_frequency_offset = sfcw->getStartFrequencyOffset();
			context.step_size = sfcw->getStepSize();
			context.step_count = static_cast<std::uint64_t>(sfcw->getStepCount());
			context.dwell_time = sfcw->getDwellTime();
			context.step_period = sfcw->getStepPeriod();
			context.sweep_period = sfcw->getSweepPeriod();
			context.sweep_count = sfcw->getSweepCount().has_value()
				? std::optional<std::uint64_t>(static_cast<std::uint64_t>(*sfcw->getSweepCount()))
				: std::nullopt;
			context.first_frequency = sfcw->firstFrequency(signal->getCarrier());
			context.last_frequency = sfcw->lastFrequency(signal->getCarrier());
			context.frequency_span = sfcw->frequencySpan();
			context.effective_bandwidth = effective_bandwidth;
			context.range_resolution = effective_bandwidth > 0.0 ? params::c() / (2.0 * effective_bandwidth) : 0.0;
			context.unambiguous_range = step_size > 0.0 ? params::c() / (2.0 * step_size) : 0.0;
		}

		[[nodiscard]] core::ReceiverStreamDescriptor::PulsedContext buildPulsedContext(const radar::Receiver* receiver)
		{
			core::ReceiverStreamDescriptor::PulsedContext context;
			if (receiver == nullptr || receiver->getMode() != radar::OperationMode::PULSED_MODE)
			{
				return context;
			}

			context.present = true;
			context.window_length = receiver->getWindowLength();
			context.window_prf = receiver->getWindowPrf();
			context.window_skip = receiver->getWindowSkip();
			context.window_count = receiver->getWindowCount();
			if (const auto* transmitter = attachedTransmitter(receiver); transmitter != nullptr)
			{
				populateWaveformIdentity(context, transmitter->getSignal());
			}
			if (context.carrier_frequency == 0.0)
			{
				if (const auto timing = receiver->getTiming(); timing)
				{
					context.carrier_frequency = timing->getFrequency();
				}
			}
			return context;
		}

		[[nodiscard]] core::ReceiverStreamDescriptor::CwContext buildCwContext(const radar::Receiver* receiver)
		{
			core::ReceiverStreamDescriptor::CwContext context;
			if (receiver == nullptr || receiver->getMode() != radar::OperationMode::CW_MODE)
			{
				return context;
			}

			context.present = true;
			if (const auto* transmitter = attachedTransmitter(receiver); transmitter != nullptr)
			{
				populateWaveformIdentity(context, transmitter->getSignal());
			}
			if (context.carrier_frequency == 0.0)
			{
				if (const auto timing = receiver->getTiming(); timing)
				{
					context.carrier_frequency = timing->getFrequency();
				}
			}
			return context;
		}

		[[nodiscard]] core::ReceiverStreamDescriptor::SfcwContext
		buildSfcwContext(const radar::Receiver* receiver,
						 const std::span<const core::ActiveStreamingSource> streaming_sources)
		{
			core::ReceiverStreamDescriptor::SfcwContext context;
			if (receiver == nullptr || receiver->getMode() != radar::OperationMode::SFCW_MODE)
			{
				return context;
			}

			context.present = true;
			if (const auto* transmitter = attachedTransmitter(receiver); transmitter != nullptr)
			{
				populateWaveformIdentity(context, transmitter->getSignal());
			}
			if (context.waveform_id == 0)
			{
				const auto found = std::ranges::find_if(streaming_sources, [](const core::ActiveStreamingSource& source)
														{ return source.is_sfcw && source.transmitter != nullptr; });
				if (found != streaming_sources.end())
				{
					populateWaveformIdentity(context, found->transmitter->getSignal());
				}
			}
			return context;
		}

		[[nodiscard]] core::ReceiverStreamDescriptor::FmcwContext
		buildFmcwContext(const radar::Receiver* receiver,
						 const std::span<const core::ActiveStreamingSource> streaming_sources)
		{
			core::ReceiverStreamDescriptor::FmcwContext context;
			if (receiver == nullptr || receiver->getMode() != radar::OperationMode::FMCW_MODE)
			{
				return context;
			}

			context.dechirp_mode = std::string(radar::dechirpModeToken(receiver->getDechirpMode()));
			const auto& reference = receiver->getDechirpReference();
			context.dechirp_reference_source = std::string(radar::dechirpReferenceSourceToken(reference.source));
			context.dechirp_reference_transmitter_id = reference.transmitter_id;
			context.dechirp_reference_transmitter_name = reference.transmitter_name;
			context.dechirp_reference_waveform_id = reference.waveform_id;
			context.dechirp_reference_waveform_name = reference.waveform_name;

			const auto* source = findFmcwContextSource(receiver, streaming_sources);
			if (source == nullptr)
			{
				return context;
			}

			const auto waveform = buildFmcwMetadata(*source);
			context.present = true;
			context.waveform_shape = waveform.waveform_shape;
			context.chirp_bandwidth = waveform.chirp_bandwidth;
			context.chirp_duration = waveform.chirp_duration;
			context.chirp_period = waveform.chirp_period;
			context.chirp_rate = waveform.chirp_rate;
			context.chirp_rate_signed = waveform.chirp_rate_signed;
			context.sweep_direction =
				source->kind == core::StreamingWaveformKind::FmcwTriangle ? "up_down" : waveform.chirp_direction;
			context.start_frequency_offset = waveform.start_frequency_offset;
			context.triangle_period = waveform.triangle_period;
			context.chirp_count = waveform.chirp_count;
			context.triangle_count = waveform.triangle_count;
			return context;
		}
	}

	core::OutputFileMetadata buildStreamingOutputMetadata(
		const radar::Receiver* receiver, const std::string& output_path, const std::size_t total_samples,
		const std::vector<core::ActiveStreamingSource>& streaming_sources, const RealType output_sample_rate)
	{
		const auto dechirp_time_spans =
			receiver->isDechirpEnabled() ? dechirpActiveTimeSpans(receiver) : std::vector<TimeSpan>{};
		return buildStreamingMetadata(receiver, output_path, total_samples, streaming_sources, output_sample_rate,
									  dechirp_time_spans);
	}

	core::ReceiverStreamDescriptor
	buildReceiverStreamDescriptor(const radar::Receiver* receiver, const RealType sample_rate,
								  const std::span<const core::ActiveStreamingSource> streaming_sources)
	{
		core::ReceiverStreamDescriptor descriptor{.receiver_id = receiver->getId(),
												  .receiver_name = receiver->getName(),
												  .mode = receiverModeToken(receiver),
												  .sample_rate = sample_rate,
												  .bandwidth = sample_rate > 0.0 ? sample_rate / 2.0 : 0.0,
												  .dechirped = receiver->isDechirpEnabled(),
												  .if_resampled = receiver->getFmcwIfResamplerPlan().has_value(),
												  .adc_bits = params::adcBits(),
												  .coordinate = buildCoordinateContext(),
												  .initial_platform_state = buildInitialPlatformState(receiver),
												  .pulsed = buildPulsedContext(receiver),
												  .cw = buildCwContext(receiver),
												  .fmcw = buildFmcwContext(receiver, streaming_sources),
												  .sfcw = buildSfcwContext(receiver, streaming_sources)};
		if (const auto timing = receiver->getTiming(); timing)
		{
			descriptor.reference_frequency = timing->getFrequency();
		}
		return descriptor;
	}

	core::ReceiverSampleBlock buildReceiverSampleBlock(const radar::Receiver* receiver,
													   const RealType first_sample_time, const RealType sample_rate,
													   const std::span<const ComplexType> samples,
													   const std::uint64_t sample_start,
													   std::shared_ptr<const core::OutputFileMetadata> file_metadata)
	{
		return buildReceiverSampleBlock(receiver, first_sample_time, sample_rate, samples, sample_start,
										std::span<const core::ActiveStreamingSource>{}, std::move(file_metadata));
	}

	core::ReceiverSampleBlock
	buildReceiverSampleBlock(const radar::Receiver* receiver, const RealType first_sample_time,
							 const RealType sample_rate, const std::span<const ComplexType> samples,
							 const std::uint64_t sample_start,
							 const std::span<const core::ActiveStreamingSource> streaming_sources,
							 std::shared_ptr<const core::OutputFileMetadata> file_metadata)
	{
		return core::ReceiverSampleBlock{.stream =
											 buildReceiverStreamDescriptor(receiver, sample_rate, streaming_sources),
										 .first_sample_time = first_sample_time,
										 .sample_rate = sample_rate,
										 .samples = samples,
										 .sample_start = sample_start,
										 .valid_data = true,
										 .calibrated_time = true,
										 .reference_lock = true,
										 .file_metadata = std::move(file_metadata)};
	}

	void runPulsedFinalizer(radar::Receiver* receiver, const std::vector<std::unique_ptr<radar::Target>>* targets,
							const std::shared_ptr<core::ProgressReporter>& reporter, const std::string& output_dir,
							const std::shared_ptr<core::OutputMetadataCollector>& metadata_collector,
							core::ReceiverOutputSink* output_sink)
	{
		(void)output_dir;
		(void)metadata_collector;
		if (output_sink == nullptr)
		{
			throw std::invalid_argument("runPulsedFinalizer requires a receiver output sink");
		}

		const auto timing_model = receiver->getTiming()->clone();
		if (!timing_model)
		{
			LOG(logging::Level::FATAL, "Failed to clone timing model for receiver '{}'", receiver->getName());
			return;
		}

		const auto sink_stream_id =
			output_sink->registerStream(buildReceiverStreamDescriptor(receiver, params::rate()));
		bool sink_stream_open = false;
		std::uint64_t sink_sample_start = 0;

		unsigned chunk_index = 0;

		LOG(logging::Level::INFO, "Finalizer thread started for receiver '{}'. Routing to output sink.",
			receiver->getName());

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

			pipeline::applyStreamingInterference(window_buffer, actual_start, dt, receiver,
												 job.active_streaming_sources, targets, streaming_tracker_cache);

			renderWindow(window_buffer, job.duration, actual_start, frac_delay, job.responses);

			if (timing_model->isEnabled())
			{
				pipeline::addPhaseNoiseToWindow(pnoise, window_buffer);
			}

			pipeline::applyDownsampling(window_buffer);
			applyThermalNoiseAtSampleRate(window_buffer,
										  receiver->getNoiseTemperature(receiver->getRotation(actual_start)),
										  receiver->getRngEngine(), params::rate());
			if (!sink_stream_open)
			{
				output_sink->openStream(sink_stream_id, actual_start);
				sink_stream_open = true;
			}
			const auto block =
				buildReceiverSampleBlock(receiver, actual_start, params::rate(), window_buffer, sink_sample_start);
			output_sink->submitBlock(block);
			sink_sample_start += static_cast<std::uint64_t>(window_buffer.size());
			++chunk_index;

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

		if (sink_stream_open)
		{
			output_sink->closeStream(sink_stream_id);
		}

		if (reporter)
		{
			reporter->report(std::format("Finished Exporting {}", receiver->getName()), 100, 100);
		}
		LOG(logging::Level::INFO, "Finalizer thread for receiver '{}' finished.", receiver->getName());
	}

}
