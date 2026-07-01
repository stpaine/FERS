// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "output_metadata.h"

#include <nlohmann/json.hpp>
#include <utility>

#include "core/parameters.h"

namespace core
{
	namespace
	{
		/// Converts one pulsed chunk metadata entry to JSON.
		nlohmann::json chunkToJson(const PulseChunkMetadata& chunk)
		{
			return {{"chunk_index", chunk.chunk_index},
					{"i_dataset", chunk.i_dataset},
					{"q_dataset", chunk.q_dataset},
					{"start_time", chunk.start_time},
					{"sample_count", chunk.sample_count},
					{"sample_start", chunk.sample_start},
					{"sample_end_exclusive", chunk.sample_end_exclusive}};
		}

		/// Converts one streaming segment metadata entry to JSON.
		nlohmann::json streamingSegmentToJson(const StreamingSegmentMetadata& segment)
		{
			nlohmann::json result = {{"start_time", segment.start_time},
									 {"end_time", segment.end_time},
									 {"sample_count", segment.sample_count},
									 {"sample_start", segment.sample_start},
									 {"sample_end_exclusive", segment.sample_end_exclusive}};
			if (segment.first_chirp_start_time.has_value())
			{
				result["first_chirp_start_time"] = *segment.first_chirp_start_time;
			}
			if (segment.emitted_chirp_count.has_value())
			{
				result["emitted_chirp_count"] = *segment.emitted_chirp_count;
			}
			if (segment.first_triangle_start_time.has_value())
			{
				result["first_triangle_start_time"] = *segment.first_triangle_start_time;
			}
			if (segment.emitted_triangle_count.has_value())
			{
				result["emitted_triangle_count"] = *segment.emitted_triangle_count;
			}
			if (segment.first_sfcw_step_start_time.has_value())
			{
				result["first_sfcw_step_start_time"] = *segment.first_sfcw_step_start_time;
			}
			if (segment.emitted_sfcw_step_count.has_value())
			{
				result["emitted_sfcw_step_count"] = *segment.emitted_sfcw_step_count;
			}
			return result;
		}

		/// Converts FMCW output metadata to JSON.
		nlohmann::json fmcwToJson(const FmcwMetadata& fmcw)
		{
			nlohmann::json result = {{"chirp_bandwidth", fmcw.chirp_bandwidth},
									 {"chirp_duration", fmcw.chirp_duration},
									 {"waveform_shape", fmcw.waveform_shape},
									 {"chirp_rate", fmcw.chirp_rate},
									 {"start_frequency_offset", fmcw.start_frequency_offset}};
			if (fmcw.waveform_shape == "linear")
			{
				result["chirp_period"] = fmcw.chirp_period;
				result["chirp_rate_signed"] = fmcw.chirp_rate_signed;
				result["chirp_direction"] = fmcw.chirp_direction;
				if (fmcw.chirp_count.has_value())
				{
					result["chirp_count"] = *fmcw.chirp_count;
				}
			}
			else if (fmcw.waveform_shape == "triangle")
			{
				if (fmcw.triangle_period.has_value())
				{
					result["triangle_period"] = *fmcw.triangle_period;
				}
				if (fmcw.triangle_count.has_value())
				{
					result["triangle_count"] = *fmcw.triangle_count;
				}
			}
			return result;
		}

		/// Converts one FMCW source segment metadata entry to JSON.
		nlohmann::json fmcwSourceSegmentToJson(const FmcwSourceSegmentMetadata& segment)
		{
			nlohmann::json result = {{"start_time", segment.start_time}, {"end_time", segment.end_time}};
			if (segment.first_chirp_start_time.has_value())
			{
				result["first_chirp_start_time"] = *segment.first_chirp_start_time;
			}
			if (segment.emitted_chirp_count.has_value())
			{
				result["emitted_chirp_count"] = *segment.emitted_chirp_count;
			}
			if (segment.first_triangle_start_time.has_value())
			{
				result["first_triangle_start_time"] = *segment.first_triangle_start_time;
			}
			if (segment.emitted_triangle_count.has_value())
			{
				result["emitted_triangle_count"] = *segment.emitted_triangle_count;
			}
			return result;
		}

		/// Converts one FMCW source metadata entry to JSON.
		nlohmann::json fmcwSourceToJson(const FmcwSourceMetadata& source)
		{
			nlohmann::json segments = nlohmann::json::array();
			for (const auto& segment : source.segments)
			{
				segments.push_back(fmcwSourceSegmentToJson(segment));
			}

			nlohmann::json result = {{"transmitter_id", source.transmitter_id},
									 {"transmitter_name", source.transmitter_name},
									 {"waveform_id", source.waveform_id},
									 {"waveform_name", source.waveform_name},
									 {"carrier_frequency", source.carrier_frequency},
									 {"segments", segments}};
			result.update(fmcwToJson(source.waveform));
			return result;
		}

		nlohmann::json sfcwToJson(const SfcwMetadata& sfcw)
		{
			nlohmann::json result = {{"carrier_frequency", sfcw.carrier_frequency},
									 {"start_frequency_offset", sfcw.start_frequency_offset},
									 {"step_size", sfcw.step_size},
									 {"step_count", sfcw.step_count},
									 {"dwell_time", sfcw.dwell_time},
									 {"step_period", sfcw.step_period},
									 {"first_frequency", sfcw.first_frequency},
									 {"last_frequency", sfcw.last_frequency},
									 {"frequency_span", sfcw.frequency_span},
									 {"effective_bandwidth", sfcw.effective_bandwidth},
									 {"range_resolution", sfcw.range_resolution},
									 {"unambiguous_range", sfcw.unambiguous_range}};
			if (sfcw.sweep_count.has_value())
			{
				result["sweep_count"] = *sfcw.sweep_count;
			}
			return result;
		}

		nlohmann::json sfcwSourceSegmentToJson(const SfcwSourceSegmentMetadata& segment)
		{
			nlohmann::json result = {{"start_time", segment.start_time}, {"end_time", segment.end_time}};
			if (segment.first_step_start_time.has_value())
			{
				result["first_step_start_time"] = *segment.first_step_start_time;
			}
			if (segment.emitted_step_count.has_value())
			{
				result["emitted_step_count"] = *segment.emitted_step_count;
			}
			return result;
		}

		nlohmann::json sfcwSourceToJson(const SfcwSourceMetadata& source)
		{
			nlohmann::json segments = nlohmann::json::array();
			for (const auto& segment : source.segments)
			{
				segments.push_back(sfcwSourceSegmentToJson(segment));
			}

			nlohmann::json result = {{"transmitter_id", source.transmitter_id},
									 {"transmitter_name", source.transmitter_name},
									 {"waveform_id", source.waveform_id},
									 {"waveform_name", source.waveform_name},
									 {"segments", segments}};
			result.update(sfcwToJson(source.waveform));
			return result;
		}

		template <typename T>
		void addOptional(nlohmann::json& result, const char* key, const std::optional<T>& value)
		{
			if (value.has_value())
			{
				result[key] = *value;
			}
		}

		void addDechirpReferenceJson(nlohmann::json& result, const OutputFileMetadata& file)
		{
			addOptional(result, "fmcw_dechirp_reference_transmitter_id", file.fmcw_dechirp_reference_transmitter_id);
			addOptional(result, "fmcw_dechirp_reference_transmitter_name",
						file.fmcw_dechirp_reference_transmitter_name);
			addOptional(result, "fmcw_dechirp_reference_waveform_id", file.fmcw_dechirp_reference_waveform_id);
			addOptional(result, "fmcw_dechirp_reference_waveform_name", file.fmcw_dechirp_reference_waveform_name);
			if (file.fmcw_dechirp_reference_waveform.has_value())
			{
				result["fmcw_dechirp_reference_waveform"] = fmcwToJson(*file.fmcw_dechirp_reference_waveform);
			}
		}

		void addFmcwIfJson(nlohmann::json& result, const OutputFileMetadata& file)
		{
			result["fmcw_if_decimation_enabled"] = file.fmcw_if_decimation_enabled;
			result["fmcw_if_legacy_full_rate"] = file.fmcw_if_legacy_full_rate;
			addOptional(result, "fmcw_if_requested_sample_rate", file.fmcw_if_requested_sample_rate);
			addOptional(result, "fmcw_if_sample_rate", file.fmcw_if_sample_rate);
			addOptional(result, "fmcw_if_input_sample_rate", file.fmcw_if_input_sample_rate);
			addOptional(result, "fmcw_if_resample_numerator", file.fmcw_if_resample_numerator);
			addOptional(result, "fmcw_if_resample_denominator", file.fmcw_if_resample_denominator);
			addOptional(result, "fmcw_if_decimation_factor", file.fmcw_if_decimation_factor);
			addOptional(result, "fmcw_if_filter_bandwidth", file.fmcw_if_filter_bandwidth);
			addOptional(result, "fmcw_if_filter_transition_width", file.fmcw_if_filter_transition_width);
			addOptional(result, "fmcw_if_filter_stopband", file.fmcw_if_filter_stopband);
			addOptional(result, "fmcw_if_filter_group_delay_seconds", file.fmcw_if_filter_group_delay_seconds);
			addOptional(result, "fmcw_if_compensated_integer_delay_samples",
						file.fmcw_if_compensated_integer_delay_samples);
			addOptional(result, "fmcw_if_compensated_fractional_delay_samples",
						file.fmcw_if_compensated_fractional_delay_samples);
			addOptional(result, "fmcw_if_warmup_discard_samples", file.fmcw_if_warmup_discard_samples);
			addOptional(result, "fmcw_if_phase_refinement", file.fmcw_if_phase_refinement);
			addOptional(result, "fmcw_if_timing_error_seconds", file.fmcw_if_timing_error_seconds);
			addOptional(result, "fmcw_if_phase_error_radians", file.fmcw_if_phase_error_radians);
			addOptional(result, "fmcw_if_noise_variance", file.fmcw_if_noise_variance);
			result["fmcw_if_group_delay_compensated"] = file.fmcw_if_group_delay_compensated;
		}

		/// Converts one output file metadata entry to JSON.
		nlohmann::json fileToJson(const OutputFileMetadata& file)
		{
			nlohmann::json chunks = nlohmann::json::array();
			for (const auto& chunk : file.chunks)
			{
				chunks.push_back(chunkToJson(chunk));
			}

			nlohmann::json streaming_segments = nlohmann::json::array();
			for (const auto& segment : file.streaming_segments)
			{
				streaming_segments.push_back(streamingSegmentToJson(segment));
			}

			nlohmann::json fmcw_sources = nlohmann::json::array();
			for (const auto& source : file.fmcw_sources)
			{
				fmcw_sources.push_back(fmcwSourceToJson(source));
			}

			nlohmann::json sfcw_sources = nlohmann::json::array();
			for (const auto& source : file.sfcw_sources)
			{
				sfcw_sources.push_back(sfcwSourceToJson(source));
			}

			nlohmann::json result = {{"receiver_id", file.receiver_id},
									 {"receiver_name", file.receiver_name},
									 {"mode", file.mode},
									 {"path", file.path},
									 {"sampling_rate", file.sampling_rate},
									 {"total_samples", file.total_samples},
									 {"sample_start", file.sample_start},
									 {"sample_end_exclusive", file.sample_end_exclusive},
									 {"pulse_count", file.pulse_count},
									 {"min_pulse_length_samples", file.min_pulse_length_samples},
									 {"max_pulse_length_samples", file.max_pulse_length_samples},
									 {"uniform_pulse_length", file.uniform_pulse_length},
									 {"chunks", chunks},
									 {"streaming_segments", streaming_segments},
									 {"fmcw_sources", fmcw_sources},
									 {"sfcw_sources", sfcw_sources},
									 {"fmcw_dechirp_mode", file.fmcw_dechirp_mode},
									 {"fmcw_dechirp_reference_source", file.fmcw_dechirp_reference_source}};
			if (file.fmcw.has_value())
			{
				result["fmcw"] = fmcwToJson(*file.fmcw);
			}
			if (file.sfcw.has_value())
			{
				result["sfcw"] = sfcwToJson(*file.sfcw);
			}
			addDechirpReferenceJson(result, file);
			addFmcwIfJson(result, file);
			return result;
		}

		/// Converts one VITA stream metadata entry to JSON.
		nlohmann::json vita49StreamToJson(const Vita49StreamMetadata& stream)
		{
			auto timestampToJson = [](const std::optional<Vita49Timestamp>& timestamp) -> nlohmann::json
			{
				if (!timestamp.has_value())
				{
					return nullptr;
				}
				return {{"integer_seconds", timestamp->integer_seconds},
						{"fractional_picoseconds", timestamp->fractional_picoseconds}};
			};

			return {{"receiver_id", stream.receiver_id},
					{"receiver_name", stream.receiver_name},
					{"stream_id", stream.stream_id},
					{"mode", stream.mode},
					{"sample_rate", stream.sample_rate},
					{"reference_frequency", stream.reference_frequency},
					{"packets_emitted", stream.packets_emitted},
					{"samples_emitted", stream.samples_emitted},
					{"packets_dropped", stream.packets_dropped},
					{"samples_dropped", stream.samples_dropped},
					{"over_range_count", stream.over_range_count},
					{"late_packet_count", stream.late_packet_count},
					{"context_packet_count", stream.context_packet_count},
					{"first_sample_time",
					 stream.first_sample_time.has_value() ? nlohmann::json(*stream.first_sample_time)
														  : nlohmann::json(nullptr)},
					{"end_sample_time",
					 stream.end_sample_time.has_value() ? nlohmann::json(*stream.end_sample_time)
														: nlohmann::json(nullptr)},
					{"first_timestamp", timestampToJson(stream.first_timestamp)},
					{"end_timestamp", timestampToJson(stream.end_timestamp)}};
		}

		/// Converts VITA output metadata to JSON.
		nlohmann::json vita49ToJson(const Vita49OutputMetadata& vita49)
		{
			nlohmann::json streams = nlohmann::json::array();
			for (const auto& stream : vita49.streams)
			{
				streams.push_back(vita49StreamToJson(stream));
			}

			nlohmann::json result = {{"endpoint", vita49.endpoint_host + ":" + std::to_string(vita49.endpoint_port)},
									 {"endpoint_host", vita49.endpoint_host},
									 {"endpoint_port", vita49.endpoint_port},
									 {"epoch_unix_nanoseconds", nullptr},
									 {"class_id", vita49.class_id},
									 {"adc_fullscale", vita49.adc_fullscale},
									 {"max_udp_payload", vita49.max_udp_payload},
									 {"queue_depth", vita49.queue_depth},
									 {"streams", streams}};
			if (vita49.epoch_unix_nanoseconds.has_value())
			{
				result["epoch_unix_nanoseconds"] = std::to_string(*vita49.epoch_unix_nanoseconds);
			}
			return result;
		}

		/// Converts a full output metadata snapshot to JSON.
		nlohmann::json metadataToJson(const OutputMetadata& metadata)
		{
			nlohmann::json files = nlohmann::json::array();
			std::vector<RealType> file_sampling_rates;
			for (const auto& file : metadata.files)
			{
				files.push_back(fileToJson(file));
				bool already_present = false;
				for (const auto rate : file_sampling_rates)
				{
					if (rate == file.sampling_rate)
					{
						already_present = true;
						break;
					}
				}
				if (!already_present)
				{
					file_sampling_rates.push_back(file.sampling_rate);
				}
			}

			nlohmann::json result = {{"schema_version", metadata.schema_version},
									 {"simulation_name", metadata.simulation_name},
									 {"output_directory", metadata.output_directory},
									 {"start_time", metadata.start_time},
									 {"end_time", metadata.end_time},
									 {"oversample_ratio", metadata.oversample_ratio},
									 {"files", files}};
			if (file_sampling_rates.empty())
			{
				result["sampling_rate"] = metadata.sampling_rate;
			}
			else if (file_sampling_rates.size() == 1)
			{
				result["sampling_rate"] = file_sampling_rates.front();
			}
			else
			{
				result["sampling_rate"] = nullptr;
				result["sampling_rates"] = file_sampling_rates;
			}
			if (metadata.vita49.has_value())
			{
				result["vita49"] = vita49ToJson(*metadata.vita49);
			}
			return result;
		}
	}

	OutputMetadataCollector::OutputMetadataCollector(std::string output_dir)
	{
		_metadata.output_directory = std::move(output_dir);
		_metadata.simulation_name = params::params.simulation_name;
		_metadata.start_time = params::startTime();
		_metadata.end_time = params::endTime();
		_metadata.sampling_rate = params::rate();
		_metadata.oversample_ratio = params::oversampleRatio();
	}

	void OutputMetadataCollector::addFile(OutputFileMetadata file_metadata)
	{
		std::scoped_lock const lock(_mutex);
		_metadata.files.push_back(std::move(file_metadata));
	}

	OutputMetadata OutputMetadataCollector::snapshot() const
	{
		std::scoped_lock const lock(_mutex);
		return _metadata;
	}

	std::string outputFileMetadataToJsonString(const OutputFileMetadata& metadata)
	{
		return fileToJson(metadata).dump(2);
	}

	std::string outputMetadataToJsonString(const OutputMetadata& metadata) { return metadataToJson(metadata).dump(2); }

	Vita49OutputMetadata vita49MetadataFromConfig(const Vita49OutputConfig& config)
	{
		return Vita49OutputMetadata{.endpoint_host = config.host,
									.endpoint_port = config.port,
									.epoch_unix_nanoseconds = config.epoch_unix_nanoseconds,
									.adc_fullscale = config.adc_fullscale,
									.max_udp_payload = config.max_udp_payload,
									.queue_depth = config.queue_depth};
	}
}
