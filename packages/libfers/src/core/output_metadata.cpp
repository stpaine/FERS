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
									 {"fmcw_dechirp_mode", file.fmcw_dechirp_mode},
									 {"fmcw_dechirp_reference_source", file.fmcw_dechirp_reference_source}};
			if (file.fmcw.has_value())
			{
				result["fmcw"] = fmcwToJson(*file.fmcw);
			}
			if (file.fmcw_dechirp_reference_transmitter_id.has_value())
			{
				result["fmcw_dechirp_reference_transmitter_id"] = *file.fmcw_dechirp_reference_transmitter_id;
			}
			if (file.fmcw_dechirp_reference_transmitter_name.has_value())
			{
				result["fmcw_dechirp_reference_transmitter_name"] = *file.fmcw_dechirp_reference_transmitter_name;
			}
			if (file.fmcw_dechirp_reference_waveform_id.has_value())
			{
				result["fmcw_dechirp_reference_waveform_id"] = *file.fmcw_dechirp_reference_waveform_id;
			}
			if (file.fmcw_dechirp_reference_waveform_name.has_value())
			{
				result["fmcw_dechirp_reference_waveform_name"] = *file.fmcw_dechirp_reference_waveform_name;
			}
			if (file.fmcw_dechirp_reference_waveform.has_value())
			{
				result["fmcw_dechirp_reference_waveform"] = fmcwToJson(*file.fmcw_dechirp_reference_waveform);
			}
			result["fmcw_if_decimation_enabled"] = file.fmcw_if_decimation_enabled;
			result["fmcw_if_legacy_full_rate"] = file.fmcw_if_legacy_full_rate;
			if (file.fmcw_if_requested_sample_rate.has_value())
			{
				result["fmcw_if_requested_sample_rate"] = *file.fmcw_if_requested_sample_rate;
			}
			if (file.fmcw_if_sample_rate.has_value())
			{
				result["fmcw_if_sample_rate"] = *file.fmcw_if_sample_rate;
			}
			if (file.fmcw_if_input_sample_rate.has_value())
			{
				result["fmcw_if_input_sample_rate"] = *file.fmcw_if_input_sample_rate;
			}
			if (file.fmcw_if_resample_numerator.has_value())
			{
				result["fmcw_if_resample_numerator"] = *file.fmcw_if_resample_numerator;
			}
			if (file.fmcw_if_resample_denominator.has_value())
			{
				result["fmcw_if_resample_denominator"] = *file.fmcw_if_resample_denominator;
			}
			if (file.fmcw_if_decimation_factor.has_value())
			{
				result["fmcw_if_decimation_factor"] = *file.fmcw_if_decimation_factor;
			}
			if (file.fmcw_if_filter_bandwidth.has_value())
			{
				result["fmcw_if_filter_bandwidth"] = *file.fmcw_if_filter_bandwidth;
			}
			if (file.fmcw_if_filter_transition_width.has_value())
			{
				result["fmcw_if_filter_transition_width"] = *file.fmcw_if_filter_transition_width;
			}
			if (file.fmcw_if_filter_stopband.has_value())
			{
				result["fmcw_if_filter_stopband"] = *file.fmcw_if_filter_stopband;
			}
			if (file.fmcw_if_filter_group_delay_seconds.has_value())
			{
				result["fmcw_if_filter_group_delay_seconds"] = *file.fmcw_if_filter_group_delay_seconds;
			}
			if (file.fmcw_if_compensated_integer_delay_samples.has_value())
			{
				result["fmcw_if_compensated_integer_delay_samples"] = *file.fmcw_if_compensated_integer_delay_samples;
			}
			if (file.fmcw_if_compensated_fractional_delay_samples.has_value())
			{
				result["fmcw_if_compensated_fractional_delay_samples"] =
					*file.fmcw_if_compensated_fractional_delay_samples;
			}
			if (file.fmcw_if_warmup_discard_samples.has_value())
			{
				result["fmcw_if_warmup_discard_samples"] = *file.fmcw_if_warmup_discard_samples;
			}
			if (file.fmcw_if_phase_refinement.has_value())
			{
				result["fmcw_if_phase_refinement"] = *file.fmcw_if_phase_refinement;
			}
			if (file.fmcw_if_timing_error_seconds.has_value())
			{
				result["fmcw_if_timing_error_seconds"] = *file.fmcw_if_timing_error_seconds;
			}
			if (file.fmcw_if_phase_error_radians.has_value())
			{
				result["fmcw_if_phase_error_radians"] = *file.fmcw_if_phase_error_radians;
			}
			if (file.fmcw_if_noise_variance.has_value())
			{
				result["fmcw_if_noise_variance"] = *file.fmcw_if_noise_variance;
			}
			result["fmcw_if_group_delay_compensated"] = file.fmcw_if_group_delay_compensated;
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
		std::scoped_lock lock(_mutex);
		_metadata.files.push_back(std::move(file_metadata));
	}

	OutputMetadata OutputMetadataCollector::snapshot() const
	{
		std::scoped_lock lock(_mutex);
		return _metadata;
	}

	std::string outputFileMetadataToJsonString(const OutputFileMetadata& metadata)
	{
		return fileToJson(metadata).dump(2);
	}

	std::string outputMetadataToJsonString(const OutputMetadata& metadata) { return metadataToJson(metadata).dump(2); }
}
