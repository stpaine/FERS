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
			return result;
		}

		/// Converts FMCW output metadata to JSON.
		nlohmann::json fmcwToJson(const FmcwMetadata& fmcw)
		{
			nlohmann::json result = {{"chirp_bandwidth", fmcw.chirp_bandwidth},
									 {"chirp_duration", fmcw.chirp_duration},
									 {"chirp_period", fmcw.chirp_period},
									 {"chirp_rate", fmcw.chirp_rate},
									 {"start_frequency_offset", fmcw.start_frequency_offset}};
			if (fmcw.chirp_count.has_value())
			{
				result["chirp_count"] = *fmcw.chirp_count;
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
									 {"total_samples", file.total_samples},
									 {"sample_start", file.sample_start},
									 {"sample_end_exclusive", file.sample_end_exclusive},
									 {"pulse_count", file.pulse_count},
									 {"min_pulse_length_samples", file.min_pulse_length_samples},
									 {"max_pulse_length_samples", file.max_pulse_length_samples},
									 {"uniform_pulse_length", file.uniform_pulse_length},
									 {"chunks", chunks},
									 {"streaming_segments", streaming_segments},
									 {"fmcw_sources", fmcw_sources}};
			if (file.fmcw.has_value())
			{
				result["fmcw"] = fmcwToJson(*file.fmcw);
			}
			return result;
		}

		/// Converts a full output metadata snapshot to JSON.
		nlohmann::json metadataToJson(const OutputMetadata& metadata)
		{
			nlohmann::json files = nlohmann::json::array();
			for (const auto& file : metadata.files)
			{
				files.push_back(fileToJson(file));
			}

			return {{"schema_version", metadata.schema_version},
					{"simulation_name", metadata.simulation_name},
					{"output_directory", metadata.output_directory},
					{"start_time", metadata.start_time},
					{"end_time", metadata.end_time},
					{"sampling_rate", metadata.sampling_rate},
					{"oversample_ratio", metadata.oversample_ratio},
					{"files", files}};
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
