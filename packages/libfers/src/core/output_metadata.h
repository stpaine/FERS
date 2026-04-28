// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#pragma once

#include <cstdint>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

#include "core/config.h"
#include "core/sim_id.h"

namespace core
{
	/// Metadata for one pulsed output chunk written to HDF5.
	struct PulseChunkMetadata
	{
		unsigned chunk_index = 0; ///< Zero-based chunk index in the receiver output sequence.
		std::string i_dataset; ///< HDF5 dataset path containing the in-phase samples.
		std::string q_dataset; ///< HDF5 dataset path containing the quadrature samples.
		RealType start_time = 0.0; ///< Simulation time of the first sample in seconds.
		std::uint64_t sample_count = 0; ///< Number of samples in the chunk.
		std::uint64_t sample_start = 0; ///< Inclusive global sample index for the chunk start.
		std::uint64_t sample_end_exclusive = 0; ///< Exclusive global sample index for the chunk end.
	};

	/// Metadata for one contiguous streaming output segment.
	struct StreamingSegmentMetadata
	{
		RealType start_time = 0.0; ///< Segment start time in seconds.
		RealType end_time = 0.0; ///< Segment end time in seconds.
		std::uint64_t sample_count = 0; ///< Number of samples emitted for the segment.
		std::uint64_t sample_start = 0; ///< Inclusive global sample index for the segment start.
		std::uint64_t sample_end_exclusive = 0; ///< Exclusive global sample index for the segment end.
		std::optional<RealType> first_chirp_start_time = std::nullopt; ///< First emitted FMCW chirp start time.
		std::optional<std::uint64_t> emitted_chirp_count = std::nullopt; ///< Number of FMCW chirps emitted.
	};

	/// FMCW waveform metadata captured for a streaming output file.
	struct FmcwMetadata
	{
		RealType chirp_bandwidth = 0.0; ///< Chirp bandwidth in hertz.
		RealType chirp_duration = 0.0; ///< Active chirp duration in seconds.
		RealType chirp_period = 0.0; ///< Chirp repetition period in seconds.
		RealType chirp_rate = 0.0; ///< Frequency sweep rate in hertz per second.
		RealType chirp_rate_signed = 0.0; ///< Signed frequency sweep rate in hertz per second.
		std::string chirp_direction = "up"; ///< Frequency sweep direction token.
		RealType start_frequency_offset = 0.0; ///< Start frequency offset relative to carrier in hertz.
		std::optional<std::uint64_t> chirp_count = std::nullopt; ///< Optional finite chirp count.
	};

	/// Metadata for one active FMCW transmitter schedule segment.
	struct FmcwSourceSegmentMetadata
	{
		RealType start_time = 0.0; ///< Transmitter segment start time in seconds.
		RealType end_time = 0.0; ///< Transmitter segment end time in seconds.
		std::optional<RealType> first_chirp_start_time = std::nullopt; ///< First emitted chirp start in the segment.
		std::optional<std::uint64_t> emitted_chirp_count = std::nullopt; ///< Number of chirps emitted in the segment.
	};

	/// Metadata for one FMCW illuminator represented in a streaming output file.
	struct FmcwSourceMetadata
	{
		SimId transmitter_id = 0; ///< FMCW transmitter SimId.
		std::string transmitter_name; ///< FMCW transmitter display name.
		SimId waveform_id = 0; ///< FMCW waveform SimId.
		std::string waveform_name; ///< FMCW waveform display name.
		RealType carrier_frequency = 0.0; ///< Waveform carrier frequency in hertz.
		FmcwMetadata waveform; ///< FMCW chirp parameters.
		std::vector<FmcwSourceSegmentMetadata> segments = {}; ///< Active transmitter segments.
	};

	/// Metadata for one receiver output file.
	struct OutputFileMetadata
	{
		SimId receiver_id = 0; ///< Receiver SimId that owns the output file.
		std::string receiver_name; ///< Receiver display name.
		std::string mode; ///< Output mode label, such as pulsed or streaming.
		std::string path; ///< Filesystem path to the generated output file.
		std::uint64_t total_samples = 0; ///< Total sample count written to the file.
		std::uint64_t sample_start = 0; ///< Inclusive global sample index for the file start.
		std::uint64_t sample_end_exclusive = 0; ///< Exclusive global sample index for the file end.
		std::uint64_t pulse_count = 0; ///< Number of pulses represented in the file.
		std::uint64_t min_pulse_length_samples = 0; ///< Minimum pulse length in samples.
		std::uint64_t max_pulse_length_samples = 0; ///< Maximum pulse length in samples.
		bool uniform_pulse_length = true; ///< True when every pulse has the same sample length.
		std::vector<PulseChunkMetadata> chunks = {}; ///< Pulsed output chunks written to the file.
		std::vector<StreamingSegmentMetadata> streaming_segments = {}; ///< Streaming segments written to the file.
		std::optional<FmcwMetadata> fmcw = std::nullopt; ///< Optional FMCW metadata for streaming outputs.
		std::vector<FmcwSourceMetadata> fmcw_sources = {}; ///< FMCW illuminators represented in the output.
	};

	/// Metadata summary for the full simulation output set.
	struct OutputMetadata
	{
		unsigned schema_version = 1; ///< Metadata schema version.
		std::string simulation_name; ///< Simulation name from the loaded scenario.
		std::string output_directory; ///< Directory containing generated output files.
		RealType start_time = 0.0; ///< Simulation start time in seconds.
		RealType end_time = 0.0; ///< Simulation end time in seconds.
		RealType sampling_rate = 0.0; ///< Output sampling rate in hertz.
		unsigned oversample_ratio = 1; ///< Oversampling ratio used during rendering.
		std::vector<OutputFileMetadata> files; ///< Metadata for each generated output file.
	};

	/// Thread-safe collector for simulation output metadata.
	class OutputMetadataCollector
	{
	public:
		/// Constructs a metadata collector for the specified output directory.
		explicit OutputMetadataCollector(std::string output_dir);

		/// Adds metadata for one generated output file.
		void addFile(OutputFileMetadata file_metadata);

		/// Returns a consistent snapshot of the collected metadata.
		[[nodiscard]] OutputMetadata snapshot() const;

	private:
		mutable std::mutex _mutex; ///< Mutex guarding the metadata aggregate.
		OutputMetadata _metadata; ///< Mutable aggregate metadata collected during simulation.
	};

	/// Serializes one output-file metadata entry to JSON.
	[[nodiscard]] std::string outputFileMetadataToJsonString(const OutputFileMetadata& metadata);

	/// Serializes a full simulation output metadata snapshot to JSON.
	[[nodiscard]] std::string outputMetadataToJsonString(const OutputMetadata& metadata);
}
