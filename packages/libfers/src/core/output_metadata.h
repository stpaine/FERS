// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#pragma once

#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

#include "core/config.h"
#include "core/sim_id.h"

namespace core
{
	struct PulseChunkMetadata
	{
		unsigned chunk_index = 0;
		std::string i_dataset;
		std::string q_dataset;
		RealType start_time = 0.0;
		std::uint64_t sample_count = 0;
		std::uint64_t sample_start = 0;
		std::uint64_t sample_end_exclusive = 0;
	};

	struct CwSegmentMetadata
	{
		RealType start_time = 0.0;
		RealType end_time = 0.0;
		std::uint64_t sample_count = 0;
		std::uint64_t sample_start = 0;
		std::uint64_t sample_end_exclusive = 0;
	};

	struct OutputFileMetadata
	{
		SimId receiver_id = 0;
		std::string receiver_name;
		std::string mode;
		std::string path;
		std::uint64_t total_samples = 0;
		std::uint64_t sample_start = 0;
		std::uint64_t sample_end_exclusive = 0;
		std::uint64_t pulse_count = 0;
		std::uint64_t min_pulse_length_samples = 0;
		std::uint64_t max_pulse_length_samples = 0;
		bool uniform_pulse_length = true;
		std::vector<PulseChunkMetadata> chunks = {};
		std::vector<CwSegmentMetadata> cw_segments = {};
	};

	struct OutputMetadata
	{
		unsigned schema_version = 1;
		std::string simulation_name;
		std::string output_directory;
		RealType start_time = 0.0;
		RealType end_time = 0.0;
		RealType sampling_rate = 0.0;
		unsigned oversample_ratio = 1;
		std::vector<OutputFileMetadata> files;
	};

	class OutputMetadataCollector
	{
	public:
		explicit OutputMetadataCollector(std::string output_dir);

		void addFile(OutputFileMetadata file_metadata);

		[[nodiscard]] OutputMetadata snapshot() const;

	private:
		mutable std::mutex _mutex;
		OutputMetadata _metadata;
	};

	[[nodiscard]] std::string outputFileMetadataToJsonString(const OutputFileMetadata& metadata);
	[[nodiscard]] std::string outputMetadataToJsonString(const OutputMetadata& metadata);
}
