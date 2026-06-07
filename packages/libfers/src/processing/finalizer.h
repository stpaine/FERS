// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file finalizer.h
 * @brief Declares the functions for the asynchronous receiver finalization pipelines.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "core/receiver_output.h"
#include "core/simulation_state.h"

namespace radar
{
	class Receiver;
	class Target;
}

namespace core
{
	struct OutputFileMetadata;
	class OutputMetadataCollector;
	class ProgressReporter;
}

namespace processing
{
	/// Builds the receiver stream descriptor used by output sinks.
	[[nodiscard]] core::ReceiverStreamDescriptor
	buildReceiverStreamDescriptor(const radar::Receiver* receiver, RealType sample_rate,
								  std::span<const core::ActiveStreamingSource> streaming_sources = {});

	/// Builds a non-owning output sample block over contiguous processed complex samples.
	[[nodiscard]] core::ReceiverSampleBlock
	buildReceiverSampleBlock(const radar::Receiver* receiver, RealType first_sample_time, RealType sample_rate,
							 std::span<const ComplexType> samples, std::uint64_t sample_start,
							 std::shared_ptr<const core::OutputFileMetadata> file_metadata = nullptr);

	/// Builds a sample block with active streaming-source context for VITA Context packets.
	[[nodiscard]] core::ReceiverSampleBlock
	buildReceiverSampleBlock(const radar::Receiver* receiver, RealType first_sample_time, RealType sample_rate,
							 std::span<const ComplexType> samples, std::uint64_t sample_start,
							 std::span<const core::ActiveStreamingSource> streaming_sources,
							 std::shared_ptr<const core::OutputFileMetadata> file_metadata);

	/// Builds HDF5 file metadata for a streaming receiver result emitted through the output sink.
	[[nodiscard]] core::OutputFileMetadata buildStreamingOutputMetadata(
		const radar::Receiver* receiver, const std::string& output_path, std::size_t total_samples,
		const std::vector<core::ActiveStreamingSource>& streaming_sources, RealType output_sample_rate);

	/**
	 * @brief The main function for a dedicated pulsed-mode receiver finalizer thread.
	 *
	 * This function runs in a loop, dequeuing and processing RenderingJobs for a
	 * specific receiver. It handles all expensive rendering, signal processing,
	 * and I/O for that receiver's data.
	 *
	 * @param receiver A pointer to the pulsed-mode receiver to process.
	 * @param targets A pointer to the world's list of targets for interference calculation.
	 * @param reporter Shared pointer to the progress reporter for status updates.
	 * @param output_dir Output directory for the simulation files.
	 */
	void runPulsedFinalizer(radar::Receiver* receiver, const std::vector<std::unique_ptr<radar::Target>>* targets,
							const std::shared_ptr<core::ProgressReporter>& reporter, const std::string& output_dir,
							const std::shared_ptr<core::OutputMetadataCollector>& metadata_collector = nullptr,
							core::ReceiverOutputSink* output_sink = nullptr);

}
