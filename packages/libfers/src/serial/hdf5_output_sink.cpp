// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "serial/hdf5_output_sink.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <format>
#include <highfive/highfive.hpp>
#include <iterator>
#include <limits>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "core/output_metadata.h"
#include "core/parameters.h"
#include "processing/finalizer_pipeline.h"
#include "processing/signal_processor.h"
#include "serial/hdf5_handler.h"

namespace serial
{
	namespace
	{
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
	}

	struct Hdf5OutputSink::Impl
	{
		struct StreamState
		{
			core::ReceiverStreamDescriptor descriptor;
			core::OutputFileMetadata metadata;
			std::vector<ComplexType> streaming_buffer;
			std::optional<core::OutputFileMetadata> streaming_metadata;
			std::unique_ptr<HighFive::File> pulsed_file;
			unsigned next_chunk_index = 0;
			bool opened = false;
			bool closed = false;
		};

		Impl(std::string out_dir, std::shared_ptr<core::OutputMetadataCollector> collector) :
			output_dir(std::move(out_dir)), metadata_collector(std::move(collector))
		{
		}

		void initializeRun(const core::OutputConfig& config, const std::string& /*simulation_name*/)
		{
			if (config.mode != core::OutputMode::Hdf5)
			{
				throw std::invalid_argument("Hdf5OutputSink requires HDF5 output mode");
			}
			std::scoped_lock const lock(mutex);
			std::filesystem::create_directories(output_dir);
			finalized = false;
		}

		std::uint32_t registerStream(const core::ReceiverStreamDescriptor& stream)
		{
			std::scoped_lock const lock(mutex);
			const auto key = streamKey(stream);
			if (const auto found = stream_ids.find(key); found != stream_ids.end())
			{
				return found->second;
			}

			const auto stream_id = next_stream_id++;
			StreamState state;
			state.descriptor = stream;
			state.metadata = core::OutputFileMetadata{.receiver_id = stream.receiver_id,
													  .receiver_name = stream.receiver_name,
													  .mode = stream.mode,
													  .path = outputPath(stream.receiver_name),
													  .sampling_rate = stream.sample_rate};
			stream_ids.emplace(key, stream_id);
			streams.emplace(stream_id, std::move(state));
			return stream_id;
		}

		void openStream(const std::uint32_t stream_id, const RealType /*first_sample_time*/)
		{
			std::scoped_lock const lock(mutex);
			auto& state = stateFor(stream_id);
			if (state.opened)
			{
				return;
			}
			std::filesystem::create_directories(output_dir);
			if (isPulsed(state))
			{
				std::scoped_lock const hdf5_lock(hdf5_global_mutex);
				state.pulsed_file = std::make_unique<HighFive::File>(state.metadata.path, HighFive::File::Truncate);
			}
			else
			{
				const auto expected_samples = expectedStreamingSamples(state.descriptor.sample_rate);
				if (expected_samples > 0u)
				{
					state.streaming_buffer.resize(expected_samples);
				}
			}
			state.opened = true;
			state.closed = false;
		}

		void submitBlock(const core::ReceiverSampleBlock& block)
		{
			std::scoped_lock const lock(mutex);
			const auto stream_id = registerStream(block.stream);
			auto& state = stateFor(stream_id);
			if (!state.opened)
			{
				openStream(stream_id, block.first_sample_time);
			}
			if (block.file_metadata)
			{
				state.streaming_metadata = *block.file_metadata;
			}

			if (isPulsed(state))
			{
				writePulsedBlock(state, block);
				return;
			}
			appendStreamingBlock(state, block);
		}

		void closeStream(const std::uint32_t stream_id)
		{
			std::scoped_lock const lock(mutex);
			auto& state = stateFor(stream_id);
			if (state.closed)
			{
				return;
			}

			if (isPulsed(state))
			{
				closePulsedStream(state);
			}
			else
			{
				closeStreamingStream(state);
			}
			state.closed = true;
		}

		core::OutputStats finalize()
		{
			std::scoped_lock const lock(mutex);
			if (!finalized)
			{
				std::vector<std::uint32_t> ids;
				ids.reserve(streams.size());
				for (const auto& [stream_id, state] : streams)
				{
					if (state.opened && !state.closed)
					{
						ids.push_back(stream_id);
					}
				}
				for (const auto stream_id : ids)
				{
					closeStream(stream_id);
				}
				finalized = true;
			}
			return core::OutputStats{.mode = core::OutputMode::Hdf5, .streams = {}};
		}

		static bool isPulsed(const StreamState& state) { return state.descriptor.mode == "pulsed"; }

		static std::size_t expectedStreamingSamples(const RealType sample_rate)
		{
			if (sample_rate <= 0.0)
			{
				return 0u;
			}
			const RealType duration = std::max<RealType>(0.0, params::endTime() - params::startTime());
			return static_cast<std::size_t>(std::ceil(duration * sample_rate));
		}

		static std::string streamKey(const core::ReceiverStreamDescriptor& stream)
		{
			return std::format("{}:{}:{}", stream.receiver_id, stream.receiver_name, stream.mode);
		}

		[[nodiscard]] std::string outputPath(const std::string& receiver_name) const
		{
			const std::filesystem::path out_path(output_dir);
			return (out_path / std::format("{}_results.h5", receiver_name)).string();
		}

		StreamState& stateFor(const std::uint32_t stream_id)
		{
			const auto found = streams.find(stream_id);
			if (found == streams.end())
			{
				throw std::out_of_range("Unknown HDF5 output stream ID");
			}
			return found->second;
		}

		void writePulsedBlock(StreamState& state, const core::ReceiverSampleBlock& block)
		{
			if (!state.pulsed_file)
			{
				throw std::logic_error("HDF5 pulsed stream is not open");
			}

			std::vector<ComplexType> chunk(block.samples.begin(), block.samples.end());
			const RealType fullscale = processing::quantizeAndScaleWindow(chunk);
			const auto chunk_index = state.next_chunk_index++;
			const auto sample_start = state.metadata.total_samples;
			core::PulseChunkMetadata chunk_metadata{.chunk_index = chunk_index,
													.i_dataset = std::format("chunk_{:06}_I", chunk_index),
													.q_dataset = std::format("chunk_{:06}_Q", chunk_index),
													.start_time = block.first_sample_time,
													.sample_count = static_cast<std::uint64_t>(chunk.size()),
													.sample_start = sample_start,
													.sample_end_exclusive =
														sample_start + static_cast<std::uint64_t>(chunk.size())};
			addChunkToFile(*state.pulsed_file, chunk, block.first_sample_time, fullscale, chunk_index, &chunk_metadata);
			state.metadata.chunks.push_back(std::move(chunk_metadata));
			state.metadata.total_samples = state.metadata.chunks.back().sample_end_exclusive;
		}

		void appendStreamingBlock(StreamState& state, const core::ReceiverSampleBlock& block)
		{
			const auto sample_start = static_cast<std::size_t>(block.sample_start);
			const auto sample_end = sample_start + block.samples.size();
			if (state.streaming_buffer.size() < sample_end)
			{
				state.streaming_buffer.resize(sample_end);
			}
			std::ranges::copy(block.samples,
							  std::next(state.streaming_buffer.begin(),
										static_cast<std::vector<ComplexType>::difference_type>(sample_start)));
		}

		void closePulsedStream(StreamState& state)
		{
			if (!state.pulsed_file)
			{
				return;
			}
			finalizePulsedMetadata(state.metadata);
			{
				std::scoped_lock const hdf5_lock(hdf5_global_mutex);
				writeOutputFileMetadataAttributes(*state.pulsed_file, state.metadata);
				state.pulsed_file.reset();
			}
			if (metadata_collector)
			{
				metadata_collector->addFile(state.metadata);
			}
		}

		void closeStreamingStream(StreamState& state)
		{
			if (state.streaming_buffer.empty())
			{
				return;
			}

			auto metadata = state.streaming_metadata.value_or(
				core::OutputFileMetadata{.receiver_id = state.descriptor.receiver_id,
										 .receiver_name = state.descriptor.receiver_name,
										 .mode = state.descriptor.mode,
										 .path = state.metadata.path,
										 .sampling_rate = state.descriptor.sample_rate});
			metadata.path = state.metadata.path;
			metadata.sampling_rate = state.descriptor.sample_rate;
			metadata.total_samples = static_cast<std::uint64_t>(state.streaming_buffer.size());
			metadata.sample_start = 0;
			metadata.sample_end_exclusive = static_cast<std::uint64_t>(state.streaming_buffer.size());

			const RealType fullscale = processing::quantizeAndScaleWindow(state.streaming_buffer);
			processing::pipeline::exportStreamingToHdf5(state.metadata.path, state.streaming_buffer, fullscale,
														state.descriptor.reference_frequency, &metadata,
														state.descriptor.sample_rate);
			if (metadata_collector)
			{
				metadata_collector->addFile(std::move(metadata));
			}
		}

		std::string output_dir;
		std::shared_ptr<core::OutputMetadataCollector> metadata_collector;
		std::recursive_mutex mutex;
		std::unordered_map<std::uint32_t, StreamState> streams;
		std::unordered_map<std::string, std::uint32_t> stream_ids;
		std::uint32_t next_stream_id = 1;
		bool finalized = false;
	};

	Hdf5OutputSink::Hdf5OutputSink(std::string output_dir,
								   std::shared_ptr<core::OutputMetadataCollector> metadata_collector) :
		_impl(std::make_unique<Impl>(std::move(output_dir), std::move(metadata_collector)))
	{
	}

	Hdf5OutputSink::~Hdf5OutputSink()
	{
		if (_impl)
		{
			try
			{
				(void)_impl->finalize();
			}
			catch (...)
			{
			}
		}
	}

	void Hdf5OutputSink::initializeRun(const core::OutputConfig& config, std::string simulation_name)
	{
		_impl->initializeRun(config, simulation_name);
	}

	std::uint32_t Hdf5OutputSink::registerStream(const core::ReceiverStreamDescriptor& stream)
	{
		return _impl->registerStream(stream);
	}

	void Hdf5OutputSink::openStream(const std::uint32_t stream_id, const RealType first_sample_time)
	{
		_impl->openStream(stream_id, first_sample_time);
	}

	void Hdf5OutputSink::submitBlock(const core::ReceiverSampleBlock& block) { _impl->submitBlock(block); }

	void Hdf5OutputSink::emitContextHeartbeat(const RealType /*simulation_time*/) {}

	void Hdf5OutputSink::closeStream(const std::uint32_t stream_id) { _impl->closeStream(stream_id); }

	core::OutputStats Hdf5OutputSink::finalize() { return _impl->finalize(); }

	std::unique_ptr<core::ReceiverOutputSink>
	makeHdf5OutputSink(std::string output_dir, std::shared_ptr<core::OutputMetadataCollector> metadata_collector)
	{
		return std::make_unique<Hdf5OutputSink>(std::move(output_dir), std::move(metadata_collector));
	}
}
