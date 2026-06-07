// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file sim_threading.h
 * @brief Header file for the main simulation runner.
 *
 * This file contains the declarations for the high-level function and engine that
 * orchestrates and manages the event-driven radar simulation.
 */

#pragma once

#include <chrono>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <span>
#include <string>
#include <thread>
#include <vector>

#include "core/config.h"
#include "core/output_config.h"
#include "core/output_metadata.h"
#include "core/parameters.h"
#include "core/receiver_output.h"
#include "core/sim_events.h"
#include "core/simulation_state.h"
#include "signal/dsp_filters.h"
#include "simulation/channel_model.h"

namespace pool
{
	class ThreadPool;
}

namespace radar
{
	class Receiver;
	class Target;
	class Transmitter;
}

namespace serial
{
	class Response;
}

namespace core
{
	class World;
	class ReceiverOutputSink;

	/**
	 * @class ProgressReporter
	 * @brief A thread-safe wrapper for the simulation progress callback.
	 *
	 * Allows multiple worker threads to report progress concurrently without race conditions.
	 */
	class ProgressReporter
	{
	public:
		/**
		 * @typedef Callback
		 * @brief Defines the signature for the progress reporting callback function.
		 */
		using Callback = std::function<void(const std::string&, int, int)>;

		/**
		 * @brief Constructs a ProgressReporter with the given callback.
		 * @param cb The callback function to wrap.
		 */
		explicit ProgressReporter(Callback cb) : _callback(std::move(cb)) {}

		/**
		 * @brief Safely reports progress to the underlying callback.
		 * @param msg The status message to report.
		 * @param current The current progress value.
		 * @param total The total progress value.
		 */
		void report(const std::string& msg, int current, int total)
		{
			if (_callback)
			{
				std::scoped_lock const lock(_mutex);
				_callback(msg, current, total);
			}
		}

	private:
		std::mutex _mutex; ///< Mutex to ensure thread-safe access to the callback.
		Callback _callback; ///< The underlying callback function.
	};

	/**
	 * @class SimulationEngine
	 * @brief Encapsulates the state and logic of the event-driven simulation loop.
	 *
	 * Breaking the simulation loop into this class allows for easily testable,
	 * focused functions with low cyclomatic complexity.
	 */
	class SimulationEngine
	{
	public:
		/**
		 * @brief Constructs the simulation engine.
		 * @param world Pointer to the simulation world containing all entities.
		 * @param pool Reference to the thread pool for asynchronous tasks.
		 * @param reporter Shared pointer to the thread-safe progress reporter.
		 * @param output_dir Output directory for the simulation files.
		 */
		SimulationEngine(World* world, pool::ThreadPool& pool, std::shared_ptr<ProgressReporter> reporter,
						 std::string output_dir, std::shared_ptr<OutputMetadataCollector> metadata_collector = nullptr,
						 ReceiverOutputSink* output_sink = nullptr, std::function<bool()> cancel_callback = nullptr,
						 bool eager_context_stream_open = false);

		/**
		 * @brief Starts and runs the main simulation loop until completion.
		 */
		void run();

		/// Returns true after cooperative cancellation has been requested.
		[[nodiscard]] bool cancelled() const noexcept { return _cancelled; }

		/**
		 * @brief Advances the time-stepped inner loop for active streaming systems.
		 * @param t_event The timestamp of the next discrete event to process up to.
		 */
		void processStreamingPhysics(RealType t_event);

		/**
		 * @brief Dispatches a discrete simulation event to its specific handler.
		 * @param event The event to process.
		 */
		void processEvent(const Event& event);

		/**
		 * @brief Handles the start of a pulsed transmission.
		 * @param tx Pointer to the transmitting radar object.
		 * @param t_event The timestamp of the transmission event.
		 */
		void handleTxPulsedStart(radar::Transmitter* tx, RealType t_event);

		/**
		 * @brief Handles the opening of a pulsed receiver's listening window.
		 * @param rx Pointer to the receiving radar object.
		 * @param t_event The timestamp of the window opening event.
		 */
		void handleRxPulsedWindowStart(radar::Receiver* rx, RealType t_event);

		/**
		 * @brief Handles the closing of a pulsed receiver's listening window, triggering finalization.
		 * @param rx Pointer to the receiving radar object.
		 * @param t_event The timestamp of the window closing event.
		 */
		void handleRxPulsedWindowEnd(radar::Receiver* rx, RealType t_event);

		/**
		 * @brief Handles a streaming transmitter turning on.
		 * @param tx Pointer to the transmitting radar object.
		 */
		void handleTxStreamingStart(const ActiveStreamingSource& source);

		/**
		 * @brief Handles a streaming transmitter turning off.
		 * @param tx Pointer to the transmitting radar object.
		 */
		void handleTxStreamingEnd(radar::Transmitter* tx);

		/**
		 * @brief Handles a streaming receiver starting to record.
		 * @param rx Pointer to the receiving radar object.
		 */
		void handleRxStreamingStart(radar::Receiver* rx);

		/**
		 * @brief Handles a streaming receiver stopping recording.
		 * @param rx Pointer to the receiving radar object.
		 */
		void handleRxStreamingEnd(radar::Receiver* rx);

	private:
		/// Calculates one streaming I/Q sample for the receiver at the specified time step.
		[[nodiscard]] ComplexType calculateStreamingSample(radar::Receiver* rx, RealType t_step,
														   const std::vector<ActiveStreamingSource>& streaming_sources,
														   ReceiverTrackerCache& tracker_cache) const;

		/// Adds tracker storage for a newly active streaming source.
		void appendStreamingTrackerSource();

		/// Removes tracker storage for a streaming source that has ended.
		void eraseStreamingTrackerSource(std::size_t source_index);

		/// Removes ended streaming sources once no future receiver sample can observe their in-flight energy.
		void cleanupInactiveStreamingSources(RealType from_time);

		/// Returns the next time at which ended streaming sources can be cleaned up.
		[[nodiscard]] std::optional<RealType> nextStreamingCleanupDeadline(RealType from_time);

		/// Returns the next streaming chunk boundary before the target event time.
		[[nodiscard]] RealType streamingChunkEnd(RealType from_time, RealType event_time);

		/// Returns true when cooperative cancellation should stop the current streaming chunk.
		[[nodiscard]] bool shouldStopStreamingChunk(std::size_t sample_index, std::size_t chunk_start_index);

		/// Processes one simulation sample for all active streaming receivers.
		void processStreamingSample(std::size_t sample_index, std::size_t first_index, std::size_t final_index,
									std::size_t progress_report_stride, RealType dt_sim);

		/// Adds one sample to every active streaming receiver.
		void appendActiveReceiverStreamingSamples(std::size_t sample_index, RealType t_step);

		/// Adds one sample to a single active streaming receiver.
		void appendReceiverStreamingSample(std::size_t receiver_index, std::size_t sample_index, RealType t_step);

		/// Builds and attaches IF resampling sinks for configured FMCW receivers.
		void initializeFmcwIfResamplers();

		/// Builds and attaches one receiver IF resampler when configured.
		void initializeFmcwIfResampler(std::size_t receiver_index);

		/// Extends resolved dechirp source spans to cover IF resampler over-render.
		void extendDechirpSourcesForIfOverrender();

		/// Flushes all pending high-rate samples buffered for IF resamplers.
		void flushFmcwIfBlocks();

		/// Flushes one receiver's pending IF-resampling high-rate block.
		void flushFmcwIfBlock(std::size_t receiver_index);

		/// Adds one high-rate sample to a receiver's IF-resampling block buffer.
		void appendFmcwIfSample(std::size_t receiver_index, RealType t_step, ComplexType sample);

		/// Adds one raw streaming sample to the live output block buffer.
		void appendStreamingOutputSample(std::size_t receiver_index, std::size_t sample_index, RealType t_step,
										 ComplexType sample);

		/// Flushes all pending live streaming output blocks.
		void flushStreamingOutputBlocks();

		/// Flushes one receiver's pending live streaming output block.
		void flushStreamingOutputBlock(std::size_t receiver_index, bool finish_downsampler = false);

		/// Emits an already processed streaming block to the selected output sink.
		void emitStreamingOutputBlock(std::size_t receiver_index, RealType first_sample_time, RealType sample_rate,
									  std::span<const ComplexType> samples, std::uint64_t sample_start);

		/// Returns the sink-visible sample rate for one live streaming receiver.
		[[nodiscard]] RealType streamingOutputSampleRate(std::size_t receiver_index) const;

		/// Registers and opens one live streaming receiver before delayed sample flushes.
		void ensureStreamingOutputStreamOpen(std::size_t receiver_index, RealType first_sample_time,
											 RealType sample_rate);

		/// Returns an initialized stateful downsampler for one streaming receiver segment.
		[[nodiscard]] fers_signal::DownsamplingSink&
		streamingDownsampler(std::size_t receiver_index, std::uint64_t input_start_index, RealType segment_start_time);

		/// Emits sink heartbeats on the continuous simulation clock up to the given time.
		void emitContextHeartbeatsThrough(RealType simulation_time);

		/// Returns the latest conservative receive time at which an ended source must still be retained.
		[[nodiscard]] std::optional<RealType> streamingSourceCleanupDeadline(const ActiveStreamingSource& source,
																			 RealType from_time) const;

		/// Returns the latest conservative receive time at which one receiver may still observe a source.
		[[nodiscard]] std::optional<RealType> receiverCleanupDeadline(const ActiveStreamingSource& source,
																	  const radar::Receiver* rx,
																	  RealType from_time) const;

		/// Creates the CW phase-noise lookup if any active timing source needs it.
		void ensureCwPhaseNoiseLookup();

		/// Returns the dechirp reference mixer for a receiver at one sample time.
		[[nodiscard]] std::optional<ComplexType> calculateDechirpMixer(radar::Receiver* rx, RealType t_step,
																	   ReceiverTrackerCache& tracker_cache) const;

		/// Adds pulsed interference into a completed high-rate IF block before resampling.
		void applyPulsedInterferenceToFmcwIfBlock(std::size_t receiver_index, std::span<ComplexType> block,
												  RealType block_start_time);

		/// Adds pulsed interference into a live streaming block before final noise/scaling.
		void applyPulsedInterferenceToStreamingBlock(std::size_t receiver_index, std::span<ComplexType> block,
													 RealType block_start_time, RealType sample_rate, bool dechirp_mix);

		/// Adds one rendered pulsed-interference slice into a streaming block.
		void addPulsedInterferenceSamples(std::span<ComplexType> block, std::span<const ComplexType> rendered_pulse,
										  long long dest_begin, long long dest_end, std::size_t crop_offset,
										  RealType block_start_time, RealType sample_rate, bool dechirp_mix,
										  radar::Receiver* receiver, ReceiverTrackerCache& tracker_cache) const;

		/// Emits summary logs for streaming receiver configuration.
		void logStreamingSummaries() const;

		/**
		 * @brief Starts dedicated finalizer threads for all pulsed receivers.
		 */
		void initializeFinalizers();

		/**
		 * @brief Routes a calculated radar response to the appropriate receiver inbox or log.
		 * @param rx Pointer to the receiving radar object.
		 * @param response The calculated response to route.
		 */
		void routeResponse(radar::Receiver* rx, std::unique_ptr<serial::Response> response) const;

		/**
		 * @brief Throttles and emits progress updates to the reporter.
		 */
		void updateProgress();

		/**
		 * @brief Throttles and emits progress updates at an explicit simulation time.
		 */
		void reportSimulationProgress(RealType t_current);

		/// Polls the host cancellation callback and latches the result.
		[[nodiscard]] bool isCancellationRequested();

		/// Collects streaming sources active anywhere within the requested time window.
		[[nodiscard]] std::vector<ActiveStreamingSource> collectStreamingSourcesForWindow(RealType start_time,
																						  RealType end_time) const;

		/**
		 * @brief Initiates the shutdown phase, waiting for all asynchronous tasks to complete.
		 */
		void shutdown();

		World* _world; ///< Pointer to the simulation world state.
		pool::ThreadPool& _pool; ///< Reference to the global thread pool.
		std::shared_ptr<ProgressReporter> _reporter; ///< Shared progress reporter instance.
		std::vector<std::jthread> _finalizer_threads; ///< Collection of dedicated pulsed finalizer threads.
		std::shared_ptr<OutputMetadataCollector> _metadata_collector; ///< Collector for generated output metadata.
		ReceiverOutputSink* _output_sink = nullptr; ///< Selected receiver output sink.
		std::function<bool()> _cancel_callback; ///< Optional cooperative cancellation callback.
		bool _eager_context_stream_open = false; ///< True when context heartbeat requires pre-data stream open.
		bool _cancelled = false; ///< Latched cancellation state.

		std::chrono::steady_clock::time_point _last_report_time; ///< Timestamp of the last progress report.
		int _last_reported_percent = -1; ///< The last reported percentage to prevent redundant updates.
		RealType _next_context_heartbeat_time = 0.0; ///< Next one-second sink heartbeat on simulation clock.

		std::string _output_dir; ///< Output directory for the simulation files.
		std::unique_ptr<simulation::CwPhaseNoiseLookup> _cw_phase_noise_lookup; ///< Cached CW phase-noise lookup.
		std::vector<ReceiverTrackerCache> _streaming_tracker_caches; ///< Per-receiver streaming tracker caches.
		std::vector<ReceiverTrackerCache> _if_pulse_tracker_caches; ///< Per-receiver dechirp trackers for pulse blocks.
		std::vector<std::vector<ComplexType>> _fmcw_if_block_buffers; ///< Pending high-rate IF blocks.
		std::vector<RealType> _fmcw_if_block_start_times; ///< Start time for each pending IF block.
		std::vector<std::vector<ComplexType>> _streaming_output_block_buffers; ///< Pending live CW/FMCW output blocks.
		std::vector<std::vector<ComplexType>> _streaming_output_processed_buffers; ///< Reused post-noise output blocks.
		std::vector<RealType> _streaming_output_block_start_times; ///< Start time for pending live output blocks.
		std::vector<std::uint64_t> _streaming_output_block_start_indices; ///< High-rate sample index for each block.
		std::vector<std::unique_ptr<fers_signal::DownsamplingSink>>
			_streaming_downsamplers; ///< Stateful downsamplers for non-dechirped streaming output.
		std::vector<std::uint64_t> _streaming_downsample_base_indices; ///< Output index at segment start.
		std::vector<RealType> _streaming_downsample_segment_start_times; ///< Segment start time at output sample 0.
		std::vector<std::uint64_t> _streaming_output_sample_cursors; ///< Output sample cursor per receiver.
		std::vector<std::uint32_t> _streaming_output_stream_ids; ///< Registered sink stream IDs for live receivers.
		std::vector<bool> _streaming_output_stream_open; ///< Live receiver stream lifecycle state.
		std::vector<std::shared_ptr<const OutputFileMetadata>>
			_streaming_output_file_metadata; ///< Per-receiver metadata attached to sink blocks.
		RealType _internal_stop_time = 0.0; ///< Physics stop time including IF over-render margin.
	};

	/**
	 * @brief Runs the unified, event-driven radar simulation.
	 *
	 * This function is the core entry point of the simulator. It advances time by
	 * processing events from a global priority queue. It handles both pulsed
	 * and continuous-wave (CW) physics, dispatching finalization tasks to
	 * worker threads for asynchronous processing.
	 *
	 * @param world A pointer to the simulation world containing all entities and state.
	 * @param pool A reference to the thread pool for executing tasks.
	 * @param progress_callback An optional callback function for reporting progress.
	 * @param output_dir Output directory for the simulation files.
	 */
	OutputMetadata runEventDrivenSim(World* world, pool::ThreadPool& pool,
									 const std::function<void(const std::string&, int, int)>& progress_callback,
									 const std::string& output_dir, const OutputConfig& output_config = OutputConfig{},
									 std::function<bool()> cancel_callback = nullptr, bool* cancelled = nullptr,
									 ReceiverOutputTelemetryCallback telemetry_callback = nullptr);
}
