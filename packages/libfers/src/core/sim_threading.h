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
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "core/config.h"
#include "core/output_metadata.h"
#include "core/parameters.h"
#include "core/sim_events.h"
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
				std::scoped_lock lock(_mutex);
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
						 std::string output_dir, std::shared_ptr<OutputMetadataCollector> metadata_collector = nullptr);

		/**
		 * @brief Starts and runs the main simulation loop until completion.
		 */
		void run();

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
		/// Per-receiver FMCW tracker state for direct and reflected streaming paths.
		struct ReceiverTrackerCache
		{
			std::vector<FmcwChirpBoundaryTracker> direct; ///< Trackers for direct paths by source index.
			std::vector<std::vector<FmcwChirpBoundaryTracker>> reflected; ///< Trackers for reflected paths.
		};

		/// Calculates one streaming I/Q sample for the receiver at the specified time step.
		[[nodiscard]] ComplexType calculateStreamingSample(radar::Receiver* rx, RealType t_step,
														   const std::vector<ActiveStreamingSource>& streaming_sources,
														   ReceiverTrackerCache& tracker_cache) const;

		/// Adds tracker storage for a newly active streaming source.
		void appendStreamingTrackerSource();

		/// Removes tracker storage for a streaming source that has ended.
		void eraseStreamingTrackerSource(std::size_t source_index);

		/// Creates the CW phase-noise lookup if any active timing source needs it.
		void ensureCwPhaseNoiseLookup();

		/// Emits summary logs for streaming receiver buffers.
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

		std::chrono::steady_clock::time_point _last_report_time; ///< Timestamp of the last progress report.
		int _last_reported_percent = -1; ///< The last reported percentage to prevent redundant updates.

		std::string _output_dir; ///< Output directory for the simulation files.
		std::unique_ptr<simulation::CwPhaseNoiseLookup> _cw_phase_noise_lookup; ///< Cached CW phase-noise lookup.
		std::vector<ReceiverTrackerCache> _streaming_tracker_caches; ///< Per-receiver streaming tracker caches.
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
									 const std::string& output_dir);
}
