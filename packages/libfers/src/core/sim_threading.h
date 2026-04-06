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
#include "core/sim_events.h"

namespace pool
{
	class ThreadPool;
}

namespace radar
{
	class Receiver;
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
						 std::string output_dir);

		/**
		 * @brief Starts and runs the main simulation loop until completion.
		 */
		void run();

		/**
		 * @brief Advances the time-stepped inner loop for active continuous-wave (CW) systems.
		 * @param t_event The timestamp of the next discrete event to process up to.
		 */
		void processCwPhysics(RealType t_event);

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
		 * @brief Handles a continuous-wave transmitter turning on.
		 * @param tx Pointer to the transmitting radar object.
		 */
		void handleTxCwStart(radar::Transmitter* tx);

		/**
		 * @brief Handles a continuous-wave transmitter turning off.
		 * @param tx Pointer to the transmitting radar object.
		 */
		void handleTxCwEnd(radar::Transmitter* tx);

		/**
		 * @brief Handles a continuous-wave receiver starting to record.
		 * @param rx Pointer to the receiving radar object.
		 */
		void handleRxCwStart(radar::Receiver* rx);

		/**
		 * @brief Handles a continuous-wave receiver stopping recording.
		 * @param rx Pointer to the receiving radar object.
		 */
		void handleRxCwEnd(radar::Receiver* rx);

		/**
		 * @brief Calculates the total complex I/Q sample for a receiver at a specific time step.
		 * @param rx Pointer to the receiving radar object.
		 * @param t_step The exact simulation time for the sample.
		 * @param cw_sources A list of currently active continuous-wave transmitters.
		 * @return The calculated complex I/Q sample combining direct and reflected paths.
		 */
		[[nodiscard]] ComplexType calculateCwSample(radar::Receiver* rx, RealType t_step,
													const std::vector<radar::Transmitter*>& cw_sources) const;

	private:
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
		 * @brief Initiates the shutdown phase, waiting for all asynchronous tasks to complete.
		 */
		void shutdown();

		World* _world; ///< Pointer to the simulation world state.
		pool::ThreadPool& _pool; ///< Reference to the global thread pool.
		std::shared_ptr<ProgressReporter> _reporter; ///< Shared progress reporter instance.
		std::vector<std::jthread> _finalizer_threads; ///< Collection of dedicated pulsed finalizer threads.

		std::chrono::steady_clock::time_point _last_report_time; ///< Timestamp of the last progress report.
		int _last_reported_percent = -1; ///< The last reported percentage to prevent redundant updates.

		std::string _output_dir; ///< Output directory for the simulation files.
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
	void runEventDrivenSim(World* world, pool::ThreadPool& pool,
						   const std::function<void(const std::string&, int, int)>& progress_callback,
						   const std::string& output_dir);
}
