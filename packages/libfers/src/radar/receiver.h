// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file receiver.h
 * @brief Radar Receiver class for managing signal reception and response handling.
 */

#pragma once

#include <condition_variable>
#include <mutex>
#include <queue>
#include <random>

#include "core/rendering_job.h"
#include "core/sim_id.h"
#include "radar_obj.h"
#include "serial/response.h"

namespace pool
{
	class ThreadPool;
}

namespace radar
{
	/**
	 * @class Receiver
	 * @brief Manages radar signal reception and response processing.
	 */
	class Receiver final : public Radar
	{
	public:
		/**
		 * @enum RecvFlag
		 * @brief Enumeration for receiver configuration flags.
		 */
		enum class RecvFlag
		{
			FLAG_NODIRECT = 1, ///< Disable direct-path reception.
			FLAG_NOPROPLOSS = 2 ///< Disable propagation-loss scaling.
		};

		/**
		 * @brief Constructs a Receiver object.
		 *
		 * @param platform The platform associated with this receiver.
		 * @param name The name of the receiver.
		 * @param seed The seed for the receiver's internal random number generator.
		 * @param mode The operational mode (PULSED_MODE or CW_MODE).
		 */
		explicit Receiver(Platform* platform, std::string name, unsigned seed, OperationMode mode,
						  const SimId id = 0) noexcept;

		~Receiver() override = default;

		Receiver(const Receiver&) = delete;

		Receiver(Receiver&&) = delete;

		Receiver& operator=(const Receiver&) = delete;

		Receiver& operator=(Receiver&&) = delete;

		/**
		 * @brief Adds a response to the receiver's pulsed-mode inbox.
		 * @param response A unique pointer to the response object.
		 */
		void addResponseToInbox(std::unique_ptr<serial::Response> response) noexcept;

		/**
		 * @brief Adds a pulsed interference response to the receiver's streaming-mode log.
		 * @param response A unique pointer to the response object.
		 */
		void addInterferenceToLog(std::unique_ptr<serial::Response> response) noexcept;

		/**
		 * @brief Checks if a specific flag is set.
		 *
		 * @param flag The flag to check.
		 * @return True if the flag is set, false otherwise.
		 */
		[[nodiscard]] bool checkFlag(RecvFlag flag) const noexcept { return (_flags & static_cast<int>(flag)) != 0; }

		/**
		 * @brief Retrieves the unique ID of the receiver.
		 *
		 * @return The receiver SimId.
		 */
		[[nodiscard]] SimId getId() const noexcept { return Radar::getId(); }

		/**
		 * @brief Retrieves the noise temperature of the receiver.
		 *
		 * @return The noise temperature.
		 */
		[[nodiscard]] RealType getNoiseTemperature() const noexcept { return _noise_temperature; }

		/**
		 * @brief Retrieves the radar window length.
		 *
		 * @return The radar window length.
		 */
		[[nodiscard]] RealType getWindowLength() const noexcept { return _window_length; }

		/**
		 * @brief Retrieves the pulse repetition frequency (PRF) of the radar window.
		 *
		 * @return The PRF of the radar window.
		 */
		[[nodiscard]] RealType getWindowPrf() const noexcept { return _window_prf; }

		/**
		 * @brief Retrieves the window skip time.
		 *
		 * @return The window skip time.
		 */
		[[nodiscard]] RealType getWindowSkip() const noexcept { return _window_skip; }

		/**
		 * @brief Gets the noise temperature for a specific angle.
		 *
		 * @param angle The angle in spherical coordinates (SVec3).
		 * @return The noise temperature at the given angle.
		 */
		[[nodiscard]] RealType getNoiseTemperature(const math::SVec3& angle) const noexcept override;

		/**
		 * @brief Retrieves the start time of a specific radar window.
		 *
		 * @param window The index of the window.
		 * @return The start time of the specified window.
		 * @throws std::logic_error If the receiver is not associated with a timing source.
		 */
		[[nodiscard]] RealType getWindowStart(unsigned window) const;

		/**
		 * @brief Gets the number of radar windows.
		 *
		 * @return The total number of radar windows.
		 */
		[[nodiscard]] unsigned getWindowCount() const noexcept;

		/**
		 * @brief Gets the receiver's internal random number generator engine.
		 * @return A mutable reference to the RNG engine.
		 */
		[[nodiscard]] std::mt19937& getRngEngine() noexcept { return _rng; }

		/**
		 * @brief Gets the operational mode of the receiver.
		 *
		 * @return The operational mode (PULSED_MODE or CW_MODE).
		 */
		[[nodiscard]] OperationMode getMode() const noexcept { return _mode; }

		/**
		 * @brief Checks if the receiver is currently active (listening).
		 * @return True if active, false otherwise.
		 */
		[[nodiscard]] bool isActive() const noexcept { return _is_active; }

		/**
		 * @brief Sets the active state of the receiver.
		 * @param active The new active state.
		 */
		void setActive(const bool active) noexcept { _is_active = active; }

		/**
		 * @brief Sets the operational mode of the receiver.
		 * @param mode The operational mode (PULSED_MODE or CW_MODE).
		 */
		void setMode(OperationMode mode) noexcept { _mode = mode; }

		/**
		 * @brief Moves all responses from the inbox into a RenderingJob.
		 * @return A vector of unique pointers to the responses.
		 */
		std::vector<std::unique_ptr<serial::Response>> drainInbox() noexcept;

		/**
		 * @brief Adds a completed RenderingJob to the finalizer queue.
		 * @param job The RenderingJob to enqueue.
		 */
		void enqueueFinalizerJob(core::RenderingJob&& job);

		/**
		 * @brief Waits for and dequeues a RenderingJob from the finalizer queue.
		 *
		 * This is a blocking call, intended for use by the dedicated finalizer thread.
		 *
		 * @param job A reference to a RenderingJob to be filled.
		 * @return `false` if a shutdown signal is received, `true` otherwise.
		 */
		bool waitAndDequeueFinalizerJob(core::RenderingJob& job);

		/**
		 * @brief Sets the properties for radar windows.
		 *
		 * @param length The length of the radar window.
		 * @param prf The pulse repetition frequency.
		 * @param skip The skip time between windows.
		 */
		void setWindowProperties(RealType length, RealType prf, RealType skip) noexcept;

		/**
		 * @brief Sets a receiver flag.
		 *
		 * @param flag The flag to set.
		 */
		void setFlag(RecvFlag flag) noexcept { _flags |= static_cast<int>(flag); }

		/**
		 * @brief Clears a receiver flag.
		 *
		 * @param flag The flag to clear.
		 */
		void clearFlag(RecvFlag flag) noexcept { _flags &= ~static_cast<int>(flag); }

		/**
		 * @brief Sets the noise temperature of the receiver.
		 *
		 * @param temp The new noise temperature.
		 * @throws std::runtime_error If the noise temperature is negative.
		 */
		void setNoiseTemperature(RealType temp);

		/**
		 * @brief Prepares the internal storage for streaming IQ data.
		 * @param numSamples The total number of samples to allocate memory for.
		 */
		void prepareStreamingData(size_t numSamples);

		/**
		 * @brief Sets a single IQ sample at a specific index for streaming simulation.
		 * @param index The index at which to store the sample.
		 * @param sample The complex IQ sample.
		 */
		void setStreamingSample(size_t index, ComplexType sample);

		/**
		 * @brief Retrieves the collected streaming IQ data.
		 * @return A constant reference to the vector of complex IQ samples.
		 */
		[[nodiscard]] const std::vector<ComplexType>& getStreamingData() const { return _streaming_iq_data; }

		/**
		 * @brief Retrieves the collected streaming IQ data for modification.
		 * @return A mutable reference to the vector of complex IQ samples.
		 */
		[[nodiscard]] std::vector<ComplexType>& getMutableStreamingData() { return _streaming_iq_data; }

		/**
		 * @brief Retrieves the log of pulsed interferences for streaming modes.
		 * @return A const reference to the vector of interference responses.
		 */
		[[nodiscard]] const std::vector<std::unique_ptr<serial::Response>>& getPulsedInterferenceLog() const
		{
			return _pulsed_interference_log;
		}

		/**
		 * @brief Sets the active schedule for the receiver.
		 * @param schedule A vector of active periods.
		 */
		void setSchedule(std::vector<SchedulePeriod> schedule);

		/**
		 * @brief Retrieves the list of active reception periods.
		 * @return A const reference to the schedule vector.
		 */
		[[nodiscard]] const std::vector<SchedulePeriod>& getSchedule() const noexcept { return _schedule; }

		/**
		 * @brief Determines the next valid window start time at or after the given time.
		 *
		 * @param time The proposed window start time.
		 * @return The actual start time, or nullopt if no valid time exists in the schedule.
		 */
		[[nodiscard]] std::optional<RealType> getNextWindowTime(RealType time) const;

	private:
		// --- Common Members ---
		bool _is_active = false; ///< True while the receiver is active.
		RealType _noise_temperature = 0; ///< The noise temperature of the receiver.
		int _flags = 0; ///< Flags for receiver configuration.
		OperationMode _mode; ///< The operational mode of the receiver.
		std::mt19937 _rng; ///< Per-object random number generator for statistical independence.
		std::vector<SchedulePeriod> _schedule; ///< The schedule of active periods.

		// --- Pulsed Mode Members ---
		RealType _window_length = 0; ///< The length of the radar window.
		RealType _window_prf = 0; ///< The pulse repetition frequency (PRF) of the radar window.
		RealType _window_skip = 0; ///< The skip time between radar windows.
		std::vector<std::unique_ptr<serial::Response>>
			_inbox; ///< Mailbox for incoming Response objects during a receive window.
		std::mutex _inbox_mutex; ///< Mutex guarding the pulsed receive inbox.
		std::queue<core::RenderingJob> _finalizer_queue; ///< Completed windows waiting for final processing.
		std::mutex _finalizer_queue_mutex; ///< Mutex guarding the finalizer queue.
		std::condition_variable _finalizer_queue_cv; ///< Condition variable for finalizer queue updates.

		// --- CW Mode Members ---
		std::vector<std::unique_ptr<serial::Response>>
			_pulsed_interference_log; ///< Log of pulsed signals that interfere with CW reception.
		std::mutex _interference_log_mutex; ///< Mutex guarding the pulsed interference log.
		std::vector<ComplexType> _streaming_iq_data; ///< Buffer for raw, simulation-long I/Q data.
		std::mutex _cw_mutex; ///< Mutex for handling CW data.
	};
}
