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
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>
#include <queue>
#include <random>
#include <span>
#include <string>
#include <string_view>

#include "core/rendering_job.h"
#include "core/sim_id.h"
#include "core/simulation_state.h"
#include "radar_obj.h"
#include "serial/response.h"
#include "signal/if_resampler.h"

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

		/// Receiver-side FMCW dechirping mode.
		enum class DechirpMode
		{
			None, ///< Output raw pre-mix streaming IQ.
			Physical, ///< Preserve timing phase-noise decorrelation in the IF output.
			Ideal ///< Ignore timing phase noise during the dechirp mix.
		};

		/// Source used to construct the receiver LO reference.
		enum class DechirpReferenceSource
		{
			None, ///< No reference configured.
			Attached, ///< Use the attached transmitter.
			Transmitter, ///< Use a named transmitter.
			Custom ///< Use a named top-level waveform with the receiver schedule.
		};

		/// Parsed and resolved dechirp reference details.
		struct DechirpReference
		{
			DechirpReferenceSource source = DechirpReferenceSource::None; ///< Reference source type.
			std::string name; ///< Parsed transmitter or waveform name when applicable.
			SimId transmitter_id = 0; ///< Resolved transmitter ID for attached/transmitter references.
			std::string transmitter_name; ///< Resolved transmitter name for attached/transmitter references.
			SimId waveform_id = 0; ///< Resolved waveform ID for custom references.
			std::string waveform_name; ///< Resolved waveform name for custom references.
		};

		/// Receiver-local FMCW IF-chain request parsed from scenario input.
		struct FmcwIfChainRequest
		{
			std::optional<RealType> sample_rate_hz = std::nullopt; ///< Requested IF ADC sample rate in hertz.
			std::optional<RealType> filter_bandwidth_hz = std::nullopt; ///< One-sided IF passband edge in hertz.
			std::optional<RealType> filter_transition_width_hz =
				std::nullopt; ///< Optional IF filter transition width in hertz.
		};

		/**
		 * @brief Constructs a Receiver object.
		 *
		 * @param platform The platform associated with this receiver.
		 * @param name The name of the receiver.
		 * @param seed The seed for the receiver's internal random number generator.
		 * @param mode The operational mode (PULSED_MODE, CW_MODE, or FMCW_MODE).
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
		 * @return The operational mode (PULSED_MODE, CW_MODE, or FMCW_MODE).
		 */
		[[nodiscard]] OperationMode getMode() const noexcept { return _mode; }

		/// Gets the configured dechirp mode.
		[[nodiscard]] DechirpMode getDechirpMode() const noexcept { return _dechirp_mode; }

		/// Returns true when the receiver emits dechirped IF data.
		[[nodiscard]] bool isDechirpEnabled() const noexcept { return _dechirp_mode != DechirpMode::None; }

		/// Gets the configured dechirp reference.
		[[nodiscard]] const DechirpReference& getDechirpReference() const noexcept { return _dechirp_reference; }

		/// Gets the optional receiver-local FMCW IF-chain request.
		[[nodiscard]] const FmcwIfChainRequest& getFmcwIfChainRequest() const noexcept { return _fmcw_if_chain; }

		/// Gets the receiver-local FMCW IF sample rate in Hz.
		[[nodiscard]] std::optional<RealType> getIfSampleRate() const noexcept { return _fmcw_if_chain.sample_rate_hz; }

		/// Gets the receiver-local FMCW IF filter bandwidth in Hz.
		[[nodiscard]] std::optional<RealType> getIfFilterBandwidth() const noexcept
		{
			return _fmcw_if_chain.filter_bandwidth_hz;
		}

		/// Gets the receiver-local FMCW IF filter transition width in Hz.
		[[nodiscard]] std::optional<RealType> getIfFilterTransitionWidth() const noexcept
		{
			return _fmcw_if_chain.filter_transition_width_hz;
		}

		/// Returns true when this receiver requests IF-rate FMCW output.
		[[nodiscard]] bool hasFmcwIfSampleRate() const noexcept { return _fmcw_if_chain.sample_rate_hz.has_value(); }

		/// Returns true when this receiver is using the online FMCW IF resampling sink.
		[[nodiscard]] bool hasFmcwIfResamplingSink() const noexcept { return _fmcw_if_sink != nullptr; }

		/// Gets the active or most recently used IF resampling plan, if any.
		[[nodiscard]] const std::optional<fers_signal::FmcwIfResamplerPlan>& getFmcwIfResamplerPlan() const noexcept
		{
			return _fmcw_if_plan;
		}

		/// Gets resolved receive-time LO source segments.
		[[nodiscard]] const std::vector<core::ActiveStreamingSource>& getDechirpSources() const noexcept
		{
			return _dechirp_sources;
		}

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
		 * @param mode The operational mode (PULSED_MODE, CW_MODE, or FMCW_MODE).
		 */
		void setMode(OperationMode mode) noexcept;

		/// Sets the receiver-side dechirp mode.
		void setDechirpMode(DechirpMode mode) noexcept;

		/// Stores the unresolved dechirp reference parsed from scenario input.
		void setDechirpReference(DechirpReference reference);

		/// Stores the receiver-local FMCW IF-chain request.
		void setFmcwIfChainRequest(FmcwIfChainRequest request) noexcept;

		/// Creates the online FMCW IF resampling sink and clears the output buffer.
		void initializeFmcwIfResampling(fers_signal::FmcwIfResamplerPlan plan);

		/// Feeds one completed high-rate dechirped block into the online IF sink.
		void consumeFmcwIfBlock(std::span<const ComplexType> block);

		/// Flushes remaining samples from the online IF sink into the output buffer.
		void flushFmcwIfResampling();

		/// Replaces resolved receive-time LO source segments.
		void setResolvedDechirpSources(std::vector<core::ActiveStreamingSource> sources);

		/// Clears resolved dechirp source segments.
		void clearResolvedDechirpSources() noexcept { _dechirp_sources.clear(); }

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
		DechirpMode _dechirp_mode{DechirpMode::None}; ///< FMCW dechirp mode.
		DechirpReference _dechirp_reference{}; ///< Configured/resolved LO reference.
		std::vector<core::ActiveStreamingSource> _dechirp_sources; ///< Resolved LO sources.
		FmcwIfChainRequest _fmcw_if_chain{}; ///< Optional receiver-local IF-chain request.
		std::optional<fers_signal::FmcwIfResamplerPlan> _fmcw_if_plan; ///< Active IF resampling plan.
		std::unique_ptr<fers_signal::FmcwIfResamplingSink> _fmcw_if_sink; ///< Online IF resampling state.
		std::uint64_t _fmcw_if_samples_to_discard = 0; ///< Remaining startup IF samples to suppress.
	};

	/// Converts a dechirp mode to its scenario token.
	[[nodiscard]] std::string_view dechirpModeToken(Receiver::DechirpMode mode) noexcept;

	/// Parses a dechirp mode scenario token.
	[[nodiscard]] Receiver::DechirpMode parseDechirpModeToken(std::string_view token);

	/// Converts a dechirp reference source to its scenario token.
	[[nodiscard]] std::string_view dechirpReferenceSourceToken(Receiver::DechirpReferenceSource source) noexcept;

	/// Parses a dechirp reference source scenario token.
	[[nodiscard]] Receiver::DechirpReferenceSource parseDechirpReferenceSourceToken(std::string_view token);
}
