// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file finalizer_pipeline.h
 * @brief Declares focused, testable pipeline steps for receiver finalization.
 *
 * This header defines the individual, single-responsibility functions that
 * constitute the signal processing pipeline for both pulsed and continuous-wave
 * receivers. By breaking the finalization process into these discrete steps,
 * each function becomes highly cohesive, easier to understand, and independently
 * testable.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <span>
#include <string>
#include <tuple>
#include <vector>

#include "core/config.h"
#include "core/simulation_state.h"

namespace timing
{
	class Timing;
}
namespace radar
{
	class Receiver;
	class Transmitter;
	class Target;
}
namespace serial
{
	class Response;
}
namespace core
{
	struct OutputFileMetadata;
}
namespace simulation
{
	struct CwPhaseNoiseLookup;
}

namespace processing::pipeline
{
	/// Half-open sample span used for dechirped LO-active finalization ranges.
	struct SampleSpan
	{
		std::size_t start = 0; ///< Inclusive first sample index in the span.
		std::size_t end_exclusive = 0; ///< Exclusive sample index at the end of the span.
	};

	/**
	 * @brief Advances the receiver's timing model to the start of the next processing window.
	 *
	 * This function handles the "dead time" between receive windows. For sync-on-pulse
	 * models, it resets the phase and skips to the window's start offset. For
	 * free-running models, it skips the number of samples corresponding to the
	 * inter-pulse period.
	 *
	 * @param timing_model The stateful timing model instance to advance.
	 * @param receiver The receiver whose properties (sync mode, PRF, etc.) determine how to advance the model.
	 * @param rate The oversampled simulation rate, used to calculate samples to skip.
	 */
	void advanceTimingModel(timing::Timing* timing_model, const radar::Receiver* receiver, RealType rate);

	/**
	 * @brief Calculates the jittered start time and fractional delay from a phase noise sample.
	 *
	 * Converts the first phase noise sample (in radians) of a window into a time
	 * jitter offset. It then decomposes the resulting "actual" start time into a
	 * component aligned with the sample clock and a fractional delay to be handled
	 * by the rendering engine's interpolation filter.
	 *
	 * @param ideal_start The perfect, jitter-free start time of the window.
	 * @param first_phase_noise The first phase noise sample (radians) for this window.
	 * @param carrier_freq The carrier frequency, needed to convert phase to time.
	 * @param rate The sampling rate, used for sample clock alignment.
	 * @return A tuple containing:
	 *         1. The sample-aligned start time (RealType).
	 *         2. The fractional sample delay (RealType).
	 */
	std::tuple<RealType, RealType> calculateJitteredStart(RealType ideal_start, RealType first_phase_noise,
														  RealType carrier_freq, RealType rate);

	/**
	 * @brief Applies streaming interference to a time window.
	 *
	 * Iterates through each sample of a processing window, calculating the combined
	 * direct and reflected path contributions from all active streaming sources at that
	 * precise moment in time. The resulting complex sample is added to the window.
	 *
	 * @param window The I/Q buffer for the receive window to which interference will be added.
	 * @param actual_start The jittered, sample-aligned start time of the window.
	 * @param dt The time step between samples (1.0 / rate).
	 * @param receiver The receiver being interfered with.
	 * @param streaming_sources A list of currently active streaming transmitters.
	 * @param targets The list of all targets for calculating reflected paths.
	 * @param tracker_cache Caller-owned reusable tracker storage for FMCW path boundary state.
	 */
	void applyStreamingInterference(std::span<ComplexType> window, RealType actual_start, RealType dt,
									const radar::Receiver* receiver,
									const std::vector<core::ActiveStreamingSource>& streaming_sources,
									const std::vector<std::unique_ptr<radar::Target>>* targets,
									core::ReceiverTrackerCache& tracker_cache,
									const simulation::CwPhaseNoiseLookup* phase_noise_lookup = nullptr);

	/**
	 * @brief Renders and applies pulsed interference to a streaming IQ buffer.
	 *
	 * Processes a log of `Response` objects that represent pulsed signals received
	 * during a streaming receiver's operation. Each response is rendered into a temporary
	 * buffer and then added to the main streaming IQ buffer at the correct time offset.
	 *
	 * @param iq_buffer The main, simulation-long IQ buffer for the streaming receiver.
	 * @param interference_log A list of `Response` objects representing the pulsed interference.
	 */
	void applyPulsedInterference(std::vector<ComplexType>& iq_buffer,
								 const std::vector<std::unique_ptr<serial::Response>>& interference_log);

	/**
	 * @brief Renders and applies pulsed interference only inside active sample spans.
	 *
	 * @param iq_buffer The main, simulation-long IQ buffer for the streaming receiver.
	 * @param interference_log A list of `Response` objects representing the pulsed interference.
	 * @param active_spans Half-open sample spans where LO-active output should receive interference.
	 * @param output_sample_rate The sample rate used by `iq_buffer`, in hertz.
	 */
	void applyPulsedInterference(std::vector<ComplexType>& iq_buffer,
								 const std::vector<std::unique_ptr<serial::Response>>& interference_log,
								 std::span<const SampleSpan> active_spans, RealType output_sample_rate);

	/**
	 * @brief Applies a pre-generated sequence of phase noise samples to an I/Q buffer.
	 *
	 * This function performs the complex multiplication `IQ_out = IQ_in * e^(j*phase_noise)`
	 * for each sample in the window, effectively modulating the phase of the signal.
	 *
	 * @param noise A span of phase noise samples in radians.
	 * @param window The I/Q buffer to be modified.
	 */
	void addPhaseNoiseToWindow(std::span<const RealType> noise, std::span<ComplexType> window);

	/**
	 * @brief Downsamples and quantizes an IQ buffer.
	 *
	 * This function performs the final processing steps. If oversampling is enabled,
	 * it first downsamples the buffer to the final output rate. It then simulates
	 * an ADC by quantizing and scaling the data.
	 *
	 * @param buffer The I/Q buffer to be processed. This is an in-out parameter; it will
	 *               be replaced by the downsampled buffer if applicable.
	 * @return The full-scale value (RealType) calculated during quantization, which is
	 *         needed for HDF5 metadata.
	 */
	RealType applyDownsamplingAndQuantization(std::vector<ComplexType>& buffer);

	/**
	 * @brief Exports a finalized streaming IQ buffer to an HDF5 file.
	 *
	 * Creates an HDF5 file, splits the complex buffer into real (I) and imaginary (Q)
	 * components, writes them as separate datasets, and adds relevant simulation
	 * parameters (sample rate, start time, etc.) as file attributes.
	 *
	 * @param filename The path to the output HDF5 file.
	 * @param iq_buffer The final, processed I/Q data to write.
	 * @param fullscale The full-scale value from quantization, saved as metadata.
	 * @param ref_freq The reference carrier frequency, saved as metadata.
	 * @param metadata Optional structured output metadata to attach to the file.
	 * @param sample_rate Output sample rate to store in the file, or `params::rate()` when unset.
	 */
	void exportStreamingToHdf5(const std::string& filename, const std::vector<ComplexType>& iq_buffer,
							   RealType fullscale, RealType ref_freq,
							   const core::OutputFileMetadata* metadata = nullptr, RealType sample_rate = 0.0);

}
