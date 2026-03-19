// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file signal_processor.h
 * @brief Header for receiver-side signal processing and rendering.
 *
 * This file contains declarations for functions that perform digital signal
 * processing on received radar signals. This includes rendering raw responses
 * into time-domain I/Q samples, adding thermal noise, and simulating ADC quantization.
 */

#pragma once

#include <memory>
#include <random>
#include <span>
#include <vector>

#include "core/config.h"

namespace serial
{
	class Response;
}

namespace pool
{
	class ThreadPool;
}

namespace processing
{
	/**
	 * @brief Renders a time-window of I/Q data from a collection of raw radar responses.
	 *
	 * This function orchestrates the process of converting abstract `Response` objects
	 * into a concrete vector of complex I/Q samples for a specific time window.
	 * It handles the superposition of multiple signals arriving at the receiver
	 * during the window and can use a thread pool for parallel processing.
	 *
	 * @param window The output vector where the rendered I/Q samples will be added.
	 * @param length The duration of the time window in seconds.
	 * @param start The start time of the window in seconds.
	 * @param fracDelay A fractional sample delay to apply for fine-grained timing.
	 * @param responses A span of unique pointers to the `Response` objects to be rendered.
	 */
	void renderWindow(std::vector<ComplexType>& window, RealType length, RealType start, RealType fracDelay,
					  std::span<const std::unique_ptr<serial::Response>> responses);

	/**
	 * @brief Applies thermal (Johnson-Nyquist) noise to a window of I/Q samples.
	 *
	 * Simulates the addition of white Gaussian noise based on the receiver's noise
	 * temperature and the simulation bandwidth.
	 *
	 * @param window A span of I/Q data to which noise will be added.
	 * @param noiseTemperature The effective noise temperature in Kelvin.
	 * @param rngEngine A random number generator engine to use for noise generation.
	 */
	void applyThermalNoise(std::span<ComplexType> window, RealType noiseTemperature, std::mt19937& rngEngine);

	/**
	 * @brief Simulates ADC quantization and scales a window of complex I/Q samples.
	 *
	 * This function first finds the maximum absolute value in the I/Q data to determine
	 * the full-scale range. It then simulates the quantization process based on the
	 * configured number of ADC bits. If no quantization is specified (adc_bits=0), it
	 * normalizes the data to a maximum amplitude of 1.0.
	 *
	 * @param window The window of complex I/Q samples to quantize and scale.
	 * @return The full-scale value used for quantization/normalization.
	 */
	RealType quantizeAndScaleWindow(std::span<ComplexType> window);
}
