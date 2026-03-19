// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file falpha_branch.h
 * @brief Implementation of the FAlphaBranch class for noise generation.
 */

#pragma once

#include <functional>
#include <memory>
#include <random>
#include <vector>

#include "core/config.h"
#include "signal/dsp_filters.h"

namespace noise
{
	/**
	 * @class FAlphaBranch
	 * @brief Class responsible for generating fractional and integer noise components.
	 *
	 * The FAlphaBranch class generates noise by applying a fractional integrator filter. It uses a series of
	 * filters and upsamplers to process and shape the noise signal.
	 */
	class FAlphaBranch
	{
	public:
		/**
		 * @brief Constructor for FAlphaBranch.
		 *
		 * @param rngEngine The random number generator engine to use.
		 * @param ffrac Fractional part of the noise generation (e.g., 0.5 for 1/f noise).
		 * @param fint Integer part of the noise generation (e.g., 1 for integration).
		 * @param pre Previous stage of the FAlphaBranch for recursive noise processing.
		 * @param last Specifies if this is the last branch in the chain of processing.
		 */
		FAlphaBranch(std::mt19937& rngEngine, RealType ffrac, unsigned fint, std::unique_ptr<FAlphaBranch> pre,
					 bool last);

		~FAlphaBranch() = default;

		FAlphaBranch(const FAlphaBranch&) = delete;

		FAlphaBranch& operator=(const FAlphaBranch&) = delete;

		FAlphaBranch(FAlphaBranch&&) = delete;

		FAlphaBranch& operator=(FAlphaBranch&&) = delete;

		/**
		 * @brief Retrieves the current noise sample.
		 *
		 * @return The current noise sample.
		 */
		RealType getSample() noexcept;

		/**
		 * @brief Flushes the branch with a new scaling factor.
		 *
		 * @param scale New scale factor to apply to the previous stage.
		 */
		void flush(RealType scale);

		/**
		 * @brief Retrieves the previous branch in the chain.
		 *
		 * @return Pointer to the previous FAlphaBranch, or nullptr if none exists.
		 */
		[[nodiscard]] FAlphaBranch* getPre() const noexcept { return _pre.get(); }

	private:
		/// Initializes the filters and sets up the initial state of the noise generator.
		void init();

		/// Refills the sample buffer with new upsampled noise values.
		void refill() noexcept;

		/// Calculates a new noise sample.
		RealType calcSample() noexcept;

		std::reference_wrapper<std::mt19937> _rng_engine_ref; ///< Reference to the RNG engine.

		std::normal_distribution<> _normal_dist; ///< Normal distribution for generating white Gaussian noise.

		std::unique_ptr<fers_signal::IirFilter> _shape_filter; ///< Filter used for shaping the noise signal.

		std::unique_ptr<fers_signal::IirFilter> _integ_filter; ///< Filter used for integrating the noise signal.

		std::unique_ptr<fers_signal::IirFilter> _highpass; ///< High-pass filter to remove low-frequency components.

		std::unique_ptr<fers_signal::DecadeUpsampler>
			_upsampler; ///< Upsampler for generating higher-frequency components.

		std::unique_ptr<FAlphaBranch> _pre; ///< Previous FAlphaBranch in the chain for recursive noise processing.

		RealType _shape_gain{1.0}; ///< Gain factor for shaping filter.

		RealType _integ_gain{1.0}; ///< Gain factor for integration filter.

		RealType _upsample_scale{}; ///< Scaling factor for the upsampled noise.

		std::vector<RealType> _buffer{}; ///< Buffer for storing upsampled noise samples.

		unsigned _buffer_samples{}; ///< Number of samples currently in the buffer.

		RealType _ffrac; ///< Fractional part of the noise generation.

		unsigned _fint; ///< Integer part of the noise generation.

		RealType _offset_sample{}; ///< Offset applied to the final noise sample.

		bool _got_offset{false}; ///< Flag indicating if the offset sample has been retrieved.

		RealType _pre_scale{1.0}; ///< Scale factor for the previous stage's noise output.

		bool _last; ///< Indicates if this is the last branch in the noise processing chain.
	};
}
