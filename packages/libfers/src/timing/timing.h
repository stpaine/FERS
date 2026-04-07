// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file timing.h
 * @brief Timing source for simulation objects.
 *
 * All objects must adhere to a common timing source, which is modeled and adjusted by the methods in this class.
 */

#pragma once

#include <cstddef>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "core/config.h"
#include "core/sim_id.h"
#include "noise/noise_generators.h"

namespace timing
{
	class PrototypeTiming;

	/**
	 * @class Timing
	 * @brief Represents a timing source for simulation.
	 */
	class Timing final
	{
	public:
		/**
		 * @brief Constructs a Timing object.
		 *
		 * @param name The name of the timing source.
		 * @param seed The seed for the timing source's internal random number generator.
		 */
		explicit Timing(std::string name, unsigned seed, const SimId id = 0) noexcept;

		~Timing() = default;

		Timing(const Timing&) = delete;

		Timing& operator=(const Timing&) = delete;

		Timing(Timing&&) = delete;

		Timing& operator=(Timing&&) = delete;

		/**
		 * @brief Gets the next sample from the timing source.
		 *
		 * @return The next sample value or 0.0 if not enabled.
		 */
		[[nodiscard]] RealType getNextSample() const noexcept { return _enabled ? _model->getSample() : 0.0; }

		/**
		 * @brief Gets the name of the timing source.
		 *
		 * @return The name of the timing source.
		 */
		[[nodiscard]] std::string getName() const noexcept { return _name; }


		/**
		 * @brief Gets the unique ID of the timing source.
		 *
		 * @return The timing source SimId.
		 */
		[[nodiscard]] SimId getId() const noexcept { return _id; }

		/**
		 * @brief Gets the initial seed used for the timing source's RNG.
		 *
		 * @return The initial seed value.
		 */
		[[nodiscard]] unsigned getSeed() const noexcept { return _seed; }

		/**
		 * @brief Checks if the timing source synchronizes on pulse.
		 *
		 * @return True if synchronized on pulse, otherwise false.
		 */
		[[nodiscard]] bool getSyncOnPulse() const noexcept { return _sync_on_pulse; }

		/**
		 * @brief Gets the frequency of the timing source.
		 *
		 * @return The frequency of the timing source.
		 */
		[[nodiscard]] RealType getFrequency() const noexcept { return _frequency; }

		/**
		 * @brief Gets the frequency offset of the timing source.
		 * @return The frequency offset.
		 */
		[[nodiscard]] RealType getFreqOffset() const noexcept { return _freq_offset; }

		/**
		 * @brief Gets the phase offset of the timing source.
		 * @return The phase offset.
		 */
		[[nodiscard]] RealType getPhaseOffset() const noexcept { return _phase_offset; }

		/**
		 * @brief Checks if the timing source is enabled.
		 *
		 * @return True if enabled, otherwise false.
		 */
		[[nodiscard]] bool isEnabled() const noexcept { return _enabled && _model && _model->enabled(); }

		/**
		 * @brief Skips a number of samples in the timing model.
		 *
		 * @param samples The number of samples to skip.
		 */
		void skipSamples(std::size_t samples) noexcept;

		/**
		 * @brief Initializes the timing model.
		 *
		 * @param timing The prototype timing configuration used for initialization.
		 */
		void initializeModel(const PrototypeTiming* timing) noexcept;

		/**
		 * @brief Resets the timing model.
		 */
		// NOLINTNEXTLINE(readability-make-member-function-const)
		void reset() noexcept
		{
			if (_model)
			{
				_model->reset();
			}
		}

		/**
		 * @brief Creates a new Timing instance based on the same prototype.
		 * @return A unique_ptr to the new Timing object.
		 * @throws std::logic_error if the timing object was not initialized from a prototype.
		 */
		[[nodiscard]] std::unique_ptr<Timing> clone() const;

	private:
		std::string _name; ///< The name of the timing source.
		SimId _id; ///< Unique ID for this timing source.
		bool _enabled{false}; ///< Flag indicating if the timing source is enabled.
		std::unique_ptr<noise::ClockModelGenerator> _model{nullptr}; ///< Noise generator model for the timing source.
		std::vector<RealType> _alphas; ///< The alpha values for the noise generator model.
		std::vector<RealType> _weights; ///< The weights for the noise generator model.
		RealType _frequency{}; ///< The frequency of the timing source.
		RealType _freq_offset{}; ///< The frequency offset of the timing source.
		RealType _phase_offset{}; ///< The phase offset of the timing source.
		bool _sync_on_pulse{false}; ///< Flag indicating if the timing source synchronizes on pulse.
		std::mt19937 _rng; ///< Per-object random number generator for statistical independence.
		const PrototypeTiming* _prototype{nullptr}; ///< Pointer to the prototype used for initialization.
		unsigned _seed; ///< The initial seed for the RNG.
	};
}
