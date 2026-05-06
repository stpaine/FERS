// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file prototype_timing.h
 * @brief Header file for the PrototypeTiming class.
 */

#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/config.h"
#include "core/sim_id.h"

namespace timing
{
	/**
	 * @class PrototypeTiming
	 * @brief Manages timing properties such as frequency, offsets, and synchronization.
	 */
	class PrototypeTiming
	{
	public:
		/**
		 * @brief Constructor for PrototypeTiming.
		 *
		 * @param name The name of the timing source.
		 */
		explicit PrototypeTiming(std::string name, const SimId id = 0) noexcept :
			_id(id == 0 ? SimIdGenerator::instance().generateId(ObjectType::Timing) : id), _name(std::move(name))
		{
		}

		~PrototypeTiming() = default;

		PrototypeTiming(const PrototypeTiming&) = default;

		PrototypeTiming(PrototypeTiming&&) = default;

		PrototypeTiming& operator=(const PrototypeTiming&) = default;

		PrototypeTiming& operator=(PrototypeTiming&&) = default;

		/**
		 * @brief Copies the alphas and weights vectors.
		 *
		 * @param alphas Reference to the vector where alpha values will be copied.
		 * @param weights Reference to the vector where weight values will be copied.
		 */
		void copyAlphas(std::vector<RealType>& alphas, std::vector<RealType>& weights) const noexcept;

		/**
		 * @brief Gets the current frequency.
		 *
		 * @return The current frequency value.
		 */
		[[nodiscard]] RealType getFrequency() const noexcept { return _frequency; }

		/**
		 * @brief Gets the name of the timing source.
		 *
		 * @return The name of the timing source.
		 */
		[[nodiscard]] std::string getName() const { return _name; }

		/**
		 * @brief Gets the unique ID of the timing source.
		 *
		 * @return The timing source SimId.
		 */
		[[nodiscard]] SimId getId() const noexcept { return _id; }

		/**
		 * @brief Checks if synchronization on pulse is enabled.
		 *
		 * @return True if synchronization on pulse is enabled, false otherwise.
		 */
		[[nodiscard]] bool getSyncOnPulse() const noexcept { return _sync_on_pulse; }

		/**
		 * @brief Gets the phase offset.
		 *
		 * @return The phase offset value.
		 */
		[[nodiscard]] std::optional<RealType> getPhaseOffset() const noexcept { return _phase_offset; }

		/**
		 * @brief Gets the frequency offset.
		 *
		 * @return The frequency offset value.
		 */
		[[nodiscard]] std::optional<RealType> getFreqOffset() const noexcept { return _freq_offset; }

		/** @brief Gets the random frequency-offset standard deviation. */
		[[nodiscard]] std::optional<RealType> getRandomFreqOffsetStdev() const noexcept { return _random_freq_stdev; }

		/** @brief Gets the random phase-offset standard deviation. */
		[[nodiscard]] std::optional<RealType> getRandomPhaseOffsetStdev() const noexcept { return _random_phase_stdev; }

		/**
		 * @brief Sets the frequency value.
		 *
		 * @param freq The frequency value to be set.
		 */
		void setFrequency(const RealType freq) noexcept { _frequency = freq; }

		/**
		 * @brief Enables synchronization on pulse.
		 */
		void setSyncOnPulse() noexcept { _sync_on_pulse = true; }

		/**
		 * @brief Disables synchronization on pulse.
		 */
		void clearSyncOnPulse() noexcept { _sync_on_pulse = false; }

		/**
		 * @brief Sets an alpha and weight value.
		 *
		 * @param alpha The alpha value to be added.
		 * @param weight The weight value to be added.
		 */
		void setAlpha(RealType alpha, RealType weight) noexcept;

		/**
		 * @brief Clears all noise entries.
		 */
		void clearNoiseEntries() noexcept;

		/**
		 * @brief Sets the frequency offset.
		 *
		 * @param offset The frequency offset to be set.
		 */
		void setFreqOffset(RealType offset) noexcept { _freq_offset = offset; }

		/**
		 * @brief Clears the frequency offset.
		 */
		void clearFreqOffset() noexcept { _freq_offset = std::nullopt; }

		/**
		 * @brief Sets the phase offset.
		 *
		 * @param offset The phase offset to be set.
		 */
		void setPhaseOffset(RealType offset) noexcept { _phase_offset = offset; }

		/**
		 * @brief Clears the phase offset.
		 */
		void clearPhaseOffset() noexcept { _phase_offset = std::nullopt; }

		/**
		 * @brief Sets a random frequency offset standard deviation.
		 *
		 * @param stdev The standard deviation for generating the random frequency offset.
		 */
		void setRandomFreqOffsetStdev(RealType stdev) noexcept { _random_freq_stdev = stdev; }

		/**
		 * @brief Sets a random phase offset standard deviation.
		 *
		 * @param stdev The standard deviation for generating the random phase offset.
		 */
		void setRandomPhaseOffsetStdev(RealType stdev) noexcept { _random_phase_stdev = stdev; }

		/**
		 * @brief Clears the random frequency offset standard deviation.
		 */
		void clearRandomFreqOffsetStdev() noexcept { _random_freq_stdev = std::nullopt; }

		/**
		 * @brief Clears the random phase offset standard deviation.
		 */
		void clearRandomPhaseOffsetStdev() noexcept { _random_phase_stdev = std::nullopt; }

		/**
		 * @brief Sets the name of the timing source.
		 *
		 * @param name The new name.
		 */
		void setName(std::string name) noexcept { _name = std::move(name); }

	private:
		SimId _id; ///< Unique ID for this timing source.
		std::string _name; ///< The name of the timing source.
		std::vector<RealType> _alphas; ///< Vector of alpha values.
		std::vector<RealType> _weights; ///< Vector of weight values.
		std::optional<RealType> _freq_offset; ///< Constant frequency offset.
		std::optional<RealType> _phase_offset; ///< Constant phase offset.
		std::optional<RealType> _random_phase_stdev; ///< Random phase offset standard deviation.
		std::optional<RealType> _random_freq_stdev; ///< Random frequency offset standard deviation.
		RealType _frequency{0}; ///< The frequency value.
		bool _sync_on_pulse{false}; ///< Flag indicating synchronization on pulse.
	};
}
