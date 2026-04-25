// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file transmitter.h
 * @brief Header file for the Transmitter class in the radar namespace.
 */

#pragma once

#include <optional>

#include "core/sim_id.h"
#include "radar_obj.h"
#include "schedule_period.h"

namespace fers_signal
{
	class FmcwChirpSignal;
	class RadarSignal;
}

namespace radar
{
	/**
	 * @class Transmitter
	 * @brief Represents a radar transmitter system.
	 */
	class Transmitter final : public Radar
	{
	public:
		/**
		 * @brief Constructor for the Transmitter class.
		 *
		 * @param platform Pointer to the platform object.
		 * @param name Name of the transmitter.
		 * @param mode The operational mode (PULSED_MODE or CW_MODE).
		 */
		Transmitter(Platform* platform, std::string name, const OperationMode mode, const SimId id = 0) noexcept :
			Radar(platform, std::move(name),
				  id == 0 ? SimIdGenerator::instance().generateId(ObjectType::Transmitter) : id),
			_mode(mode)
		{
		}

		~Transmitter() override = default;

		Transmitter(const Transmitter&) = delete;

		Transmitter& operator=(const Transmitter&) = delete;

		Transmitter(Transmitter&&) = delete;

		Transmitter& operator=(Transmitter&&) = delete;

		/**
		 * @brief Retrieves the pulse repetition frequency (PRF).
		 *
		 * @return The current PRF of the transmitter.
		 */
		[[nodiscard]] RealType getPrf() const noexcept { return _prf; }

		/**
		 * @brief Retrieves the radar signal currently being transmitted.
		 *
		 * @return Pointer to the RadarSignal object being transmitted.
		 */
		[[nodiscard]] fers_signal::RadarSignal* getSignal() const noexcept { return _signal; }

		[[nodiscard]] const fers_signal::FmcwChirpSignal* getFmcwSignal() const noexcept;

		[[nodiscard]] bool isStreamingMode() const noexcept
		{
			return _mode == OperationMode::CW_MODE || _mode == OperationMode::FMCW_MODE;
		}

		/**
		 * @brief Retrieves the unique ID of the transmitter.
		 *
		 * @return The transmitter SimId.
		 */
		[[nodiscard]] SimId getId() const noexcept { return Radar::getId(); }

		/**
		 * @brief Gets the operational mode of the transmitter.
		 *
		 * @return The operational mode (PULSED_MODE or CW_MODE).
		 */
		[[nodiscard]] OperationMode getMode() const noexcept { return _mode; }

		/**
		 * @brief Sets the operational mode of the transmitter.
		 *
		 * @param mode The operational mode (PULSED_MODE or CW_MODE).
		 */
		void setMode(OperationMode mode) noexcept { _mode = mode; }

		/**
		 * @brief Sets the radar signal wave to be transmitted.
		 *
		 * @param pulse Pointer to the RadarSignal object representing the wave.
		 */
		void setWave(fers_signal::RadarSignal* pulse) noexcept { _signal = pulse; }

		/**
		 * @brief Sets the radar signal wave to be transmitted.
		 *
		 * @param signal Pointer to the RadarSignal object to be transmitted.
		 */
		void setSignal(fers_signal::RadarSignal* signal) noexcept { _signal = signal; }

		/**
		 * @brief Sets the pulse repetition frequency (PRF) of the transmitter.
		 *
		 * @param mprf Desired PRF to be set.
		 */
		void setPrf(RealType mprf) noexcept;

		/**
		 * @brief Sets the active schedule for the transmitter.
		 *
		 * The provided schedule should be pre-validated and sorted.
		 * @param schedule A vector of active periods.
		 */
		void setSchedule(std::vector<SchedulePeriod> schedule);

		/**
		 * @brief Retrieves the list of active transmission periods.
		 * @return A const reference to the schedule vector.
		 */
		[[nodiscard]] const std::vector<SchedulePeriod>& getSchedule() const noexcept { return _schedule; }

		/**
		 * @brief Determines the valid simulation time for a pulse at or after the given time.
		 *
		 * If the requested time falls within an active period, it is returned.
		 * If it falls in a gap between periods, the start of the next period is returned.
		 * If it falls after all periods, std::nullopt is returned.
		 * If no schedule is defined, the transmitter is considered "always on".
		 *
		 * @param time The proposed pulse time.
		 * @return The actual pulse time, or nullopt if no valid time exists.
		 */
		[[nodiscard]] std::optional<RealType> getNextPulseTime(RealType time) const;

	private:
		fers_signal::RadarSignal* _signal = nullptr; ///< Pointer to the radar signal being transmitted.

		RealType _prf = {}; ///< The pulse repetition frequency (PRF) of the transmitter.

		OperationMode _mode; ///< The operational mode of the transmitter.
		std::vector<SchedulePeriod> _schedule; ///< The schedule of active periods.
	};

}
