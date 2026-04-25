// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file radar_obj.h
 * @brief Defines the Radar class and associated functionality.
 */

#pragma once

#include "core/sim_id.h"
#include "object.h"

namespace timing
{
	class Timing;
}

namespace antenna
{
	class Antenna;
}

namespace radar
{
	class Platform;

	/**
	 * @enum OperationMode
	 * @brief Defines the operational mode of a radar component.
	 */
	enum class OperationMode
	{
		PULSED_MODE, ///< The component operates in a pulsed mode.
		CW_MODE, ///< The component operates in a continuous-wave mode.
		FMCW_MODE ///< The component operates in an FMCW streaming mode.
	};

	/**
	 * @class Radar
	 * @brief Represents a radar system on a platform.
	 */
	class Radar : public Object
	{
	public:
		/**
		 * @brief Constructs a Radar object.
		 *
		 * @param platform Pointer to the platform on which the radar is mounted.
		 * @param name Name of the radar object.
		 */
		Radar(Platform* platform, std::string name, const SimId id = 0) noexcept :
			Object(platform, std::move(name), ObjectType::Unknown,
				   id == 0 ? SimIdGenerator::instance().generateId(ObjectType::Unknown) : id)
		{
		}

		~Radar() override = default;

		Radar(const Radar&) = delete;

		Radar& operator=(const Radar&) = delete;

		Radar(Radar&&) = delete;

		Radar& operator=(Radar&&) = delete;

		/**
		 * @brief Retrieves the attached radar object.
		 *
		 * @return Pointer to the attached radar object.
		 */
		[[nodiscard]] const Radar* getAttached() const noexcept { return _attached; }

		/**
		 * @brief Retrieves the unique ID of the radar object.
		 *
		 * @return The radar SimId.
		 */
		[[nodiscard]] SimId getId() const noexcept { return Object::getId(); }

		/**
		 * @brief Gets the antenna associated with this radar.
		 *
		 * @return Pointer to the associated antenna.
		 */
		[[nodiscard]] const antenna::Antenna* getAntenna() const noexcept { return _antenna; }

		/**
		 * @brief Calculates the radar gain based on input angles and wavelength.
		 *
		 * @param angle The radar's pointing angle.
		 * @param refangle The reference angle for comparison.
		 * @param wavelength The wavelength of the radar signal.
		 * @return The calculated radar gain.
		 */
		[[nodiscard]] RealType getGain(const math::SVec3& angle, const math::SVec3& refangle,
									   RealType wavelength) const;

		/**
		 * @brief Gets the noise temperature of the radar.
		 *
		 * @param angle The angle at which the noise temperature is calculated.
		 * @return The calculated noise temperature.
		 */
		[[nodiscard]] virtual RealType getNoiseTemperature(const math::SVec3& angle) const noexcept;

		/**
		 * @brief Retrieves the timing source for the radar.
		 *
		 * @return Shared pointer to the timing source.
		 */
		[[nodiscard]] std::shared_ptr<timing::Timing> getTiming() const;

		/**
		 * @brief Sets the timing source for the radar.
		 *
		 * @param tim Shared pointer to the timing source to set.
		 */
		void setTiming(const std::shared_ptr<timing::Timing>& tim);

		/**
		 * @brief Sets the antenna for the radar.
		 *
		 * @param ant Pointer to the antenna to set.
		 */
		void setAntenna(const antenna::Antenna* ant);

		/**
		 * @brief Attaches another radar object to this radar.
		 *
		 * @param obj Pointer to the radar object to attach.
		 * @throws std::runtime_error If another object is already attached.
		 */
		void setAttached(const Radar* obj);

	protected:
		std::shared_ptr<timing::Timing> _timing; ///< Timing source for the radar.

	private:
		const antenna::Antenna* _antenna{nullptr}; ///< Antenna object associated with the radar.
		const Radar* _attached{nullptr}; ///< Attached radar object.
	};
}
