// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file object.h
 * @brief Base class for all physical objects in the radar system.
 */

#pragma once

#include "core/sim_id.h"
#include "platform.h"

namespace radar
{
	/**
	 * @class Object
	 * @brief Represents a physical object in the radar system.
	 */
	class Object
	{
	public:
		/**
		 * @brief Constructor for Object.
		 *
		 * @param platform Pointer to the Platform object associated with this Object.
		 * @param name The name of the Object.
		 */
		Object(Platform* platform, std::string name, const ObjectType type, const SimId id = 0) noexcept :
			_platform(platform), _id(id == 0 ? SimIdGenerator::instance().generateId(type) : id), _name(std::move(name))
		{
		}

		virtual ~Object() = default;
		Object(const Object&) = delete;
		Object& operator=(const Object&) = delete;
		Object(Object&&) noexcept = default;
		Object& operator=(Object&&) noexcept = default;

		/**
		 * @brief Retrieves the position of the object.
		 *
		 * @param time The time at which to get the position of the object.
		 * @return A math::Vec3 representing the position of the object.
		 */
		[[nodiscard]] math::Vec3 getPosition(const RealType time) const { return _platform->getPosition(time); }

		/**
		 * @brief Retrieves the rotation of the object.
		 *
		 * @param time The time at which to get the rotation of the object.
		 * @return A math::SVec3 representing the rotation of the object.
		 */
		[[nodiscard]] math::SVec3 getRotation(const RealType time) const { return _platform->getRotation(time); }

		/**
		 * @brief Retrieves the associated platform of the object.
		 *
		 * @return A pointer to the Platform object associated with this object.
		 */
		[[nodiscard]] Platform* getPlatform() const noexcept { return _platform; }

		/**
		 * @brief Retrieves the unique ID of the object.
		 *
		 * @return The unique SimId of the object.
		 */
		[[nodiscard]] SimId getId() const noexcept { return _id; }

		/**
		 * @brief Retrieves the name of the object.
		 *
		 * @return A const reference to the string representing the object's name.
		 */
		[[nodiscard]] const std::string& getName() const noexcept { return _name; }

		/**
		 * @brief Sets the name of the object.
		 *
		 * @param name The new name for the object.
		 */
		void setName(std::string name) noexcept { _name = std::move(name); }

	private:
		Platform* _platform; ///< Pointer to the Platform object associated with this Object.
		SimId _id; ///< Unique ID for this object.
		std::string _name; ///< The name of the Object.
	};
}
