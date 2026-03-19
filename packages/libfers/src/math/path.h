// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file path.h
 * @brief Provides the definition and functionality of the Path class
 * for handling coordinate-based paths with different interpolation types.
 */

#pragma once

#include <vector>

#include "coord.h"
#include "core/config.h"

namespace math
{
	class Vec3;

	/**
	 * @class Path
	 * @brief Represents a path with coordinates and allows for various interpolation methods.
	 */
	class Path
	{
	public:
		/**
		 * @brief Types of interpolation supported by the Path class.
		 */
		enum class InterpType
		{
			INTERP_STATIC,
			INTERP_LINEAR,
			INTERP_CUBIC
		};

		/**
		 * @brief Constructs a Path object with a specified interpolation type.
		 *
		 * @param type The interpolation type (default is INTERP_STATIC).
		 */
		explicit Path(const InterpType type = InterpType::INTERP_STATIC) noexcept : _type(type) {}

		~Path() = default;

		Path(const Path&) = delete;

		Path(Path&&) = delete;

		Path& operator=(const Path&) = delete;

		Path& operator=(Path&&) = delete;

		/**
		 * @brief Adds a coordinate to the path.
		 *
		 * @param coord The coordinate to be added.
		 */
		void addCoord(const Coord& coord) noexcept;

		/**
		 * @brief Finalizes the path, preparing it for interpolation.
		 */
		void finalize();

		/**
		 * @brief Retrieves the current interpolation type of the path.
		 *
		 * @return The interpolation type of the path.
		 */
		[[nodiscard]] InterpType getType() const noexcept { return _type; }

		/**
		 * @brief Gets the list of coordinates in the path.
		 *
		 * @return A constant reference to the vector of coordinates.
		 */
		[[nodiscard]] const std::vector<Coord>& getCoords() const noexcept { return _coords; }

		/**
		 * @brief Retrieves the position at a given time along the path.
		 *
		 * @param t The time parameter at which to get the position.
		 * @return The interpolated position as a Vec3 object.
		 * @throws PathException If finalize() has not been called before this method.
		 */
		[[nodiscard]] Vec3 getPosition(RealType t) const;

		/**
		 * @brief Retrieves the velocity at a given time along the path.
		 *
		 * @param t The time parameter at which to get the velocity.
		 * @return The interpolated velocity as a Vec3 object.
		 * @throws PathException If finalize() has not been called before this method.
		 */
		[[nodiscard]] Vec3 getVelocity(RealType t) const;

		/**
		 * @brief Changes the interpolation type.
		 *
		 * @param settype The new interpolation type to be used.
		 */
		void setInterp(InterpType settype) noexcept;

	private:
		std::vector<Coord> _coords; ///< The list of coordinates in the path.
		std::vector<Coord> _dd; ///< The list of second derivatives for cubic interpolation.
		bool _final{false}; ///< Flag indicating whether the path has been finalized.
		InterpType _type; ///< The current interpolation type of the path.
	};
}
