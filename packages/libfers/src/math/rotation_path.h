// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file rotation_path.h
 * @brief Defines the RotationPath class for handling rotational paths with different interpolation types.
 */

#pragma once

#include <vector>

#include "coord.h"
#include "core/config.h"
#include "geometry_ops.h"

namespace math
{
	/**
	 * @class RotationPath
	 * @brief Manages rotational paths with different interpolation techniques.
	 */
	class RotationPath
	{
	public:
		/**
		 * @enum InterpType
		 * @brief Enumeration for types of interpolation.
		 */
		enum class InterpType
		{
			INTERP_STATIC,
			INTERP_CONSTANT,
			INTERP_LINEAR,
			INTERP_CUBIC
		};

		/**
		 * @brief Constructor for RotationPath.
		 *
		 * @param type The type of interpolation (default is static).
		 */
		explicit RotationPath(const InterpType type = InterpType::INTERP_STATIC) noexcept : _type(type) {}

		~RotationPath() = default;

		RotationPath(const RotationPath&) = delete;

		RotationPath& operator=(const RotationPath&) = delete;

		RotationPath(RotationPath&&) = delete;

		RotationPath& operator=(RotationPath&&) = delete;

		/**
		 * @brief Adds a rotation coordinate to the path.
		 *
		 * @param coord The rotation coordinate to be added.
		 */
		void addCoord(const RotationCoord& coord) noexcept;

		/**
		 * @brief Finalizes the rotation path for interpolation.
		 */
		void finalize();

		/**
		 * @brief Gets the list of rotation coordinates.
		 *
		 * @return A const reference to the vector of rotation coordinates.
		 */
		[[nodiscard]] const std::vector<RotationCoord>& getCoords() const noexcept { return _coords; }

		/**
		 * @brief Gets the starting rotation coordinate.
		 *
		 * @return The starting rotation coordinate.
		 */
		[[nodiscard]] RotationCoord getStart() const noexcept { return _start; }

		/**
		 * @brief Gets the rate of change for the rotation.
		 *
		 * @return The rotation rate.
		 */
		[[nodiscard]] RotationCoord getRate() const noexcept { return _rate; }

		/**
		 * @brief Gets the interpolation type of the path.
		 *
		 * @return The interpolation type.
		 */
		[[nodiscard]] InterpType getType() const noexcept { return _type; }

		/**
		 * @brief Gets the rotational position at a given time.
		 *
		 * @param t The time value for which to calculate the position.
		 * @return The calculated position as an SVec3.
		 * @throws PathException if the path has not been finalized.
		 */
		[[nodiscard]] SVec3 getPosition(RealType t) const;

		/**
		 * @brief Sets the starting rotation coordinate.
		 *
		 * @param start The new starting rotation coordinate.
		 */
		void setStart(const RotationCoord& start) noexcept { _start = start; }

		/**
		 * @brief Sets the rate of change for the rotation.
		 *
		 * @param rate The new rate of change for the rotation.
		 */
		void setRate(const RotationCoord& rate) noexcept { _rate = rate; }

		/**
		 * @brief Sets the interpolation type for the path.
		 *
		 * @param setinterp The new interpolation type to be used.
		 */
		void setInterp(InterpType setinterp) noexcept;

		/**
		 * @brief Sets constant rate interpolation.
		 *
		 * @param setstart The starting rotation coordinate.
		 * @param setrate The rate of change for the rotation.
		 */
		void setConstantRate(const RotationCoord& setstart, const RotationCoord& setrate) noexcept;

	private:
		std::vector<RotationCoord> _coords; ///< Vector of rotation coordinates in the path.
		std::vector<RotationCoord> _dd; ///< Vector of second derivatives for cubic interpolation.
		bool _final{false}; ///< Boolean flag to indicate whether the path has been finalized.
		RotationCoord _start{}; ///< Starting rotation coordinate.
		RotationCoord _rate{}; ///< Rate of change for constant interpolation.
		InterpType _type{InterpType::INTERP_STATIC}; ///< Interpolation type used by the path.
	};
}
