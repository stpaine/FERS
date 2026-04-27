// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#pragma once

#include "core/config.h"
#include "core/parameters.h"
#include "math/coord.h"

namespace serial::rotation_angle_utils
{
	/// Converts a value from the given external unit to radians.
	inline RealType unit_to_radians(const RealType value, const params::RotationAngleUnit unit) noexcept
	{
		return unit == params::RotationAngleUnit::Degrees ? value * (PI / 180.0) : value;
	}

	/// Converts a radian value to the requested external unit.
	inline RealType radians_to_unit(const RealType value, const params::RotationAngleUnit unit) noexcept
	{
		return unit == params::RotationAngleUnit::Degrees ? value * 180.0 / PI : value;
	}

	/// Converts external compass azimuth/elevation into internal rotation coordinates.
	inline math::RotationCoord external_rotation_to_internal(const RealType azimuth, const RealType elevation,
															 const RealType time,
															 const params::RotationAngleUnit unit) noexcept
	{
		return {(PI / 2.0) - unit_to_radians(azimuth, unit), unit_to_radians(elevation, unit), time};
	}

	/// Converts external compass azimuth/elevation rates into internal rotation rates.
	inline math::RotationCoord external_rotation_rate_to_internal(const RealType azimuth_rate,
																  const RealType elevation_rate, const RealType time,
																  const params::RotationAngleUnit unit) noexcept
	{
		return {-unit_to_radians(azimuth_rate, unit), unit_to_radians(elevation_rate, unit), time};
	}

	/// Converts an internal azimuth angle to the external compass convention.
	inline RealType internal_azimuth_to_external(const RealType azimuth, const params::RotationAngleUnit unit) noexcept
	{
		return radians_to_unit((PI / 2.0) - azimuth, unit);
	}

	/// Converts an internal elevation angle to the external unit.
	inline RealType internal_elevation_to_external(const RealType elevation,
												   const params::RotationAngleUnit unit) noexcept
	{
		return radians_to_unit(elevation, unit);
	}

	/// Converts an internal azimuth rate to the external compass convention.
	inline RealType internal_azimuth_rate_to_external(const RealType azimuth_rate,
													  const params::RotationAngleUnit unit) noexcept
	{
		return radians_to_unit(-azimuth_rate, unit);
	}

	/// Converts an internal elevation rate to the external unit.
	inline RealType internal_elevation_rate_to_external(const RealType elevation_rate,
														const params::RotationAngleUnit unit) noexcept
	{
		return radians_to_unit(elevation_rate, unit);
	}
}
