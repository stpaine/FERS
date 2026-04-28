// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file geometry_ops.cpp
 * @brief Implementation of geometry classes.
 */

#include "geometry_ops.h"

namespace
{
	/**
	 * @brief Normalizes an angle to the signed principal range.
	 *
	 * Wraps the input angle modulo 2π and returns the equivalent signed angle in
	 * the range (-π, π]. This is useful for angular differences where the shortest
	 * signed rotation is required.
	 *
	 * @param angle The input angle in radians.
	 * @return The equivalent wrapped angle in radians.
	 */
	[[nodiscard]] RealType wrapSignedAngle(const RealType angle) noexcept
	{
		RealType wrapped = std::remainder(angle, 2 * PI);

		// Preserve the previous operator- behavior:
		// old code mapped -PI to PI, giving the range (-PI, PI].
		if (wrapped <= -PI)
		{
			wrapped += 2 * PI;
		}

		return wrapped;
	}
}

namespace math
{
	Vec3::Vec3(const SVec3& svec) noexcept :
		x(svec.length * std::cos(svec.azimuth) * std::cos(svec.elevation)),
		y(svec.length * std::sin(svec.azimuth) * std::cos(svec.elevation)), z(svec.length * std::sin(svec.elevation))
	{
	}

	Vec3& Vec3::operator+=(const Vec3& b) noexcept
	{
		x += b.x;
		y += b.y;
		z += b.z;
		return *this;
	}

	Vec3& Vec3::operator-=(const Vec3& b) noexcept
	{
		x -= b.x;
		y -= b.y;
		z -= b.z;
		return *this;
	}

	Vec3& Vec3::operator*=(const Vec3& b) noexcept
	{
		x *= b.x;
		y *= b.y;
		z *= b.z;
		return *this;
	}

	Vec3& Vec3::operator*=(const Matrix3& m) noexcept
	{
		const RealType* mat = m.getData();
		const Vec3 v(x, y, z);
		x = mat[0] * v.x + mat[1] * v.y + mat[2] * v.z;
		y = mat[3] * v.x + mat[4] * v.y + mat[5] * v.z;
		z = mat[6] * v.x + mat[7] * v.y + mat[8] * v.z;
		return *this;
	}

	Vec3& Vec3::operator*=(const RealType b) noexcept
	{
		x *= b;
		y *= b;
		z *= b;
		return *this;
	}

	Vec3& Vec3::operator/=(const RealType b) noexcept
	{
		x /= b;
		y /= b;
		z /= b;
		return *this;
	}

	SVec3::SVec3(const Vec3& vec) noexcept : length(vec.length())
	{
		elevation = std::asin(vec.z / length);
		azimuth = std::atan2(vec.y, vec.x);
	}

	SVec3& SVec3::operator*=(const RealType b) noexcept
	{
		length *= b;
		return *this;
	}

	SVec3& SVec3::operator/=(const RealType b) noexcept
	{
		length /= b;
		return *this;
	}

	SVec3 operator+(const SVec3& a, const SVec3& b) noexcept
	{
		const RealType new_azimuth = wrapSignedAngle(a.azimuth + b.azimuth);
		const RealType new_elevation = std::fmod(a.elevation + b.elevation, PI);

		return {a.length + b.length, new_azimuth, new_elevation};
	}

	SVec3 operator-(const SVec3& a, const SVec3& b) noexcept
	{
		const RealType new_azimuth = wrapSignedAngle(a.azimuth - b.azimuth);

		// Elevation difference is typically simpler and doesn't need wrapping
		// unless 180 degree elevation sweeps are expected.
		const RealType new_elevation = std::fmod(a.elevation - b.elevation, PI);

		return {a.length - b.length, new_azimuth, new_elevation};
	}
}
