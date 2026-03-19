// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file geometry_ops.h
 * @brief Classes and operations for 3D geometry.
 */

#pragma once

#include <array>
#include <cmath>

#include "core/config.h"

namespace math
{
	class Vec3;

	/**
	 * @class Matrix3
	 * @brief A class representing a 3x3 matrix.
	 */
	class Matrix3
	{
	public:
		std::array<RealType, 9> elements{}; ///< The 3x3 matrix elements

		/**
		 * @brief Get the matrix data as a constant pointer.
		 *
		 * @return A constant pointer to the matrix elements.
		 */
		[[nodiscard]] const RealType* getData() const noexcept { return elements.data(); }

		/**
		 * @brief Get the matrix data as a modifiable pointer.
		 *
		 * @return A pointer to the matrix elements.
		 */
		RealType* getData() noexcept { return elements.data(); }
	};

	/**
	 * @class SVec3
	 * @brief A class representing a vector in spherical coordinates.
	 */
	class SVec3
	{
	public:
		RealType length{}; ///< The length of the vector
		RealType azimuth{}; ///< The azimuth angle of the vector
		RealType elevation{}; ///< The elevation angle of the vector

		SVec3() noexcept = default;
		SVec3(const SVec3&) noexcept = default;
		SVec3(SVec3&&) noexcept = default;
		SVec3& operator=(const SVec3&) noexcept = default;
		SVec3& operator=(SVec3&&) noexcept = default;

		/**
		 * @brief Parameterized constructor for SVec3.
		 *
		 * @param length The length or magnitude of the vector.
		 * @param azimuth The azimuthal angle of the vector.
		 * @param elevation The elevation angle of the vector.
		 */
		constexpr SVec3(const RealType length, const RealType azimuth, const RealType elevation) noexcept :
			length(length), azimuth(azimuth), elevation(elevation)
		{
		}

		/**
		 * @brief Constructs a spherical vector from a rectangular Vec3.
		 *
		 * @param vec A rectangular vector (Vec3) to convert.
		 */
		explicit SVec3(const Vec3& vec) noexcept;

		/**
		 * @brief Scalar multiplication assignment for SVec3.
		 *
		 * @param b The scalar value.
		 * @return A reference to the updated SVec3.
		 */
		SVec3& operator*=(RealType b) noexcept;

		/**
		 * @brief Scalar division assignment for SVec3.
		 *
		 * @param b The scalar value.
		 * @return A reference to the updated SVec3.
		 */
		SVec3& operator/=(RealType b) noexcept;
	};

	/**
	 * @class Vec3
	 * @brief A class representing a vector in rectangular coordinates.
	 */
	class Vec3
	{
	public:
		RealType x{}; ///< The x component of the vector
		RealType y{}; ///< The y component of the vector
		RealType z{}; ///< The z component of the vector

		Vec3() noexcept = default;
		Vec3(const Vec3&) noexcept = default;
		Vec3(Vec3&&) noexcept = default;
		Vec3& operator=(const Vec3&) noexcept = default;
		Vec3& operator=(Vec3&&) noexcept = default;

		/**
		 * @brief Parameterized constructor for Vec3.
		 *
		 * @param x The x component.
		 * @param y The y component.
		 * @param z The z component.
		 */
		constexpr Vec3(const RealType x, const RealType y, const RealType z) noexcept : x(x), y(y), z(z) {}

		/**
		 * @brief Constructs a rectangular vector from a spherical SVec3.
		 *
		 * @param svec A spherical vector (SVec3) to convert.
		 */
		explicit Vec3(const SVec3& svec) noexcept;

		/**
		 * @brief Addition assignment operator for Vec3.
		 *
		 * @param b The vector to add.
		 * @return A reference to the updated Vec3.
		 */
		Vec3& operator+=(const Vec3& b) noexcept;

		/**
		 * @brief Subtraction assignment operator for Vec3.
		 *
		 * @param b The vector to subtract.
		 * @return A reference to the updated Vec3.
		 */
		Vec3& operator-=(const Vec3& b) noexcept;

		/**
		 * @brief Multiplication assignment operator for Vec3.
		 *
		 * @param b The vector to multiply by.
		 * @return A reference to the updated Vec3.
		 */
		Vec3& operator*=(const Vec3& b) noexcept;

		/**
		 * @brief Matrix multiplication assignment for Vec3.
		 *
		 * @param m The matrix to multiply by.
		 * @return A reference to the updated Vec3.
		 */
		Vec3& operator*=(const Matrix3& m) noexcept;

		/**
		 * @brief Scalar multiplication assignment for Vec3.
		 *
		 * @param b The scalar value.
		 * @return A reference to the updated Vec3.
		 */
		Vec3& operator*=(RealType b) noexcept;

		/**
		 * @brief Scalar division assignment for Vec3.
		 *
		 * @param b The scalar value.
		 * @return A reference to the updated Vec3.
		 */
		Vec3& operator/=(RealType b) noexcept;

		/**
		 * @brief Addition operator for Vec3.
		 *
		 * @param value The scalar value to add.
		 * @return The resulting vector.
		 */
		Vec3 operator+(const RealType value) const { return {x + value, y + value, z + value}; }

		Vec3 operator-() const { return {-x, -y, -z}; }

		/**
		 * @brief Calculates the length (magnitude) of the vector.
		 *
		 * @return The length of the vector.
		 */
		[[nodiscard]] RealType length() const noexcept { return std::sqrt(x * x + y * y + z * z); }
	};

	/**
	 * @brief Computes the dot product of two Vec3 vectors.
	 *
	 * @param a The first vector.
	 * @param b The second vector.
	 * @return The dot product of the vectors.
	 */
	inline RealType dotProduct(const Vec3& a, const Vec3& b) noexcept { return a.x * b.x + a.y * b.y + a.z * b.z; }

	// Cross product
	// Note: This function is not used in the codebase
	/*inline Vec3 crossProduct(const Vec3& a, const Vec3& b) noexcept
	{
		return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
	}*/

	/// Multiplies two Vec3 vectors component-wise
	inline Vec3 operator*(const Vec3& a, const Vec3& b) noexcept { return {a.x * b.x, a.y * b.y, a.z * b.z}; }

	/// Adds two Vec3 vectors component-wise
	inline Vec3 operator+(const Vec3& a, const Vec3& b) noexcept { return {a.x + b.x, a.y + b.y, a.z + b.z}; }

	/// Subtracts two Vec3 vectors component-wise
	inline Vec3 operator-(const Vec3& a, const Vec3& b) noexcept { return {a.x - b.x, a.y - b.y, a.z - b.z}; }

	/// Divides two Vec3 vectors component-wise
	inline Vec3 operator/(const Vec3& a, const Vec3& b) { return {a.x / b.x, a.y / b.y, a.z / b.z}; }

	/// Multiplies a Vec3 vector by a scalar value
	inline Vec3 operator*(const Vec3& a, const RealType b) noexcept { return {a.x * b, a.y * b, a.z * b}; }

	/// Divides a Vec3 vector by a scalar value
	inline Vec3 operator/(const Vec3& a, const RealType b) noexcept { return {a.x / b, a.y / b, a.z / b}; }

	/// Divides a scalar value by a Vec3 vector
	inline Vec3 operator/(const RealType a, const Vec3& b) noexcept { return {a / b.x, a / b.y, a / b.z}; }

	/// Adds two SVec3 vectors
	SVec3 operator+(const SVec3& a, const SVec3& b) noexcept;

	/// Subtracts two SVec3 vectors
	SVec3 operator-(const SVec3& a, const SVec3& b) noexcept;
}
