// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file path_utils.h
 * @brief Utility functions for path interpolation and exception handling.
 *
 * The cubic interpolation functions are based on methods described in "Numerical Recipes
 * in C, Second Edition" by Press et al., but the code here is distinct from the original.
 */

#pragma once

#include <algorithm>
#include <concepts>
#include <stdexcept>
#include <vector>

#include "core/config.h"

namespace math
{
	/**
	 * @class PathException
	 * @brief Exception class for handling path-related errors.
	 */
	class PathException final : public std::runtime_error
	{
	public:
		/**
		 * @brief Constructor for PathException.
		 *
		 * @param description A detailed description of the error.
		 */
		explicit PathException(const std::string& description) :
			std::runtime_error("Error While Executing Path Code: " + description)
		{
		}
	};
}

/**
 * @brief Concept for types that can be interpolated.
 *
 * @tparam T The type to be checked for interpolation.
 * @param a The first value to be interpolated.
 * @param b The second value to be interpolated.
 * @param t The interpolation factor.
 */
template <typename T>
concept Interpolatable = requires(T a, T b, RealType t) {
	{ a - b } -> std::same_as<T>;
	{ a * t } -> std::same_as<T>;
	{ a + b } -> std::same_as<T>;
	{ a.t } -> std::convertible_to<RealType>;
};

/**
 * @brief Interpolates a static position from a list of coordinates.
 *
 * @tparam T The type of the coordinate, which must satisfy the Interpolatable concept.
 * @param coord The output coordinate to be set.
 * @param coords A vector of coordinates from which the first will be selected.
 * @throws PathException if the list of coordinates is empty.
 */
template <Interpolatable T>
void getPositionStatic(T& coord, const std::vector<T>& coords)
{
	if (coords.empty())
	{
		throw math::PathException("coord list empty during GetPositionStatic");
	}
	coord = coords[0];
}

/**
 * @brief Performs linear interpolation between coordinate points.
 *
 * @tparam T The type of the coordinate, which must satisfy the Interpolatable concept.
 * @param t The interpolation factor (usually time) to determine the position.
 * @param coord The output coordinate that will be interpolated.
 * @param coords A vector of coordinates to interpolate between.
 * @throws PathException if the list of coordinates is empty.
 */
template <Interpolatable T>
void getPositionLinear(RealType t, T& coord, const std::vector<T>& coords)
{
	if (coords.empty())
	{
		throw math::PathException("coord list empty during GetPositionLinear");
	}

	auto xrp = std::ranges::upper_bound(coords, t, {}, &T::t);

	if (xrp == coords.begin())
	{
		coord = *xrp;
	}
	else if (xrp == coords.end())
	{
		coord = *(xrp - 1);
	}
	else
	{
		auto xri = std::distance(coords.begin(), xrp);
		auto xli = xri - 1;

		const RealType iw = coords[xri].t - coords[xli].t;
		const RealType rw = (coords[xri].t - t) / iw;
		const RealType lw = 1 - rw;

		coord = coords[xri] * lw + coords[xli] * rw;
	}
	coord.t = t;
}

/**
 * @brief Performs cubic spline interpolation between coordinate points.
 * The method used for calculating the spline is from "Numerical Recipes in C."
 *
 * @tparam T The type of the coordinate, which must satisfy the Interpolatable concept.
 * @param t The interpolation factor (usually time) to determine the position.
 * @param coord The output coordinate that will be interpolated.
 * @param coords A vector of coordinates to interpolate between.
 * @param dd A vector of second derivatives used in the cubic interpolation.
 * @throws PathException if the list of coordinates is empty.
 */
template <Interpolatable T>
void getPositionCubic(RealType t, T& coord, const std::vector<T>& coords, const std::vector<T>& dd)
{
	if (coords.empty())
	{
		throw math::PathException("coord list empty during GetPositionCubic");
	}

	auto xrp = std::ranges::upper_bound(coords, t, {}, &T::t);

	if (xrp == coords.begin())
	{
		coord = *xrp;
	}
	else if (xrp == coords.end())
	{
		coord = *(xrp - 1);
	}
	else
	{
		auto xri = std::distance(coords.begin(), xrp);
		auto xli = xri - 1;

		const RealType xrd = coords[xri].t - t;
		const RealType xld = t - coords[xli].t;
		const RealType iw = coords[xri].t - coords[xli].t;
		const RealType iws = iw * iw / 6.0;
		const RealType a = xrd / iw;
		const RealType b = xld / iw;
		const RealType c = (a * a * a - a) * iws;
		const RealType d = (b * b * b - b) * iws;

		coord = coords[xli] * a + coords[xri] * b + dd[xli] * c + dd[xri] * d;
	}
	coord.t = t;
}

/**
 * @brief Finalizes cubic spline interpolation by calculating second derivatives.
 *
 * @tparam T The type of the coordinate, which must satisfy the Interpolatable concept.
 * @param coords A vector of coordinates for which second derivatives will be calculated.
 * @param dd The output vector that will store the calculated second derivatives.
 * @throws PathException if there are fewer than two points for interpolation.
 */
template <Interpolatable T>
void finalizeCubic(const std::vector<T>& coords, std::vector<T>& dd)
{
	const int size = static_cast<int>(coords.size());
	if (size < 2)
	{
		throw math::PathException("Not enough points for cubic interpolation");
	}

	std::vector<T> tmp(size);
	dd.resize(size);

	dd.front() = 0;
	dd.back() = 0;

	for (int i = 1; i < size - 1; ++i)
	{
		const T yrd = coords[i + 1] - coords[i];
		const T yld = coords[i] - coords[i - 1];
		const RealType xrd = coords[i + 1].t - coords[i].t;
		const RealType xld = coords[i].t - coords[i - 1].t;
		const RealType iw = coords[i + 1].t - coords[i - 1].t;

		if (iw <= EPSILON)
		{
			dd[i] = 0.0;
			tmp[i] = 0.0;
			continue;
		}

		const T yrd_xrd = (xrd <= EPSILON) ? (yrd * 0.0) : (yrd / xrd);
		const T yld_xld = (xld <= EPSILON) ? (yld * 0.0) : (yld / xld);

		const RealType si = xld / iw;
		const T p = dd[i - 1] * si + 2.0;
		dd[i] = (si - 1.0) / p;
		tmp[i] = ((yrd_xrd - yld_xld) * 6.0 / iw - tmp[i - 1] * si) / p;
	}

	for (int i = size - 2; i >= 0; --i)
	{
		dd[i] = dd[i] * dd[i + 1] + tmp[i];
	}
}
