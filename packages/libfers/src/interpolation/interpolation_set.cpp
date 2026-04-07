// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2007-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file interpolation_set.cpp
 * @brief Implementation file for interpolation of sets of data.
 */

#include "interpolation_set.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <ranges>
#include <stdexcept>

namespace interp
{
	template <RealConcept T>
	std::optional<T> InterpSetData::value(T x) const noexcept
	{
		if (_data.empty())
		{
			return std::nullopt;
		}

		// TODO: RealConcept is broader than the implementation really supports?
		// std::is_arithmetic_v<T> admits integral types
		// That may be intended, but for value(T x) it means interpolation with
		// integer queries returns std::optional<int>, and the final interpolated
		// double gets truncated by:
		// return static_cast<T>(...)
		// If the class is conceptually for real-valued interpolation,
		// std::floating_point<T> would be a better constraint
		const RealType x_value = static_cast<RealType>(x);
		const auto iter = _data.lower_bound(x_value);

		if (iter == _data.begin())
		{
			return static_cast<T>(iter->second);
		}
		if (iter == _data.end())
		{
			const auto prev = std::prev(iter);
			return static_cast<T>(prev->second);
		}
		if (iter->first == x_value)
		{
			return static_cast<T>(iter->second);
		}

		auto prev = std::prev(iter);
		const auto [x1, y1] = *prev;
		const auto [x2, y2] = *iter;

		return static_cast<T>(y2 * (x_value - x1) / (x2 - x1) + y1 * (x2 - x_value) / (x2 - x1));
	}

	double InterpSetData::max() const noexcept
	{
		auto values = _data | std::views::values;

		const auto max_element =
			std::ranges::max_element(values, [](const double a, const double b) { return std::abs(a) < std::abs(b); });

		return max_element != values.end() ? std::abs(*max_element) : 0.0;
	}

	template <RealConcept T>
	void InterpSetData::divide(T a)
	{
		if (a == 0)
		{
			throw std::invalid_argument("Division by zero is not allowed.");
		}

		std::ranges::for_each(_data | std::views::values, [a](auto& value) { value /= static_cast<double>(a); });
	}

	// Explicit instantiations for double and float
	template std::optional<double> InterpSetData::value<double>(double) const noexcept;

	template void InterpSetData::divide<double>(double);

	template std::optional<float> InterpSetData::value<float>(float) const noexcept;

	template void InterpSetData::divide<float>(float);
}
