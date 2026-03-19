// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file interpolation_filter.h
 * @brief Interpolation filter implementation using Kaiser windowing.
 */

#pragma once

#include <expected>
#include <span>
#include <vector>

#include "core/config.h"

namespace interp
{
	/**
	 * @class InterpFilter
	 * @brief Provides methods to generate interpolation filters using Kaiser windows.
	 */
	class InterpFilter
	{
	public:
		InterpFilter(const InterpFilter&) = delete;
		InterpFilter(InterpFilter&&) = delete;
		InterpFilter& operator=(const InterpFilter&) = delete;
		InterpFilter& operator=(InterpFilter&&) = delete;
		~InterpFilter() = default;

		/**
		 * @brief Computes the sinc function for a given input.
		 *
		 * @param x The input value for which the sinc function is computed.
		 * @return The computed sinc value.
		 */
		static constexpr RealType sinc(const RealType x) noexcept
		{
			return x == 0.0 ? 1.0 : std::sin(x * PI) / (x * PI);
		}

		/**
		 * @brief Computes the Kaiser window function for a given input.
		 *
		 * @param x The input value for the Kaiser window calculation.
		 * @return The computed window value, or an error message if computation fails.
		 */
		[[nodiscard]] std::expected<RealType, std::string> kaiserWinCompute(RealType x) const noexcept;

		/**
		 * @brief Computes the interpolation filter value for a given input.
		 *
		 * @param x The input value for which the interpolation filter is computed.
		 * @return The computed filter value, or an error message if computation fails.
		 */
		[[nodiscard]] std::expected<RealType, std::string> interpFilter(RealType x) const noexcept;

		/**
		 * @brief Retrieves a span of precomputed filter values for a given delay.
		 *
		 * @param delay The delay value for which the filter is retrieved.
		 * @return A span of precomputed filter values.
		 * @throws std::runtime_error If the delay value is out of range.
		 */
		[[nodiscard]] std::span<const RealType> getFilter(RealType delay) const;

		/**
		 * @brief Retrieves the singleton instance of the `InterpFilter` class.
		 *
		 * @return The singleton instance of the `InterpFilter` class.
		 */
		static InterpFilter& getInstance() noexcept;

	private:
		InterpFilter(); ///< Private constructor to prevent instantiation.

		RealType _alpha; ///< The alpha value for the Kaiser window.
		RealType _beta = 5; ///< The beta value for the Kaiser window.
		RealType _bessel_beta; ///< The Bessel function value for the Kaiser window.
		int _length; ///< The length of the filter.
		int _table_filters; ///< The number of filters in the table.
		std::vector<RealType> _filter_table; ///< The table of precomputed filters.
	};
}
