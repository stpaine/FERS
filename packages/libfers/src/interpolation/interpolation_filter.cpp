// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file interpolation_filter.cpp
 * @brief Implementation of the InterpFilter class.
 */

#include "interpolation_filter.h"

#include <stdexcept>

#include "core/logging.h"
#include "core/parameters.h"

using logging::Level;

namespace
{
	/**
	 * @brief Computes the modified Bessel function of the first kind for x.
	 *
	 * @param x The input value for which the Bessel function is computed.
	 * @return The computed Bessel function value, or an error message if computation fails.
	 */
	std::expected<RealType, std::string> besselI0(const RealType x)
	{
		// Use the polynomial approximation from section 9.8 of
		// "Handbook of Mathematical Functions" by Abramowitz and Stegun
		if (x < 0.0)
		{
			return std::unexpected("Modified Bessel approximation only valid for x > 0");
		}
		if (RealType t = x / 3.75; t <= 1.0)
		{
			t *= t;
			return 1.0 +
				t * (3.5156229 + t * (3.0899424 + t * (1.2067492 + t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))));
		}
		else
		{
			const RealType i0 = 0.39894228 +
				t *
					(0.01328592 +
					 t *
						 (0.00225319 +
						  t *
							  (-0.00157565 +
							   t *
								   (0.00916281 +
									t * (-0.02057706 + t * (0.02635537 + t * (-0.01647633 + t * 0.00392377)))))));
			return i0 * std::exp(x) / std::sqrt(x);
		}
	}
}

namespace interp
{
	InterpFilter& InterpFilter::getInstance() noexcept
	{
		static InterpFilter instance;
		return instance;
	}

	std::expected<RealType, std::string> InterpFilter::kaiserWinCompute(const RealType x) const noexcept
	{
		if (x < 0 || x > _alpha * 2)
		{
			return 0;
		}
		if (auto bessel = besselI0(_beta * std::sqrt(1 - std::pow((x - _alpha) / _alpha, 2))); bessel)
		{
			return *bessel / _bessel_beta;
		}
		else
		{
			return std::unexpected(bessel.error());
		}
	}

	std::expected<RealType, std::string> InterpFilter::interpFilter(const RealType x) const noexcept
	{
		if (auto kaiser = kaiserWinCompute(x + _alpha); kaiser)
		{
			return *kaiser * sinc(x);
		}
		else
		{
			return std::unexpected(kaiser.error());
		}
	}

	InterpFilter::InterpFilter()
	{
		_length = static_cast<int>(params::renderFilterLength());
		_table_filters = 1000;
		_filter_table = std::vector<RealType>(_table_filters * _length);

		_alpha = std::floor(params::renderFilterLength() / 2.0);
		if (auto bessel = besselI0(_beta); bessel)
		{
			_bessel_beta = *bessel;
		}
		else
		{
			LOG(Level::FATAL, "Bessel function calculation failed: {}", bessel.error());
			throw std::runtime_error("Bessel function calculation failed");
		}

		const int hfilt = _table_filters / 2;

		LOG(Level::DEBUG, "Building table of {} filters", _table_filters);

		for (int i = -hfilt; i < hfilt; ++i)
		{
			const RealType delay = i / static_cast<RealType>(hfilt);
			for (int j = static_cast<int>(-_alpha); j < _alpha; ++j)
			{
				if (auto interp = interpFilter(j - delay); interp)
				{
					_filter_table[(i + hfilt) * _length + j + static_cast<int>(_alpha)] = *interp;
				}
				else
				{
					LOG(Level::FATAL, "Interpolation filter calculation failed: {}", interp.error());
					throw std::runtime_error("Interpolation filter calculation failed");
				}
			}
		}

		LOG(Level::DEBUG, "Filter table complete");
	}

	std::span<const RealType> InterpFilter::getFilter(RealType delay) const
	{
		if (delay < -1 || delay > 1)
		{
			LOG(Level::FATAL, "Requested delay filter value out of range: {}", delay);
			throw std::runtime_error("Requested delay filter value out of range");
		}

		const auto filt = static_cast<unsigned>((delay + 1) * (_table_filters / 2.0));
		return std::span{&_filter_table[filt * _length], static_cast<size_t>(_length)};
	}
}
