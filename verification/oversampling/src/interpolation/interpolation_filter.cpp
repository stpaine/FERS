// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// Upstream reference:
//   packages/libfers/src/interpolation/interpolation_filter.cpp

#include "interpolation_filter.h"

#include <cmath>
#include <stdexcept>

#include "core/logging.h"
#include "core/parameters.h"

using logging::Level;

namespace
{
	std::expected<RealType, std::string> besselI0(const RealType x)
	{
		if (x < 0.0)
		{
			return std::unexpected("Modified Bessel approximation only valid for x > 0");
		}
		RealType t = x / 3.75;
		if (t <= 1.0)
		{
			t *= t;
			return 1.0 +
				t * (3.5156229 +
					 t * (3.0899424 + t * (1.2067492 + t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))));
		}

		const RealType i0 =
			0.39894228 +
			t *
				(0.01328592 +
				 t *
					 (0.00225319 +
					  t *
						  (-0.00157565 +
						   t * (0.00916281 + t * (-0.02057706 + t * (0.02635537 + t * (-0.01647633 + t * 0.00392377)))))));
		return i0 * std::exp(x) / std::sqrt(x);
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
		auto bessel = besselI0(_beta * std::sqrt(1 - std::pow((x - _alpha) / _alpha, 2)));
		if (bessel)
		{
			return *bessel / _bessel_beta;
		}

		return std::unexpected(bessel.error());
	}

	std::expected<RealType, std::string> InterpFilter::interpFilter(const RealType x) const noexcept
	{
		auto kaiser = kaiserWinCompute(x + _alpha);
		if (kaiser)
		{
			return *kaiser * sinc(x);
		}

		return std::unexpected(kaiser.error());
	}

	InterpFilter::InterpFilter()
	{
		_length = static_cast<int>(params::renderFilterLength());
		_table_filters = 1000;
		_filter_table = std::vector<RealType>(static_cast<size_t>(_table_filters * _length));

		_alpha = std::floor(params::renderFilterLength() / 2.0);
		if (auto bessel = besselI0(_beta); bessel)
		{
			_bessel_beta = *bessel;
		}
		else
		{
			LOG(Level::FATAL, "Bessel function calculation failed: ", bessel.error());
			throw std::runtime_error("Bessel function calculation failed");
		}

		const int hfilt = _table_filters / 2;

		for (int i = -hfilt; i < hfilt; ++i)
		{
			const RealType delay = i / static_cast<RealType>(hfilt);
			for (int j = static_cast<int>(-_alpha); j < _alpha; ++j)
			{
				if (auto interp = interpFilter(j - delay); interp)
				{
					_filter_table[static_cast<size_t>((i + hfilt) * _length + j + static_cast<int>(_alpha))] = *interp;
				}
				else
				{
					LOG(Level::FATAL, "Interpolation filter calculation failed: ", interp.error());
					throw std::runtime_error("Interpolation filter calculation failed");
				}
			}
		}
	}

	std::span<const RealType> InterpFilter::getFilter(const RealType delay) const
	{
		if (delay < -1 || delay > 1)
		{
			LOG(Level::FATAL, "Requested delay filter value out of range: ", delay);
			throw std::runtime_error("Requested delay filter value out of range");
		}

		const auto filt = static_cast<unsigned>((delay + 1) * (_table_filters / 2.0));
		return std::span{&_filter_table[static_cast<size_t>(filt * _length)], static_cast<size_t>(_length)};
	}
}

