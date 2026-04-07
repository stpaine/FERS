// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file path.cpp
 * @brief Implementation of the Path class.
 */

#include "path.h"

#include <algorithm>
#include <cstddef>

#include "coord.h"
#include "core/logging.h"
#include "geometry_ops.h"
#include "path_utils.h"

using logging::Level;

namespace math
{
	void Path::addCoord(const Coord& coord) noexcept
	{
		auto comp = [](const Coord& a, const Coord& b) { return a.t < b.t; };

		const auto iter = std::ranges::lower_bound(_coords, coord, comp);
		_coords.insert(iter, coord);
		_final = false;
	}

	Vec3 Path::getPosition(const RealType t) const
	{
		if (!_final)
		{
			LOG(Level::FATAL, "Finalize not called before GetPosition");
			throw PathException("Finalize not called before GetPosition");
		}

		Coord coord{};
		switch (_type)
		{
		case InterpType::INTERP_STATIC:
			getPositionStatic(coord, _coords);
			break;
		case InterpType::INTERP_LINEAR:
			getPositionLinear(t, coord, _coords);
			break;
		case InterpType::INTERP_CUBIC:
			getPositionCubic(t, coord, _coords, _dd);
			break;
		}
		return coord.pos;
	}

	Vec3 Path::getVelocity(const RealType t) const
	{
		if (!_final)
		{
			LOG(Level::FATAL, "Finalize not called before GetVelocity");
			throw PathException("Finalize not called before GetVelocity");
		}

		if (_coords.empty())
		{
			return {0, 0, 0};
		}

		switch (_type)
		{
		case InterpType::INTERP_STATIC:
			return {0, 0, 0};

		case InterpType::INTERP_LINEAR:
			{
				auto xrp = std::ranges::upper_bound(_coords, t, {}, &Coord::t);
				auto idx = std::distance(_coords.begin(), xrp);

				// Clamp to valid segments
				if (idx <= 0)
					idx = 1;
				if (static_cast<std::size_t>(idx) >= _coords.size())
					idx = static_cast<decltype(idx)>(_coords.size() - 1);

				if (idx < 1 || static_cast<std::size_t>(idx) >= _coords.size())
					return {0, 0, 0}; // Should not happen if size >= 1 and clamp works

				const auto right_idx = static_cast<std::size_t>(idx);
				const auto left_idx = right_idx - 1;

				const auto& p1 = _coords[left_idx];
				const auto& p2 = _coords[right_idx];
				const RealType dt = p2.t - p1.t;

				if (dt <= EPSILON)
					return {0, 0, 0};

				return (p2.pos - p1.pos) / dt;
			}

		case InterpType::INTERP_CUBIC:
			{
				auto xrp = std::ranges::upper_bound(_coords, t, {}, &Coord::t);
				std::ptrdiff_t xri;
				if (xrp == _coords.begin())
					xri = 1;
				else if (xrp == _coords.end())
					xri = static_cast<std::ptrdiff_t>(_coords.size() - 1);
				else
					xri = std::distance(_coords.begin(), xrp);

				if (xri < 1 || static_cast<std::size_t>(xri) >= _coords.size())
					return {0, 0, 0};

				const auto right_idx = static_cast<std::size_t>(xri);
				const auto left_idx = right_idx - 1;

				const RealType h = _coords[right_idx].t - _coords[left_idx].t;
				if (h <= EPSILON)
					return {0, 0, 0};

				const RealType a = (_coords[right_idx].t - t) / h;
				const RealType b = (t - _coords[left_idx].t) / h;

				// Derivative coefficients
				// da/dt = -1/h
				// db/dt = 1/h
				// dc/dt = -h/6 * (3a^2 - 1)
				// dd/dt = h/6 * (3b^2 - 1)

				const RealType da = -1.0 / h;
				const RealType db = 1.0 / h;
				const RealType dc = -h / 6.0 * (3.0 * a * a - 1.0);
				const RealType dd_coeff = h / 6.0 * (3.0 * b * b - 1.0);

				return _coords[left_idx].pos * da + _coords[right_idx].pos * db + _dd[left_idx].pos * dc +
					_dd[right_idx].pos * dd_coeff;
			}
		}
		return {0, 0, 0};
	}

	void Path::finalize()
	{
		if (!_final)
		{
			switch (_type)
			{
			case InterpType::INTERP_STATIC:
			case InterpType::INTERP_LINEAR:
				break;
			case InterpType::INTERP_CUBIC:
				finalizeCubic<Coord>(_coords, _dd);
				break;
			}
			_final = true;
		}
	}

	void Path::setInterp(const InterpType settype) noexcept
	{
		_type = settype;
		_final = false;
	}
}
