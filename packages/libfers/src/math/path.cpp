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
				size_t idx = std::distance(_coords.begin(), xrp);

				// Clamp to valid segments
				if (idx == 0)
					idx = 1;
				if (idx >= _coords.size())
					idx = _coords.size() - 1;

				if (idx < 1)
					return {0, 0, 0}; // Should not happen if size >= 1 and clamp works

				const auto& p1 = _coords[idx - 1];
				const auto& p2 = _coords[idx];
				const RealType dt = p2.t - p1.t;

				if (dt <= EPSILON)
					return {0, 0, 0};

				return (p2.pos - p1.pos) / dt;
			}

		case InterpType::INTERP_CUBIC:
			{
				auto xrp = std::ranges::upper_bound(_coords, t, {}, &Coord::t);
				size_t xri;
				if (xrp == _coords.begin())
					xri = 1;
				else if (xrp == _coords.end())
					xri = _coords.size() - 1;
				else
					xri = std::distance(_coords.begin(), xrp);

				if (xri < 1 || xri >= _coords.size())
					return {0, 0, 0};

				size_t xli = xri - 1;

				const RealType h = _coords[xri].t - _coords[xli].t;
				if (h <= EPSILON)
					return {0, 0, 0};

				const RealType a = (_coords[xri].t - t) / h;
				const RealType b = (t - _coords[xli].t) / h;

				// Derivative coefficients
				// da/dt = -1/h
				// db/dt = 1/h
				// dc/dt = -h/6 * (3a^2 - 1)
				// dd/dt = h/6 * (3b^2 - 1)

				const RealType da = -1.0 / h;
				const RealType db = 1.0 / h;
				const RealType dc = -h / 6.0 * (3.0 * a * a - 1.0);
				const RealType dd_coeff = h / 6.0 * (3.0 * b * b - 1.0);

				return _coords[xli].pos * da + _coords[xri].pos * db + _dd[xli].pos * dc + _dd[xri].pos * dd_coeff;
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
