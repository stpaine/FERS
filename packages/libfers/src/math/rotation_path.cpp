// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file rotation_path.cpp
 * @brief Implementation of the RotationPath class.
 */

#include "rotation_path.h"

#include <algorithm>
#include <cmath>

#include "coord.h"
#include "geometry_ops.h"
#include "path_utils.h"

namespace math
{
	void RotationPath::addCoord(const RotationCoord& coord) noexcept
	{
		const auto iter = std::lower_bound(_coords.begin(), _coords.end(), coord);
		_coords.insert(iter, coord);
		_final = false;
	}

	SVec3 RotationPath::getPosition(const RealType t) const
	{
		if (!_final)
		{
			throw PathException("Finalize not called before getPosition in RotationPath.");
		}
		RotationCoord coord{};

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
		case InterpType::INTERP_CONSTANT:
			coord.azimuth = std::fmod(t * _rate.azimuth + _start.azimuth, 2 * PI);
			coord.elevation = std::fmod(t * _rate.elevation + _start.elevation, 2 * PI);
			break;
		default:
			throw PathException("Unknown interpolation type.");
		}

		return {1, coord.azimuth, coord.elevation};
	}

	void RotationPath::finalize()
	{
		if (!_final)
		{
			if (_type == InterpType::INTERP_CUBIC)
			{
				finalizeCubic(_coords, _dd);
			}
			_final = true;
		}
	}

	void RotationPath::setInterp(const InterpType setinterp) noexcept
	{
		_type = setinterp;
		_final = false;
	}

	void RotationPath::setConstantRate(const RotationCoord& setstart, const RotationCoord& setrate) noexcept
	{
		_start = setstart;
		_rate = setrate;
		_type = InterpType::INTERP_CONSTANT;
		_final = true;
	}
}
