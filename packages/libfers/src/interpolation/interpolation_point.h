// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file interpolation_point.h
 * @brief Defines a structure to store interpolation point data for signal processing.
 */

#pragma once

namespace interp
{
	/**
	 * @struct InterpPoint
	 * @brief Stores data for an interpolation point.
	 */
	struct InterpPoint
	{
		RealType power{}; ///< Power level of the signal at the interpolation point.
		RealType time{}; ///< Time at which the interpolation point is recorded.
		RealType delay{}; ///< Delay associated with the signal at the interpolation point.
		RealType phase{}; ///< Phase of the signal at the interpolation point.
	};
}
