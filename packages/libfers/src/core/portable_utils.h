// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file portable_utils.h
 * @brief Utility functions for mathematical and system operations.
 */

#pragma once

#include <cmath>
#include <thread>

#include "config.h"

namespace core
{
	/**
	 * @brief Computes the Bessel function of the first kind (order 1) for a given value.
	 *
	 * @param x The value for which the Bessel function is to be computed.
	 * @return The computed value of the Bessel function of the first kind (order 1).
	 */
	inline RealType besselJ1(const RealType x) noexcept { return j1(x); }

	/**
	 * @brief Detects the number of CPUs on the machine.
	 *
	 * @return The number of CPUs detected, or 1 if detection fails.
	 */
	inline unsigned countProcessors() noexcept
	{
		if (const unsigned hardware_threads = std::thread::hardware_concurrency(); hardware_threads > 0)
		{
			return hardware_threads;
		}
		LOG(logging::Level::ERROR, "Unable to get CPU count, assuming 1.");
		return 1;
	}
}
