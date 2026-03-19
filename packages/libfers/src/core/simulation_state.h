// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file simulation_state.h
 * @brief Defines the global state for the event-driven simulation engine.
 */

#pragma once

#include <vector>

#include "config.h"
#include "radar/transmitter.h"

namespace core
{
	/**
	 * @struct SimulationState
	 * @brief Holds the dynamic global state of the simulation.
	 *
	 * This includes the master simulation clock and lists of active objects
	 * that are needed for calculations across different event types.
	 */
	struct SimulationState
	{
		/// The master simulation clock, advanced by the event loop.
		RealType t_current = 0.0;

		/// A global list of all currently active continuous-wave transmitters.
		std::vector<radar::Transmitter*> active_cw_transmitters;
	};
}
