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

#include <cstddef>
#include <optional>
#include <vector>

#include "config.h"
#include "radar/transmitter.h"

namespace fers_signal
{
	class FmcwChirpSignal;
}

namespace core
{
	struct ActiveStreamingSource
	{
		const radar::Transmitter* transmitter = nullptr;
		RealType segment_start = 0.0;
		RealType segment_end = 0.0;

		// Cached for one TX_STREAMING_START/TX_STREAMING_END segment. The world owns
		// the waveform, so the raw FMCW pointer is stable while this source is active.
		RealType carrier_freq = 0.0;
		RealType amplitude = 0.0;
		bool is_fmcw = false;

		const fers_signal::FmcwChirpSignal* fmcw = nullptr;
		RealType chirp_duration = 0.0;
		RealType chirp_period = 0.0;
		RealType chirp_rate = 0.0;
		RealType start_freq_off = 0.0;
		RealType two_pi_f0 = 0.0;
		RealType pi_alpha = 0.0;
		std::optional<std::size_t> chirp_count;
	};

	[[nodiscard]] ActiveStreamingSource makeActiveSource(const radar::Transmitter* tx, RealType segment_start,
														 RealType segment_end);

	struct FmcwChirpBoundaryTracker
	{
		bool initialized = false;
		RealType t_n = 0.0;
		std::size_t n_current = 0;
	};

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

		/// A global list of all currently active streaming transmitters.
		std::vector<ActiveStreamingSource> active_streaming_transmitters;
	};
}
