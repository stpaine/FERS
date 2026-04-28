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
#include <cstdint>
#include <optional>
#include <vector>

#include "config.h"
#include "radar/transmitter.h"

namespace fers_signal
{
	class FmcwChirpSignal;
	class FmcwTriangleSignal;
}

namespace core
{
	/// Streaming waveform shape cached for a currently active source.
	enum class StreamingWaveformKind
	{
		Cw,
		FmcwLinear,
		FmcwTriangle
	};

	/// Cached description of an active streaming transmitter segment.
	struct ActiveStreamingSource
	{
		const radar::Transmitter* transmitter = nullptr; ///< Transmitter active during this segment.
		RealType segment_start = 0.0; ///< Segment start time in seconds.
		RealType segment_end = 0.0; ///< Segment end time in seconds.

		RealType carrier_freq = 0.0; ///< Cached carrier frequency in hertz.
		RealType amplitude = 0.0; ///< Cached emitted signal amplitude.
		StreamingWaveformKind kind = StreamingWaveformKind::Cw; ///< Cached streaming waveform shape.
		bool is_fmcw = false; ///< Compatibility flag for any FMCW source.

		const fers_signal::FmcwChirpSignal* fmcw = nullptr; ///< Stable pointer to the linear FMCW waveform, if any.
		const fers_signal::FmcwTriangleSignal* triangle = nullptr; ///< Stable pointer to the triangle waveform, if any.
		RealType chirp_duration = 0.0; ///< Cached FMCW chirp duration in seconds.
		RealType chirp_period = 0.0; ///< Cached FMCW chirp period in seconds.
		RealType chirp_rate = 0.0; ///< Cached FMCW chirp rate in hertz per second.
		RealType signed_chirp_rate = 0.0; ///< Cached signed FMCW chirp rate in hertz per second.
		RealType start_freq_off = 0.0; ///< Cached FMCW start frequency offset in hertz.
		RealType two_pi_f0 = 0.0; ///< Cached two-pi carrier angular frequency factor.
		RealType s_pi_alpha = 0.0; ///< Cached signed pi-scaled FMCW chirp-rate factor.
		std::optional<std::size_t> chirp_count; ///< Optional finite chirp count for the segment.

		RealType two_pi_f0_plus_B = 0.0; ///< Triangle down-leg linear coefficient.
		RealType pi_alpha = 0.0; ///< Triangle up-leg quadratic coefficient.
		RealType neg_pi_alpha = 0.0; ///< Triangle down-leg quadratic coefficient.
		RealType mod_phi_up = 0.0; ///< Triangle leg phase increment modulo 2*pi.
		RealType mod_phi_tri = 0.0; ///< Triangle period phase increment modulo 2*pi.
		RealType triangle_period = 0.0; ///< Cached full triangle period in seconds.
		std::optional<std::size_t> triangle_count; ///< Optional finite triangle count for the segment.
	};

	/// Builds an active-source cache from a streaming transmitter and segment bounds.
	[[nodiscard]] ActiveStreamingSource makeActiveSource(const radar::Transmitter* tx, RealType segment_start,
														 RealType segment_end);

	/// Returns the first FMCW chirp start inside the absolute interval, if one exists.
	[[nodiscard]] std::optional<RealType> firstFmcwChirpStart(const ActiveStreamingSource& source,
															  RealType active_start, RealType active_end);

	/// Counts FMCW chirps that start inside the absolute interval.
	[[nodiscard]] std::uint64_t countFmcwChirpStarts(const ActiveStreamingSource& source, RealType active_start,
													 RealType active_end);

	/// Returns the first FMCW triangle start inside the absolute interval, if one exists.
	[[nodiscard]] std::optional<RealType> firstFmcwTriangleStart(const ActiveStreamingSource& source,
																 RealType active_start, RealType active_end);

	/// Counts FMCW triangles that start inside the absolute interval.
	[[nodiscard]] std::uint64_t countFmcwTriangleStarts(const ActiveStreamingSource& source, RealType active_start,
														RealType active_end);

	/// Tracks the current FMCW chirp boundary for a streaming path.
	struct FmcwChirpBoundaryTracker
	{
		bool initialized = false; ///< True after the tracker has been initialized for a path.
		RealType t_n = 0.0; ///< Current chirp boundary time in seconds.
		std::size_t n_current = 0; ///< Current zero-based chirp index.
		bool triangle_initialized = false; ///< True after the triangle tracker has been initialized.
		std::size_t triangle_leg = 0; ///< Current triangle leg index: 0 for up-leg, 1 for down-leg.
		std::size_t triangle_index = 0; ///< Current zero-based triangle index.
		RealType triangle_t_leg = 0.0; ///< Current triangle leg boundary time in seconds.
		RealType triangle_phi_base = 0.0; ///< Current modular triangle base phase in radians.
	};

	/// Tracks the current FMCW triangle leg boundary for a streaming path.
	struct FmcwTriangleBoundaryTracker
	{
		bool initialized = false; ///< True after the tracker has been initialized for a path.
		std::size_t m_current = 0; ///< Current leg index: 0 for up-leg, 1 for down-leg.
		std::size_t M_current = 0; ///< Current zero-based triangle index.
		RealType t_leg = 0.0; ///< Current leg boundary time in seconds.
		RealType phi_base = 0.0; ///< Current modular base phase in radians.
	};

	/// Tracks the current streaming waveform boundary for a path.
	struct StreamingWaveformTracker
	{
		FmcwChirpBoundaryTracker linear;
		FmcwTriangleBoundaryTracker triangle;
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
