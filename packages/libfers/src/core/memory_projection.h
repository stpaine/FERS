// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file memory_projection.h
 * @brief Startup memory and output-size projection helpers for simulations.
 */

#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "core/config.h"

namespace timing
{
	class Timing;
}

namespace core
{
	class World;

	/**
	 * @brief Describes a projected byte count and whether it saturated during arithmetic.
	 */
	struct ByteProjection
	{
		std::uint64_t bytes = 0; ///< Projected byte count, clamped to `uint64_t` max on overflow.
		bool overflowed = false; ///< True if any arithmetic used to produce `bytes` overflowed.
	};

	/**
	 * @brief Captures startup memory and rendered-output projections for a simulation.
	 *
	 * The projection separates memory that is expected to be allocated by the simulation
	 * run itself from the resident baseline already present at startup. Byte totals are
	 * estimates intended for logging, diagnostics, and API reporting.
	 */
	struct SimulationMemoryProjection
	{
		RealType duration_seconds = 0.0; ///< Simulated duration covered by the projection.
		RealType simulation_sample_rate_hz = 0.0; ///< Oversampled simulation rate used for sample counts.
		unsigned oversample_ratio = 1; ///< Oversampling ratio applied to the configured output rate.

		std::uint64_t streaming_sample_count = 0; ///< Samples held by each full-duration streaming receiver.
		std::uint64_t phase_noise_sample_count = 0; ///< Samples held by each enabled phase-noise lookup.
		std::uint64_t phase_noise_timing_count = 0; ///< Unique streaming timing sources considered.
		std::uint64_t enabled_phase_noise_timing_count = 0; ///< Unique timing sources with phase noise enabled.
		std::uint64_t streaming_receiver_count = 0; ///< Receivers that keep full-duration streaming IQ buffers.
		std::uint64_t pulsed_receiver_count = 0; ///< Receivers that render finite pulsed receive windows.
		std::uint64_t pulsed_window_count = 0; ///< Projected count of pulsed receive windows.
		std::uint64_t rendered_hdf5_sample_count = 0; ///< Projected IQ samples written to HDF5 datasets.

		ByteProjection phase_noise_lookup; ///< Projected memory for all enabled timing-source lookup tables.
		ByteProjection streaming_iq_buffers; ///< Projected full-duration streaming receiver IQ memory.
		ByteProjection allocated_streaming_iq_buffers; ///< Streaming IQ capacity already allocated at startup.
		ByteProjection rendered_hdf5_payload; ///< Projected raw I/Q dataset payload bytes written to HDF5.
		std::optional<ByteProjection> resident_baseline; ///< Startup RSS excluding allocated streaming IQ buffers.
		std::optional<ByteProjection> projected_total_footprint; ///< Projected total resident footprint at peak.
		std::optional<std::uint64_t> current_resident_set; ///< Process RSS at projection time, when available.
	};

	/**
	 * @brief Collects unique timing sources used by CW/FMCW transmitters and receivers.
	 * @param world The simulation world to inspect.
	 * @return A vector of unique timing sources that may need phase-noise lookup tables.
	 */
	[[nodiscard]] std::vector<std::shared_ptr<timing::Timing>> collectCwPhaseNoiseTimings(const World& world);

	/**
	 * @brief Projects startup memory and rendered-output sizes for a simulation world.
	 * @param world The simulation world to inspect.
	 * @return A populated memory projection for the current simulation parameters.
	 */
	[[nodiscard]] SimulationMemoryProjection projectSimulationMemory(const World& world);

	/**
	 * @brief Formats a byte count using binary units.
	 * @param bytes The byte count to format.
	 * @return A human-readable string such as `512 B` or `1.50 MiB`.
	 */
	[[nodiscard]] std::string formatByteSize(std::uint64_t bytes);

	/**
	 * @brief Serializes a simulation memory projection as JSON.
	 * @param projection The projection to serialize.
	 * @return A formatted JSON string suitable for API responses and diagnostics.
	 */
	[[nodiscard]] std::string memoryProjectionToJsonString(const SimulationMemoryProjection& projection);

	/**
	 * @brief Logs the projected simulation memory footprint for the provided world.
	 * @param world The simulation world to inspect and report.
	 */
	void logSimulationMemoryProjection(const World& world);
}
