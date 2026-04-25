// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file channel_model.h
 * @brief Header for radar channel propagation and interaction models.
 *
 * This file contains the declarations for functions that model the radar channel,
 * including direct and reflected signal path calculations,  and
 * power scaling according to the radar range equation. These functions form the
 * core physics engine for determining the properties of a received signal based
 * on the positions and velocities of transmitters, receivers, and targets.
 */

#pragma once

#include <chrono>
#include <exception>
#include <memory>
#include <span>
#include <unordered_map>
#include <vector>

#include "core/config.h"
#include "core/sim_id.h"
#include "core/simulation_state.h"
#include "math/geometry_ops.h"

namespace core
{
	class World;
}
namespace timing
{
	class Timing;
}
namespace radar
{
	class Receiver;

	class Transmitter;

	class Target;
}

namespace serial
{
	class Response;
}

namespace fers_signal
{
	class RadarSignal;
}

namespace simulation
{
	struct CwPhaseNoiseBuffer
	{
		RealType start_time{};
		RealType dt{};
		std::vector<RealType> samples;

		[[nodiscard]] RealType sampleAt(RealType time) const noexcept;
	};

	struct CwPhaseNoiseLookup
	{
		RealType start_time{};
		RealType end_time{};
		RealType dt{};
		std::unordered_map<SimId, CwPhaseNoiseBuffer> buffers;

		[[nodiscard]] static CwPhaseNoiseLookup build(std::span<const std::shared_ptr<timing::Timing>> timings,
													  RealType start_time, RealType end_time);
		[[nodiscard]] RealType sample(const timing::Timing* timing, RealType time) const noexcept;
		[[nodiscard]] RealType phaseDifference(const timing::Timing* rx_timing, RealType rx_time,
											   const timing::Timing* tx_timing, RealType tx_time) const noexcept;
	};

	/**
	 * @struct ReResults
	 * @brief Stores the intermediate results of a radar equation calculation for a single time point.
	 */
	struct ReResults
	{
		RealType power; /**< Power scaling factor (dimensionless, relative to transmitted power). */
		RealType delay; /**< Signal propagation delay in seconds. */
		RealType phase; /**< Phase shift in radians due to propagation delay. */
	};

	/**
	 * @class RangeError
	 * @brief Exception thrown when a range calculation fails, typically due to objects being too close.
	 */
	class RangeError final : public std::exception
	{
	public:
		/**
		 * @brief Provides the error message for the exception.
		 * @return A C-style string describing the error.
		 */
		[[nodiscard]] const char* what() const noexcept override
		{
			return "Range error in radar equation calculations";
		}
	};

	/**
	 * @brief Solves the bistatic radar equation for a reflected path (Tx -> Tgt -> Rx).
	 *
	 * This function calculates the signal properties (power, delay, phase)
	 * for a signal traveling from a transmitter, reflecting off a target, and arriving at a receiver.
	 * It accounts for antenna gains, target RCS, and propagation loss.
	 *
	 * @param trans Pointer to the transmitter.
	 * @param recv Pointer to the receiver.
	 * @param targ Pointer to the target.
	 * @param time The time at which the pulse is transmitted.
	 * @param wave Pointer to the transmitted radar signal.
	 * @param results Output struct to store the calculation results.
	 * @throws RangeError If the target is too close to the transmitter or receiver.
	 */
	void solveRe(const radar::Transmitter* trans, const radar::Receiver* recv, const radar::Target* targ,
				 const std::chrono::duration<RealType>& time, const fers_signal::RadarSignal* wave, ReResults& results);

	/**
	 * @brief Solves the radar equation for a direct path (Tx -> Rx).
	 *
	 * This function calculates the signal properties for a direct line-of-sight signal
	 * traveling from a transmitter to a receiver.
	 *
	 * @param trans Pointer to the transmitter.
	 * @param recv Pointer to the receiver.
	 * @param time The time at which the pulse is transmitted.
	 * @param wave Pointer to the transmitted radar signal.
	 * @param results Output struct to store the calculation results.
	 * @throws RangeError If the transmitter and receiver are too close.
	 */
	void solveReDirect(const radar::Transmitter* trans, const radar::Receiver* recv,
					   const std::chrono::duration<RealType>& time, const fers_signal::RadarSignal* wave,
					   ReResults& results);

	/**
	 * @brief Calculates the complex envelope contribution for a direct propagation path (Tx -> Rx) at a specific time.
	 * This function is used for Continuous Wave (CW) simulations.
	 *
	 * @param trans The transmitter.
	 * @param recv The receiver.
	 * @param timeK The current simulation time.
	 * @return The complex I/Q sample contribution for this path.
	 */
	ComplexType calculateDirectPathContribution(const radar::Transmitter* trans, const radar::Receiver* recv,
												RealType timeK, const CwPhaseNoiseLookup* phase_noise_lookup = nullptr);

	ComplexType calculateStreamingDirectPathContribution(const core::ActiveStreamingSource& source,
														 const radar::Receiver* recv, RealType timeK,
														 const CwPhaseNoiseLookup* phase_noise_lookup = nullptr,
														 core::FmcwChirpBoundaryTracker* chirp_tracker = nullptr);

	/**
	 * @brief Calculates the complex envelope contribution for a reflected path (Tx -> Tgt -> Rx) at a specific time.
	 * This function is used for Continuous Wave (CW) simulations.
	 *
	 * @param trans The transmitter.
	 * @param recv The receiver.
	 * @param targ The target.
	 * @param timeK The current simulation time.
	 * @return The complex I/Q sample contribution for this path.
	 */
	ComplexType calculateReflectedPathContribution(const radar::Transmitter* trans, const radar::Receiver* recv,
												   const radar::Target* targ, RealType timeK,
												   const CwPhaseNoiseLookup* phase_noise_lookup = nullptr);

	ComplexType calculateStreamingReflectedPathContribution(const core::ActiveStreamingSource& source,
															const radar::Receiver* recv, const radar::Target* targ,
															RealType timeK,
															const CwPhaseNoiseLookup* phase_noise_lookup = nullptr,
															core::FmcwChirpBoundaryTracker* chirp_tracker = nullptr);

	/**
	 * @brief Creates a Response object by simulating a signal's interaction over its duration.
	 *
	 * This function iterates over the duration of a transmitted pulse, calling the
	 * appropriate channel model function (`solveRe` or `solveReDirect`) at discrete
	 * time steps to generate a series of `InterpPoint`s. These points capture the
	 * time-varying properties of the received signal and are collected into a `Response` object.
	 *
	 * @param trans Pointer to the transmitter.
	 * @param recv Pointer to the receiver.
	 * @param signal Pointer to the transmitted pulse signal.
	 * @param startTime The absolute simulation time when the pulse transmission starts.
	 * @param targ Optional pointer to a target. If null, a direct path is simulated.
	 * @return A unique pointer to the generated Response object.
	 * @throws RangeError If the channel model reports an invalid geometry.
	 * @throws std::runtime_error If the simulation parameters result in zero time steps.
	 */
	std::unique_ptr<serial::Response> calculateResponse(const radar::Transmitter* trans, const radar::Receiver* recv,
														const fers_signal::RadarSignal* signal, RealType startTime,
														const radar::Target* targ = nullptr);

	/**
	 * @enum LinkType
	 * @brief Categorizes the visual link for rendering.
	 */
	enum class LinkType
	{
		Monostatic, ///< Combined Tx/Rx path
		BistaticTxTgt, ///< Illuminator path
		BistaticTgtRx, ///< Scattered path
		DirectTxRx ///< Interference path
	};

	/**
	 * @enum LinkQuality
	 * @brief Describes the radiometric quality of the link.
	 */
	enum class LinkQuality
	{
		Strong, ///< SNR > 0 dB
		Weak ///< SNR < 0 dB (Geometric line of sight, but below noise floor)
	};

	/**
	 * @struct PreviewLink
	 * @brief A calculated link segment for 3D visualization.
	 */
	struct PreviewLink
	{
		LinkType type;
		LinkQuality quality;
		std::string label;
		SimId source_id; // The start of this specific link segment
		SimId dest_id; // The end of this specific link segment
		SimId origin_id; // The original source of energy (Transmitter)
		double rcs{-1.0}; // RCS in m^2 for this path; -1.0 if not applicable
		double actual_power_dbm{-999.0}; // Received power with actual RCS applied; -999 if not applicable
	};

	/**
	 * @brief Calculates all visual links for the current world state at a specific time.
	 *
	 * This function utilizes the core radar equation helpers to determine visibility,
	 * power levels, and SNR for all Tx/Rx/Target combinations. It is lightweight
	 * and does not update simulation state.
	 *
	 * @param world The simulation world containing radar components.
	 * @param time The time at which to calculate geometry.
	 * @return A vector of renderable links.
	 */
	std::vector<PreviewLink> calculatePreviewLinks(const core::World& world, RealType time);
}
