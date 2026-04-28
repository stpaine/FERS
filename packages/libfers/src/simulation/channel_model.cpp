// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file channel_model.cpp
 * @brief Implementation of radar channel propagation and interaction models.
 *
 * This file provides the implementations for the functions that model the radar
 * channel, as declared in channel_model.h. It contains the core physics calculations
 * that determine signal properties based on geometry, velocity, and object characteristics.
 */

#include "channel_model.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string_view>
#include <unordered_map>

#include "core/logging.h"
#include "core/parameters.h"
#include "core/sim_id.h"
#include "core/world.h"
#include "interpolation/interpolation_point.h"
#include "math/geometry_ops.h"
#include "radar/radar_obj.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/response.h"
#include "signal/radar_signal.h"
#include "timing/timing.h"

using fers_signal::RadarSignal;
using logging::Level;
using math::SVec3;
using math::Vec3;
using radar::Receiver;
using radar::Target;
using radar::Transmitter;

namespace
{
	/// Initializes an FMCW chirp tracker from a retarded time.
	void initializeFmcwTracker(const core::ActiveStreamingSource& source, const RealType t_ret,
							   core::FmcwChirpBoundaryTracker& tracker)
	{
		tracker.initialized = true;
		if (t_ret >= source.segment_start)
		{
			const RealType time_since_segment_start = t_ret - source.segment_start;
			tracker.n_current = static_cast<std::size_t>(std::floor(time_since_segment_start / source.chirp_period));
			tracker.t_n = source.segment_start + static_cast<RealType>(tracker.n_current) * source.chirp_period;
			return;
		}

		tracker.n_current = 0;
		tracker.t_n = source.segment_start;
	}

	/**
	 * @struct LinkGeometry
	 * @brief Holds geometric properties of a path segment between two points.
	 */
	struct LinkGeometry
	{
		Vec3 u_vec; ///< Unit vector pointing from Source to Destination.
		RealType dist{}; ///< Distance between Source and Destination.
	};

	/**
	 * @brief Computes the geometry (distance and direction) between two points.
	 * @param p_from Starting position.
	 * @param p_to Ending position.
	 * @return LinkGeometry containing distance and unit vector.
	 * @throws RangeError if the distance is too small (<= EPSILON).
	 */
	LinkGeometry computeLink(const Vec3& p_from, const Vec3& p_to)
	{
		const Vec3 vec = p_to - p_from;
		const RealType dist = vec.length();

		if (dist <= EPSILON)
		{
			// LOG(Level::FATAL) is handled by the caller or generic exception handler if needed,
			// but for RangeError strictly we just throw here to keep it pure.
			// However, existing code logged FATAL before throwing.
			// We'll throw RangeError, and let callers decide if they want to log or return 0.
			throw simulation::RangeError();
		}

		return {vec / dist, dist};
	}

	/**
	 * @brief Calculates the antenna gain for a specific direction and time.
	 * @param radar The radar object (Transmitter or Receiver).
	 * @param direction_vec The unit vector pointing AWAY from the antenna towards the target/receiver.
	 * @param time The simulation time for rotation lookup.
	 * @param lambda The signal wavelength.
	 * @return The linear gain value.
	 */
	RealType computeAntennaGain(const radar::Radar* radar, const Vec3& direction_vec, RealType time, RealType lambda)
	{
		return radar->getGain(SVec3(direction_vec), radar->getRotation(time), lambda);
	}

	/**
	 * @brief Computes the power scaling factor for a direct path (Friis Transmission Equation).
	 * @param tx_gain Transmitter gain (linear).
	 * @param rx_gain Receiver gain (linear).
	 * @param lambda Wavelength (meters).
	 * @param dist Distance (meters).
	 * @param no_prop_loss If true, distance-based attenuation is ignored.
	 * @return The power scaling factor (Pr / Pt).
	 */
	RealType computeDirectPathPower(RealType tx_gain, RealType rx_gain, RealType lambda, RealType dist,
									bool no_prop_loss)
	{
		const RealType numerator = tx_gain * rx_gain * lambda * lambda;
		RealType denominator = 16.0 * PI * PI; // (4 * PI)^2

		if (!no_prop_loss)
		{
			denominator *= dist * dist;
		}

		return numerator / denominator;
	}

	/**
	 * @brief Computes the power scaling factor for a reflected path (Bistatic Radar Range Equation).
	 * @param tx_gain Transmitter gain (linear).
	 * @param rx_gain Receiver gain (linear).
	 * @param rcs Target Radar Cross Section (m^2).
	 * @param lambda Wavelength (meters).
	 * @param r_tx Distance from Transmitter to Target.
	 * @param r_rx Distance from Target to Receiver.
	 * @param no_prop_loss If true, distance-based attenuation is ignored.
	 * @return The power scaling factor (Pr / Pt).
	 */
	RealType computeReflectedPathPower(RealType tx_gain, RealType rx_gain, RealType rcs, RealType lambda, RealType r_tx,
									   RealType r_rx, bool no_prop_loss)
	{
		const RealType numerator = tx_gain * rx_gain * rcs * lambda * lambda;
		RealType denominator = 64.0 * PI * PI * PI; // (4 * PI)^3

		if (!no_prop_loss)
		{
			denominator *= r_tx * r_tx * r_rx * r_rx;
		}

		return numerator / denominator;
	}

	/**
	 * @brief Computes the non-coherent phase shift due to timing offsets.
	 *
	 * Used for CW simulation where LO effects are modeled analytically.
	 *
	 * @param tx The transmitter.
	 * @param rx The receiver.
	 * @param time The current simulation time.
	 * @return The phase shift in radians.
	 */
	RealType computeTimingPhase(const Transmitter* tx, const Receiver* rx, RealType time)
	{
		const auto tx_timing = tx->getTiming();
		const auto rx_timing = rx->getTiming();
		const RealType delta_f = tx_timing->getFreqOffset() - rx_timing->getFreqOffset();
		const RealType delta_phi = tx_timing->getPhaseOffset() - rx_timing->getPhaseOffset();
		return 2 * PI * delta_f * time + delta_phi;
	}

	/// Computes timing phase with an optional precomputed CW phase-noise lookup.
	RealType computeTimingPhase(const Transmitter* tx, const Receiver* rx, const RealType rx_time,
								const RealType tx_time, const simulation::CwPhaseNoiseLookup* phase_noise_lookup)
	{
		if (phase_noise_lookup == nullptr)
		{
			return computeTimingPhase(tx, rx, rx_time);
		}
		return phase_noise_lookup->phaseDifference(rx->getTiming().get(), rx_time, tx->getTiming().get(), tx_time);
	}

	/// Checks whether received power is above the thermal noise floor.
	bool isSignalStrong(RealType power_watts, RealType temp_kelvin)
	{
		// Use configured rate or default to 1Hz if unconfigured to prevent divide-by-zero or silly values
		const RealType bw = params::rate() > 0 ? params::rate() : 1.0;
		const RealType noise_floor = params::boltzmannK() * (temp_kelvin > 0 ? temp_kelvin : 290.0) * bw;
		return power_watts > noise_floor;
	}

	/**
	 * @brief Converts power in watts to decibels milliwatts (dBm).
	 *
	 * @param watts Power in watts.
	 * @return Power in dBm. Returns -999.0 dBm for non-positive input.
	 */
	RealType wattsToDbm(RealType watts)
	{
		if (watts <= 0)
		{
			return -999.0;
		}
		return 10.0 * std::log10(watts * 1000.0);
	}

	/**
	 * @brief Converts power in watts to decibels (dB).
	 *
	 * @param watts Power in watts.
	 * @return Power in decibels (dB). Returns -999.0 dB for non-positive input.
	 */
	RealType wattsToDb(RealType watts)
	{
		if (watts <= 0)
		{
			return -999.0;
		}
		return 10.0 * std::log10(watts);
	}

	/**
	 * @brief Checks if a component is active at the given time based on its schedule.
	 *
	 * @param schedule The component's operating schedule.
	 * @param time The current simulation time.
	 * @return true If the schedule is empty (implied always on) or if time is within a period.
	 * @return false If the schedule is populated but the time is outside all periods.
	 */
	bool isComponentActive(const std::vector<radar::SchedulePeriod>& schedule, RealType time)
	{
		if (schedule.empty())
		{
			return true;
		}
		for (const auto& period : schedule)
		{
			if (time >= period.start && time <= period.end)
			{
				return true;
			}
		}
		return false;
	}

	/// Builds a compatibility streaming-source cache for classic CW paths.
	core::ActiveStreamingSource makeClassicStreamingSource(const Transmitter* trans)
	{
		auto source = core::makeActiveSource(trans, params::startTime(), std::numeric_limits<RealType>::max());
		if (!source.is_fmcw)
		{
			source.segment_start = std::numeric_limits<RealType>::lowest();
		}
		return source;
	}

	/// Computes streaming phase for CW or FMCW sources at a receiver time.
	bool computeStreamingPhase(const core::ActiveStreamingSource& source, const RealType rx_time, const RealType tau,
							   core::FmcwChirpBoundaryTracker* const chirp_tracker, RealType& phase_out)
	{
		if (source.carrier_freq <= 0.0)
		{
			return false;
		}

		const RealType t_ret = rx_time - tau;
		if (t_ret < source.segment_start || t_ret >= source.segment_end)
		{
			return false;
		}

		if (source.is_fmcw)
		{
			if (chirp_tracker == nullptr)
			{
				const RealType time_since_segment_start = t_ret - source.segment_start;
				const auto chirp_index =
					static_cast<std::size_t>(std::floor(time_since_segment_start / source.chirp_period));
				if (source.chirp_count.has_value() && chirp_index >= *source.chirp_count)
				{
					return false;
				}
				const RealType chirp_time =
					time_since_segment_start - static_cast<RealType>(chirp_index) * source.chirp_period;
				if (chirp_time < 0.0 || chirp_time >= source.chirp_duration)
				{
					return false;
				}
				phase_out = -2.0 * PI * source.carrier_freq * tau + source.two_pi_f0 * chirp_time +
					source.s_pi_alpha * chirp_time * chirp_time;
				return true;
			}

			if (!chirp_tracker->initialized)
			{
				initializeFmcwTracker(source, t_ret, *chirp_tracker);
			}

			while (t_ret >= chirp_tracker->t_n + source.chirp_period)
			{
				chirp_tracker->t_n += source.chirp_period;
				++chirp_tracker->n_current;
			}

			if (source.chirp_count.has_value() && chirp_tracker->n_current >= *source.chirp_count)
			{
				return false;
			}

			const RealType u_ret = t_ret - chirp_tracker->t_n;
			// This lower-bound guard is required for causality: cold-path initialization pins the
			// tracker to the segment start when the retarded time precedes the first chirp, and the
			// advance loop can only move forward. Without the u_ret < 0 check, we'd evaluate the chirp
			// polynomial at negative local time and extrapolate a non-causal signal before turn-on.
			if (u_ret < 0.0 || u_ret >= source.chirp_duration)
			{
				return false;
			}

			phase_out =
				-2.0 * PI * source.carrier_freq * tau + source.two_pi_f0 * u_ret + source.s_pi_alpha * u_ret * u_ret;
			return true;
		}

		phase_out = -2.0 * PI * source.carrier_freq * tau;
		return true;
	}

	/// Returns average radiated power for preview visualization.
	RealType previewRadiatedPower(const RadarSignal* waveform)
	{
		if (waveform == nullptr)
		{
			return 0.0;
		}
		if (const auto* fmcw = waveform->getFmcwChirpSignal(); fmcw != nullptr)
		{
			return waveform->getPower() * (fmcw->getChirpDuration() / fmcw->getChirpPeriod());
		}
		return waveform->getPower();
	}

	/// Formats received or radiated power as a dBm preview label.
	std::string formatPreviewDbmLabel(const RealType watts, const std::string_view prefix = {})
	{
		return std::format("{}{:.1f} dBm", prefix, wattsToDbm(watts));
	}

	/// Formats target illumination density as a dBW-per-square-meter preview label.
	std::string formatPreviewDbwPerSquareMeterLabel(const RealType watts)
	{
		return std::format("{:.1f} dBW/m\u00B2", wattsToDb(watts));
	}
}

namespace simulation
{
	RealType CwPhaseNoiseBuffer::sampleAt(const RealType time) const noexcept
	{
		if (samples.empty())
		{
			return 0.0;
		}
		if ((time <= start_time) || (samples.size() == 1) || (dt <= 0.0))
		{
			return samples.front();
		}

		const RealType position = (time - start_time) / dt;
		const auto last_index = static_cast<RealType>(samples.size() - 1);
		if (position >= last_index)
		{
			return samples.back();
		}

		const auto lower_index = static_cast<std::size_t>(position);
		const RealType fraction = position - static_cast<RealType>(lower_index);
		return samples[lower_index] + fraction * (samples[lower_index + 1] - samples[lower_index]);
	}

	CwPhaseNoiseLookup CwPhaseNoiseLookup::build(const std::span<const std::shared_ptr<timing::Timing>> timings,
												 const RealType start_time, const RealType end_time)
	{
		CwPhaseNoiseLookup lookup{};
		lookup.start_time = start_time;
		lookup.end_time = std::max(start_time, end_time);

		const RealType sample_rate = params::rate() * params::oversampleRatio();
		lookup.dt = sample_rate > 0.0 ? (1.0 / sample_rate) : 1.0;

		const auto sample_count =
			static_cast<std::size_t>(std::ceil((lookup.end_time - lookup.start_time) / lookup.dt) + 1.0);
		constexpr std::size_t phase_noise_warning_threshold_bytes = 500ULL * 1024ULL * 1024ULL;
		const auto bytes_per_buffer = sample_count * sizeof(RealType);

		for (const auto& timing : timings)
		{
			if (!timing || lookup.buffers.contains(timing->getId()))
			{
				continue;
			}

			CwPhaseNoiseBuffer buffer{};
			buffer.start_time = lookup.start_time;
			buffer.dt = lookup.dt;
			if (timing->isEnabled())
			{
				auto timing_clone = timing->clone();
				if (lookup.start_time > 0.0)
				{
					const auto skip_count = static_cast<std::size_t>(std::llround(lookup.start_time / lookup.dt));
					timing_clone->skipSamples(skip_count);
				}
				// TODO: Replace whole-simulation CW lookup generation with chunked/streaming generation.
				if (bytes_per_buffer > phase_noise_warning_threshold_bytes)
				{
					LOG(Level::WARNING,
						"CW phase-noise lookup for timing '{}' allocates {} bytes; large scenarios need chunked "
						"streaming.",
						timing->getName(), bytes_per_buffer);
				}
				buffer.samples.resize(sample_count);
				std::ranges::generate(buffer.samples, [&] { return timing_clone->getNextSample(); });
			}

			lookup.buffers.emplace(timing->getId(), std::move(buffer));
		}

		return lookup;
	}

	RealType CwPhaseNoiseLookup::sample(const timing::Timing* const timing, const RealType time) const noexcept
	{
		if (timing == nullptr)
		{
			return 0.0;
		}
		const auto it = buffers.find(timing->getId());
		if (it == buffers.end())
		{
			return 0.0;
		}
		return it->second.sampleAt(time);
	}

	RealType CwPhaseNoiseLookup::phaseDifference(const timing::Timing* const rx_timing, const RealType rx_time,
												 const timing::Timing* const tx_timing,
												 const RealType tx_time) const noexcept
	{
		return sample(tx_timing, tx_time) - sample(rx_timing, rx_time);
	}

	void solveRe(const Transmitter* trans, const Receiver* recv, const Target* targ,
				 const std::chrono::duration<RealType>& time, const RadarSignal* wave, ReResults& results)
	{
		// Note: RangeError log messages are handled by the original catch block in calculateResponse
		// or explicitly here if strict adherence to original logging is required.
		// Using the helper logic which throws RangeError on epsilon check.

		const RealType t_val = time.count();
		const auto p_tx = trans->getPosition(t_val);
		const auto p_rx = recv->getPosition(t_val);
		const auto p_tgt = targ->getPosition(t_val);

		// Link 1: Tx -> Target
		LinkGeometry link_tx_tgt;
		// Link 2: Target -> Rx (Note: Vector for calculation is Tgt->Rx)
		LinkGeometry link_tgt_rx;

		try
		{
			link_tx_tgt = computeLink(p_tx, p_tgt);
			link_tgt_rx = computeLink(p_tgt, p_rx); // Vector Tgt -> Rx
		}
		catch (const RangeError&)
		{
			LOG(Level::INFO, "Transmitter or Receiver too close to Target for accurate simulation");
			throw;
		}

		results.delay = (link_tx_tgt.dist + link_tgt_rx.dist) / params::c();

		// Calculate RCS
		// Note: getRcs expects (InAngle, OutAngle).
		// InAngle: Tx -> Tgt (link_tx_tgt.u_vec)
		// OutAngle: Rx -> Tgt (Opposite of Tgt->Rx, so -link_tgt_rx.u_vec)
		// This matches existing logic.
		SVec3 in_angle(link_tx_tgt.u_vec);
		SVec3 out_angle(-link_tgt_rx.u_vec);
		const auto rcs = targ->getRcs(in_angle, out_angle, t_val);

		const auto wavelength = params::c() / wave->getCarrier();

		// Tx Gain: Direction Tx -> Tgt
		const auto tx_gain = computeAntennaGain(trans, link_tx_tgt.u_vec, t_val, wavelength);
		// Rx Gain: Direction Rx -> Tgt (Opposite of Tgt->Rx).
		// Time is time + delay.
		const auto rx_gain = computeAntennaGain(recv, -link_tgt_rx.u_vec, results.delay + t_val, wavelength);

		const bool no_loss = recv->checkFlag(Receiver::RecvFlag::FLAG_NOPROPLOSS);
		results.power =
			computeReflectedPathPower(tx_gain, rx_gain, rcs, wavelength, link_tx_tgt.dist, link_tgt_rx.dist, no_loss);

		results.phase = -results.delay * 2 * PI * wave->getCarrier();
	}

	void solveReDirect(const Transmitter* trans, const Receiver* recv, const std::chrono::duration<RealType>& time,
					   const RadarSignal* wave, ReResults& results)
	{
		const RealType t_val = time.count();
		const auto p_tx = trans->getPosition(t_val);
		const auto p_rx = recv->getPosition(t_val);

		LinkGeometry link;
		try
		{
			link = computeLink(p_tx, p_rx); // Vector Tx -> Rx
		}
		catch (const RangeError&)
		{
			LOG(Level::INFO, "Transmitter or Receiver too close for accurate simulation");
			throw;
		}

		results.delay = link.dist / params::c();
		const RealType wavelength = params::c() / wave->getCarrier();

		// Discrepancy Fix: Original code used (Rx - Tx) for Receiver Gain but (Tx - Rx) logic for Transmitter gain
		// was ambiguous/incorrect (using `tpos - rpos` which is Rx->Tx).
		// Per `calculateDirectPathContribution` preference:
		// Tx Gain uses Vector Tx -> Rx.
		// Rx Gain uses Vector Rx -> Tx.

		const auto tx_gain = computeAntennaGain(trans, link.u_vec, t_val, wavelength);
		const auto rx_gain = computeAntennaGain(recv, -link.u_vec, t_val + results.delay, wavelength);

		const bool no_loss = recv->checkFlag(Receiver::RecvFlag::FLAG_NOPROPLOSS);
		results.power = computeDirectPathPower(tx_gain, rx_gain, wavelength, link.dist, no_loss);

		results.phase = -results.delay * 2 * PI * wave->getCarrier();
	}

	ComplexType calculateDirectPathContribution(const Transmitter* trans, const Receiver* recv, const RealType timeK,
												const CwPhaseNoiseLookup* const phase_noise_lookup)
	{
		return calculateStreamingDirectPathContribution(makeClassicStreamingSource(trans), recv, timeK,
														phase_noise_lookup);
	}

	ComplexType calculateStreamingDirectPathContribution(const core::ActiveStreamingSource& source,
														 const Receiver* recv, const RealType timeK,
														 const CwPhaseNoiseLookup* const phase_noise_lookup,
														 core::FmcwChirpBoundaryTracker* const chirp_tracker)
	{
		const auto* const trans = source.transmitter;
		if (trans == nullptr)
		{
			return {0.0, 0.0};
		}
		// Check for co-location to prevent singularities.
		// If they share the same platform, we assume they are isolated (no leakage) or explicit
		// monostatic handling is required (which is not modeled via the far-field path).
		if (trans->getPlatform() == recv->getPlatform())
		{
			return {0.0, 0.0};
		}

		const auto p_tx = trans->getPlatform()->getPosition(timeK);
		const auto p_rx = recv->getPlatform()->getPosition(timeK);

		LinkGeometry link;
		try
		{
			link = computeLink(p_tx, p_rx);
		}
		catch (const RangeError&)
		{
			return {0.0, 0.0};
		}

		const RealType tau = link.dist / params::c();
		if (source.carrier_freq <= 0.0)
		{
			return {0.0, 0.0};
		}
		const RealType lambda = params::c() / source.carrier_freq;

		// Tx Gain: Direction Tx -> Rx
		const RealType tx_gain = computeAntennaGain(trans, link.u_vec, timeK, lambda);
		// Rx Gain: Direction Rx -> Tx (-u_vec)
		const RealType rx_gain = computeAntennaGain(recv, -link.u_vec, timeK + tau, lambda);

		const bool no_loss = recv->checkFlag(Receiver::RecvFlag::FLAG_NOPROPLOSS);
		const RealType scaling_factor = computeDirectPathPower(tx_gain, rx_gain, lambda, link.dist, no_loss);

		// Include Signal Power
		const RealType amplitude = source.amplitude * std::sqrt(scaling_factor);

		// Carrier Phase
		RealType phase = 0.0;
		if (!computeStreamingPhase(source, timeK, tau, chirp_tracker, phase))
		{
			return {0.0, 0.0};
		}
		ComplexType contribution = std::polar(amplitude, phase);

		// Non-coherent Local Oscillator Effects
		const RealType non_coherent_phase = computeTimingPhase(trans, recv, timeK, timeK - tau, phase_noise_lookup);
		contribution *= std::polar(1.0, non_coherent_phase);

		return contribution;
	}

	ComplexType calculateReflectedPathContribution(const Transmitter* trans, const Receiver* recv, const Target* targ,
												   const RealType timeK,
												   const CwPhaseNoiseLookup* const phase_noise_lookup)
	{
		return calculateStreamingReflectedPathContribution(makeClassicStreamingSource(trans), recv, targ, timeK,
														   phase_noise_lookup);
	}

	ComplexType calculateStreamingReflectedPathContribution(const core::ActiveStreamingSource& source,
															const Receiver* recv, const Target* targ,
															const RealType timeK,
															const CwPhaseNoiseLookup* const phase_noise_lookup,
															core::FmcwChirpBoundaryTracker* const chirp_tracker)
	{
		const auto* const trans = source.transmitter;
		if (trans == nullptr)
		{
			return {0.0, 0.0};
		}
		// Check for co-location involving the target.
		// We do not model a platform tracking itself (R=0) or illuminating itself (R=0).
		if (trans->getPlatform() == targ->getPlatform() || recv->getPlatform() == targ->getPlatform())
		{
			return {0.0, 0.0};
		}

		const auto p_tx = trans->getPlatform()->getPosition(timeK);
		const auto p_rx = recv->getPlatform()->getPosition(timeK);
		const auto p_tgt = targ->getPlatform()->getPosition(timeK);

		LinkGeometry link_tx_tgt;
		LinkGeometry link_tgt_rx;

		try
		{
			link_tx_tgt = computeLink(p_tx, p_tgt);
			link_tgt_rx = computeLink(p_tgt, p_rx);
		}
		catch (const RangeError&)
		{
			return {0.0, 0.0};
		}

		const RealType tau = (link_tx_tgt.dist + link_tgt_rx.dist) / params::c();
		if (source.carrier_freq <= 0.0)
		{
			return {0.0, 0.0};
		}
		const RealType lambda = params::c() / source.carrier_freq;

		// RCS Lookups: In (Tx->Tgt), Out (Rx->Tgt = - (Tgt->Rx))
		SVec3 in_angle(link_tx_tgt.u_vec);
		SVec3 out_angle(-link_tgt_rx.u_vec);
		const RealType rcs = targ->getRcs(in_angle, out_angle, timeK);

		// Tx Gain: Direction Tx -> Tgt
		const RealType tx_gain = computeAntennaGain(trans, link_tx_tgt.u_vec, timeK, lambda);
		// Rx Gain: Direction Rx -> Tgt (- (Tgt->Rx)). Time: timeK + tau.
		const RealType rx_gain = computeAntennaGain(recv, -link_tgt_rx.u_vec, timeK + tau, lambda);

		const bool no_loss = recv->checkFlag(Receiver::RecvFlag::FLAG_NOPROPLOSS);
		const RealType scaling_factor =
			computeReflectedPathPower(tx_gain, rx_gain, rcs, lambda, link_tx_tgt.dist, link_tgt_rx.dist, no_loss);

		// Include Signal Power
		const RealType amplitude = source.amplitude * std::sqrt(scaling_factor);

		RealType phase = 0.0;
		if (!computeStreamingPhase(source, timeK, tau, chirp_tracker, phase))
		{
			return {0.0, 0.0};
		}
		ComplexType contribution = std::polar(amplitude, phase);

		// Non-coherent Local Oscillator Effects
		const RealType non_coherent_phase = computeTimingPhase(trans, recv, timeK, timeK - tau, phase_noise_lookup);
		contribution *= std::polar(1.0, non_coherent_phase);

		return contribution;
	}

	std::unique_ptr<serial::Response> calculateResponse(const Transmitter* trans, const Receiver* recv,
														const RadarSignal* signal, const RealType startTime,
														const Target* targ)
	{
		// If calculating direct path (no target) and components are co-located:
		// 1. If explicitly attached (monostatic), skip (internal leakage handled elsewhere).
		// 2. If independent but on the same platform, distance is 0. Far-field logic (1/R^2)
		//    diverges. We skip calculation to avoid RangeError crashes, assuming
		//    no direct coupling/interference for co-located far-field antennas.
		if (targ == nullptr && (trans->getAttached() == recv || trans->getPlatform() == recv->getPlatform()))
		{
			return nullptr;
		}

		// If calculating reflected path and target is co-located with either Tx or Rx:
		// Skip to avoid singularity. Simulating a radar tracking its own platform
		// requires near-field clutter models, not point-target RCS models.
		if (targ != nullptr &&
			(targ->getPlatform() == trans->getPlatform() || targ->getPlatform() == recv->getPlatform()))
		{
			LOG(Level::TRACE,
				"Skipping reflected path calculation for Target {} co-located with Transmitter {} or Receiver {}",
				targ->getName(), trans->getName(), recv->getName());
			return nullptr;
		}

		const auto start_time_chrono = std::chrono::duration<RealType>(startTime);
		const auto end_time_chrono = start_time_chrono + std::chrono::duration<RealType>(signal->getLength());
		const auto sample_time_chrono = std::chrono::duration<RealType>(1.0 / params::simSamplingRate());
		const int point_count = static_cast<int>(std::ceil(signal->getLength() / sample_time_chrono.count()));

		if ((targ != nullptr) && point_count == 0)
		{
			LOG(Level::FATAL, "No time points are available for execution!");
			throw std::runtime_error("No time points are available for execution!");
		}

		auto response = std::make_unique<serial::Response>(signal, trans);

		try
		{
			for (int i = 0; i <= point_count; ++i)
			{
				const auto current_time =
					i < point_count ? start_time_chrono + i * sample_time_chrono : end_time_chrono;

				ReResults results{};
				if (targ != nullptr)
				{
					solveRe(trans, recv, targ, current_time, signal, results);
				}
				else
				{
					solveReDirect(trans, recv, current_time, signal, results);
				}

				interp::InterpPoint point{.power = results.power,
										  .time = current_time.count() + results.delay,
										  .delay = results.delay,
										  .phase = results.phase};
				response->addInterpPoint(point);
			}
		}
		catch (const RangeError&)
		{
			LOG(Level::INFO, "Receiver or Transmitter too close for accurate simulation");
			throw; // Re-throw to be caught by the runner
		}

		return response;
	}

	std::vector<PreviewLink> calculatePreviewLinks(const core::World& world, const RealType time)
	{
		std::vector<PreviewLink> links;

		// Default wavelength (1GHz) if no waveform is attached, to allow geometric visualization
		const RealType lambda_default = 0.3;

		for (const auto& tx : world.getTransmitters())
		{
			// 1. Check Transmitter Schedule
			if (!isComponentActive(tx->getSchedule(), time))
			{
				continue;
			}

			const auto p_tx = tx->getPosition(time);
			const auto* waveform = tx->getSignal();
			const RealType pt = previewRadiatedPower(waveform);
			const RealType lambda = (waveform != nullptr) ? (params::c() / waveform->getCarrier()) : lambda_default;

			// --- PRE-CALCULATE ILLUMINATOR PATHS (Tx -> Tgt) ---
			// These depend only on the Transmitter and Targets. We calculate them once per Tx
			// to avoid duplication when multiple receivers are present.
			// We do NOT filter out Monostatic transmitters here. Even in a monostatic setup,
			// the "Illuminator" link provides distinct info (Power Density at Target) compared
			// to the "Monostatic" link (Received Power at Rx). By being outside the Rx loop,
			// it is guaranteed to be unique per Target.
			for (const auto& tgt : world.getTargets())
			{
				const auto p_tgt = tgt->getPosition(time);
				const Vec3 vec_tx_tgt = p_tgt - p_tx;
				const RealType r1 = vec_tx_tgt.length();

				if (r1 <= EPSILON)
					continue;

				const Vec3 u_tx_tgt = vec_tx_tgt / r1;
				// Tx Gain: Tx -> Tgt
				const RealType gt = computeAntennaGain(tx.get(), u_tx_tgt, time, lambda);

				// Power Density at Target: S = (Pt * Gt) / (4 * pi * R1^2)
				const RealType p_density = (pt * gt) / (4.0 * PI * r1 * r1);

				links.push_back({.type = LinkType::BistaticTxTgt,
								 .quality = LinkQuality::Strong,
								 .label = formatPreviewDbwPerSquareMeterLabel(p_density),
								 .source_id = tx->getId(),
								 .dest_id = tgt->getId(),
								 .origin_id = tx->getId()});
			}

			for (const auto& rx : world.getReceivers())
			{
				// 2. Check Receiver Schedule
				if (!isComponentActive(rx->getSchedule(), time))
				{
					continue;
				}

				const auto p_rx = rx->getPosition(time);
				const bool is_monostatic = (tx->getAttached() == rx.get());
				const bool no_loss = rx->checkFlag(Receiver::RecvFlag::FLAG_NOPROPLOSS);

				if (is_monostatic)
				{
					// --- Monostatic Case ---
					// Calculates the round-trip Received Power (dBm)
					for (const auto& tgt : world.getTargets())
					{
						const auto p_tgt = tgt->getPosition(time);

						// Use computeLink helper logic, but catch exception manually by checking dist
						// to avoid try/catch overhead in render loop
						const Vec3 vec_tx_tgt = p_tgt - p_tx;
						const RealType dist = vec_tx_tgt.length();

						if (dist <= EPSILON)
							continue;

						const Vec3 u_tx_tgt = vec_tx_tgt / dist; // Unit vec Tx -> Tgt

						// Reusing internal helpers
						// Tx Gain: Direction Tx -> Tgt
						const RealType gt = computeAntennaGain(tx.get(), u_tx_tgt, time, lambda);
						// Rx Gain: Direction Rx -> Tgt (which is -u_tx_tgt because Rx is at Tx)
						const RealType gr = computeAntennaGain(rx.get(), u_tx_tgt, time, lambda);

						// RCS: In = Tx->Tgt, Out = Tgt->Rx (which is -(Tx->Tgt))
						SVec3 in_angle(u_tx_tgt);
						SVec3 out_angle(-u_tx_tgt);
						const RealType rcs = tgt->getRcs(in_angle, out_angle, time);

						// Reusing Reflected Path Helper
						// r_tx = dist, r_rx = dist
						const RealType power_ratio =
							computeReflectedPathPower(gt, gr, rcs, lambda, dist, dist, no_loss);

						const RealType pr_watts = pt * power_ratio;

						// Label shows power with unit RCS (sigma=1) so it varies only with geometry,
						// not with the fluctuation model. Actual RCS is carried separately for the tooltip.
						const RealType pr_unit_watts =
							pt * computeReflectedPathPower(gt, gr, 1.0, lambda, dist, dist, no_loss);

						links.push_back({.type = LinkType::Monostatic,
										 .quality = isSignalStrong(pr_unit_watts, rx->getNoiseTemperature())
											 ? LinkQuality::Strong
											 : LinkQuality::Weak,
										 .label = formatPreviewDbmLabel(pr_unit_watts),
										 .source_id = tx->getId(),
										 .dest_id = tgt->getId(),
										 .origin_id = tx->getId(),
										 .rcs = rcs,
										 .actual_power_dbm = wattsToDbm(pr_watts)});
					}
				}
				else
				{
					// --- Bistatic Case ---

					// 1. Direct Path (Interference)
					if (!rx->checkFlag(Receiver::RecvFlag::FLAG_NODIRECT))
					{
						const Vec3 vec_direct = p_rx - p_tx;
						const RealType dist = vec_direct.length();

						if (dist > EPSILON)
						{
							const Vec3 u_tx_rx = vec_direct / dist;

							// Tx Gain: Tx -> Rx
							const RealType gt = computeAntennaGain(tx.get(), u_tx_rx, time, lambda);
							// Rx Gain: Rx -> Tx (which is -u_tx_rx)
							const RealType gr = computeAntennaGain(rx.get(), -u_tx_rx, time, lambda);

							const RealType power_ratio = computeDirectPathPower(gt, gr, lambda, dist, no_loss);
							const RealType pr_watts = pt * power_ratio;

							links.push_back({.type = LinkType::DirectTxRx,
											 .quality = LinkQuality::Strong,
											 .label = formatPreviewDbmLabel(pr_watts, "Direct: "),
											 .source_id = tx->getId(),
											 .dest_id = rx->getId(),
											 .origin_id = tx->getId()});
						}
					}

					// 2. Reflected Path (Scattered Leg: Tgt -> Rx)
					for (const auto& tgt : world.getTargets())
					{
						const auto p_tgt = tgt->getPosition(time);
						const Vec3 vec_tx_tgt = p_tgt - p_tx;
						const Vec3 vec_tgt_rx = p_rx - p_tgt; // Vector Tgt -> Rx

						const RealType r1 = vec_tx_tgt.length();
						const RealType r2 = vec_tgt_rx.length();

						if (r1 <= EPSILON || r2 <= EPSILON)
							continue;

						const Vec3 u_tx_tgt = vec_tx_tgt / r1;
						const Vec3 u_tgt_rx = vec_tgt_rx / r2;

						// Tx Gain: Tx -> Tgt
						const RealType gt = computeAntennaGain(tx.get(), u_tx_tgt, time, lambda);
						// Rx Gain: Rx -> Tgt (Vector is Tgt->Rx, so look vector is -u_tgt_rx)
						const RealType gr = computeAntennaGain(rx.get(), -u_tgt_rx, time, lambda);

						// RCS
						SVec3 in_angle(u_tx_tgt);
						SVec3 out_angle(-u_tgt_rx); // Out angle is usually defined from Target OUTWARDS
						const RealType rcs = tgt->getRcs(in_angle, out_angle, time);

						const RealType power_ratio = computeReflectedPathPower(gt, gr, rcs, lambda, r1, r2, no_loss);
						const RealType pr_watts = pt * power_ratio;

						// Label shows power with unit RCS (sigma=1) so it varies only with geometry.
						const RealType pr_unit_watts =
							pt * computeReflectedPathPower(gt, gr, 1.0, lambda, r1, r2, no_loss);

						// Note: Illuminator leg (Tx->Tgt) was handled in the outer loop.

						// Leg 2: Scattered
						links.push_back({.type = LinkType::BistaticTgtRx,
										 .quality = isSignalStrong(pr_unit_watts, rx->getNoiseTemperature())
											 ? LinkQuality::Strong
											 : LinkQuality::Weak,
										 .label = formatPreviewDbmLabel(pr_unit_watts),
										 .source_id = tgt->getId(),
										 .dest_id = rx->getId(),
										 .origin_id = tx->getId(),
										 .rcs = rcs,
										 .actual_power_dbm = wattsToDbm(pr_watts)});
					}
				}
			}
		}
		return links;
	}
}
