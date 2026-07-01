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

	/// Initializes an FMCW triangle tracker from a retarded time.
	void initializeFmcwTriangleTracker(const core::ActiveStreamingSource& source, const RealType t_ret,
									   core::FmcwChirpBoundaryTracker& tracker)
	{
		tracker.triangle_initialized = true;
		if (t_ret < source.segment_start)
		{
			tracker.triangle_index = 0;
			tracker.triangle_leg = 0;
			tracker.triangle_t_leg = source.segment_start;
			tracker.triangle_phi_base = 0.0;
			return;
		}

		const RealType delta = t_ret - source.segment_start;
		tracker.triangle_index = static_cast<std::size_t>(std::floor(delta / source.triangle_period));
		const RealType local_triangle_time =
			delta - static_cast<RealType>(tracker.triangle_index) * source.triangle_period;
		tracker.triangle_leg = local_triangle_time < source.chirp_duration ? 0U : 1U;
		tracker.triangle_t_leg = source.segment_start +
			static_cast<RealType>(tracker.triangle_index) * source.triangle_period +
			(tracker.triangle_leg == 1U ? source.chirp_duration : 0.0);
		tracker.triangle_phi_base =
			std::fmod(static_cast<RealType>(tracker.triangle_index) * source.mod_phi_tri, 2.0 * PI) +
			(tracker.triangle_leg == 1U ? source.mod_phi_up : 0.0);
		if (tracker.triangle_phi_base >= 2.0 * PI)
		{
			tracker.triangle_phi_base -= 2.0 * PI;
		}
		if (tracker.triangle_phi_base < 0.0)
		{
			tracker.triangle_phi_base += 2.0 * PI;
		}
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

	struct StreamingWaveformEvaluation
	{
		RealType phase = 0.0;
		RealType rf_frequency = 0.0;
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

	/// Computes deterministic timing phase for one timing source when no lookup is available.
	RealType computeSingleTimingPhase(const timing::Timing* const timing, const RealType time)
	{
		if (timing == nullptr)
		{
			return 0.0;
		}
		return 2.0 * PI * timing->getFreqOffset() * time + timing->getPhaseOffset();
	}

	/// Computes timing phase with an optional lookup and streaming phase-application mode.
	RealType computeTimingPhase(const Transmitter* tx, const Receiver* rx, const RealType rx_time,
								const RealType tx_time, const simulation::CwPhaseNoiseLookup* phase_noise_lookup,
								const simulation::StreamingTimingPhaseMode mode)
	{
		if (mode == simulation::StreamingTimingPhaseMode::None)
		{
			return 0.0;
		}
		if (mode == simulation::StreamingTimingPhaseMode::TransmitterOnly)
		{
			if (phase_noise_lookup == nullptr)
			{
				return computeSingleTimingPhase(tx->getTiming().get(), tx_time);
			}
			return phase_noise_lookup->sample(tx->getTiming().get(), tx_time);
		}
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
		if (source.kind == core::StreamingWaveformKind::Cw)
		{
			source.segment_start = std::numeric_limits<RealType>::lowest();
		}
		return source;
	}

	bool computeLinearFmcwPhaseWithoutTracker(const core::ActiveStreamingSource& source, const RealType t_ret,
											  const RealType tau, RealType& phase_out)
	{
		const RealType time_since_segment_start = t_ret - source.segment_start;
		const auto chirp_index = static_cast<std::size_t>(std::floor(time_since_segment_start / source.chirp_period));
		if (source.chirp_count.has_value() && chirp_index >= *source.chirp_count)
		{
			return false;
		}
		const RealType chirp_time = time_since_segment_start - static_cast<RealType>(chirp_index) * source.chirp_period;
		if (chirp_time < 0.0 || chirp_time >= source.chirp_duration)
		{
			return false;
		}
		phase_out = -2.0 * PI * source.carrier_freq * tau + source.two_pi_f0 * chirp_time +
			source.s_pi_alpha * chirp_time * chirp_time;
		return true;
	}

	bool computeLinearFmcwPhaseWithTracker(const core::ActiveStreamingSource& source, const RealType t_ret,
										   const RealType tau, core::FmcwChirpBoundaryTracker& chirp_tracker,
										   RealType& phase_out)
	{
		if (!chirp_tracker.initialized)
		{
			initializeFmcwTracker(source, t_ret, chirp_tracker);
		}

		while (t_ret >= chirp_tracker.t_n + source.chirp_period)
		{
			chirp_tracker.t_n += source.chirp_period;
			++chirp_tracker.n_current;
		}

		if (source.chirp_count.has_value() && chirp_tracker.n_current >= *source.chirp_count)
		{
			return false;
		}

		const RealType u_ret = t_ret - chirp_tracker.t_n;
		if (u_ret < 0.0 || u_ret >= source.chirp_duration)
		{
			return false;
		}

		phase_out =
			-2.0 * PI * source.carrier_freq * tau + source.two_pi_f0 * u_ret + source.s_pi_alpha * u_ret * u_ret;
		return true;
	}

	bool computeTriangleFmcwPhaseWithoutTracker(const core::ActiveStreamingSource& source, const RealType t_ret,
												const RealType tau, RealType& phase_out)
	{
		const RealType delta = t_ret - source.segment_start;
		const auto triangle_index = static_cast<std::size_t>(std::floor(delta / source.triangle_period));
		if (source.triangle_count.has_value() && triangle_index >= *source.triangle_count)
		{
			return false;
		}
		const RealType local_triangle_time = delta - static_cast<RealType>(triangle_index) * source.triangle_period;
		const bool down_leg = local_triangle_time >= source.chirp_duration;
		const RealType u_ret = down_leg ? local_triangle_time - source.chirp_duration : local_triangle_time;
		if (u_ret < 0.0 || u_ret >= source.chirp_duration)
		{
			return false;
		}
		const RealType phi_base = std::fmod(static_cast<RealType>(triangle_index) * source.mod_phi_tri, 2.0 * PI) +
			(down_leg ? source.mod_phi_up : 0.0);
		const RealType modular_phi_base = phi_base >= 2.0 * PI ? phi_base - 2.0 * PI : phi_base;
		const RealType linear_coeff = down_leg ? source.two_pi_f0_plus_B : source.two_pi_f0;
		const RealType quad_coeff = down_leg ? source.neg_pi_alpha : source.pi_alpha;
		phase_out = -2.0 * PI * source.carrier_freq * tau + modular_phi_base + linear_coeff * u_ret +
			quad_coeff * u_ret * u_ret;
		return true;
	}

	void advanceTriangleTracker(const core::ActiveStreamingSource& source, const RealType t_ret,
								core::FmcwChirpBoundaryTracker& chirp_tracker)
	{
		while (t_ret >= chirp_tracker.triangle_t_leg + source.chirp_duration)
		{
			chirp_tracker.triangle_t_leg += source.chirp_duration;
			chirp_tracker.triangle_leg = 1U - chirp_tracker.triangle_leg;
			if (chirp_tracker.triangle_leg == 0U)
			{
				++chirp_tracker.triangle_index;
			}
			chirp_tracker.triangle_phi_base += source.mod_phi_up;
			if (chirp_tracker.triangle_phi_base >= 2.0 * PI)
			{
				chirp_tracker.triangle_phi_base -= 2.0 * PI;
			}
		}
	}

	bool computeTriangleFmcwPhaseWithTracker(const core::ActiveStreamingSource& source, const RealType t_ret,
											 const RealType tau, core::FmcwChirpBoundaryTracker& chirp_tracker,
											 RealType& phase_out)
	{
		if (!chirp_tracker.triangle_initialized)
		{
			initializeFmcwTriangleTracker(source, t_ret, chirp_tracker);
		}

		advanceTriangleTracker(source, t_ret, chirp_tracker);
		if (source.triangle_count.has_value() && chirp_tracker.triangle_index >= *source.triangle_count)
		{
			return false;
		}

		const RealType u_ret = t_ret - chirp_tracker.triangle_t_leg;
		if (u_ret < 0.0 || u_ret >= source.chirp_duration)
		{
			return false;
		}

		const bool down_leg = chirp_tracker.triangle_leg == 1U;
		const RealType linear_coeff = down_leg ? source.two_pi_f0_plus_B : source.two_pi_f0;
		const RealType quad_coeff = down_leg ? source.neg_pi_alpha : source.pi_alpha;
		phase_out = -2.0 * PI * source.carrier_freq * tau + chirp_tracker.triangle_phi_base + linear_coeff * u_ret +
			quad_coeff * u_ret * u_ret;
		return true;
	}

	/// Computes streaming waveform phase and active RF at a receiver time.
	bool computeStreamingEvaluation(const core::ActiveStreamingSource& source, const RealType rx_time,
									const RealType tau, core::FmcwChirpBoundaryTracker* const chirp_tracker,
									StreamingWaveformEvaluation& eval)
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

		if (source.kind == core::StreamingWaveformKind::FmcwLinear)
		{
			eval.rf_frequency = source.carrier_freq;
			if (chirp_tracker == nullptr)
			{
				return computeLinearFmcwPhaseWithoutTracker(source, t_ret, tau, eval.phase);
			}
			return computeLinearFmcwPhaseWithTracker(source, t_ret, tau, *chirp_tracker, eval.phase);
		}

		if (source.kind == core::StreamingWaveformKind::FmcwTriangle)
		{
			eval.rf_frequency = source.carrier_freq;
			if (chirp_tracker == nullptr)
			{
				return computeTriangleFmcwPhaseWithoutTracker(source, t_ret, tau, eval.phase);
			}
			return computeTriangleFmcwPhaseWithTracker(source, t_ret, tau, *chirp_tracker, eval.phase);
		}

		if (source.kind == core::StreamingWaveformKind::Sfcw)
		{
			if (source.sfcw == nullptr)
			{
				return false;
			}
			const auto step = source.sfcw->activeStepAt(t_ret - source.segment_start, source.carrier_freq);
			if (!step.has_value() || step->rf_frequency <= 0.0)
			{
				return false;
			}
			eval.rf_frequency = step->rf_frequency;
			eval.phase = -2.0 * PI * step->rf_frequency * tau;
			return true;
		}

		eval.rf_frequency = source.carrier_freq;
		eval.phase = -2.0 * PI * source.carrier_freq * tau;
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
		if (waveform->isFmcwTriangle())
		{
			return waveform->getPower();
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
														 core::FmcwChirpBoundaryTracker* const chirp_tracker,
														 const StreamingTimingPhaseMode timing_phase_mode)
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
		StreamingWaveformEvaluation eval;
		if (!computeStreamingEvaluation(source, timeK, tau, chirp_tracker, eval))
		{
			return {0.0, 0.0};
		}
		const RealType lambda = params::c() / eval.rf_frequency;

		// Tx Gain: Direction Tx -> Rx
		const RealType tx_gain = computeAntennaGain(trans, link.u_vec, timeK, lambda);
		// Rx Gain: Direction Rx -> Tx (-u_vec)
		const RealType rx_gain = computeAntennaGain(recv, -link.u_vec, timeK + tau, lambda);

		const bool no_loss = recv->checkFlag(Receiver::RecvFlag::FLAG_NOPROPLOSS);
		const RealType scaling_factor = computeDirectPathPower(tx_gain, rx_gain, lambda, link.dist, no_loss);

		// Include Signal Power
		const RealType amplitude = source.amplitude * std::sqrt(scaling_factor);

		// Carrier Phase
		ComplexType contribution = std::polar(amplitude, eval.phase);

		// Non-coherent Local Oscillator Effects
		const RealType non_coherent_phase =
			computeTimingPhase(trans, recv, timeK, timeK - tau, phase_noise_lookup, timing_phase_mode);
		contribution *= std::polar(1.0, non_coherent_phase);

		return contribution;
	}

	bool calculateStreamingReferencePhase(const core::ActiveStreamingSource& source, const RealType timeK,
										  core::FmcwChirpBoundaryTracker* const chirp_tracker, RealType& phase_out)
	{
		StreamingWaveformEvaluation eval;
		if (!computeStreamingEvaluation(source, timeK, 0.0, chirp_tracker, eval))
		{
			return false;
		}
		phase_out = eval.phase;
		return true;
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
															core::FmcwChirpBoundaryTracker* const chirp_tracker,
															const StreamingTimingPhaseMode timing_phase_mode)
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
		StreamingWaveformEvaluation eval;
		if (!computeStreamingEvaluation(source, timeK, tau, chirp_tracker, eval))
		{
			return {0.0, 0.0};
		}
		const RealType lambda = params::c() / eval.rf_frequency;

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

		ComplexType contribution = std::polar(amplitude, eval.phase);

		// Non-coherent Local Oscillator Effects
		const RealType non_coherent_phase =
			computeTimingPhase(trans, recv, timeK, timeK - tau, phase_noise_lookup, timing_phase_mode);
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

				interp::InterpPoint const point{.power = results.power,
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

	namespace
	{
		struct PreviewTransmitterContext
		{
			const Transmitter& transmitter;
			Vec3 position;
			RealType radiated_power;
			RealType lambda;
		};

		struct PreviewReceiverContext
		{
			const Receiver& receiver;
			Vec3 position;
			bool no_loss;
		};

		PreviewTransmitterContext makePreviewTransmitterContext(const Transmitter& transmitter, const RealType time)
		{
			const auto* waveform = transmitter.getSignal();
			return PreviewTransmitterContext{.transmitter = transmitter,
											 .position = transmitter.getPosition(time),
											 .radiated_power = previewRadiatedPower(waveform),
											 .lambda =
												 (waveform != nullptr) ? (params::c() / waveform->getCarrier()) : 0.3};
		}

		PreviewReceiverContext makePreviewReceiverContext(const Receiver& receiver, const RealType time)
		{
			return PreviewReceiverContext{.receiver = receiver,
										  .position = receiver.getPosition(time),
										  .no_loss = receiver.checkFlag(Receiver::RecvFlag::FLAG_NOPROPLOSS)};
		}

		void addIlluminatorLinks(std::vector<PreviewLink>& links, const PreviewTransmitterContext& tx_ctx,
								 const core::World& world, const RealType time)
		{
			for (const auto& target : world.getTargets())
			{
				const auto target_position = target->getPosition(time);
				const Vec3 vec_tx_tgt = target_position - tx_ctx.position;
				const RealType range = vec_tx_tgt.length();
				if (range <= EPSILON)
				{
					continue;
				}

				const Vec3 u_tx_tgt = vec_tx_tgt / range;
				const RealType gain = computeAntennaGain(&tx_ctx.transmitter, u_tx_tgt, time, tx_ctx.lambda);
				const RealType power_density = (tx_ctx.radiated_power * gain) / (4.0 * PI * range * range);
				links.push_back({.type = LinkType::BistaticTxTgt,
								 .quality = LinkQuality::Strong,
								 .label = formatPreviewDbwPerSquareMeterLabel(power_density),
								 .display_value = wattsToDb(power_density),
								 .source_id = tx_ctx.transmitter.getId(),
								 .dest_id = target->getId(),
								 .origin_id = tx_ctx.transmitter.getId()});
			}
		}

		void addMonostaticLinks(std::vector<PreviewLink>& links, const PreviewTransmitterContext& tx_ctx,
								const PreviewReceiverContext& rx_ctx, const core::World& world, const RealType time)
		{
			for (const auto& target : world.getTargets())
			{
				const auto target_position = target->getPosition(time);
				const Vec3 vec_tx_tgt = target_position - tx_ctx.position;
				const RealType range = vec_tx_tgt.length();
				if (range <= EPSILON)
				{
					continue;
				}

				const Vec3 u_tx_tgt = vec_tx_tgt / range;
				const RealType gt = computeAntennaGain(&tx_ctx.transmitter, u_tx_tgt, time, tx_ctx.lambda);
				const RealType gr = computeAntennaGain(&rx_ctx.receiver, u_tx_tgt, time, tx_ctx.lambda);
				SVec3 in_angle(u_tx_tgt);
				SVec3 out_angle(-u_tx_tgt);
				const RealType rcs = target->getRcs(in_angle, out_angle, time);
				const RealType power_ratio =
					computeReflectedPathPower(gt, gr, rcs, tx_ctx.lambda, range, range, rx_ctx.no_loss);
				const RealType pr_watts = tx_ctx.radiated_power * power_ratio;
				const RealType pr_unit_watts = tx_ctx.radiated_power *
					computeReflectedPathPower(gt, gr, 1.0, tx_ctx.lambda, range, range, rx_ctx.no_loss);

				links.push_back({.type = LinkType::Monostatic,
								 .quality = isSignalStrong(pr_unit_watts, rx_ctx.receiver.getNoiseTemperature())
									 ? LinkQuality::Strong
									 : LinkQuality::Weak,
								 .label = formatPreviewDbmLabel(pr_unit_watts),
								 .display_value = wattsToDbm(pr_unit_watts),
								 .source_id = tx_ctx.transmitter.getId(),
								 .dest_id = target->getId(),
								 .origin_id = tx_ctx.transmitter.getId(),
								 .rcs = rcs,
								 .actual_power_dbm = wattsToDbm(pr_watts)});
			}
		}

		void addDirectLink(std::vector<PreviewLink>& links, const PreviewTransmitterContext& tx_ctx,
						   const PreviewReceiverContext& rx_ctx, const RealType time)
		{
			if (rx_ctx.receiver.checkFlag(Receiver::RecvFlag::FLAG_NODIRECT))
			{
				return;
			}

			const Vec3 vec_direct = rx_ctx.position - tx_ctx.position;
			const RealType range = vec_direct.length();
			if (range <= EPSILON)
			{
				return;
			}

			const Vec3 u_tx_rx = vec_direct / range;
			const RealType gt = computeAntennaGain(&tx_ctx.transmitter, u_tx_rx, time, tx_ctx.lambda);
			const RealType gr = computeAntennaGain(&rx_ctx.receiver, -u_tx_rx, time, tx_ctx.lambda);
			const RealType power_ratio = computeDirectPathPower(gt, gr, tx_ctx.lambda, range, rx_ctx.no_loss);
			const RealType pr_watts = tx_ctx.radiated_power * power_ratio;

			links.push_back({.type = LinkType::DirectTxRx,
							 .quality = LinkQuality::Strong,
							 .label = formatPreviewDbmLabel(pr_watts, "Direct: "),
							 .display_value = wattsToDbm(pr_watts),
							 .source_id = tx_ctx.transmitter.getId(),
							 .dest_id = rx_ctx.receiver.getId(),
							 .origin_id = tx_ctx.transmitter.getId()});
		}

		void addBistaticTargetReceiverLinks(std::vector<PreviewLink>& links, const PreviewTransmitterContext& tx_ctx,
											const PreviewReceiverContext& rx_ctx, const core::World& world,
											const RealType time)
		{
			for (const auto& target : world.getTargets())
			{
				const auto target_position = target->getPosition(time);
				const Vec3 vec_tx_tgt = target_position - tx_ctx.position;
				const Vec3 vec_tgt_rx = rx_ctx.position - target_position;
				const RealType r1 = vec_tx_tgt.length();
				const RealType r2 = vec_tgt_rx.length();
				if (r1 <= EPSILON || r2 <= EPSILON)
				{
					continue;
				}

				const Vec3 u_tx_tgt = vec_tx_tgt / r1;
				const Vec3 u_tgt_rx = vec_tgt_rx / r2;
				const RealType gt = computeAntennaGain(&tx_ctx.transmitter, u_tx_tgt, time, tx_ctx.lambda);
				const RealType gr = computeAntennaGain(&rx_ctx.receiver, -u_tgt_rx, time, tx_ctx.lambda);
				SVec3 in_angle(u_tx_tgt);
				SVec3 out_angle(-u_tgt_rx);
				const RealType rcs = target->getRcs(in_angle, out_angle, time);
				const RealType power_ratio =
					computeReflectedPathPower(gt, gr, rcs, tx_ctx.lambda, r1, r2, rx_ctx.no_loss);
				const RealType pr_watts = tx_ctx.radiated_power * power_ratio;
				const RealType pr_unit_watts = tx_ctx.radiated_power *
					computeReflectedPathPower(gt, gr, 1.0, tx_ctx.lambda, r1, r2, rx_ctx.no_loss);

				links.push_back({.type = LinkType::BistaticTgtRx,
								 .quality = isSignalStrong(pr_unit_watts, rx_ctx.receiver.getNoiseTemperature())
									 ? LinkQuality::Strong
									 : LinkQuality::Weak,
								 .label = formatPreviewDbmLabel(pr_unit_watts),
								 .display_value = wattsToDbm(pr_unit_watts),
								 .source_id = target->getId(),
								 .dest_id = rx_ctx.receiver.getId(),
								 .origin_id = tx_ctx.transmitter.getId(),
								 .rcs = rcs,
								 .actual_power_dbm = wattsToDbm(pr_watts)});
			}
		}

		void addReceiverLinks(std::vector<PreviewLink>& links, const PreviewTransmitterContext& tx_ctx,
							  const PreviewReceiverContext& rx_ctx, const core::World& world, const RealType time)
		{
			if (tx_ctx.transmitter.getAttached() == &rx_ctx.receiver)
			{
				addMonostaticLinks(links, tx_ctx, rx_ctx, world, time);
				return;
			}

			addDirectLink(links, tx_ctx, rx_ctx, time);
			addBistaticTargetReceiverLinks(links, tx_ctx, rx_ctx, world, time);
		}
	}

	std::vector<PreviewLink> calculatePreviewLinks(const core::World& world, const RealType time)
	{
		std::vector<PreviewLink> links;

		for (const auto& tx : world.getTransmitters())
		{
			if (!isComponentActive(tx->getSchedule(), time))
			{
				continue;
			}

			const auto tx_ctx = makePreviewTransmitterContext(*tx, time);
			addIlluminatorLinks(links, tx_ctx, world, time);

			for (const auto& rx : world.getReceivers())
			{
				if (!isComponentActive(rx->getSchedule(), time))
				{
					continue;
				}
				const auto rx_ctx = makePreviewReceiverContext(*rx, time);
				addReceiverLinks(links, tx_ctx, rx_ctx, world, time);
			}
		}
		return links;
	}
}
