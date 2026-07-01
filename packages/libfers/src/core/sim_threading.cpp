// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file sim_threading.cpp
 * @brief Implements the core event-driven simulation engine.
 *
 * This file contains the primary simulation loop, which orchestrates the entire
 * simulation process. It operates on a unified, event-driven model capable of
 * handling both pulsed and continuous-wave (CW) radar systems concurrently.
 */

#include "sim_threading.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <optional>
#include <utility>

#include "logging.h"
#include "math/path_utils.h"
#include "memory_projection.h"
#include "parameters.h"
#include "processing/finalizer.h"
#include "processing/finalizer_pipeline.h"
#include "processing/signal_processor.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/hdf5_output_sink.h"
#include "serial/response.h"
#include "serial/vita49/vita49_output_sink.h"
#include "signal/if_resampler.h"
#include "signal/radar_signal.h"
#include "sim_events.h"
#include "simulation/channel_model.h"
#include "thread_pool.h"
#include "timing/timing.h"
#include "world.h"

using logging::Level;
using radar::OperationMode;
using radar::Receiver;
using radar::Transmitter;

namespace core
{
	namespace
	{
		constexpr std::size_t fmcw_if_block_size = 1024;
		constexpr std::size_t streaming_output_block_size = 4096;

		[[nodiscard]] std::size_t expectedStreamingOutputSamples(const RealType sample_rate)
		{
			return static_cast<std::size_t>(
				std::ceil(std::max<RealType>(0.0, params::endTime() - params::startTime()) * sample_rate));
		}

		[[nodiscard]] Vita49StreamMetadata streamStatsToMetadata(const ReceiverStreamStats& stats)
		{
			return Vita49StreamMetadata{.receiver_id = stats.receiver_id,
										.receiver_name = stats.receiver_name,
										.stream_id = stats.stream_id,
										.mode = stats.mode,
										.sample_rate = stats.sample_rate,
										.reference_frequency = stats.reference_frequency,
										.packets_emitted = stats.packets_emitted,
										.samples_emitted = stats.samples_emitted,
										.packets_dropped = stats.packets_dropped,
										.samples_dropped = stats.samples_dropped,
										.over_range_count = stats.over_range_count,
										.late_packet_count = stats.late_packet_count,
										.context_packet_count = stats.context_packets,
										.first_sample_time = stats.first_sample_time,
										.end_sample_time = stats.end_sample_time,
										.first_timestamp = stats.first_timestamp,
										.end_timestamp = stats.end_timestamp};
		}

		[[nodiscard]] std::string fmcwCountToken(const std::optional<std::size_t>& count)
		{
			return count.has_value() ? std::format("{}", *count) : std::string("unbounded");
		}

		void logUnscheduledFmcwChirpSummary(const Transmitter& transmitter, const fers_signal::FmcwChirpSignal& fmcw,
											const std::string& direction, const std::string& configured_count,
											const RealType duty_cycle, const RealType average_power)
		{
			const RealType active_start = params::startTime();
			const auto source = makeActiveSource(&transmitter, active_start, params::endTime());
			const auto total_chirp_count = countFmcwChirpStarts(source, active_start, source.segment_end);
			LOG(Level::INFO,
				"FMCW transmitter '{}' shape=linear {} B={} Hz T_c={} s T_rep={} s f_0={} Hz alpha={} Hz/s "
				"duty_cycle={} chirp_count={} total_chirp_count={} average_power={} W",
				transmitter.getName(), direction, fmcw.getChirpBandwidth(), fmcw.getChirpDuration(),
				fmcw.getChirpPeriod(), fmcw.getStartFrequencyOffset(), fmcw.getChirpRate(), duty_cycle,
				configured_count, total_chirp_count, average_power);
		}

		void logScheduledFmcwChirpSummary(const Transmitter& transmitter, const fers_signal::FmcwChirpSignal& fmcw,
										  const std::string& direction, const std::string& configured_count,
										  const RealType duty_cycle, const RealType average_power)
		{
			std::uint64_t total_chirp_count = 0;
			for (const auto& period : transmitter.getSchedule())
			{
				const RealType active_start = std::max(params::startTime(), period.start);
				const auto source =
					makeActiveSource(&transmitter, period.start, std::min(params::endTime(), period.end));
				const auto segment_chirp_count = countFmcwChirpStarts(source, active_start, source.segment_end);
				total_chirp_count += segment_chirp_count;
				LOG(Level::INFO,
					"FMCW transmitter '{}' segment [{}, {}] shape=linear {} B={} Hz T_c={} s T_rep={} s f_0={} "
					"Hz alpha={} Hz/s duty_cycle={} chirp_count={} segment_chirp_count={} total_chirp_count={} "
					"average_power={} W",
					transmitter.getName(), period.start, source.segment_end, direction, fmcw.getChirpBandwidth(),
					fmcw.getChirpDuration(), fmcw.getChirpPeriod(), fmcw.getStartFrequencyOffset(), fmcw.getChirpRate(),
					duty_cycle, configured_count, segment_chirp_count, total_chirp_count, average_power);
			}
		}

		void logFmcwChirpSummary(const Transmitter& transmitter, const fers_signal::RadarSignal& waveform,
								 const fers_signal::FmcwChirpSignal& fmcw)
		{
			const RealType duty_cycle = fmcw.getChirpDuration() / fmcw.getChirpPeriod();
			const RealType average_power = waveform.getPower() * duty_cycle;
			const std::string direction(fers_signal::fmcwChirpDirectionToken(fmcw.getDirection()));
			const auto configured_count = fmcwCountToken(fmcw.getChirpCount());
			if (transmitter.getSchedule().empty())
			{
				logUnscheduledFmcwChirpSummary(transmitter, fmcw, direction, configured_count, duty_cycle,
											   average_power);
				return;
			}
			logScheduledFmcwChirpSummary(transmitter, fmcw, direction, configured_count, duty_cycle, average_power);
		}

		void logUnscheduledFmcwTriangleSummary(const Transmitter& transmitter,
											   const fers_signal::FmcwTriangleSignal& triangle,
											   const std::string& configured_count, const RealType average_power)
		{
			const RealType active_start = params::startTime();
			const auto source = makeActiveSource(&transmitter, active_start, params::endTime());
			const auto total_triangle_count = countFmcwTriangleStarts(source, active_start, source.segment_end);
			LOG(Level::INFO,
				"FMCW transmitter '{}' shape=triangle B={} Hz T_c={} s T_tri={} s f_0={} Hz alpha={} Hz/s "
				"duty_cycle=1 triangle_count={} total_triangle_count={} average_power={} W",
				transmitter.getName(), triangle.getChirpBandwidth(), triangle.getChirpDuration(),
				triangle.getTrianglePeriod(), triangle.getStartFrequencyOffset(), triangle.getChirpRate(),
				configured_count, total_triangle_count, average_power);
		}

		void logScheduledFmcwTriangleSummary(const Transmitter& transmitter,
											 const fers_signal::FmcwTriangleSignal& triangle,
											 const std::string& configured_count, const RealType average_power)
		{
			std::uint64_t total_triangle_count = 0;
			for (const auto& period : transmitter.getSchedule())
			{
				const RealType active_start = std::max(params::startTime(), period.start);
				const auto source =
					makeActiveSource(&transmitter, period.start, std::min(params::endTime(), period.end));
				const auto segment_triangle_count = countFmcwTriangleStarts(source, active_start, source.segment_end);
				total_triangle_count += segment_triangle_count;
				LOG(Level::INFO,
					"FMCW transmitter '{}' segment [{}, {}] shape=triangle B={} Hz T_c={} s T_tri={} s f_0={} "
					"Hz alpha={} Hz/s duty_cycle=1 triangle_count={} segment_triangle_count={} "
					"total_triangle_count={} average_power={} W",
					transmitter.getName(), period.start, source.segment_end, triangle.getChirpBandwidth(),
					triangle.getChirpDuration(), triangle.getTrianglePeriod(), triangle.getStartFrequencyOffset(),
					triangle.getChirpRate(), configured_count, segment_triangle_count, total_triangle_count,
					average_power);
			}
		}

		void logFmcwTriangleSummary(const Transmitter& transmitter, const fers_signal::RadarSignal& waveform,
									const fers_signal::FmcwTriangleSignal& triangle)
		{
			const RealType average_power = waveform.getPower();
			const auto configured_count = fmcwCountToken(triangle.getTriangleCount());
			if (transmitter.getSchedule().empty())
			{
				logUnscheduledFmcwTriangleSummary(transmitter, triangle, configured_count, average_power);
				return;
			}
			logScheduledFmcwTriangleSummary(transmitter, triangle, configured_count, average_power);
		}

		void logUnscheduledSfcwSummary(const Transmitter& transmitter, const fers_signal::RadarSignal& waveform,
									   const fers_signal::SteppedFrequencySignal& sfcw,
									   const std::string& configured_count, const RealType duty_cycle,
									   const RealType average_power)
		{
			const RealType active_start = params::startTime();
			const auto source = makeActiveSource(&transmitter, active_start, params::endTime());
			const auto total_step_count = countSfcwStepStarts(source, active_start, source.segment_end);
			LOG(Level::INFO,
				"SFCW transmitter '{}' steps={} df={} Hz dwell={} s step_period={} s sweep_period={} s "
				"f_first={} Hz f_last={} Hz B_eff={} Hz range_resolution={} m unambiguous_range={} m duty_cycle={} "
				"sweep_count={} total_step_count={} average_power={} W",
				transmitter.getName(), sfcw.getStepCount(), sfcw.getStepSize(), sfcw.getDwellTime(),
				sfcw.getStepPeriod(), sfcw.getSweepPeriod(), sfcw.firstFrequency(waveform.getCarrier()),
				sfcw.lastFrequency(waveform.getCarrier()), sfcw.effectiveBandwidth(),
				params::c() / (2.0 * sfcw.effectiveBandwidth()), params::c() / (2.0 * std::abs(sfcw.getStepSize())),
				duty_cycle, configured_count, total_step_count, average_power);
		}

		void logScheduledSfcwSummary(const Transmitter& transmitter, const fers_signal::RadarSignal& waveform,
									 const fers_signal::SteppedFrequencySignal& sfcw,
									 const std::string& configured_count, const RealType duty_cycle,
									 const RealType average_power)
		{
			std::uint64_t total_step_count = 0;
			for (const auto& period : transmitter.getSchedule())
			{
				const RealType active_start = std::max(params::startTime(), period.start);
				const auto source =
					makeActiveSource(&transmitter, period.start, std::min(params::endTime(), period.end));
				const auto segment_step_count = countSfcwStepStarts(source, active_start, source.segment_end);
				total_step_count += segment_step_count;
				LOG(Level::INFO,
					"SFCW transmitter '{}' segment [{}, {}] steps={} df={} Hz dwell={} s step_period={} s "
					"sweep_period={} s f_first={} Hz f_last={} Hz B_eff={} Hz range_resolution={} m "
					"unambiguous_range={} m duty_cycle={} sweep_count={} segment_step_count={} total_step_count={} "
					"average_power={} W",
					transmitter.getName(), period.start, source.segment_end, sfcw.getStepCount(), sfcw.getStepSize(),
					sfcw.getDwellTime(), sfcw.getStepPeriod(), sfcw.getSweepPeriod(),
					sfcw.firstFrequency(waveform.getCarrier()), sfcw.lastFrequency(waveform.getCarrier()),
					sfcw.effectiveBandwidth(), params::c() / (2.0 * sfcw.effectiveBandwidth()),
					params::c() / (2.0 * std::abs(sfcw.getStepSize())), duty_cycle, configured_count,
					segment_step_count, total_step_count, average_power);
			}
		}

		void logSfcwSummary(const Transmitter& transmitter, const fers_signal::RadarSignal& waveform,
							const fers_signal::SteppedFrequencySignal& sfcw)
		{
			const RealType duty_cycle = sfcw.getDwellTime() / sfcw.getStepPeriod();
			const RealType average_power = waveform.getPower() * duty_cycle;
			const auto configured_count = fmcwCountToken(sfcw.getSweepCount());
			if (transmitter.getSchedule().empty())
			{
				logUnscheduledSfcwSummary(transmitter, waveform, sfcw, configured_count, duty_cycle, average_power);
				return;
			}
			logScheduledSfcwSummary(transmitter, waveform, sfcw, configured_count, duty_cycle, average_power);
		}

		[[nodiscard]] bool isStreamingReceiver(const Receiver* const receiver) noexcept
		{
			return receiver != nullptr &&
				(receiver->getMode() == OperationMode::CW_MODE || receiver->getMode() == OperationMode::FMCW_MODE ||
				 receiver->getMode() == OperationMode::SFCW_MODE);
		}

		[[nodiscard]] bool activePastUserEnd(const Receiver* const receiver) noexcept
		{
			if (receiver == nullptr)
			{
				return false;
			}
			if (receiver->getSchedule().empty())
			{
				return true;
			}
			return std::ranges::any_of(receiver->getSchedule(),
									   [](const auto& period) { return period.end > params::endTime(); });
		}

		[[nodiscard]] std::size_t streamingSampleIndexAtOrAfter(const RealType time, const RealType dt_sim)
		{
			if (dt_sim <= 0.0 || time <= params::startTime())
			{
				return 0;
			}
			return static_cast<std::size_t>(std::ceil((time - params::startTime()) / dt_sim));
		}

		struct PositionBounds
		{
			math::Vec3 min{};
			math::Vec3 max{};
			bool valid{false};
			bool unbounded{false};
		};

		[[nodiscard]] bool isFinite(const math::Vec3& point) noexcept
		{
			return std::isfinite(point.x) && std::isfinite(point.y) && std::isfinite(point.z);
		}

		void includePoint(PositionBounds& bounds, const math::Vec3& point) noexcept
		{
			if (!isFinite(point))
			{
				bounds.unbounded = true;
				return;
			}
			if (!bounds.valid)
			{
				bounds.min = point;
				bounds.max = point;
				bounds.valid = true;
				return;
			}
			bounds.min.x = std::min(bounds.min.x, point.x);
			bounds.min.y = std::min(bounds.min.y, point.y);
			bounds.min.z = std::min(bounds.min.z, point.z);
			bounds.max.x = std::max(bounds.max.x, point.x);
			bounds.max.y = std::max(bounds.max.y, point.y);
			bounds.max.z = std::max(bounds.max.z, point.z);
		}

		[[nodiscard]] RealType axisValue(const math::Vec3& point, const std::size_t axis) noexcept
		{
			switch (axis)
			{
			case 0:
				return point.x;
			case 1:
				return point.y;
			default:
				return point.z;
			}
		}

		[[nodiscard]] RealType axisValue(const std::array<RealType, 3>& values, const std::size_t axis) noexcept
		{
			switch (axis)
			{
			case 0:
				return values[0];
			case 1:
				return values[1];
			default:
				return values[2];
			}
		}

		[[nodiscard]] RealType& axisValue(std::array<RealType, 3>& values, const std::size_t axis) noexcept
		{
			switch (axis)
			{
			case 0:
				return values[0];
			case 1:
				return values[1];
			default:
				return values[2];
			}
		}

		[[nodiscard]] RealType axisDistanceBound(const PositionBounds& lhs, const PositionBounds& rhs,
												 const std::size_t axis) noexcept
		{
			const RealType lhs_min = axisValue(lhs.min, axis);
			const RealType lhs_max = axisValue(lhs.max, axis);
			const RealType rhs_min = axisValue(rhs.min, axis);
			const RealType rhs_max = axisValue(rhs.max, axis);
			return std::max(std::abs(lhs_max - rhs_min), std::abs(rhs_max - lhs_min));
		}

		[[nodiscard]] RealType maxDistanceBetweenBounds(const PositionBounds& lhs, const PositionBounds& rhs) noexcept
		{
			if (lhs.unbounded || rhs.unbounded || !lhs.valid || !rhs.valid)
			{
				return std::numeric_limits<RealType>::infinity();
			}
			const RealType dx = axisDistanceBound(lhs, rhs, 0);
			const RealType dy = axisDistanceBound(lhs, rhs, 1);
			const RealType dz = axisDistanceBound(lhs, rhs, 2);
			return std::sqrt(dx * dx + dy * dy + dz * dz);
		}

		[[nodiscard]] std::array<RealType, 3> coordinateAxes(const math::Coord& coord) noexcept
		{
			return {coord.pos.x, coord.pos.y, coord.pos.z};
		}

		void includeCubicVelocityRoot(PositionBounds& bounds, const math::Path& path, const RealType segment_start,
									  const RealType segment_length, const RealType root_u, const RealType lower_u,
									  const RealType upper_u)
		{
			if (root_u < lower_u || root_u > upper_u)
			{
				return;
			}
			includePoint(bounds, path.getPosition(segment_start + root_u * segment_length));
		}

		void includeCubicPositionExtrema(PositionBounds& bounds, const math::Path& path,
										 const std::vector<math::Coord>& coords,
										 const std::vector<math::Coord>& second_derivatives, const std::size_t index,
										 const RealType lower_u, const RealType upper_u)
		{
			const RealType segment_length = coords[index + 1].t - coords[index].t;
			if (segment_length <= EPSILON)
			{
				return;
			}
			const auto left = coordinateAxes(coords[index]);
			const auto right = coordinateAxes(coords[index + 1]);
			const auto dd_left = coordinateAxes(second_derivatives[index]);
			const auto dd_right = coordinateAxes(second_derivatives[index + 1]);
			const RealType h2 = segment_length * segment_length;

			for (std::size_t axis = 0; axis < 3; ++axis)
			{
				const RealType dd_right_axis = axisValue(dd_right, axis);
				const RealType dd_left_axis = axisValue(dd_left, axis);
				const RealType a = 0.5 * h2 * (dd_right_axis - dd_left_axis);
				const RealType b = h2 * dd_left_axis;
				const RealType c = (axisValue(right, axis) - axisValue(left, axis)) +
					(h2 / 6.0) * (-2.0 * dd_left_axis - dd_right_axis);

				if (std::abs(a) <= EPSILON)
				{
					if (std::abs(b) > EPSILON)
					{
						includeCubicVelocityRoot(bounds, path, coords[index].t, segment_length, -c / b, lower_u,
												 upper_u);
					}
					continue;
				}

				const RealType discriminant = b * b - 4.0 * a * c;
				if (discriminant < -EPSILON)
				{
					continue;
				}
				const RealType sqrt_discriminant = std::sqrt(std::max(0.0, discriminant));
				includeCubicVelocityRoot(bounds, path, coords[index].t, segment_length,
										 (-b - sqrt_discriminant) / (2.0 * a), lower_u, upper_u);
				includeCubicVelocityRoot(bounds, path, coords[index].t, segment_length,
										 (-b + sqrt_discriminant) / (2.0 * a), lower_u, upper_u);
			}
		}

		[[nodiscard]] PositionBounds pathPositionBounds(const math::Path& path, const RealType start,
														const RealType end)
		{
			PositionBounds bounds;
			if (start >= end)
			{
				return bounds;
			}

			try
			{
				includePoint(bounds, path.getPosition(start));
				includePoint(bounds, path.getPosition(end));
			}
			catch (const math::PathException&)
			{
				bounds.unbounded = true;
				return bounds;
			}

			const auto& coords = path.getCoords();
			if (coords.empty() || path.getType() == math::Path::InterpType::INTERP_STATIC)
			{
				return bounds;
			}

			for (const auto& coord : coords)
			{
				if (coord.t >= start && coord.t <= end)
				{
					includePoint(bounds, coord.pos);
				}
			}

			if (path.getType() != math::Path::InterpType::INTERP_CUBIC || coords.size() < 2)
			{
				return bounds;
			}

			std::vector<math::Coord> second_derivatives;
			try
			{
				finalizeCubic<math::Coord>(coords, second_derivatives);
			}
			catch (const math::PathException&)
			{
				bounds.unbounded = true;
				return bounds;
			}

			for (std::size_t index = 0; index + 1 < coords.size(); ++index)
			{
				const RealType segment_start = coords[index].t;
				const RealType segment_end = coords[index + 1].t;
				const RealType segment_length = segment_end - segment_start;
				if (segment_length <= EPSILON || end < segment_start || start > segment_end)
				{
					continue;
				}

				const RealType lower_u =
					std::clamp((std::max(start, segment_start) - segment_start) / segment_length, 0.0, 1.0);
				const RealType upper_u =
					std::clamp((std::min(end, segment_end) - segment_start) / segment_length, 0.0, 1.0);
				if (lower_u <= upper_u)
				{
					includeCubicPositionExtrema(bounds, path, coords, second_derivatives, index, lower_u, upper_u);
				}
			}
			return bounds;
		}

		struct QuadraticVelocityExtremum
		{
			RealType a;
			RealType b;
			RealType c;
			RealType segment_length;
			RealType root_u;
			RealType lower_u;
			RealType upper_u;
		};

		void includeQuadraticVelocityExtremum(std::array<RealType, 3>& max_abs_velocity, const std::size_t axis,
											  const QuadraticVelocityExtremum& extremum) noexcept
		{
			if (extremum.root_u < extremum.lower_u || extremum.root_u > extremum.upper_u ||
				extremum.segment_length <= EPSILON)
			{
				return;
			}
			const RealType velocity =
				(extremum.a * extremum.root_u * extremum.root_u + extremum.b * extremum.root_u + extremum.c) /
				extremum.segment_length;
			if (std::isfinite(velocity))
			{
				RealType& axis_max_velocity = axisValue(max_abs_velocity, axis);
				axis_max_velocity = std::max(axis_max_velocity, std::abs(velocity));
			}
			else
			{
				axisValue(max_abs_velocity, axis) = std::numeric_limits<RealType>::infinity();
			}
		}

		void includeCubicVelocityBounds(std::array<RealType, 3>& max_abs_velocity,
										const std::vector<math::Coord>& coords,
										const std::vector<math::Coord>& second_derivatives, const std::size_t index,
										const RealType lower_u, const RealType upper_u)
		{
			const RealType segment_length = coords[index + 1].t - coords[index].t;
			if (segment_length <= EPSILON)
			{
				return;
			}
			const auto left = coordinateAxes(coords[index]);
			const auto right = coordinateAxes(coords[index + 1]);
			const auto dd_left = coordinateAxes(second_derivatives[index]);
			const auto dd_right = coordinateAxes(second_derivatives[index + 1]);
			const RealType h2 = segment_length * segment_length;

			for (std::size_t axis = 0; axis < 3; ++axis)
			{
				const RealType dd_right_axis = axisValue(dd_right, axis);
				const RealType dd_left_axis = axisValue(dd_left, axis);
				const RealType a = 0.5 * h2 * (dd_right_axis - dd_left_axis);
				const RealType b = h2 * dd_left_axis;
				const RealType c = (axisValue(right, axis) - axisValue(left, axis)) +
					(h2 / 6.0) * (-2.0 * dd_left_axis - dd_right_axis);
				includeQuadraticVelocityExtremum(max_abs_velocity, axis,
												 QuadraticVelocityExtremum{.a = a,
																		   .b = b,
																		   .c = c,
																		   .segment_length = segment_length,
																		   .root_u = lower_u,
																		   .lower_u = lower_u,
																		   .upper_u = upper_u});
				includeQuadraticVelocityExtremum(max_abs_velocity, axis,
												 QuadraticVelocityExtremum{.a = a,
																		   .b = b,
																		   .c = c,
																		   .segment_length = segment_length,
																		   .root_u = upper_u,
																		   .lower_u = lower_u,
																		   .upper_u = upper_u});

				if (std::abs(a) > EPSILON)
				{
					includeQuadraticVelocityExtremum(max_abs_velocity, axis,
													 QuadraticVelocityExtremum{.a = a,
																			   .b = b,
																			   .c = c,
																			   .segment_length = segment_length,
																			   .root_u = -b / (2.0 * a),
																			   .lower_u = lower_u,
																			   .upper_u = upper_u});
				}
			}
		}

		[[nodiscard]] RealType pathSpeedBound(const math::Path& path, const RealType start, const RealType end)
		{
			if (start >= end)
			{
				return 0.0;
			}

			const auto& coords = path.getCoords();
			if (coords.empty() || path.getType() == math::Path::InterpType::INTERP_STATIC || coords.size() < 2)
			{
				return 0.0;
			}

			if (path.getType() == math::Path::InterpType::INTERP_LINEAR)
			{
				RealType max_speed = 0.0;
				for (std::size_t index = 0; index + 1 < coords.size(); ++index)
				{
					const RealType segment_start = coords[index].t;
					const RealType segment_end = coords[index + 1].t;
					const RealType segment_length = segment_end - segment_start;
					if (segment_length <= EPSILON || end < segment_start || start > segment_end)
					{
						continue;
					}
					max_speed =
						std::max(max_speed, (coords[index + 1].pos - coords[index].pos).length() / segment_length);
				}
				return max_speed;
			}

			std::vector<math::Coord> second_derivatives;
			try
			{
				finalizeCubic<math::Coord>(coords, second_derivatives);
			}
			catch (const math::PathException&)
			{
				return std::numeric_limits<RealType>::infinity();
			}

			std::array<RealType, 3> max_abs_velocity{0.0, 0.0, 0.0};
			for (std::size_t index = 0; index + 1 < coords.size(); ++index)
			{
				const RealType segment_start = coords[index].t;
				const RealType segment_end = coords[index + 1].t;
				const RealType segment_length = segment_end - segment_start;
				if (segment_length <= EPSILON || end < segment_start || start > segment_end)
				{
					continue;
				}
				const RealType lower_u =
					std::clamp((std::max(start, segment_start) - segment_start) / segment_length, 0.0, 1.0);
				const RealType upper_u =
					std::clamp((std::min(end, segment_end) - segment_start) / segment_length, 0.0, 1.0);
				if (lower_u <= upper_u)
				{
					includeCubicVelocityBounds(max_abs_velocity, coords, second_derivatives, index, lower_u, upper_u);
				}
			}
			return std::sqrt(max_abs_velocity[0] * max_abs_velocity[0] + max_abs_velocity[1] * max_abs_velocity[1] +
							 max_abs_velocity[2] * max_abs_velocity[2]);
		}

		[[nodiscard]] std::optional<RealType>
		deadlineFromTailKinematics(const RealType tail_end, const RealType interval_start, const RealType interval_end,
								   const RealType delay_at_start, const RealType distance_rate_bound,
								   const RealType max_delay_bound)
		{
			const RealType propagation_speed = params::c();
			if (propagation_speed <= 0.0 || interval_start >= interval_end)
			{
				return std::nullopt;
			}

			if (std::isfinite(distance_rate_bound) && distance_rate_bound < propagation_speed)
			{
				const RealType retarded_start = interval_start - delay_at_start;
				if (retarded_start >= tail_end)
				{
					return std::nullopt;
				}

				const RealType min_retarded_slope = 1.0 - distance_rate_bound / propagation_speed;
				const RealType deadline = interval_start + (tail_end - retarded_start) / min_retarded_slope;
				if (deadline <= interval_start)
				{
					return std::nullopt;
				}
				return std::min(interval_end, deadline);
			}

			if (!std::isfinite(max_delay_bound))
			{
				return interval_end;
			}
			const RealType deadline = std::min(interval_end, tail_end + max_delay_bound);
			if (deadline <= interval_start)
			{
				return std::nullopt;
			}
			return deadline;
		}

		[[nodiscard]] std::optional<RealType> directPathCleanupDeadline(const ActiveStreamingSource& source,
																		const Receiver* const rx,
																		const RealType interval_start,
																		const RealType interval_end)
		{
			const auto* const tx = source.transmitter;
			if (tx == nullptr || rx == nullptr || tx->getPlatform() == rx->getPlatform() || params::c() <= 0.0)
			{
				return std::nullopt;
			}

			const auto* const tx_path = tx->getPlatform()->getMotionPath();
			const auto* const rx_path = rx->getPlatform()->getMotionPath();
			const RealType distance_at_start =
				(rx_path->getPosition(interval_start) - tx_path->getPosition(interval_start)).length();
			const RealType max_delay_bound =
				maxDistanceBetweenBounds(pathPositionBounds(*tx_path, interval_start, interval_end),
										 pathPositionBounds(*rx_path, interval_start, interval_end)) /
				params::c();
			const RealType distance_rate_bound = pathSpeedBound(*tx_path, interval_start, interval_end) +
				pathSpeedBound(*rx_path, interval_start, interval_end);
			return deadlineFromTailKinematics(source.segment_end, interval_start, interval_end,
											  distance_at_start / params::c(), distance_rate_bound, max_delay_bound);
		}

		[[nodiscard]] std::optional<RealType> reflectedPathCleanupDeadline(const ActiveStreamingSource& source,
																		   const Receiver* const rx,
																		   const radar::Target* const target,
																		   const RealType interval_start,
																		   const RealType interval_end)
		{
			const auto* const tx = source.transmitter;
			if (tx == nullptr || rx == nullptr || target == nullptr || params::c() <= 0.0 ||
				tx->getPlatform() == target->getPlatform() || rx->getPlatform() == target->getPlatform())
			{
				return std::nullopt;
			}

			const auto* const tx_path = tx->getPlatform()->getMotionPath();
			const auto* const rx_path = rx->getPlatform()->getMotionPath();
			const auto* const target_path = target->getPlatform()->getMotionPath();
			const auto tx_position = tx_path->getPosition(interval_start);
			const auto rx_position = rx_path->getPosition(interval_start);
			const auto target_position = target_path->getPosition(interval_start);
			const RealType distance_at_start =
				(target_position - tx_position).length() + (rx_position - target_position).length();

			const PositionBounds tx_bounds = pathPositionBounds(*tx_path, interval_start, interval_end);
			const PositionBounds rx_bounds = pathPositionBounds(*rx_path, interval_start, interval_end);
			const PositionBounds target_bounds = pathPositionBounds(*target_path, interval_start, interval_end);
			const RealType max_delay_bound = (maxDistanceBetweenBounds(tx_bounds, target_bounds) +
											  maxDistanceBetweenBounds(target_bounds, rx_bounds)) /
				params::c();
			const RealType tx_speed = pathSpeedBound(*tx_path, interval_start, interval_end);
			const RealType rx_speed = pathSpeedBound(*rx_path, interval_start, interval_end);
			const RealType target_speed = pathSpeedBound(*target_path, interval_start, interval_end);
			const RealType distance_rate_bound = tx_speed + rx_speed + 2.0 * target_speed;

			return deadlineFromTailKinematics(source.segment_end, interval_start, interval_end,
											  distance_at_start / params::c(), distance_rate_bound, max_delay_bound);
		}

		/// Builds an active streaming source for a transmitter at an event timestamp.
		std::optional<ActiveStreamingSource> streamingSourceAtEvent(const Transmitter* const transmitter,
																	const RealType timestamp,
																	const RealType internal_stop_time)
		{
			if (transmitter == nullptr || !transmitter->isStreamingMode())
			{
				return std::nullopt;
			}

			const auto& schedule = transmitter->getSchedule();
			if (schedule.empty())
			{
				const RealType segment_start = params::startTime();
				auto source = makeActiveSource(transmitter, segment_start, internal_stop_time);
				if (timestamp >= segment_start && timestamp < source.segment_end)
				{
					return source;
				}
				return std::nullopt;
			}

			// TODO: O(N) Schedule Lookups - Since the schedule is guaranteed to be sorted (enforced by
			// `processRawSchedule`), should be using `std::lower_bound` or binary search to find the relevant period in
			// $O(\log N)$ time.
			for (const auto& period : schedule)
			{
				const RealType active_start = std::max(params::startTime(), period.start);
				auto source = makeActiveSource(transmitter, period.start, std::min(internal_stop_time, period.end));
				if (timestamp >= active_start && timestamp < source.segment_end)
				{
					return source;
				}
			}
			return std::nullopt;
		}
	}

	SimulationEngine::SimulationEngine(World* world, pool::ThreadPool& pool, std::shared_ptr<ProgressReporter> reporter,
									   std::string output_dir,
									   std::shared_ptr<OutputMetadataCollector> metadata_collector,
									   ReceiverOutputSink* output_sink, std::function<bool()> cancel_callback,
									   const bool eager_context_stream_open) :
		_world(world), _pool(pool), _reporter(std::move(reporter)), _metadata_collector(std::move(metadata_collector)),
		_output_sink(output_sink), _cancel_callback(std::move(cancel_callback)),
		_eager_context_stream_open(eager_context_stream_open), _last_report_time(std::chrono::steady_clock::now()),
		_next_context_heartbeat_time(params::startTime() + 1.0), _output_dir(std::move(output_dir)),
		_internal_stop_time(params::endTime())
	{
		_streaming_tracker_caches.resize(_world->getReceivers().size());
		_if_pulse_tracker_caches.resize(_world->getReceivers().size());
		_fmcw_if_block_buffers.resize(_world->getReceivers().size());
		_fmcw_if_block_start_times.resize(_world->getReceivers().size(), params::startTime());
		_streaming_output_block_buffers.resize(_world->getReceivers().size());
		_streaming_output_processed_buffers.resize(_world->getReceivers().size());
		_streaming_output_block_start_times.resize(_world->getReceivers().size(), params::startTime());
		_streaming_output_block_start_indices.resize(_world->getReceivers().size(), 0);
		_streaming_downsamplers.resize(_world->getReceivers().size());
		_streaming_downsample_base_indices.resize(_world->getReceivers().size(), 0);
		_streaming_downsample_segment_start_times.resize(_world->getReceivers().size(), params::startTime());
		_streaming_output_sample_cursors.resize(_world->getReceivers().size(), 0);
		_streaming_output_stream_ids.resize(_world->getReceivers().size(), 0);
		_streaming_output_stream_open.resize(_world->getReceivers().size(), false);
		_streaming_output_file_metadata.resize(_world->getReceivers().size());
		for (auto& block : _streaming_output_block_buffers)
		{
			block.reserve(streaming_output_block_size);
		}
		for (auto& block : _streaming_output_processed_buffers)
		{
			block.reserve(streaming_output_block_size);
		}
	}

	void SimulationEngine::run()
	{
		if (_reporter)
		{
			_reporter->report("Initializing event-driven simulation...", 0, 100);
		}

		initializeFmcwIfResamplers();

		logSimulationMemoryProjection(*_world);

		initializeFinalizers();

		LOG(Level::INFO, "Starting unified event-driven simulation loop.");
		logStreamingSummaries();

		auto& event_queue = _world->getEventQueue();
		auto& state = _world->getSimulationState();
		const RealType end_time = _internal_stop_time;

		while (!event_queue.empty() && state.t_current <= end_time)
		{
			if (isCancellationRequested())
			{
				break;
			}
			const Event event = event_queue.top();
			event_queue.pop();

			processStreamingPhysics(event.timestamp);
			if (isCancellationRequested())
			{
				break;
			}
			flushFmcwIfBlocks();
			flushStreamingOutputBlocks();

			state.t_current = event.timestamp;

			processEvent(event);
			updateProgress();
		}

		const bool has_active_if_overrender =
			std::ranges::any_of(_world->getReceivers(), [](const auto& receiver)
								{ return receiver->isActive() && receiver->hasFmcwIfResamplingSink(); });
		if (!isCancellationRequested() && has_active_if_overrender)
		{
			processStreamingPhysics(end_time);
		}
		flushFmcwIfBlocks();
		flushStreamingOutputBlocks();

		shutdown();
	}

	void SimulationEngine::logStreamingSummaries() const
	{
		for (const auto& transmitter_ptr : _world->getTransmitters())
		{
			const auto* waveform = transmitter_ptr->getSignal();
			if (waveform == nullptr)
			{
				continue;
			}

			if (const auto* fmcw = waveform->getFmcwChirpSignal(); fmcw != nullptr)
			{
				logFmcwChirpSummary(*transmitter_ptr, *waveform, *fmcw);
			}
			else if (const auto* triangle = waveform->getFmcwTriangleSignal(); triangle != nullptr)
			{
				logFmcwTriangleSummary(*transmitter_ptr, *waveform, *triangle);
			}
			else if (const auto* sfcw = waveform->getSteppedFrequencySignal(); sfcw != nullptr)
			{
				logSfcwSummary(*transmitter_ptr, *waveform, *sfcw);
			}
		}
	}

	void SimulationEngine::initializeFinalizers()
	{
		if (_output_sink == nullptr)
		{
			return;
		}
		for (const auto& receiver_ptr : _world->getReceivers())
		{
			if (receiver_ptr->getMode() == OperationMode::PULSED_MODE)
			{
				_finalizer_threads.emplace_back(processing::runPulsedFinalizer, receiver_ptr.get(),
												&_world->getTargets(), _reporter, _output_dir, _metadata_collector,
												_output_sink);
			}
		}
	}

	void SimulationEngine::initializeFmcwIfResamplers()
	{
		_internal_stop_time = params::endTime();
		for (std::size_t receiver_index = 0; receiver_index < _world->getReceivers().size(); ++receiver_index)
		{
			initializeFmcwIfResampler(receiver_index);
		}

		if (_internal_stop_time > params::endTime())
		{
			extendDechirpSourcesForIfOverrender();
		}
	}

	void SimulationEngine::initializeFmcwIfResampler(const std::size_t receiver_index)
	{
		const auto& receiver_ptr = _world->getReceivers()[receiver_index];
		if (!receiver_ptr->isDechirpEnabled() || !receiver_ptr->hasFmcwIfSampleRate())
		{
			return;
		}

		const auto& request = receiver_ptr->getFmcwIfChainRequest();
		const RealType output_rate = request.sample_rate_hz.value_or(0.0);
		const RealType bandwidth = request.filter_bandwidth_hz.value_or(0.40 * output_rate);
		const fers_signal::FmcwIfResamplerRequest resampler_request{
			.input_sample_rate_hz = params::rate() * static_cast<RealType>(params::oversampleRatio()),
			.output_sample_rate_hz = output_rate,
			.filter_bandwidth_hz = bandwidth,
			.filter_transition_width_hz = request.filter_transition_width_hz};
		auto plan = fers_signal::planFmcwIfResampler(resampler_request);
		const RealType block_time = static_cast<RealType>(fmcw_if_block_size) / resampler_request.input_sample_rate_hz;
		const RealType over_render = plan.group_delay_seconds + 1.0 / plan.actual_output_sample_rate_hz + block_time;
		_internal_stop_time = std::max(_internal_stop_time, params::endTime() + over_render);
		if (_output_sink != nullptr)
		{
			receiver_ptr->setFmcwIfOutputCallback(
				[this, receiver_index](const std::span<const ComplexType> samples, const std::uint64_t sample_start)
				{
					const auto& receiver = _world->getReceivers()[receiver_index];
					const auto& if_plan = receiver->getFmcwIfResamplerPlan();
					if (!if_plan.has_value())
					{
						return;
					}
					const RealType first_sample_time = params::startTime() +
						static_cast<RealType>(sample_start) / if_plan->actual_output_sample_rate_hz;
					emitStreamingOutputBlock(receiver_index, first_sample_time, if_plan->actual_output_sample_rate_hz,
											 samples, sample_start);
				});
		}
		const RealType actual_output_sample_rate_hz = plan.actual_output_sample_rate_hz;
		const auto overall_ratio = plan.overall_ratio;
		const RealType filter_bandwidth_hz = plan.filter_bandwidth_hz;
		const RealType filter_transition_width_hz = plan.filter_transition_width_hz;
		receiver_ptr->initializeFmcwIfResampling(std::move(plan));
		LOG(Level::INFO,
			"Receiver '{}' enabled FMCW IF resampling: input_rate={} Hz requested_output_rate={} Hz "
			"actual_output_rate={} Hz ratio={}/{} passband={} Hz transition={} Hz.",
			receiver_ptr->getName(), resampler_request.input_sample_rate_hz, output_rate, actual_output_sample_rate_hz,
			overall_ratio.numerator, overall_ratio.denominator, filter_bandwidth_hz, filter_transition_width_hz);
	}

	void SimulationEngine::extendDechirpSourcesForIfOverrender()
	{
		for (const auto& receiver_ptr : _world->getReceivers())
		{
			if (!receiver_ptr->hasFmcwIfSampleRate())
			{
				continue;
			}

			auto dechirp_sources = receiver_ptr->getDechirpSources();
			for (auto& source : dechirp_sources)
			{
				if (std::abs(source.segment_end - params::endTime()) > 1.0e-12)
				{
					continue;
				}
				if (source.transmitter == nullptr || source.transmitter->getSchedule().empty())
				{
					source.segment_end = _internal_stop_time;
					continue;
				}
				for (const auto& period : source.transmitter->getSchedule())
				{
					if (period.start <= params::endTime() && period.end > params::endTime())
					{
						source.segment_end = std::min(_internal_stop_time, period.end);
						break;
					}
				}
			}
			receiver_ptr->setResolvedDechirpSources(std::move(dechirp_sources));
		}
	}

	void SimulationEngine::ensureCwPhaseNoiseLookup()
	{
		if (_cw_phase_noise_lookup)
		{
			return;
		}

		const auto timings = collectCwPhaseNoiseTimings(*_world);
		RealType lookup_start = _world->earliestPhaseNoiseLookupStart();
		for (const auto& source : _world->getSimulationState().active_streaming_transmitters)
		{
			lookup_start = std::min(lookup_start, source.segment_start);
		}
		_cw_phase_noise_lookup = std::make_unique<simulation::CwPhaseNoiseLookup>(
			simulation::CwPhaseNoiseLookup::build(timings, lookup_start, _internal_stop_time));
	}

	void SimulationEngine::processStreamingPhysics(const RealType t_event)
	{
		auto& state = _world->getSimulationState();
		auto& t_current = state.t_current;

		if (t_event <= t_current)
		{
			return;
		}

		const RealType dt_sim = 1.0 / (params::rate() * params::oversampleRatio());
		const auto first_index = streamingSampleIndexAtOrAfter(t_current, dt_sim);
		const auto final_index = streamingSampleIndexAtOrAfter(t_event, dt_sim);
		const auto sample_count = final_index - first_index;
		const auto progress_report_stride = std::max<std::size_t>(1, sample_count / 1000);

		ensureCwPhaseNoiseLookup();

		while (t_current < t_event && !isCancellationRequested())
		{
			cleanupInactiveStreamingSources(t_current);

			const RealType chunk_end = streamingChunkEnd(t_current, t_event);
			if (chunk_end <= t_current)
			{
				break;
			}

			const auto start_index = streamingSampleIndexAtOrAfter(t_current, dt_sim);
			const auto end_index = streamingSampleIndexAtOrAfter(chunk_end, dt_sim);
			for (size_t sample_index = start_index; sample_index < end_index; ++sample_index)
			{
				if (shouldStopStreamingChunk(sample_index, start_index))
				{
					break;
				}
				processStreamingSample(sample_index, first_index, final_index, progress_report_stride, dt_sim);
			}

			t_current = chunk_end;
			emitContextHeartbeatsThrough(t_current);
		}
		cleanupInactiveStreamingSources(t_current);
	}

	std::optional<RealType> SimulationEngine::nextStreamingCleanupDeadline(const RealType from_time)
	{
		const auto& active_streaming_transmitters = _world->getSimulationState().active_streaming_transmitters;
		std::optional<RealType> next_deadline;
		for (const auto& source : active_streaming_transmitters)
		{
			if (source.segment_end > from_time)
			{
				continue;
			}
			const auto cleanup_deadline = streamingSourceCleanupDeadline(source, from_time);
			if (cleanup_deadline.has_value() && *cleanup_deadline > from_time &&
				(!next_deadline.has_value() || *cleanup_deadline < *next_deadline))
			{
				next_deadline = cleanup_deadline;
			}
		}
		return next_deadline;
	}

	RealType SimulationEngine::streamingChunkEnd(const RealType from_time, const RealType event_time)
	{
		if (const auto cleanup_deadline = nextStreamingCleanupDeadline(from_time);
			cleanup_deadline.has_value() && *cleanup_deadline < event_time)
		{
			return *cleanup_deadline;
		}
		return event_time;
	}

	bool SimulationEngine::shouldStopStreamingChunk(const std::size_t sample_index, const std::size_t chunk_start_index)
	{
		return ((sample_index - chunk_start_index) % 1024) == 0 && isCancellationRequested();
	}

	void SimulationEngine::processStreamingSample(const std::size_t sample_index, const std::size_t first_index,
												  const std::size_t final_index,
												  const std::size_t progress_report_stride, const RealType dt_sim)
	{
		const RealType t_step = params::startTime() + static_cast<RealType>(sample_index) * dt_sim;
		appendActiveReceiverStreamingSamples(sample_index, t_step);

		if (_output_sink != nullptr && t_step >= _next_context_heartbeat_time)
		{
			emitContextHeartbeatsThrough(t_step);
		}
		if (((sample_index - first_index) % progress_report_stride) == 0 || sample_index + 1 == final_index)
		{
			reportSimulationProgress(t_step);
		}
	}

	void SimulationEngine::appendActiveReceiverStreamingSamples(const std::size_t sample_index, const RealType t_step)
	{
		for (std::size_t receiver_index = 0; receiver_index < _world->getReceivers().size(); ++receiver_index)
		{
			appendReceiverStreamingSample(receiver_index, sample_index, t_step);
		}
	}

	void SimulationEngine::appendReceiverStreamingSample(const std::size_t receiver_index,
														 const std::size_t sample_index, const RealType t_step)
	{
		const auto& receiver_ptr = _world->getReceivers()[receiver_index];
		if (!isStreamingReceiver(receiver_ptr.get()) || !receiver_ptr->isActive())
		{
			return;
		}

		const auto& active_streaming_transmitters = _world->getSimulationState().active_streaming_transmitters;
		ComplexType const sample = calculateStreamingSample(receiver_ptr.get(), t_step, active_streaming_transmitters,
															_streaming_tracker_caches[receiver_index]);
		if (receiver_ptr->hasFmcwIfResamplingSink())
		{
			appendFmcwIfSample(receiver_index, t_step, sample);
		}
		else if (_output_sink != nullptr)
		{
			appendStreamingOutputSample(receiver_index, sample_index, t_step, sample);
		}
	}

	void SimulationEngine::appendFmcwIfSample(const std::size_t receiver_index, const RealType t_step,
											  const ComplexType sample)
	{
		auto& block = _fmcw_if_block_buffers[receiver_index];
		if (block.empty())
		{
			_fmcw_if_block_start_times[receiver_index] = t_step;
		}
		block.push_back(sample);
		if (block.size() >= fmcw_if_block_size)
		{
			flushFmcwIfBlock(receiver_index);
		}
	}

	void SimulationEngine::appendStreamingOutputSample(const std::size_t receiver_index, const std::size_t sample_index,
													   const RealType t_step, const ComplexType sample)
	{
		if (_eager_context_stream_open)
		{
			ensureStreamingOutputStreamOpen(receiver_index, t_step, streamingOutputSampleRate(receiver_index));
		}
		auto& block = _streaming_output_block_buffers[receiver_index];
		if (block.empty())
		{
			_streaming_output_block_start_times[receiver_index] = t_step;
			_streaming_output_block_start_indices[receiver_index] = static_cast<std::uint64_t>(sample_index);
		}
		block.push_back(sample);
		if (block.size() >= streaming_output_block_size)
		{
			flushStreamingOutputBlock(receiver_index);
		}
	}

	void SimulationEngine::flushStreamingOutputBlocks()
	{
		for (std::size_t receiver_index = 0; receiver_index < _streaming_output_block_buffers.size(); ++receiver_index)
		{
			flushStreamingOutputBlock(receiver_index);
		}
	}

	void SimulationEngine::flushStreamingOutputBlock(const std::size_t receiver_index, const bool finish_downsampler)
	{
		if (_output_sink == nullptr || receiver_index >= _world->getReceivers().size())
		{
			return;
		}

		auto& block = _streaming_output_block_buffers[receiver_index];
		if (block.empty())
		{
			if (finish_downsampler && _streaming_downsamplers[receiver_index])
			{
				auto& downsampler = *_streaming_downsamplers[receiver_index];
				const auto output_start_index = downsampler.outputSampleCount();
				downsampler.finish();
				auto output = downsampler.takeOutput();
				if (!output.empty())
				{
					const RealType output_sample_rate = params::rate();
					const RealType output_start_time = _streaming_downsample_segment_start_times[receiver_index] +
						static_cast<RealType>(output_start_index) / output_sample_rate;
					emitStreamingOutputBlock(receiver_index, output_start_time, output_sample_rate, output,
											 _streaming_downsample_base_indices[receiver_index] + output_start_index);
				}
				_streaming_downsamplers[receiver_index].reset();
			}
			return;
		}

		const auto& receiver = _world->getReceivers()[receiver_index];
		const bool dechirped = receiver->isDechirpEnabled();
		const RealType input_sample_rate = params::rate() * static_cast<RealType>(params::oversampleRatio());
		const RealType block_start_time = _streaming_output_block_start_times[receiver_index];
		const auto input_start_index = _streaming_output_block_start_indices[receiver_index];

		applyPulsedInterferenceToStreamingBlock(receiver_index, block, block_start_time, input_sample_rate, dechirped);

		RealType output_sample_rate = input_sample_rate;
		RealType output_start_time = block_start_time;
		std::uint64_t output_sample_start = input_start_index;
		std::vector<ComplexType> downsampled_block;
		if (!dechirped && params::oversampleRatio() > 1)
		{
			auto& downsampler = streamingDownsampler(receiver_index, input_start_index, block_start_time);
			const auto output_start_index = downsampler.outputSampleCount();
			downsampler.consume(block);
			if (finish_downsampler)
			{
				downsampler.finish();
			}
			downsampled_block = downsampler.takeOutput();
			output_sample_rate = params::rate();
			output_sample_start = _streaming_downsample_base_indices[receiver_index] + output_start_index;
			output_start_time = _streaming_downsample_segment_start_times[receiver_index] +
				static_cast<RealType>(output_start_index) / output_sample_rate;
		}
		else if (!dechirped)
		{
			output_sample_rate = params::rate();
			output_sample_start = input_start_index;
			output_start_time = block_start_time;
		}

		const auto output_samples = !downsampled_block.empty()
			? std::span<const ComplexType>(downsampled_block.data(), downsampled_block.size())
			: std::span<const ComplexType>(block.data(), block.size());
		if (!output_samples.empty() && (!downsampled_block.empty() || dechirped || params::oversampleRatio() <= 1))
		{
			emitStreamingOutputBlock(receiver_index, output_start_time, output_sample_rate, output_samples,
									 output_sample_start);
		}
		block.clear();
		if (finish_downsampler && _streaming_downsamplers[receiver_index])
		{
			_streaming_downsamplers[receiver_index].reset();
		}
	}

	fers_signal::DownsamplingSink& SimulationEngine::streamingDownsampler(const std::size_t receiver_index,
																		  const std::uint64_t input_start_index,
																		  const RealType segment_start_time)
	{
		if (!_streaming_downsamplers[receiver_index])
		{
			_streaming_downsamplers[receiver_index] = std::make_unique<fers_signal::DownsamplingSink>();
			_streaming_downsample_base_indices[receiver_index] =
				input_start_index / std::max<unsigned>(1, _streaming_downsamplers[receiver_index]->ratio());
			_streaming_downsample_segment_start_times[receiver_index] = segment_start_time;
		}
		return *_streaming_downsamplers[receiver_index];
	}

	RealType SimulationEngine::streamingOutputSampleRate(const std::size_t receiver_index) const
	{
		if (receiver_index >= _world->getReceivers().size())
		{
			return 0.0;
		}

		const auto& receiver = _world->getReceivers()[receiver_index];
		if (receiver->hasFmcwIfResamplingSink())
		{
			const auto& if_plan = receiver->getFmcwIfResamplerPlan();
			return if_plan.has_value() ? if_plan->actual_output_sample_rate_hz : 0.0;
		}
		if (receiver->isDechirpEnabled())
		{
			return params::rate() * static_cast<RealType>(params::oversampleRatio());
		}
		return params::rate();
	}

	void SimulationEngine::ensureStreamingOutputStreamOpen(const std::size_t receiver_index,
														   const RealType first_sample_time, const RealType sample_rate)
	{
		if (_output_sink == nullptr || receiver_index >= _world->getReceivers().size() || sample_rate <= 0.0)
		{
			return;
		}
		if (_streaming_output_stream_ids[receiver_index] != 0 && _streaming_output_stream_open[receiver_index] &&
			_streaming_output_file_metadata[receiver_index])
		{
			return;
		}

		const auto& receiver = _world->getReceivers()[receiver_index];
		auto streaming_sources = collectStreamingSourcesForWindow(params::startTime(), params::endTime());
		if (_streaming_output_stream_ids[receiver_index] == 0)
		{
			_streaming_output_stream_ids[receiver_index] = _output_sink->registerStream(
				processing::buildReceiverStreamDescriptor(receiver.get(), sample_rate, streaming_sources));
		}
		if (!_streaming_output_file_metadata[receiver_index])
		{
			_streaming_output_file_metadata[receiver_index] =
				std::make_shared<OutputFileMetadata>(processing::buildStreamingOutputMetadata(
					receiver.get(), "", expectedStreamingOutputSamples(sample_rate), streaming_sources, sample_rate));
		}
		if (!_streaming_output_stream_open[receiver_index])
		{
			_output_sink->openStream(_streaming_output_stream_ids[receiver_index], first_sample_time);
			_streaming_output_stream_open[receiver_index] = true;
		}
	}

	void SimulationEngine::emitStreamingOutputBlock(const std::size_t receiver_index, const RealType first_sample_time,
													const RealType sample_rate,
													const std::span<const ComplexType> samples,
													const std::uint64_t sample_start)
	{
		if (_output_sink == nullptr || samples.empty() || receiver_index >= _world->getReceivers().size())
		{
			return;
		}

		const auto& receiver = _world->getReceivers()[receiver_index];
		auto& processed = _streaming_output_processed_buffers[receiver_index];
		processed.assign(samples.begin(), samples.end());
		processing::applyThermalNoiseAtSampleRate(processed, receiver->getNoiseTemperature(), receiver->getRngEngine(),
												  sample_rate);

		auto streaming_sources = collectStreamingSourcesForWindow(params::startTime(), params::endTime());
		ensureStreamingOutputStreamOpen(receiver_index, first_sample_time, sample_rate);

		const auto block = processing::buildReceiverSampleBlock(receiver.get(), first_sample_time, sample_rate,
																processed, sample_start, streaming_sources,
																_streaming_output_file_metadata[receiver_index]);
		_output_sink->submitBlock(block);
		_streaming_output_sample_cursors[receiver_index] = sample_start + static_cast<std::uint64_t>(processed.size());
	}

	void SimulationEngine::emitContextHeartbeatsThrough(const RealType simulation_time)
	{
		if (_output_sink == nullptr)
		{
			return;
		}
		if (_next_context_heartbeat_time > simulation_time)
		{
			return;
		}

		if (simulation_time - _next_context_heartbeat_time < 1.0)
		{
			_output_sink->emitContextHeartbeat(_next_context_heartbeat_time);
			_next_context_heartbeat_time += 1.0;
			return;
		}

		_output_sink->emitContextHeartbeat(simulation_time);
		_next_context_heartbeat_time = simulation_time + 1.0;
	}

	void SimulationEngine::flushFmcwIfBlocks()
	{
		for (std::size_t receiver_index = 0; receiver_index < _fmcw_if_block_buffers.size(); ++receiver_index)
		{
			flushFmcwIfBlock(receiver_index);
		}
	}

	void SimulationEngine::flushFmcwIfBlock(const std::size_t receiver_index)
	{
		if (receiver_index >= _world->getReceivers().size())
		{
			return;
		}
		auto& block = _fmcw_if_block_buffers[receiver_index];
		if (block.empty())
		{
			return;
		}
		const auto& receiver = _world->getReceivers()[receiver_index];
		if (!receiver->hasFmcwIfResamplingSink())
		{
			block.clear();
			return;
		}

		applyPulsedInterferenceToFmcwIfBlock(receiver_index, block, _fmcw_if_block_start_times[receiver_index]);
		receiver->consumeFmcwIfBlock(block, _fmcw_if_block_start_times[receiver_index]);
		block.clear();
	}

	void SimulationEngine::applyPulsedInterferenceToFmcwIfBlock(const std::size_t receiver_index,
																std::span<ComplexType> block,
																const RealType block_start_time)
	{
		applyPulsedInterferenceToStreamingBlock(receiver_index, block, block_start_time,
												params::rate() * static_cast<RealType>(params::oversampleRatio()),
												true);
	}

	void SimulationEngine::addPulsedInterferenceSamples(std::span<ComplexType> block,
														std::span<const ComplexType> rendered_pulse,
														const long long dest_begin, const long long dest_end,
														const std::size_t crop_offset, const RealType block_start_time,
														const RealType sample_rate, const bool dechirp_mix,
														Receiver* receiver, ReceiverTrackerCache& tracker_cache) const
	{
		for (long long dest = dest_begin; dest < dest_end; ++dest)
		{
			const RealType t_sample = block_start_time + static_cast<RealType>(dest) / sample_rate;
			const auto source_index = crop_offset + static_cast<std::size_t>(dest - dest_begin);
			if (source_index >= rendered_pulse.size())
			{
				continue;
			}
			if (dechirp_mix)
			{
				const auto mixer = calculateDechirpMixer(receiver, t_sample, tracker_cache);
				if (!mixer.has_value())
				{
					continue;
				}
				block[static_cast<std::size_t>(dest)] += *mixer * std::conj(rendered_pulse[source_index]);
			}
			else
			{
				block[static_cast<std::size_t>(dest)] += rendered_pulse[source_index];
			}
		}
	}

	void SimulationEngine::applyPulsedInterferenceToStreamingBlock(const std::size_t receiver_index,
																   std::span<ComplexType> block,
																   const RealType block_start_time,
																   const RealType sample_rate, const bool dechirp_mix)
	{
		if (block.empty() || receiver_index >= _world->getReceivers().size())
		{
			return;
		}

		const auto& receiver = _world->getReceivers()[receiver_index];
		if (!std::isfinite(sample_rate) || sample_rate <= 0.0)
		{
			return;
		}
		const RealType block_end_time = block_start_time + static_cast<RealType>(block.size()) / sample_rate;
		auto& tracker_cache = _if_pulse_tracker_caches[receiver_index];

		receiver->prunePulsedInterferenceEndingBefore(block_start_time);
		for (const auto& response : receiver->getPulsedInterferenceLog())
		{
			const RealType pulse_rate = response->sampleRate();
			const unsigned pulse_size = response->sampleCount();
			if (pulse_rate <= 0.0 || pulse_size == 0)
			{
				continue;
			}

			const RealType pulse_start_time = response->startTime();
			const RealType pulse_end_time = pulse_start_time + static_cast<RealType>(pulse_size) / pulse_rate;
			if (pulse_end_time <= block_start_time || pulse_start_time >= block_end_time)
			{
				continue;
			}

			const RealType overlap_start = std::max(pulse_start_time, block_start_time);
			const RealType overlap_end = std::min(pulse_end_time, block_end_time);
			const auto dest_begin = static_cast<long long>(
				std::max<RealType>(0.0, std::ceil((overlap_start - block_start_time) * sample_rate)));
			const auto dest_end = static_cast<long long>(std::min<RealType>(
				static_cast<RealType>(block.size()), std::ceil((overlap_end - block_start_time) * sample_rate)));
			if (dest_begin >= dest_end)
			{
				continue;
			}

			const auto interp_padding = static_cast<long long>(params::renderFilterLength()) / 2 + 1;
			const long long padded_begin = dest_begin - interp_padding;
			const long long padded_end = dest_end + interp_padding;
			const RealType render_start = block_start_time + static_cast<RealType>(padded_begin) / sample_rate;
			const auto render_count = static_cast<std::size_t>(padded_end - padded_begin);
			const auto rendered_pulse = response->renderSlice(sample_rate, render_start, render_count, 0.0);
			const auto crop_offset = static_cast<std::size_t>(dest_begin - padded_begin);
			addPulsedInterferenceSamples(block, rendered_pulse, dest_begin, dest_end, crop_offset, block_start_time,
										 sample_rate, dechirp_mix, receiver.get(), tracker_cache);
		}
	}

	std::optional<ComplexType> SimulationEngine::calculateDechirpMixer(Receiver* rx, const RealType t_step,
																	   ReceiverTrackerCache& tracker_cache) const
	{
		RealType reference_phase = 0.0;
		const auto& dechirp_sources = rx->getDechirpSources();
		if (tracker_cache.dechirp_reference.size() < dechirp_sources.size())
		{
			tracker_cache.dechirp_reference.resize(dechirp_sources.size());
		}

		if (!tracker_cache.last_dechirp_time.has_value() || t_step < *tracker_cache.last_dechirp_time)
		{
			tracker_cache.active_dechirp_source_index = 0;
			std::ranges::fill(tracker_cache.dechirp_reference, FmcwChirpBoundaryTracker{});
		}
		tracker_cache.last_dechirp_time = t_step;

		bool reference_active = false;
		auto& source_index = tracker_cache.active_dechirp_source_index;
		while (source_index < dechirp_sources.size() && t_step >= dechirp_sources[source_index].segment_end)
		{
			++source_index;
		}
		if (source_index < dechirp_sources.size())
		{
			const auto& reference_source = dechirp_sources[source_index];
			if (t_step >= reference_source.segment_start && t_step < reference_source.segment_end &&
				simulation::calculateStreamingReferencePhase(
					reference_source, t_step, &tracker_cache.dechirp_reference[source_index], reference_phase))
			{
				reference_active = true;
			}
		}

		if (!reference_active)
		{
			return std::nullopt;
		}

		RealType receiver_phase = 0.0;
		if (rx->getDechirpMode() == Receiver::DechirpMode::Physical && _cw_phase_noise_lookup)
		{
			receiver_phase = _cw_phase_noise_lookup->sample(rx->getTiming().get(), t_step);
		}
		return std::polar(1.0, reference_phase + receiver_phase);
	}

	ComplexType SimulationEngine::calculateStreamingSample(Receiver* rx, const RealType t_step,
														   const std::vector<ActiveStreamingSource>& streaming_sources,
														   ReceiverTrackerCache& tracker_cache) const
	{
		const bool dechirping = rx->isDechirpEnabled();
		std::optional<ComplexType> dechirp_mixer;
		if (dechirping)
		{
			dechirp_mixer = calculateDechirpMixer(rx, t_step, tracker_cache);
			if (!dechirp_mixer.has_value())
			{
				return {0.0, 0.0};
			}
		}

		const auto timing_phase_mode = !dechirping ? simulation::StreamingTimingPhaseMode::ReceiverRelative
												   : (rx->getDechirpMode() == Receiver::DechirpMode::Ideal
														  ? simulation::StreamingTimingPhaseMode::None
														  : simulation::StreamingTimingPhaseMode::TransmitterOnly);

		ComplexType total_sample{0.0, 0.0};
		for (std::size_t source_index = 0; source_index < streaming_sources.size(); ++source_index)
		{
			const auto& streaming_source = streaming_sources[source_index];
			if (!rx->checkFlag(Receiver::RecvFlag::FLAG_NODIRECT))
			{
				total_sample += simulation::calculateStreamingDirectPathContribution(
					streaming_source, rx, t_step, _cw_phase_noise_lookup.get(), &tracker_cache.direct[source_index],
					timing_phase_mode);
			}
			for (std::size_t target_index = 0; target_index < _world->getTargets().size(); ++target_index)
			{
				const auto& target_ptr = _world->getTargets()[target_index];
				total_sample += simulation::calculateStreamingReflectedPathContribution(
					streaming_source, rx, target_ptr.get(), t_step, _cw_phase_noise_lookup.get(),
					&tracker_cache.reflected[source_index][target_index], timing_phase_mode);
			}
		}

		if (!dechirping)
		{
			return total_sample;
		}

		// Mixing Convention: s_IF = s_ref * conj(s_rx)
		// This convention is chosen to ensure that:
		// 1. Stationary targets (positive delay tau) result in a POSITIVE beat frequency (f_b = alpha * tau).
		// 2. In physical dechirp mode, phase noise from the same LO source partially cancels
		//    at short ranges (Range Correlation Effect).
		// 3. For an up-chirp, a receding target (negative RF Doppler) results in a
		//    higher IF frequency (f_IF = f_b + |f_d|).
		return *dechirp_mixer * std::conj(total_sample);
	}

	void SimulationEngine::appendStreamingTrackerSource()
	{
		const std::size_t target_count = _world->getTargets().size();

		for (auto& cache : _streaming_tracker_caches)
		{
			cache.direct.emplace_back();
			cache.reflected.emplace_back(target_count);
		}
	}

	void SimulationEngine::eraseStreamingTrackerSource(const std::size_t source_index)
	{
		for (auto& cache : _streaming_tracker_caches)
		{
			if (source_index < cache.direct.size())
			{
				cache.direct.erase(cache.direct.begin() + static_cast<std::ptrdiff_t>(source_index));
			}
			if (source_index < cache.reflected.size())
			{
				cache.reflected.erase(cache.reflected.begin() + static_cast<std::ptrdiff_t>(source_index));
			}
		}
	}

	void SimulationEngine::cleanupInactiveStreamingSources(const RealType from_time)
	{
		auto& sources = _world->getSimulationState().active_streaming_transmitters;
		for (std::size_t source_index = sources.size(); source_index > 0; --source_index)
		{
			const std::size_t index = source_index - 1;
			if (sources[index].segment_end > from_time)
			{
				continue;
			}
			const auto cleanup_deadline = streamingSourceCleanupDeadline(sources[index], from_time);
			if (cleanup_deadline.has_value() && from_time < *cleanup_deadline)
			{
				continue;
			}

			sources.erase(sources.begin() + static_cast<std::ptrdiff_t>(index));
			eraseStreamingTrackerSource(index);
		}
	}

	std::optional<RealType> SimulationEngine::streamingSourceCleanupDeadline(const ActiveStreamingSource& source,
																			 const RealType from_time) const
	{
		if (source.transmitter == nullptr || source.carrier_freq <= 0.0)
		{
			return std::nullopt;
		}

		std::optional<RealType> latest_deadline;
		for (const auto& receiver_ptr : _world->getReceivers())
		{
			const auto receiver_deadline = receiverCleanupDeadline(source, receiver_ptr.get(), from_time);
			if (receiver_deadline.has_value() &&
				(!latest_deadline.has_value() || *receiver_deadline > *latest_deadline))
			{
				latest_deadline = receiver_deadline;
			}
		}
		return latest_deadline;
	}

	std::optional<RealType> SimulationEngine::receiverCleanupDeadline(const ActiveStreamingSource& source,
																	  const Receiver* const rx,
																	  const RealType from_time) const
	{
		if (!isStreamingReceiver(rx))
		{
			return std::nullopt;
		}

		const auto update_latest = [](std::optional<RealType>& latest, const std::optional<RealType> candidate)
		{
			if (candidate.has_value() && (!latest.has_value() || *candidate > *latest))
			{
				latest = candidate;
			}
		};

		const auto interval_deadline = [&](const RealType interval_start,
										   const RealType interval_end) -> std::optional<RealType>
		{
			const RealType start = std::max({params::startTime(), from_time, interval_start});
			const RealType end = std::min(params::endTime(), interval_end);
			if (start >= end)
			{
				return std::nullopt;
			}

			std::optional<RealType> latest;
			if (!rx->checkFlag(Receiver::RecvFlag::FLAG_NODIRECT))
			{
				update_latest(latest, directPathCleanupDeadline(source, rx, start, end));
			}
			for (const auto& target_ptr : _world->getTargets())
			{
				update_latest(latest, reflectedPathCleanupDeadline(source, rx, target_ptr.get(), start, end));
			}
			return latest;
		};

		std::optional<RealType> latest_deadline;
		const auto& schedule = rx->getSchedule();
		if (schedule.empty())
		{
			update_latest(latest_deadline, interval_deadline(params::startTime(), params::endTime()));
			return latest_deadline;
		}

		for (const auto& period : schedule)
		{
			update_latest(latest_deadline, interval_deadline(period.start, period.end));
		}
		return latest_deadline;
	}

	void SimulationEngine::processEvent(const Event& event)
	{
		// NOLINTBEGIN(cppcoreguidelines-pro-type-static-cast-downcast)
		switch (event.type)
		{
		case EventType::TX_PULSED_START:
			handleTxPulsedStart(static_cast<Transmitter*>(event.source_object), event.timestamp);
			break;
		case EventType::RX_PULSED_WINDOW_START:
			handleRxPulsedWindowStart(static_cast<Receiver*>(event.source_object), event.timestamp);
			break;
		case EventType::RX_PULSED_WINDOW_END:
			handleRxPulsedWindowEnd(static_cast<Receiver*>(event.source_object), event.timestamp);
			break;
		case EventType::TX_STREAMING_START:
			if (const auto source = streamingSourceAtEvent(static_cast<Transmitter*>(event.source_object),
														   event.timestamp, _internal_stop_time);
				source.has_value())
			{
				handleTxStreamingStart(*source);
			}
			break;
		case EventType::TX_STREAMING_END:
			handleTxStreamingEnd(static_cast<Transmitter*>(event.source_object));
			break;
		case EventType::RX_STREAMING_START:
			handleRxStreamingStart(static_cast<Receiver*>(event.source_object));
			break;
		case EventType::RX_STREAMING_END:
			handleRxStreamingEnd(static_cast<Receiver*>(event.source_object));
			break;
		}
		// NOLINTEND(cppcoreguidelines-pro-type-static-cast-downcast)
	}

	void SimulationEngine::routeResponse(Receiver* rx, std::unique_ptr<serial::Response> response) const
	{
		if (!response)
		{
			return;
		}
		if (rx->getMode() == OperationMode::PULSED_MODE)
		{
			rx->addResponseToInbox(std::move(response));
		}
		else
		{
			rx->addInterferenceToLog(std::move(response));
		}
	}

	void SimulationEngine::handleTxPulsedStart(Transmitter* tx, const RealType t_event)
	{
		for (const auto& rx_ptr : _world->getReceivers())
		{
			if (!rx_ptr->checkFlag(Receiver::RecvFlag::FLAG_NODIRECT))
			{
				routeResponse(rx_ptr.get(), simulation::calculateResponse(tx, rx_ptr.get(), tx->getSignal(), t_event));
			}
			for (const auto& target_ptr : _world->getTargets())
			{
				routeResponse(
					rx_ptr.get(),
					simulation::calculateResponse(tx, rx_ptr.get(), tx->getSignal(), t_event, target_ptr.get()));
			}
		}

		const RealType next_theoretical_time = t_event + 1.0 / tx->getPrf();
		if (const auto next_pulse_opt = tx->getNextPulseTime(next_theoretical_time);
			next_pulse_opt && *next_pulse_opt <= params::endTime())
		{
			_world->getEventQueue().push({*next_pulse_opt, EventType::TX_PULSED_START, tx});
		}
	}

	void SimulationEngine::handleRxPulsedWindowStart(Receiver* rx, const RealType t_event)
	{
		rx->setActive(true);
		_world->getEventQueue().push({t_event + rx->getWindowLength(), EventType::RX_PULSED_WINDOW_END, rx});
	}

	void SimulationEngine::handleRxPulsedWindowEnd(Receiver* rx, const RealType t_event)
	{
		rx->setActive(false);
		const auto active_streaming_sources =
			collectStreamingSourcesForWindow(t_event - rx->getWindowLength(), t_event);

		RenderingJob job{.ideal_start_time = t_event - rx->getWindowLength(),
						 .duration = rx->getWindowLength(),
						 .responses = rx->drainInbox(),
						 .active_streaming_sources = active_streaming_sources};

		rx->enqueueFinalizerJob(std::move(job));

		const RealType next_theoretical = t_event - rx->getWindowLength() + 1.0 / rx->getWindowPrf();
		if (const auto next_start = rx->getNextWindowTime(next_theoretical);
			next_start && *next_start <= params::endTime())
		{
			_world->getEventQueue().push({*next_start, EventType::RX_PULSED_WINDOW_START, rx});
		}
	}

	void SimulationEngine::handleTxStreamingStart(const ActiveStreamingSource& source)
	{
		_world->getSimulationState().active_streaming_transmitters.push_back(source);
		appendStreamingTrackerSource();
	}

	void SimulationEngine::handleTxStreamingEnd(Transmitter* tx)
	{
		(void)tx;
		// A transmitter stop is a transmit-time boundary, not an instantaneous receive-time cutoff.
		// Ended sources are removed only after all future receive-time samples fail the retarded-time gate.
		cleanupInactiveStreamingSources(_world->getSimulationState().t_current);
	}

	void SimulationEngine::handleRxStreamingStart(Receiver* rx)
	{
		rx->setActive(true);
		const auto receiver_it = std::ranges::find_if(_world->getReceivers(), [rx](const auto& receiver_ptr)
													  { return receiver_ptr.get() == rx; });
		if (receiver_it != _world->getReceivers().end())
		{
			const auto receiver_index = static_cast<std::size_t>(receiver_it - _world->getReceivers().begin());
			_streaming_downsamplers[receiver_index].reset();
			if (_eager_context_stream_open)
			{
				ensureStreamingOutputStreamOpen(receiver_index, _world->getSimulationState().t_current,
												streamingOutputSampleRate(receiver_index));
			}
		}
		if (rx->hasFmcwIfResamplingSink())
		{
			rx->beginFmcwIfResamplingSegment(_world->getSimulationState().t_current);
		}
	}

	void SimulationEngine::handleRxStreamingEnd(Receiver* rx)
	{
		const auto receiver_it = std::ranges::find_if(_world->getReceivers(), [rx](const auto& receiver_ptr)
													  { return receiver_ptr.get() == rx; });
		if (receiver_it != _world->getReceivers().end())
		{
			const auto receiver_index = static_cast<std::size_t>(receiver_it - _world->getReceivers().begin());
			flushFmcwIfBlock(receiver_index);
			flushStreamingOutputBlock(receiver_index, true);
		}
		if (rx->hasFmcwIfResamplingSink() && _world->getSimulationState().t_current >= params::endTime() &&
			_world->getSimulationState().t_current < _internal_stop_time && activePastUserEnd(rx))
		{
			return;
		}
		if (rx->hasFmcwIfResamplingSink())
		{
			rx->endFmcwIfResamplingSegment();
		}
		if (_output_sink != nullptr && receiver_it != _world->getReceivers().end())
		{
			const auto receiver_index = static_cast<std::size_t>(receiver_it - _world->getReceivers().begin());
			if (_streaming_output_stream_open[receiver_index])
			{
				_output_sink->closeStream(_streaming_output_stream_ids[receiver_index]);
				_streaming_output_stream_open[receiver_index] = false;
			}
		}
		rx->setActive(false);
	}

	void SimulationEngine::updateProgress() { reportSimulationProgress(_world->getSimulationState().t_current); }

	bool SimulationEngine::isCancellationRequested()
	{
		if (_cancelled)
		{
			return true;
		}
		if (_cancel_callback && _cancel_callback())
		{
			_cancelled = true;
			LOG(Level::INFO, "Simulation cancellation requested.");
			if (_reporter)
			{
				_reporter->report("Simulation cancelled", 100, 100);
			}
			return true;
		}
		return false;
	}

	void SimulationEngine::reportSimulationProgress(const RealType t_current)
	{
		if (!_reporter)
		{
			return;
		}

		const RealType start_time = params::startTime();
		const RealType end_time = params::endTime();
		const RealType duration = end_time - start_time;
		const RealType progress_fraction = duration > 0.0 ? (t_current - start_time) / duration : 1.0;
		const int progress = static_cast<int>(
			std::clamp(progress_fraction * 100.0, static_cast<RealType>(0.0), static_cast<RealType>(100.0)));

		if (const auto now = std::chrono::steady_clock::now();
			progress != _last_reported_percent || now - _last_report_time >= std::chrono::milliseconds(100))
		{
			_reporter->report(std::format("Simulating... {:.2f}s / {:.2f}s", t_current, end_time), progress, 100);
			_last_reported_percent = progress;
			_last_report_time = now;
		}
	}

	std::vector<ActiveStreamingSource> SimulationEngine::collectStreamingSourcesForWindow(const RealType start_time,
																						  const RealType end_time) const
	{
		// A segment that ended before this window can still be in flight at the receiver.
		(void)start_time;
		std::vector<ActiveStreamingSource> sources;
		for (const auto& transmitter_ptr : _world->getTransmitters())
		{
			if (!transmitter_ptr->isStreamingMode())
			{
				continue;
			}

			const auto append_candidate = [&](const RealType segment_start, const RealType segment_end)
			{
				auto source = makeActiveSource(transmitter_ptr.get(), segment_start, segment_end);
				if (source.segment_start < source.segment_end && source.segment_start < end_time)
				{
					sources.push_back(source);
				}
			};

			if (transmitter_ptr->getSchedule().empty())
			{
				append_candidate(params::startTime(), params::endTime());
				continue;
			}

			for (const auto& period : transmitter_ptr->getSchedule())
			{
				append_candidate(period.start, std::min(params::endTime(), period.end));
			}
		}
		return sources;
	}

	void SimulationEngine::shutdown()
	{
		LOG(Level::INFO, "Simulation compute loop finished. Waiting for receiver finalization tasks...");
		if (_reporter)
		{
			_reporter->report("Simulation compute finished. Waiting for receiver finalization...", 100, 100);
		}

		for (std::size_t receiver_index = 0; receiver_index < _world->getReceivers().size(); ++receiver_index)
		{
			const auto& receiver_ptr = _world->getReceivers()[receiver_index];
			if (isStreamingReceiver(receiver_ptr.get()))
			{
				if (_output_sink != nullptr)
				{
					flushFmcwIfBlock(receiver_index);
					receiver_ptr->flushFmcwIfResampling();
					flushStreamingOutputBlock(receiver_index, true);
					if (_streaming_output_stream_open[receiver_index])
					{
						_output_sink->closeStream(_streaming_output_stream_ids[receiver_index]);
						_streaming_output_stream_open[receiver_index] = false;
					}
				}
			}
			else if (receiver_ptr->getMode() == OperationMode::PULSED_MODE)
			{
				RenderingJob shutdown_job{};
				shutdown_job.duration = -1.0;
				receiver_ptr->enqueueFinalizerJob(std::move(shutdown_job));
			}
		}

		_pool.wait();
		for (auto& finalizer_thread : _finalizer_threads)
		{
			if (finalizer_thread.joinable())
			{
				finalizer_thread.join();
			}
		}

		LOG(Level::INFO, "All finalization tasks complete.");
	}

	OutputMetadata runEventDrivenSim(World* world, pool::ThreadPool& pool,
									 const std::function<void(const std::string&, int, int)>& progress_callback,
									 const std::string& output_dir, const OutputConfig& output_config,
									 std::function<bool()> cancel_callback, bool* cancelled,
									 ReceiverOutputTelemetryCallback telemetry_callback)
	{
		if (cancelled != nullptr)
		{
			*cancelled = false;
		}
		auto reporter = std::make_shared<ProgressReporter>(progress_callback);
		auto metadata_collector = std::make_shared<OutputMetadataCollector>(output_dir);
		std::unique_ptr<ReceiverOutputSink> output_sink;
		if (isVita49Enabled(output_config))
		{
			output_sink = serial::vita49::makeVita49OutputSink(std::move(telemetry_callback));
			output_sink->initializeRun(output_config, params::params.simulation_name);
		}
		else
		{
			output_sink = serial::makeHdf5OutputSink(output_dir, metadata_collector);
			output_sink->initializeRun(output_config, params::params.simulation_name);
		}

		SimulationEngine engine(world, pool, reporter, output_dir, metadata_collector, output_sink.get(),
								std::move(cancel_callback), isVita49Enabled(output_config));
		engine.run();
		if (cancelled != nullptr)
		{
			*cancelled = engine.cancelled();
		}
		if (isVita49Enabled(output_config))
		{
			LOG(Level::INFO, "Waiting for VITA output stream drain...");
			reporter->report("Waiting for VITA output stream drain...", 100, 100);
		}
		const auto stats = output_sink->finalize();
		reporter->report(engine.cancelled() ? "Simulation cancelled" : "Simulation complete", 100, 100);
		LOG(Level::INFO, "Event-driven simulation loop finished.");
		auto metadata = metadata_collector->snapshot();
		if (output_sink)
		{
			if (isVita49Enabled(output_config))
			{
				auto vita49_metadata = vita49MetadataFromConfig(output_config.vita49);
				if (stats.epoch_unix_nanoseconds.has_value())
				{
					vita49_metadata.epoch_unix_nanoseconds = stats.epoch_unix_nanoseconds;
				}
				for (const auto& stream : stats.streams)
				{
					vita49_metadata.streams.push_back(streamStatsToMetadata(stream));
				}
				metadata.vita49 = std::move(vita49_metadata);
			}
		}
		return metadata;
	}
}
