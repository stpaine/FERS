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
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
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
		[[nodiscard]] bool isStreamingReceiver(const Receiver* const receiver) noexcept
		{
			return receiver != nullptr &&
				(receiver->getMode() == OperationMode::CW_MODE || receiver->getMode() == OperationMode::FMCW_MODE);
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
				const RealType a = 0.5 * h2 * (dd_right[axis] - dd_left[axis]);
				const RealType b = h2 * dd_left[axis];
				const RealType c = (right[axis] - left[axis]) + (h2 / 6.0) * (-2.0 * dd_left[axis] - dd_right[axis]);

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

		void includeQuadraticVelocityExtremum(std::array<RealType, 3>& max_abs_velocity, const std::size_t axis,
											  const RealType a, const RealType b, const RealType c,
											  const RealType segment_length, const RealType root_u,
											  const RealType lower_u, const RealType upper_u) noexcept
		{
			if (root_u < lower_u || root_u > upper_u || segment_length <= EPSILON)
			{
				return;
			}
			const RealType velocity = (a * root_u * root_u + b * root_u + c) / segment_length;
			if (std::isfinite(velocity))
			{
				max_abs_velocity[axis] = std::max(max_abs_velocity[axis], std::abs(velocity));
			}
			else
			{
				max_abs_velocity[axis] = std::numeric_limits<RealType>::infinity();
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
				const RealType a = 0.5 * h2 * (dd_right[axis] - dd_left[axis]);
				const RealType b = h2 * dd_left[axis];
				const RealType c = (right[axis] - left[axis]) + (h2 / 6.0) * (-2.0 * dd_left[axis] - dd_right[axis]);
				includeQuadraticVelocityExtremum(max_abs_velocity, axis, a, b, c, segment_length, lower_u, lower_u,
												 upper_u);
				includeQuadraticVelocityExtremum(max_abs_velocity, axis, a, b, c, segment_length, upper_u, lower_u,
												 upper_u);

				if (std::abs(a) > EPSILON)
				{
					includeQuadraticVelocityExtremum(max_abs_velocity, axis, a, b, c, segment_length, -b / (2.0 * a),
													 lower_u, upper_u);
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
																	const RealType timestamp)
		{
			if (transmitter == nullptr || !transmitter->isStreamingMode())
			{
				return std::nullopt;
			}

			const auto& schedule = transmitter->getSchedule();
			if (schedule.empty())
			{
				const RealType segment_start = params::startTime();
				auto source = makeActiveSource(transmitter, segment_start, params::endTime());
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
				auto source = makeActiveSource(transmitter, period.start, std::min(params::endTime(), period.end));
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
									   std::shared_ptr<OutputMetadataCollector> metadata_collector) :
		_world(world), _pool(pool), _reporter(std::move(reporter)), _metadata_collector(std::move(metadata_collector)),
		_last_report_time(std::chrono::steady_clock::now()), _output_dir(std::move(output_dir))
	{
		_streaming_tracker_caches.resize(_world->getReceivers().size());
	}

	void SimulationEngine::run()
	{
		if (_reporter)
		{
			_reporter->report("Initializing event-driven simulation...", 0, 100);
		}

		logSimulationMemoryProjection(*_world);

		initializeFinalizers();

		LOG(Level::INFO, "Starting unified event-driven simulation loop.");
		logStreamingSummaries();

		auto& event_queue = _world->getEventQueue();
		auto& state = _world->getSimulationState();
		const RealType end_time = params::endTime();

		while (!event_queue.empty() && state.t_current <= end_time)
		{
			const Event event = event_queue.top();
			event_queue.pop();

			processStreamingPhysics(event.timestamp);

			state.t_current = event.timestamp;

			processEvent(event);
			updateProgress();
		}

		shutdown();
	}

	void SimulationEngine::logStreamingSummaries() const
	{
		for (const auto& transmitter_ptr : _world->getTransmitters())
		{
			const auto* waveform = transmitter_ptr->getSignal();
			if (waveform == nullptr || !waveform->isFmcwFamily())
			{
				continue;
			}

			if (const auto* fmcw = waveform->getFmcwChirpSignal(); fmcw != nullptr)
			{
				const RealType duty_cycle = fmcw->getChirpDuration() / fmcw->getChirpPeriod();
				const RealType average_power = waveform->getPower() * duty_cycle;
				const auto direction = fers_signal::fmcwChirpDirectionToken(fmcw->getDirection());
				const auto configured_count = fmcw->getChirpCount().has_value()
					? std::format("{}", *fmcw->getChirpCount())
					: std::string("unbounded");
				if (transmitter_ptr->getSchedule().empty())
				{
					const RealType active_start = params::startTime();
					const auto source = makeActiveSource(transmitter_ptr.get(), active_start, params::endTime());
					const auto total_chirp_count = countFmcwChirpStarts(source, active_start, source.segment_end);
					LOG(Level::INFO,
						"FMCW transmitter '{}' shape=linear {} B={} Hz T_c={} s T_rep={} s f_0={} Hz alpha={} Hz/s "
						"duty_cycle={} chirp_count={} total_chirp_count={} average_power={} W",
						transmitter_ptr->getName(), direction, fmcw->getChirpBandwidth(), fmcw->getChirpDuration(),
						fmcw->getChirpPeriod(), fmcw->getStartFrequencyOffset(), fmcw->getChirpRate(), duty_cycle,
						configured_count, total_chirp_count, average_power);
				}
				else
				{
					std::uint64_t total_chirp_count = 0;
					for (const auto& period : transmitter_ptr->getSchedule())
					{
						const RealType active_start = std::max(params::startTime(), period.start);
						const auto source = makeActiveSource(transmitter_ptr.get(), period.start,
															 std::min(params::endTime(), period.end));
						const auto segment_chirp_count = countFmcwChirpStarts(source, active_start, source.segment_end);
						total_chirp_count += segment_chirp_count;
						LOG(Level::INFO,
							"FMCW transmitter '{}' segment [{}, {}] shape=linear {} B={} Hz T_c={} s T_rep={} s f_0={} "
							"Hz alpha={} Hz/s duty_cycle={} chirp_count={} segment_chirp_count={} total_chirp_count={} "
							"average_power={} W",
							transmitter_ptr->getName(), period.start, source.segment_end, direction,
							fmcw->getChirpBandwidth(), fmcw->getChirpDuration(), fmcw->getChirpPeriod(),
							fmcw->getStartFrequencyOffset(), fmcw->getChirpRate(), duty_cycle, configured_count,
							segment_chirp_count, total_chirp_count, average_power);
					}
				}
			}
			else if (const auto* triangle = waveform->getFmcwTriangleSignal(); triangle != nullptr)
			{
				const RealType average_power = waveform->getPower();
				const auto configured_count = triangle->getTriangleCount().has_value()
					? std::format("{}", *triangle->getTriangleCount())
					: std::string("unbounded");
				if (transmitter_ptr->getSchedule().empty())
				{
					const RealType active_start = params::startTime();
					const auto source = makeActiveSource(transmitter_ptr.get(), active_start, params::endTime());
					const auto total_triangle_count = countFmcwTriangleStarts(source, active_start, source.segment_end);
					LOG(Level::INFO,
						"FMCW transmitter '{}' shape=triangle B={} Hz T_c={} s T_tri={} s f_0={} Hz alpha={} Hz/s "
						"duty_cycle=1 triangle_count={} total_triangle_count={} average_power={} W",
						transmitter_ptr->getName(), triangle->getChirpBandwidth(), triangle->getChirpDuration(),
						triangle->getTrianglePeriod(), triangle->getStartFrequencyOffset(), triangle->getChirpRate(),
						configured_count, total_triangle_count, average_power);
				}
				else
				{
					std::uint64_t total_triangle_count = 0;
					for (const auto& period : transmitter_ptr->getSchedule())
					{
						const RealType active_start = std::max(params::startTime(), period.start);
						const auto source = makeActiveSource(transmitter_ptr.get(), period.start,
															 std::min(params::endTime(), period.end));
						const auto segment_triangle_count =
							countFmcwTriangleStarts(source, active_start, source.segment_end);
						total_triangle_count += segment_triangle_count;
						LOG(Level::INFO,
							"FMCW transmitter '{}' segment [{}, {}] shape=triangle B={} Hz T_c={} s T_tri={} s f_0={} "
							"Hz alpha={} Hz/s duty_cycle=1 triangle_count={} segment_triangle_count={} "
							"total_triangle_count={} average_power={} W",
							transmitter_ptr->getName(), period.start, source.segment_end, triangle->getChirpBandwidth(),
							triangle->getChirpDuration(), triangle->getTrianglePeriod(),
							triangle->getStartFrequencyOffset(), triangle->getChirpRate(), configured_count,
							segment_triangle_count, total_triangle_count, average_power);
					}
				}
			}
		}
	}

	void SimulationEngine::initializeFinalizers()
	{
		for (const auto& receiver_ptr : _world->getReceivers())
		{
			if (receiver_ptr->getMode() == OperationMode::PULSED_MODE)
			{
				_finalizer_threads.emplace_back(processing::runPulsedFinalizer, receiver_ptr.get(),
												&_world->getTargets(), _reporter, _output_dir, _metadata_collector);
			}
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
			simulation::CwPhaseNoiseLookup::build(timings, lookup_start, params::endTime()));
	}

	void SimulationEngine::processStreamingPhysics(const RealType t_event)
	{
		auto& state = _world->getSimulationState();
		auto& t_current = state.t_current;
		auto& active_streaming_transmitters = state.active_streaming_transmitters;

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

		const auto next_cleanup_deadline = [&](const RealType from_time) -> std::optional<RealType>
		{
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
		};

		while (t_current < t_event)
		{
			cleanupInactiveStreamingSources(t_current);

			RealType chunk_end = t_event;
			if (const auto cleanup_deadline = next_cleanup_deadline(t_current);
				cleanup_deadline.has_value() && *cleanup_deadline < chunk_end)
			{
				chunk_end = *cleanup_deadline;
			}
			if (chunk_end <= t_current)
			{
				break;
			}

			const auto start_index = streamingSampleIndexAtOrAfter(t_current, dt_sim);
			const auto end_index = streamingSampleIndexAtOrAfter(chunk_end, dt_sim);
			for (size_t sample_index = start_index; sample_index < end_index; ++sample_index)
			{
				const RealType t_step = params::startTime() + static_cast<RealType>(sample_index) * dt_sim;

				for (std::size_t receiver_index = 0; receiver_index < _world->getReceivers().size(); ++receiver_index)
				{
					const auto& receiver_ptr = _world->getReceivers()[receiver_index];
					if ((receiver_ptr->getMode() == OperationMode::CW_MODE ||
						 receiver_ptr->getMode() == OperationMode::FMCW_MODE) &&
						receiver_ptr->isActive())
					{
						ComplexType sample =
							calculateStreamingSample(receiver_ptr.get(), t_step, active_streaming_transmitters,
													 _streaming_tracker_caches[receiver_index]);
						receiver_ptr->setStreamingSample(sample_index, sample);
					}
				}
				if (((sample_index - first_index) % progress_report_stride) == 0 || sample_index + 1 == final_index)
				{
					reportSimulationProgress(t_step);
				}
			}

			t_current = chunk_end;
		}
		cleanupInactiveStreamingSources(t_current);
	}

	ComplexType SimulationEngine::calculateStreamingSample(Receiver* rx, const RealType t_step,
														   const std::vector<ActiveStreamingSource>& streaming_sources,
														   ReceiverTrackerCache& tracker_cache) const
	{
		const bool dechirping = rx->isDechirpEnabled();
		RealType reference_phase = 0.0;
		if (dechirping)
		{
			const auto& dechirp_sources = rx->getDechirpSources();
			if (tracker_cache.dechirp_reference.size() < dechirp_sources.size())
			{
				tracker_cache.dechirp_reference.resize(dechirp_sources.size());
			}

			if (!tracker_cache.last_dechirp_time.has_value() || t_step < *tracker_cache.last_dechirp_time)
			{
				tracker_cache.active_dechirp_source_index = 0;
				std::fill(tracker_cache.dechirp_reference.begin(), tracker_cache.dechirp_reference.end(),
						  FmcwChirpBoundaryTracker{});
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

		RealType receiver_phase = 0.0;
		if (rx->getDechirpMode() == Receiver::DechirpMode::Physical && _cw_phase_noise_lookup)
		{
			receiver_phase = _cw_phase_noise_lookup->sample(rx->getTiming().get(), t_step);
		}

		// Mixing Convention: s_IF = s_ref * conj(s_rx)
		// This convention is chosen to ensure that:
		// 1. Stationary targets (positive delay tau) result in a POSITIVE beat frequency (f_b = alpha * tau).
		// 2. In physical dechirp mode, phase noise from the same LO source partially cancels
		//    at short ranges (Range Correlation Effect).
		// 3. For an up-chirp, a receding target (negative RF Doppler) results in a
		//    higher IF frequency (f_IF = f_b + |f_d|).
		return std::polar(1.0, reference_phase + receiver_phase) * std::conj(total_sample);
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
			if (const auto source =
					streamingSourceAtEvent(static_cast<Transmitter*>(event.source_object), event.timestamp);
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

	void SimulationEngine::handleRxStreamingStart(Receiver* rx) { rx->setActive(true); }

	void SimulationEngine::handleRxStreamingEnd(Receiver* rx) { rx->setActive(false); }

	void SimulationEngine::updateProgress() { reportSimulationProgress(_world->getSimulationState().t_current); }

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
		LOG(Level::INFO, "Main simulation loop finished. Waiting for finalization tasks...");
		if (_reporter)
		{
			_reporter->report("Main simulation finished. Waiting for data export...", 100, 100);
		}

		const auto streaming_sources = collectStreamingSourcesForWindow(params::startTime(), params::endTime());
		for (const auto& receiver_ptr : _world->getReceivers())
		{
			if (receiver_ptr->getMode() == OperationMode::CW_MODE ||
				receiver_ptr->getMode() == OperationMode::FMCW_MODE)
			{
				_pool.enqueue(processing::finalizeStreamingReceiver, receiver_ptr.get(), &_pool, _reporter, _output_dir,
							  _metadata_collector, streaming_sources);
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
		if (_reporter)
		{
			_reporter->report("Simulation complete", 100, 100);
		}
		LOG(Level::INFO, "Event-driven simulation loop finished.");
	}

	OutputMetadata runEventDrivenSim(World* world, pool::ThreadPool& pool,
									 const std::function<void(const std::string&, int, int)>& progress_callback,
									 const std::string& output_dir)
	{
		auto reporter = std::make_shared<ProgressReporter>(progress_callback);
		auto metadata_collector = std::make_shared<OutputMetadataCollector>(output_dir);
		SimulationEngine engine(world, pool, reporter, output_dir, metadata_collector);
		engine.run();
		return metadata_collector->snapshot();
	}
}
