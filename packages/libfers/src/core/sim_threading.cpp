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
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <optional>
#include <unordered_map>
#include <utility>

#include "logging.h"
#include "parameters.h"
#include "processing/finalizer.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
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
		RealType clippedStreamingSegmentEnd(const Transmitter& transmitter, const RealType segment_start,
											const RealType segment_end)
		{
			RealType clipped_end = std::min(params::endTime(), segment_end);
			if (const auto* fmcw = transmitter.getFmcwSignal(); (fmcw != nullptr) && fmcw->getChirpCount().has_value())
			{
				clipped_end =
					std::min(clipped_end,
							 segment_start + static_cast<RealType>(*fmcw->getChirpCount()) * fmcw->getChirpPeriod());
			}
			return clipped_end;
		}

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
				const RealType segment_end = clippedStreamingSegmentEnd(*transmitter, segment_start, params::endTime());
				if (timestamp >= segment_start && timestamp < segment_end)
				{
					return ActiveStreamingSource{
						.transmitter = transmitter, .segment_start = segment_start, .segment_end = segment_end};
				}
				return std::nullopt;
			}

			for (const auto& period : schedule)
			{
				const RealType active_start = std::max(params::startTime(), period.start);
				const RealType active_end = clippedStreamingSegmentEnd(*transmitter, period.start, period.end);
				if (timestamp >= active_start && timestamp < active_end)
				{
					return ActiveStreamingSource{
						.transmitter = transmitter, .segment_start = period.start, .segment_end = active_end};
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
			if (const auto* fmcw = transmitter_ptr->getFmcwSignal(); fmcw != nullptr)
			{
				const RealType duty_cycle = fmcw->getChirpDuration() / fmcw->getChirpPeriod();
				const RealType average_power = transmitter_ptr->getSignal()->getPower() * duty_cycle;
				const auto configured_count = fmcw->getChirpCount().has_value()
					? std::format("{}", *fmcw->getChirpCount())
					: std::string("unbounded");
				if (transmitter_ptr->getSchedule().empty())
				{
					LOG(Level::INFO,
						"FMCW transmitter '{}' B={} Hz T_c={} s T_rep={} s f_0={} Hz alpha={} Hz/s duty_cycle={} "
						"chirp_count={} average_power={} W",
						transmitter_ptr->getName(), fmcw->getChirpBandwidth(), fmcw->getChirpDuration(),
						fmcw->getChirpPeriod(), fmcw->getStartFrequencyOffset(), fmcw->getChirpRate(), duty_cycle,
						configured_count, average_power);
				}
				else
				{
					for (const auto& period : transmitter_ptr->getSchedule())
					{
						const RealType active_end =
							std::min(period.end,
									 fmcw->getChirpCount().has_value()
										 ? (period.start +
											static_cast<RealType>(*fmcw->getChirpCount()) * fmcw->getChirpPeriod())
										 : period.end);
						LOG(Level::INFO,
							"FMCW transmitter '{}' segment [{}, {}] B={} Hz T_c={} s T_rep={} s f_0={} Hz alpha={} "
							"Hz/s duty_cycle={} chirp_count={} average_power={} W",
							transmitter_ptr->getName(), period.start, active_end, fmcw->getChirpBandwidth(),
							fmcw->getChirpDuration(), fmcw->getChirpPeriod(), fmcw->getStartFrequencyOffset(),
							fmcw->getChirpRate(), duty_cycle, configured_count, average_power);
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

		std::vector<std::shared_ptr<timing::Timing>> timings;
		std::unordered_map<SimId, std::shared_ptr<timing::Timing>> unique_timings;

		for (const auto& transmitter_ptr : _world->getTransmitters())
		{
			if (!transmitter_ptr->isStreamingMode())
			{
				continue;
			}
			unique_timings.try_emplace(transmitter_ptr->getTiming()->getId(), transmitter_ptr->getTiming());
		}

		for (const auto& receiver_ptr : _world->getReceivers())
		{
			if ((receiver_ptr->getMode() != OperationMode::CW_MODE) &&
				(receiver_ptr->getMode() != OperationMode::FMCW_MODE))
			{
				continue;
			}
			unique_timings.try_emplace(receiver_ptr->getTiming()->getId(), receiver_ptr->getTiming());
		}

		timings.reserve(unique_timings.size());
		for (const auto& entry : unique_timings)
		{
			timings.push_back(entry.second);
		}

		_cw_phase_noise_lookup = std::make_unique<simulation::CwPhaseNoiseLookup>(
			simulation::CwPhaseNoiseLookup::build(timings, params::startTime(), params::endTime()));
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
		const auto start_index = static_cast<size_t>(std::ceil((t_current - params::startTime()) / dt_sim));
		const auto end_index = static_cast<size_t>(std::ceil((t_event - params::startTime()) / dt_sim));

		ensureCwPhaseNoiseLookup();

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
		}
	}

	ComplexType SimulationEngine::calculateStreamingSample(Receiver* rx, const RealType t_step,
														   const std::vector<ActiveStreamingSource>& streaming_sources,
														   ReceiverTrackerCache& tracker_cache) const
	{
		ComplexType total_sample{0.0, 0.0};
		for (std::size_t source_index = 0; source_index < streaming_sources.size(); ++source_index)
		{
			const auto& streaming_source = streaming_sources[source_index];
			if (!rx->checkFlag(Receiver::RecvFlag::FLAG_NODIRECT))
			{
				total_sample += simulation::calculateStreamingDirectPathContribution(
					streaming_source, rx, t_step, _cw_phase_noise_lookup.get(), &tracker_cache.direct[source_index]);
			}
			for (std::size_t target_index = 0; target_index < _world->getTargets().size(); ++target_index)
			{
				const auto& target_ptr = _world->getTargets()[target_index];
				total_sample += simulation::calculateStreamingReflectedPathContribution(
					streaming_source, rx, target_ptr.get(), t_step, _cw_phase_noise_lookup.get(),
					&tracker_cache.reflected[source_index][target_index]);
			}
		}
		return total_sample;
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
		// Keep the source as an in-flight candidate; per-path retarded-time gating zeros it once
		// rx_time - tau reaches segment_end.
		// TODO(perf): tx source now remains in active_streaming_transmitters forever. Find a clean solution for this.
	}

	void SimulationEngine::handleRxStreamingStart(Receiver* rx) { rx->setActive(true); }

	void SimulationEngine::handleRxStreamingEnd(Receiver* rx) { rx->setActive(false); }

	void SimulationEngine::updateProgress()
	{
		if (!_reporter)
		{
			return;
		}

		const RealType t_current = _world->getSimulationState().t_current;
		const RealType end_time = params::endTime();
		const int progress = static_cast<int>(t_current / end_time * 100.0);

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
				if (segment_start < segment_end && segment_start < end_time)
				{
					sources.push_back({.transmitter = transmitter_ptr.get(),
									   .segment_start = segment_start,
									   .segment_end = segment_end});
				}
			};

			if (transmitter_ptr->getSchedule().empty())
			{
				RealType segment_end = params::endTime();
				if (const auto* fmcw = transmitter_ptr->getFmcwSignal();
					(fmcw != nullptr) && fmcw->getChirpCount().has_value())
				{
					segment_end = std::min(segment_end,
										   params::startTime() +
											   static_cast<RealType>(*fmcw->getChirpCount()) * fmcw->getChirpPeriod());
				}
				append_candidate(params::startTime(), segment_end);
				continue;
			}

			for (const auto& period : transmitter_ptr->getSchedule())
			{
				RealType segment_end = std::min(params::endTime(), period.end);
				if (const auto* fmcw = transmitter_ptr->getFmcwSignal();
					(fmcw != nullptr) && fmcw->getChirpCount().has_value())
				{
					segment_end =
						std::min(segment_end,
								 period.start + static_cast<RealType>(*fmcw->getChirpCount()) * fmcw->getChirpPeriod());
				}
				append_candidate(period.start, segment_end);
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

		for (const auto& receiver_ptr : _world->getReceivers())
		{
			if (receiver_ptr->getMode() == OperationMode::CW_MODE ||
				receiver_ptr->getMode() == OperationMode::FMCW_MODE)
			{
				_pool.enqueue(processing::finalizeStreamingReceiver, receiver_ptr.get(), &_pool, _reporter, _output_dir,
							  _metadata_collector);
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
