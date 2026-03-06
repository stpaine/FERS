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
#include <format>
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
#include "world.h"

using logging::Level;
using radar::OperationMode;
using radar::Receiver;
using radar::Transmitter;

namespace core
{
	SimulationEngine::SimulationEngine(World* world, pool::ThreadPool& pool,
									   std::shared_ptr<ProgressReporter> reporter) :
		_world(world), _pool(pool), _reporter(std::move(reporter)), _last_report_time(std::chrono::steady_clock::now())
	{
	}

	void SimulationEngine::run()
	{
		if (_reporter)
		{
			_reporter->report("Initializing event-driven simulation...", 0, 100);
		}

		initializeFinalizers();

		LOG(Level::INFO, "Starting unified event-driven simulation loop.");

		auto& event_queue = _world->getEventQueue();
		auto& state = _world->getSimulationState();
		const RealType end_time = params::endTime();

		while (!event_queue.empty() && state.t_current <= end_time)
		{
			const Event event = event_queue.top();
			event_queue.pop();

			processCwPhysics(event.timestamp);

			state.t_current = event.timestamp;

			processEvent(event);
			updateProgress();
		}

		shutdown();
	}

	void SimulationEngine::initializeFinalizers()
	{
		for (const auto& receiver_ptr : _world->getReceivers())
		{
			if (receiver_ptr->getMode() == OperationMode::PULSED_MODE)
			{
				_finalizer_threads.emplace_back(processing::runPulsedFinalizer, receiver_ptr.get(),
												&_world->getTargets(), _reporter);
			}
		}
	}

	void SimulationEngine::processCwPhysics(const RealType t_event)
	{
		auto& [t_current, active_cw_transmitters] = _world->getSimulationState();

		if (t_event <= t_current)
		{
			return;
		}

		const RealType dt_sim = 1.0 / (params::rate() * params::oversampleRatio());
		const auto start_index = static_cast<size_t>(std::ceil((t_current - params::startTime()) / dt_sim));
		const auto end_index = static_cast<size_t>(std::ceil((t_event - params::startTime()) / dt_sim));

		for (size_t sample_index = start_index; sample_index < end_index; ++sample_index)
		{
			const RealType t_step = params::startTime() + static_cast<RealType>(sample_index) * dt_sim;

			for (const auto& receiver_ptr : _world->getReceivers())
			{
				if (receiver_ptr->getMode() == OperationMode::CW_MODE && receiver_ptr->isActive())
				{
					ComplexType sample = calculateCwSample(receiver_ptr.get(), t_step, active_cw_transmitters);
					receiver_ptr->setCwSample(sample_index, sample);
				}
			}
		}
	}

	ComplexType SimulationEngine::calculateCwSample(Receiver* rx, const RealType t_step,
													const std::vector<Transmitter*>& cw_sources) const
	{
		ComplexType total_sample{0.0, 0.0};
		for (const auto& cw_source : cw_sources)
		{
			if (!rx->checkFlag(Receiver::RecvFlag::FLAG_NODIRECT))
			{
				total_sample += simulation::calculateDirectPathContribution(cw_source, rx, t_step);
			}
			for (const auto& target_ptr : _world->getTargets())
			{
				total_sample += simulation::calculateReflectedPathContribution(cw_source, rx, target_ptr.get(), t_step);
			}
		}
		return total_sample;
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
		case EventType::TX_CW_START:
			handleTxCwStart(static_cast<Transmitter*>(event.source_object));
			break;
		case EventType::TX_CW_END:
			handleTxCwEnd(static_cast<Transmitter*>(event.source_object));
			break;
		case EventType::RX_CW_START:
			handleRxCwStart(static_cast<Receiver*>(event.source_object));
			break;
		case EventType::RX_CW_END:
			handleRxCwEnd(static_cast<Receiver*>(event.source_object));
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
		const auto& active_cw_transmitters = _world->getSimulationState().active_cw_transmitters;

		RenderingJob job{.ideal_start_time = t_event - rx->getWindowLength(),
						 .duration = rx->getWindowLength(),
						 .responses = rx->drainInbox(),
						 .active_cw_sources = active_cw_transmitters};

		rx->enqueueFinalizerJob(std::move(job));

		const RealType next_theoretical = t_event - rx->getWindowLength() + 1.0 / rx->getWindowPrf();
		if (const auto next_start = rx->getNextWindowTime(next_theoretical);
			next_start && *next_start <= params::endTime())
		{
			_world->getEventQueue().push({*next_start, EventType::RX_PULSED_WINDOW_START, rx});
		}
	}

	void SimulationEngine::handleTxCwStart(Transmitter* tx)
	{
		_world->getSimulationState().active_cw_transmitters.push_back(tx);
	}

	void SimulationEngine::handleTxCwEnd(Transmitter* tx)
	{
		auto& cw_txs = _world->getSimulationState().active_cw_transmitters;
		std::erase(cw_txs, tx);
	}

	void SimulationEngine::handleRxCwStart(Receiver* rx) { rx->setActive(true); }

	void SimulationEngine::handleRxCwEnd(Receiver* rx) { rx->setActive(false); }

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

	void SimulationEngine::shutdown()
	{
		LOG(Level::INFO, "Main simulation loop finished. Waiting for finalization tasks...");
		if (_reporter)
		{
			_reporter->report("Main simulation finished. Waiting for data export...", 100, 100);
		}

		for (const auto& receiver_ptr : _world->getReceivers())
		{
			if (receiver_ptr->getMode() == OperationMode::CW_MODE)
			{
				_pool.enqueue(processing::finalizeCwReceiver, receiver_ptr.get(), &_pool, _reporter);
			}
			else if (receiver_ptr->getMode() == OperationMode::PULSED_MODE)
			{
				RenderingJob shutdown_job{.duration = -1.0};
				receiver_ptr->enqueueFinalizerJob(std::move(shutdown_job));
			}
		}

		_pool.wait();

		LOG(Level::INFO, "All finalization tasks complete.");
		if (_reporter)
		{
			_reporter->report("Simulation complete", 100, 100);
		}
		LOG(Level::INFO, "Event-driven simulation loop finished.");
	}

	void runEventDrivenSim(World* world, pool::ThreadPool& pool,
						   const std::function<void(const std::string&, int, int)>& progress_callback)
	{
		auto reporter = std::make_shared<ProgressReporter>(progress_callback);
		SimulationEngine engine(world, pool, reporter);
		engine.run();
	}
}
