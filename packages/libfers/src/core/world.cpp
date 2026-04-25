// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file world.cpp
 * @brief Implementation of the World class for the radar simulation environment.
 */

#include "world.h"

#include <iomanip>
#include <limits>
#include <sstream>
#include <unordered_map>

#include "antenna/antenna_factory.h"
#include "core/sim_events.h"
#include "core/sim_id.h"
#include "parameters.h"
#include "radar/radar_obj.h"
#include "signal/radar_signal.h"
#include "timing/prototype_timing.h"
#include "timing/timing.h"

using antenna::Antenna;
using fers_signal::RadarSignal;
using radar::Platform;
using radar::Receiver;
using radar::Target;
using radar::Transmitter;
using timing::PrototypeTiming;

namespace core
{
	void World::add(std::unique_ptr<Platform> plat) noexcept { _platforms.push_back(std::move(plat)); }

	void World::add(std::unique_ptr<Transmitter> trans) noexcept { _transmitters.push_back(std::move(trans)); }

	void World::add(std::unique_ptr<Receiver> recv) noexcept { _receivers.push_back(std::move(recv)); }

	void World::add(std::unique_ptr<Target> target) noexcept { _targets.push_back(std::move(target)); }

	void World::add(std::unique_ptr<RadarSignal> waveform)
	{
		const SimId id = waveform->getId();
		if (_waveforms.contains(id))
		{
			throw std::runtime_error("A waveform with the ID " + std::to_string(id) + " already exists.");
		}
		_waveforms[id] = std::move(waveform);
	}

	void World::add(std::unique_ptr<Antenna> antenna)
	{
		const SimId id = antenna->getId();
		if (_antennas.contains(id))
		{
			throw std::runtime_error("An antenna with the ID " + std::to_string(id) + " already exists.");
		}
		_antennas[id] = std::move(antenna);
	}

	void World::add(std::unique_ptr<PrototypeTiming> timing)
	{
		const SimId id = timing->getId();
		if (_timings.contains(id))
		{
			throw std::runtime_error("A timing source with the ID " + std::to_string(id) + " already exists.");
		}
		_timings[id] = std::move(timing);
	}

	RadarSignal* World::findWaveform(const SimId id)
	{
		const auto it = _waveforms.find(id);
		return it != _waveforms.end() ? it->second.get() : nullptr;
	}

	Antenna* World::findAntenna(const SimId id)
	{
		const auto it = _antennas.find(id);
		return it != _antennas.end() ? it->second.get() : nullptr;
	}

	PrototypeTiming* World::findTiming(const SimId id)
	{
		const auto it = _timings.find(id);
		return it != _timings.end() ? it->second.get() : nullptr;
	}

	Platform* World::findPlatform(const SimId id)
	{
		for (auto& p : _platforms)
		{
			if (p->getId() == id)
				return p.get();
		}
		return nullptr;
	}

	Transmitter* World::findTransmitter(const SimId id)
	{
		for (auto& tx : _transmitters)
			if (tx->getId() == id)
				return tx.get();
		return nullptr;
	}

	Receiver* World::findReceiver(const SimId id)
	{
		for (auto& rx : _receivers)
			if (rx->getId() == id)
				return rx.get();
		return nullptr;
	}

	Target* World::findTarget(const SimId id)
	{
		for (auto& tgt : _targets)
			if (tgt->getId() == id)
				return tgt.get();
		return nullptr;
	}

	void World::replace(std::unique_ptr<Target> target)
	{
		const SimId id = target->getId();
		for (auto& t : _targets)
		{
			if (t->getId() == id)
			{
				t = std::move(target);
				return;
			}
		}
		_targets.push_back(std::move(target));
	}

	void World::replace(std::unique_ptr<Antenna> antenna)
	{
		const SimId id = antenna->getId();
		const Antenna* new_ptr = antenna.get();

		std::unique_ptr<Antenna> old_owned;
		const Antenna* old_ptr = nullptr;

		if (auto it = _antennas.find(id); it != _antennas.end())
		{
			old_owned = std::move(it->second);
			old_ptr = old_owned.get();
			it->second = std::move(antenna);
		}
		else
		{
			_antennas[id] = std::move(antenna);
		}

		if ((old_ptr != nullptr) && old_ptr != new_ptr)
		{
			for (auto& tx : _transmitters)
				if (tx->getAntenna() == old_ptr)
					tx->setAntenna(new_ptr);

			for (auto& rx : _receivers)
				if (rx->getAntenna() == old_ptr)
					rx->setAntenna(new_ptr);
		}
	}

	void World::replace(std::unique_ptr<RadarSignal> waveform)
	{
		const SimId id = waveform->getId();
		RadarSignal* new_ptr = waveform.get();

		std::unique_ptr<RadarSignal> old_owned;
		const RadarSignal* old_ptr = nullptr;

		if (auto it = _waveforms.find(id); it != _waveforms.end())
		{
			old_owned = std::move(it->second);
			old_ptr = old_owned.get();
			it->second = std::move(waveform);
		}
		else
		{
			_waveforms[id] = std::move(waveform);
		}

		if ((old_ptr != nullptr) && old_ptr != new_ptr)
		{
			for (auto& tx : _transmitters)
				if (tx->getSignal() == old_ptr)
					tx->setSignal(new_ptr);
		}
	}

	void World::replace(std::unique_ptr<PrototypeTiming> timing)
	{
		const SimId id = timing->getId();
		const PrototypeTiming* new_ptr = timing.get();

		std::unique_ptr<PrototypeTiming> old_owned;
		const PrototypeTiming* old_ptr = nullptr;

		if (auto it = _timings.find(id); it != _timings.end())
		{
			old_owned = std::move(it->second);
			old_ptr = old_owned.get();
			it->second = std::move(timing);
		}
		else
		{
			_timings[id] = std::move(timing);
		}

		std::unordered_map<const timing::Timing*, std::shared_ptr<timing::Timing>> refreshed_instances;
		auto refresh_timing = [id, new_ptr, &refreshed_instances](auto& radar_obj)
		{
			const auto current_timing = radar_obj->getTiming();
			if (!current_timing || (current_timing->getId() != id))
			{
				return;
			}

			const timing::Timing* const timing_key = current_timing.get();
			const auto [it, inserted] = refreshed_instances.try_emplace(timing_key);
			if (inserted)
			{
				auto refreshed =
					std::make_shared<timing::Timing>(new_ptr->getName(), current_timing->getSeed(), new_ptr->getId());
				refreshed->initializeModel(new_ptr);
				it->second = std::move(refreshed);
			}
			radar_obj->setTiming(it->second);
		};

		if ((old_ptr != nullptr) && old_ptr != new_ptr)
		{
			for (auto& tx : _transmitters)
				refresh_timing(tx);

			for (auto& rx : _receivers)
				refresh_timing(rx);
		}
	}

	void World::clear() noexcept
	{
		_platforms.clear();
		_transmitters.clear();
		_receivers.clear();
		_targets.clear();
		_waveforms.clear();
		_antennas.clear();
		_timings.clear();
		_event_queue = {};
		_simulation_state = {};
	}

	void World::scheduleInitialEvents()
	{
		const RealType sim_start = params::startTime();
		const RealType sim_end = params::endTime();

		for (const auto& transmitter : _transmitters)
		{
			if (transmitter->getMode() == radar::OperationMode::PULSED_MODE)
			{
				// Find the first valid pulse time starting from the simulation start time.
				if (auto start_time = transmitter->getNextPulseTime(sim_start); start_time)
				{
					if (*start_time <= sim_end)
					{
						_event_queue.push({*start_time, EventType::TX_PULSED_START, transmitter.get()});
					}
				}
			}
			else
			{
				const auto& schedule = transmitter->getSchedule();
				const auto clip_streaming_end = [&](const RealType start, const RealType end)
				{
					RealType clipped_end = std::min(sim_end, end);
					// FMCW chirp_count is a per-schedule-segment cap. Each scheduled period starts its own count.
					if (const auto* fmcw = transmitter->getFmcwSignal();
						(fmcw != nullptr) && fmcw->getChirpCount().has_value())
					{
						clipped_end =
							std::min(clipped_end,
									 start + static_cast<RealType>(*fmcw->getChirpCount()) * fmcw->getChirpPeriod());
					}
					return clipped_end;
				};
				if (schedule.empty())
				{
					const RealType end = clip_streaming_end(sim_start, sim_end);
					if (sim_start < end)
					{
						_event_queue.push({sim_start, EventType::TX_STREAMING_START, transmitter.get()});
						_event_queue.push({end, EventType::TX_STREAMING_END, transmitter.get()});
					}
				}
				else
				{
					for (const auto& period : schedule)
					{
						const RealType start = std::max(sim_start, period.start);
						const RealType end = clip_streaming_end(period.start, period.end);

						if (start < end)
						{
							_event_queue.push({start, EventType::TX_STREAMING_START, transmitter.get()});
							_event_queue.push({end, EventType::TX_STREAMING_END, transmitter.get()});
						}
					}
				}
			}
		}

		for (const auto& receiver : _receivers)
		{
			if (receiver->getMode() == radar::OperationMode::PULSED_MODE)
			{
				// Schedule the first receive window checking against schedule
				const RealType nominal_start = receiver->getWindowStart(0);
				if (auto start = receiver->getNextWindowTime(nominal_start); start && *start < params::endTime())
				{
					_event_queue.push({*start, EventType::RX_PULSED_WINDOW_START, receiver.get()});
				}
			}
			else
			{
				const auto& schedule = receiver->getSchedule();
				if (schedule.empty())
				{
					_event_queue.push({params::startTime(), EventType::RX_STREAMING_START, receiver.get()});
					_event_queue.push({params::endTime(), EventType::RX_STREAMING_END, receiver.get()});
				}
				else
				{
					for (const auto& period : schedule)
					{
						const RealType start = std::max(params::startTime(), period.start);
						const RealType end = std::min(params::endTime(), period.end);
						if (start < end)
						{
							_event_queue.push({start, EventType::RX_STREAMING_START, receiver.get()});
							_event_queue.push({end, EventType::RX_STREAMING_END, receiver.get()});
						}
					}
				}
			}
		}
	}

	std::string World::dumpEventQueue() const
	{
		if (_event_queue.empty())
		{
			return "Event Queue is empty.\n";
		}

		std::stringstream ss;
		ss << std::fixed << std::setprecision(6);

		const std::string separator = "--------------------------------------------------------------------";
		const std::string title = "| Event Queue Contents (" + std::to_string(_event_queue.size()) + " events)";
		if (separator.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
		{
			throw std::runtime_error("Separator width exceeds stream formatting limits.");
		}
		const int title_width = static_cast<int>(separator.size()) - 1;

		ss << separator << "\n"
		   << std::left << std::setw(title_width) << title << "|\n"
		   << separator << "\n"
		   << "| " << std::left << std::setw(12) << "Timestamp" << " | " << std::setw(21) << "Event Type" << " | "
		   << std::setw(25) << "Source Object" << " |\n"
		   << separator << "\n";

		auto queue_copy = _event_queue;

		while (!queue_copy.empty())
		{
			const auto [timestamp, event_type, source_object] = queue_copy.top();
			queue_copy.pop();

			ss << "| " << std::right << std::setw(12) << timestamp << " | " << std::left << std::setw(21)
			   << toString(event_type) << " | " << std::left << std::setw(25) << source_object->getName() << " |\n";
		}
		ss << separator << "\n";

		return ss.str();
	}
}
