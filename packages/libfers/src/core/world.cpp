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

#include <algorithm>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

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

	void World::add(std::unique_ptr<Transmitter> trans) noexcept
	{
		_transmitters_by_name[trans->getName()] = trans.get();
		_transmitters.push_back(std::move(trans));
	}

	void World::add(std::unique_ptr<Receiver> recv) noexcept { _receivers.push_back(std::move(recv)); }

	void World::add(std::unique_ptr<Target> target) noexcept { _targets.push_back(std::move(target)); }

	void World::add(std::unique_ptr<RadarSignal> waveform)
	{
		const SimId id = waveform->getId();
		if (_waveforms.contains(id))
		{
			throw std::runtime_error("A waveform with the ID " + std::to_string(id) + " already exists.");
		}
		_waveform_ids_by_name[waveform->getName()] = id;
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

	Transmitter* World::findTransmitterByName(const std::string& name)
	{
		const auto it = _transmitters_by_name.find(name);
		return it != _transmitters_by_name.end() ? it->second : nullptr;
	}

	Receiver* World::findReceiver(const SimId id)
	{
		for (auto& rx : _receivers)
			if (rx->getId() == id)
				return rx.get();
		return nullptr;
	}

	RadarSignal* World::findWaveformByName(const std::string& name)
	{
		const auto it = _waveform_ids_by_name.find(name);
		return it != _waveform_ids_by_name.end() ? findWaveform(it->second) : nullptr;
	}

	RealType World::earliestPhaseNoiseLookupStart() const
	{
		const auto include_streaming_interval_start = [](std::optional<RealType>& earliest,
														 const RealType segment_start, const RealType segment_end,
														 const bool allow_pre_start)
		{
			const RealType sim_start = params::startTime();
			const RealType sim_end = params::endTime();
			if (segment_end <= sim_start || segment_start >= sim_end)
			{
				return;
			}

			const RealType required_start =
				allow_pre_start && segment_start < sim_start ? segment_start : std::max(sim_start, segment_start);
			earliest = earliest.has_value() ? std::min(*earliest, required_start) : required_start;
		};

		std::optional<RealType> earliest;
		for (const auto& transmitter_ptr : _transmitters)
		{
			if (transmitter_ptr == nullptr || !transmitter_ptr->isStreamingMode())
			{
				continue;
			}

			const auto& schedule = transmitter_ptr->getSchedule();
			if (schedule.empty())
			{
				include_streaming_interval_start(earliest, params::startTime(), params::endTime(), false);
				continue;
			}

			for (const auto& period : schedule)
			{
				include_streaming_interval_start(earliest, period.start, period.end, true);
			}
		}

		for (const auto& receiver_ptr : _receivers)
		{
			if (receiver_ptr == nullptr ||
				(receiver_ptr->getMode() != radar::OperationMode::CW_MODE &&
				 receiver_ptr->getMode() != radar::OperationMode::FMCW_MODE))
			{
				continue;
			}

			const auto& schedule = receiver_ptr->getSchedule();
			if (schedule.empty())
			{
				include_streaming_interval_start(earliest, params::startTime(), params::endTime(), false);
				continue;
			}

			for (const auto& period : schedule)
			{
				include_streaming_interval_start(earliest, period.start, period.end, false);
			}
		}

		return earliest.value_or(params::startTime());
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
		const std::string new_name = waveform->getName();

		std::unique_ptr<RadarSignal> old_owned;
		const RadarSignal* old_ptr = nullptr;

		if (auto it = _waveforms.find(id); it != _waveforms.end())
		{
			old_owned = std::move(it->second);
			old_ptr = old_owned.get();
			_waveform_ids_by_name.erase(old_owned->getName());
			_waveform_ids_by_name[new_name] = id;
			it->second = std::move(waveform);
		}
		else
		{
			_waveform_ids_by_name[new_name] = id;
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
		_transmitters_by_name.clear();
		_receivers.clear();
		_targets.clear();
		_waveforms.clear();
		_waveform_ids_by_name.clear();
		_antennas.clear();
		_timings.clear();
		_event_queue = {};
		_simulation_state = {};
	}

	void World::swap(World& other) noexcept
	{
		using std::swap;

		_platforms.swap(other._platforms);
		_transmitters.swap(other._transmitters);
		_receivers.swap(other._receivers);
		_targets.swap(other._targets);
		_waveforms.swap(other._waveforms);
		_waveform_ids_by_name.swap(other._waveform_ids_by_name);
		_transmitters_by_name.swap(other._transmitters_by_name);
		_antennas.swap(other._antennas);
		_timings.swap(other._timings);
		_event_queue.swap(other._event_queue);
		swap(_simulation_state, other._simulation_state);
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
				if (schedule.empty())
				{
					const RealType end = makeActiveSource(transmitter.get(), sim_start, sim_end).segment_end;
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
						const RealType end =
							makeActiveSource(transmitter.get(), period.start, std::min(sim_end, period.end))
								.segment_end;

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

	void World::resolveReceiverDechirpReferences()
	{
		const auto append_transmitter_sources = [](const Transmitter* const tx)
		{
			std::vector<ActiveStreamingSource> sources;
			if (tx->getSchedule().empty())
			{
				auto source = makeActiveSource(tx, params::startTime(), params::endTime());
				if (source.segment_start < source.segment_end)
				{
					sources.push_back(source);
				}
				return sources;
			}

			for (const auto& period : tx->getSchedule())
			{
				auto source = makeActiveSource(tx, period.start, std::min(params::endTime(), period.end));
				if (source.segment_start < source.segment_end && source.segment_end > params::startTime())
				{
					sources.push_back(source);
				}
			}
			return sources;
		};

		const auto append_waveform_sources = [](const RadarSignal* const waveform, const Receiver* const rx)
		{
			std::vector<ActiveStreamingSource> sources;
			if (rx->getSchedule().empty())
			{
				auto source = makeActiveSourceFromWaveform(waveform, params::startTime(), params::endTime());
				if (source.segment_start < source.segment_end)
				{
					sources.push_back(source);
				}
				return sources;
			}

			for (const auto& period : rx->getSchedule())
			{
				auto source =
					makeActiveSourceFromWaveform(waveform, period.start, std::min(params::endTime(), period.end));
				if (source.segment_start < source.segment_end && source.segment_end > params::startTime())
				{
					sources.push_back(source);
				}
			}
			return sources;
		};

		const auto validate_transmitter = [](const Transmitter* const tx, const std::string& owner)
		{
			if (tx == nullptr)
			{
				throw std::runtime_error(owner + " references a missing dechirp transmitter.");
			}
			if (tx->getMode() != radar::OperationMode::FMCW_MODE || tx->getSignal() == nullptr ||
				!tx->getSignal()->isFmcwFamily())
			{
				throw std::runtime_error(owner + " dechirp reference transmitter '" + tx->getName() +
										 "' must be an FMCW transmitter with an FMCW waveform.");
			}
		};

		for (const auto& rx_ptr : _receivers)
		{
			auto& rx = *rx_ptr;
			rx.clearResolvedDechirpSources();
			if (!rx.isDechirpEnabled())
			{
				continue;
			}
			if (rx.getMode() != radar::OperationMode::FMCW_MODE)
			{
				throw std::runtime_error("Receiver '" + rx.getName() + "' enables dechirping outside FMCW mode.");
			}

			auto reference = rx.getDechirpReference();
			std::vector<ActiveStreamingSource> sources;
			const std::string owner = "Receiver '" + rx.getName() + "'";
			switch (reference.source)
			{
			case Receiver::DechirpReferenceSource::Attached:
				{
					const auto* const tx = dynamic_cast<const Transmitter*>(rx.getAttached());
					validate_transmitter(tx, owner);
					reference.transmitter_id = tx->getId();
					reference.transmitter_name = tx->getName();
					sources = append_transmitter_sources(tx);
					break;
				}
			case Receiver::DechirpReferenceSource::Transmitter:
				{
					auto* const tx = findTransmitterByName(reference.name);
					validate_transmitter(tx, owner);
					reference.transmitter_id = tx->getId();
					reference.transmitter_name = tx->getName();
					sources = append_transmitter_sources(tx);
					break;
				}
			case Receiver::DechirpReferenceSource::Custom:
				{
					auto* const waveform = findWaveformByName(reference.name);
					if (waveform == nullptr || !waveform->isFmcwFamily())
					{
						throw std::runtime_error(owner + " custom dechirp reference waveform '" + reference.name +
												 "' must be a top-level FMCW waveform.");
					}
					reference.waveform_id = waveform->getId();
					reference.waveform_name = waveform->getName();
					sources = append_waveform_sources(waveform, &rx);
					break;
				}
			case Receiver::DechirpReferenceSource::None:
				throw std::runtime_error(owner + " enables dechirping without a dechirp reference.");
			}

			if (sources.empty())
			{
				throw std::runtime_error(owner + " dechirp reference has no active LO segments in the simulation.");
			}
			std::sort(sources.begin(), sources.end(),
					  [](const ActiveStreamingSource& lhs, const ActiveStreamingSource& rhs)
					  { return lhs.segment_start < rhs.segment_start; });
			rx.setDechirpReference(std::move(reference));
			rx.setResolvedDechirpSources(std::move(sources));
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
