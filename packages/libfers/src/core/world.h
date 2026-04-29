// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file world.h
 * @brief Header file for the World class in the simulator.
 */

#pragma once

#include <memory>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/sim_events.h"
#include "core/sim_id.h"
#include "core/simulation_state.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "signal/radar_signal.h"
#include "timing/prototype_timing.h"

namespace core
{
	/**
	 * @class World
	 * @brief The World class manages the simulator environment.
	 */
	class World
	{
	public:
		World() = default;

		~World() noexcept = default;

		World(const World&) = delete;

		World& operator=(const World&) = delete;

		World(World&&) = delete;

		World& operator=(World&&) = delete;

		/**
		 * @brief Adds a radar platform to the simulation world.
		 *
		 * @param plat A unique pointer to a Platform object.
		 */
		void add(std::unique_ptr<radar::Platform> plat) noexcept;

		/**
		 * @brief Adds a radar transmitter to the simulation world.
		 *
		 * @param trans A unique pointer to a Transmitter object.
		 */
		void add(std::unique_ptr<radar::Transmitter> trans) noexcept;

		/**
		 * @brief Adds a radar receiver to the simulation world.
		 *
		 * @param recv A unique pointer to a Receiver object.
		 */
		void add(std::unique_ptr<radar::Receiver> recv) noexcept;

		/**
		 * @brief Adds a radar target to the simulation world.
		 *
		 * @param target A unique pointer to a Target object.
		 */
		void add(std::unique_ptr<radar::Target> target) noexcept;

		/**
		 * @brief Adds a radar signal (waveform) to the simulation world.
		 *
		 * @param waveform A unique pointer to a RadarSignal object.
		 * @throws std::runtime_error if a waveform with the same ID already exists.
		 */
		void add(std::unique_ptr<fers_signal::RadarSignal> waveform);

		/**
		 * @brief Adds an antenna to the simulation world.
		 *
		 * @param antenna A unique pointer to an Antenna object.
		 * @throws std::runtime_error if an antenna with the same ID already exists.
		 */
		void add(std::unique_ptr<antenna::Antenna> antenna);

		/**
		 * @brief Adds a timing source to the simulation world.
		 *
		 * @param timing A unique pointer to a PrototypeTiming object.
		 * @throws std::runtime_error if a timing source with the same ID already exists.
		 */
		void add(std::unique_ptr<timing::PrototypeTiming> timing);

		/**
		 * @brief Finds a radar signal by ID.
		 *
		 * @param id The ID of the radar signal to find.
		 * @return A pointer to the RadarSignal if found, or nullptr if not found.
		 */
		[[nodiscard]] fers_signal::RadarSignal* findWaveform(const SimId id);

		/**
		 * @brief Finds an antenna by ID.
		 *
		 * @param id The ID of the antenna to find.
		 * @return A pointer to the Antenna if found, or nullptr if not found.
		 */
		[[nodiscard]] antenna::Antenna* findAntenna(const SimId id);

		/**
		 * @brief Finds a timing source by ID.
		 *
		 * @param id The ID of the timing source to find.
		 * @return A pointer to the PrototypeTiming if found, or nullptr if not found.
		 */
		[[nodiscard]] timing::PrototypeTiming* findTiming(const SimId id);

		/**
		 * @brief Finds a platform by ID.
		 *
		 * @param id The ID of the platform to find.
		 * @return A pointer to the Platform if found, or nullptr if not found.
		 */
		[[nodiscard]] radar::Platform* findPlatform(const SimId id);

		/**
		 * @brief Finds a transmitter by ID.
		 *
		 * @param id The ID of the transmitter to find.
		 * @return A pointer to the Transmitter if found, or nullptr if not found.
		 */
		[[nodiscard]] radar::Transmitter* findTransmitter(const SimId id);

		/**
		 * @brief Finds a transmitter by name.
		 *
		 * @param name The transmitter name to find.
		 * @return A pointer to the Transmitter if found, or nullptr if not found.
		 */
		[[nodiscard]] radar::Transmitter* findTransmitterByName(const std::string& name);

		/**
		 * @brief Finds a receiver by ID.
		 *
		 * @param id The ID of the receiver to find.
		 * @return A pointer to the Receiver if found, or nullptr if not found.
		 */
		[[nodiscard]] radar::Receiver* findReceiver(const SimId id);

		/**
		 * @brief Finds a waveform by name.
		 *
		 * @param name The waveform name to find.
		 * @return A pointer to the RadarSignal if found, or nullptr if not found.
		 */
		[[nodiscard]] fers_signal::RadarSignal* findWaveformByName(const std::string& name);

		/**
		 * @brief Finds a target by ID.
		 *
		 * @param id The ID of the target to find.
		 * @return A pointer to the Target if found, or nullptr if not found.
		 */
		[[nodiscard]] radar::Target* findTarget(const SimId id);

		/**
		 * @brief Replaces an existing target, updating internal pointers.
		 * @param target Unique pointer to the new target.
		 */
		void replace(std::unique_ptr<radar::Target> target);

		/**
		 * @brief Replaces an existing antenna, updating internal pointers.
		 * @param antenna Unique pointer to the new antenna.
		 */
		void replace(std::unique_ptr<antenna::Antenna> antenna);

		/**
		 * @brief Replaces an existing waveform, updating internal pointers.
		 * @param waveform Unique pointer to the new waveform.
		 */
		void replace(std::unique_ptr<fers_signal::RadarSignal> waveform);

		/**
		 * @brief Replaces an existing timing prototype and refreshes dependent radar timing models.
		 * @param timing Unique pointer to the new timing prototype.
		 */
		void replace(std::unique_ptr<timing::PrototypeTiming> timing);

		/**
		 * @brief Retrieves the list of platforms.
		 *
		 * @return A const reference to a vector of unique pointers to Platform objects.
		 */
		[[nodiscard]] const std::vector<std::unique_ptr<radar::Platform>>& getPlatforms() const noexcept
		{
			return _platforms;
		}

		/**
		 * @brief Retrieves the list of radar targets.
		 *
		 * @return A const reference to a vector of unique pointers to Target objects.
		 */
		[[nodiscard]] const std::vector<std::unique_ptr<radar::Target>>& getTargets() const noexcept
		{
			return _targets;
		}

		/**
		 * @brief Retrieves the list of radar receivers.
		 *
		 * @return A const reference to a vector of unique pointers to Receiver objects.
		 */
		[[nodiscard]] const std::vector<std::unique_ptr<radar::Receiver>>& getReceivers() const noexcept
		{
			return _receivers;
		}

		/**
		 * @brief Retrieves the list of radar transmitters.
		 *
		 * @return A const reference to a vector of unique pointers to Transmitter objects.
		 */
		[[nodiscard]] const std::vector<std::unique_ptr<radar::Transmitter>>& getTransmitters() const noexcept
		{
			return _transmitters;
		}

		/**
		 * @brief Retrieves the map of radar signals (waveforms).
		 * @return A const reference to the map of signal names to RadarSignal objects.
		 */
		[[nodiscard]] const std::unordered_map<SimId, std::unique_ptr<fers_signal::RadarSignal>>&
		getWaveforms() const noexcept
		{
			return _waveforms;
		}

		/**
		 * @brief Retrieves the map of antennas.
		 * @return A const reference to the map of antenna names to Antenna objects.
		 */
		[[nodiscard]] const std::unordered_map<SimId, std::unique_ptr<antenna::Antenna>>& getAntennas() const noexcept
		{
			return _antennas;
		}

		/**
		 * @brief Retrieves the map of timing prototypes.
		 * @return A const reference to the map of timing names to PrototypeTiming objects.
		 */
		[[nodiscard]] const std::unordered_map<SimId, std::unique_ptr<timing::PrototypeTiming>>&
		getTimings() const noexcept
		{
			return _timings;
		}

		/**
		 * @brief Clears all objects and assets from the simulation world.
		 */
		void clear() noexcept;

		/**
		 * @brief Exchanges all owned world state with another world.
		 */
		void swap(World& other) noexcept;

		/**
		 * @brief Populates the event queue with the initial events for the simulation.
		 * This method should be called after all simulation objects have been parsed and added to the world.
		 */
		void scheduleInitialEvents();

		/**
		 * @brief Resolves and validates receiver FMCW dechirp references after all components are loaded.
		 *
		 * @throws std::runtime_error if a dechirped receiver has an invalid or inactive LO reference.
		 */
		void resolveReceiverDechirpReferences();

		/**
		 * @brief Dumps the current state of the event queue to a string for debugging.
		 * @return A formatted string representing the contents of the event queue.
		 */
		[[nodiscard]] std::string dumpEventQueue() const;

		/**
		 * @brief Gets a mutable reference to the global event queue.
		 * @return A reference to the priority queue of events.
		 */
		[[nodiscard]] std::priority_queue<Event, std::vector<Event>, EventComparator>& getEventQueue() noexcept
		{
			return _event_queue;
		}

		/**
		 * @brief Gets a mutable reference to the global simulation state.
		 * @return A reference to the SimulationState object.
		 */
		[[nodiscard]] SimulationState& getSimulationState() noexcept { return _simulation_state; }

	private:
		std::vector<std::unique_ptr<radar::Platform>> _platforms; ///< Owned radar platforms.

		std::vector<std::unique_ptr<radar::Transmitter>> _transmitters; ///< Owned transmitters.

		std::vector<std::unique_ptr<radar::Receiver>> _receivers; ///< Owned receivers.

		std::vector<std::unique_ptr<radar::Target>> _targets; ///< Owned targets.

		std::unordered_map<SimId, std::unique_ptr<fers_signal::RadarSignal>> _waveforms; ///< Owned waveform assets.

		std::unordered_map<SimId, std::unique_ptr<antenna::Antenna>> _antennas; ///< Owned antenna assets.

		std::unordered_map<SimId, std::unique_ptr<timing::PrototypeTiming>> _timings; ///< Owned timing prototypes.

		std::priority_queue<Event, std::vector<Event>, EventComparator> _event_queue; ///< Pending simulation events.

		SimulationState _simulation_state; ///< Mutable runtime simulation state.
	};
}
