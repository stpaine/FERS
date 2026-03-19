// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

/**
 * @file json_serializer.h
 * @brief Provides functions to serialize and deserialize the simulation world to/from JSON.
 *
 * This module is the primary data interchange layer between the C++ core engine and
 * the user interface. JSON was chosen as the format for its native support in web
 * technologies (like the React/TypeScript frontend), its human-readability, and its
 * lightweight nature. This serializer defines the contract for how C++ simulation
 * objects are represented in JSON, enabling the UI to read, modify, and write
 * back the entire simulation state.
 */

#pragma once

#include <nlohmann/json.hpp>
#include <random>

namespace core
{
	class World;
}

namespace serial
{
	/**
	 * @brief Serializes the entire simulation world into a nlohmann::json object.
	 *
	 * This function traverses the `core::World` object model and constructs a JSON
	 * representation. It is designed to produce a format that is convenient for the
	 * frontend to consume. This involves translating internal data formats (e.g.,
	 * angles in radians) to a more UI-friendly format (e.g., compass degrees) and
	 * restructuring complex object relationships (like monostatic radars) into simpler
	 * representations.
	 *
	 * @param world The world object to serialize.
	 * @return A nlohmann::json object representing the world.
	 */
	nlohmann::json world_to_json(const core::World& world);

	/**
	 * @brief Deserializes a nlohmann::json object and reconstructs the simulation world.
	 *
	 * This function is the counterpart to `world_to_json`. It performs a full state
	 * replacement by clearing the existing world and rebuilding it from the provided
	 * JSON. This "replace" strategy simplifies state management, guaranteeing that the
	 * C++ core is always perfectly synchronized with the state provided by the UI without
	 * requiring complex diffing or patching logic. It also handles re-seeding the master
	 * random number generator to ensure that loading a state also restores its
	 * deterministic behavior.
	 *
	 * @param j The json object to deserialize.
	 * @param world The world object to populate.
	 * @param masterSeeder A reference to the master random number generator, which will be re-seeded.
	 */
	void json_to_world(const nlohmann::json& j, core::World& world, std::mt19937& masterSeeder);
}
