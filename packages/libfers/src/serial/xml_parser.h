// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file xml_parser.h
 * @brief High-level facade for parsing XML configuration files into the FERS simulation environment.
 *
 * This header provides the main entry points for loading FERS scenario definitions
 * from XML files or strings. It handles validation, file inclusion, asset loading,
 * and pushing the parsed configuration into the global simulation state.
 */

#pragma once

#include <random>
#include <string>

namespace core
{
	class World;
}

namespace serial
{
	/**
	 * @brief Parses a simulation configuration from an XML file.
	 *
	 * This function acts as the primary facade for the simulator's XML loading pipeline.
	 * It performs the following steps:
	 * 1. Resets the target `World` and global parameters (`params::params`).
	 * 2. Loads the main XML file.
	 * 3. Recursively finds and merges any `<include>` files into the main document.
	 * 4. Optionally validates the combined document against the built-in DTD and XSD schemas.
	 * 5. Uses the internal parser utilities to instantiate simulation objects.
	 * 6. Updates the global `params::params` with the parsed context parameters.
	 *
	 * @param filename The filesystem path to the main XML simulation script.
	 * @param world A pointer to the `World` object to be populated with parsed components.
	 * @param validate A boolean indicating whether to perform strict XML schema validation.
	 * @param masterSeeder A reference to the master random number generator used for assigning independent seeds to
	 * components.
	 *
	 * @throws XmlException if the XML is malformed, fails schema validation, or contains invalid scenario logic.
	 * @throws std::runtime_error for file I/O errors or other critical setup issues.
	 */
	void parseSimulation(const std::string& filename, core::World* world, bool validate, std::mt19937& masterSeeder);

	/**
	 * @brief Parses a simulation configuration directly from an XML string in memory.
	 *
	 * Similar to `parseSimulation`, but operates on a raw string instead of a file.
	 * Because it does not load from the filesystem, `<include>` tags are ignored,
	 * and any relative paths for file-backed assets (like waveforms or antennas)
	 * will be resolved against the current working directory (`.`).
	 *
	 * @param xmlContent The raw XML string containing the scenario definition.
	 * @param world A pointer to the `World` object to be populated with parsed components.
	 * @param validate A boolean indicating whether to perform strict XML schema validation.
	 * @param masterSeeder A reference to the master random number generator used for assigning independent seeds to
	 * components.
	 *
	 * @throws XmlException if the XML string is malformed, fails schema validation, or contains invalid scenario logic.
	 */
	void parseSimulationFromString(const std::string& xmlContent, core::World* world, bool validate,
								   std::mt19937& masterSeeder);
}
