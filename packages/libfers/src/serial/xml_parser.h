// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file xml_parser.h
 * @brief Header file for parsing XML configuration files for simulation.
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
	 * This function loads an XML file, merges any included files, validates it against
	 * the FERS DTD and XSD schemas, and then populates the simulation `World` object
	 * with the parsed parameters and components.
	 *
	 * @param filename The path to the main XML simulation script.
	 * @param world A pointer to the `World` object to be populated.
	 * @param validate A boolean indicating whether to perform XML validation.
	 * @param masterSeeder A reference to the master random number generator used for seeding components.
	 * @throws XmlException if the XML is malformed, fails validation, or contains invalid data.
	 * @throws std::runtime_error for file I/O errors or other critical issues.
	 */
	void parseSimulation(const std::string& filename, core::World* world, bool validate, std::mt19937& masterSeeder);

	void parseSimulationFromString(const std::string& xmlContent, core::World* world, bool validate,
								   std::mt19937& masterSeeder);
}
