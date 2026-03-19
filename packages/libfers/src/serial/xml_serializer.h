// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file xml_serializer.h
 * @brief Provides functions to serialize the simulation world back into the FERS XML format.
 */

#pragma once

#include <string>

namespace core
{
	class World;
}

namespace serial
{
	/**
	 * @brief Serializes the entire simulation world into an XML formatted string.
	 *
	 * This function serves as the reverse of the XML parser. It is essential for allowing
	 * users to modify a scenario in a UI and then export their changes back into a valid
	 * FERS XML file that can be used by the CLI or shared. It iterates through the
	 * in-memory `core::World` object and reconstructs the corresponding XML structure.
	 *
	 * @param world The world object to serialize.
	 * @return A string containing the XML representation of the world.
	 */
	std::string world_to_xml_string(const core::World& world);
}
