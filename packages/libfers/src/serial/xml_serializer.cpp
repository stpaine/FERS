// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file xml_serializer.cpp
 * @brief Implementation for serializing the simulation world to XML.
 *
 * This file acts as a facade, reading the global simulation parameters
 * and delegating the object-by-object XML reconstruction to the utils layer.
 */

#include "xml_serializer.h"

#include <ranges>

#include "core/parameters.h"
#include "core/world.h"
#include "serial/libxml_wrapper.h"
#include "xml_serializer_utils.h"

namespace serial
{
	std::string world_to_xml_string(const core::World& world)
	{
		XmlDocument doc;
		xmlNodePtr sim_node = xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("simulation"));
		XmlElement root(sim_node);
		doc.setRootElement(root);

		const auto& p = params::params;

		if (!p.simulation_name.empty())
		{
			root.setAttribute("name", p.simulation_name);
		}
		else
		{
			root.setAttribute("name", "FERS Scenario");
		}

		const XmlElement params_elem = root.addChild("parameters");
		xml_serializer_utils::serializeParameters(params_elem, p);

		for (const auto& waveform : world.getWaveforms() | std::views::values)
		{
			XmlElement waveform_elem = root.addChild("waveform");
			xml_serializer_utils::serializeWaveform(*waveform, waveform_elem);
		}
		for (const auto& timing : world.getTimings() | std::views::values)
		{
			XmlElement timing_elem = root.addChild("timing");
			xml_serializer_utils::serializeTiming(*timing, timing_elem);
		}
		for (const auto& antenna : world.getAntennas() | std::views::values)
		{
			XmlElement antenna_elem = root.addChild("antenna");
			xml_serializer_utils::serializeAntenna(*antenna, antenna_elem);
		}
		for (const auto& platform : world.getPlatforms())
		{
			XmlElement plat_elem = root.addChild("platform");
			xml_serializer_utils::serializePlatform(*platform, world, plat_elem);
		}

		return doc.dumpToString();
	}
}
