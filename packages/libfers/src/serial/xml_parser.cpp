// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file xml_parser.cpp
 * @brief Implementation file for parsing XML configuration files for simulation.
 */

#include "xml_parser.h"

#include <filesystem>

#include "core/logging.h"
#include "core/parameters.h"
#include "core/world.h"
#include "libxml_wrapper.h"
#include "xml_parser_utils.h"

namespace serial
{
	void parseSimulation(const std::string& filename, core::World* world, const bool validate,
						 std::mt19937& masterSeeder)
	{
		world->clear();
		params::params.reset();

		XmlDocument main_doc;
		if (!main_doc.loadFile(filename))
		{
			throw XmlException("Failed to load main XML file: " + filename);
		}

		const std::filesystem::path main_dir = std::filesystem::path(filename).parent_path();
		const bool did_combine = xml_parser_utils::addIncludeFilesToMainDocument(main_doc, main_dir);

		if (validate)
		{
			xml_parser_utils::validateXml(did_combine, main_doc);
		}
		else
		{
			LOG(logging::Level::DEBUG, "Skipping XML validation.");
		}

		xml_parser_utils::ParserContext ctx;
		ctx.world = world;
		ctx.master_seeder = &masterSeeder;
		ctx.base_dir = main_dir;
		ctx.loaders = xml_parser_utils::createDefaultAssetLoaders();

		xml_parser_utils::processParsedDocument(main_doc, ctx);

		// Push the isolated context parameters into global application parameters
		params::params = ctx.parameters;
	}

	void parseSimulationFromString(const std::string& xmlContent, core::World* world, const bool validate,
								   std::mt19937& masterSeeder)
	{
		world->clear();
		params::params.reset();

		XmlDocument doc;
		if (!doc.loadString(xmlContent))
		{
			throw XmlException("Failed to parse XML from memory string.");
		}

		if (validate)
		{
			// Note: <include> tags are not processed when loading from a string.
			xml_parser_utils::validateXml(false, doc);
		}
		else
		{
			LOG(logging::Level::DEBUG, "Skipping XML validation.");
		}

		// When loading from a string, there's no base directory for relative asset paths.
		// The UI/caller is responsible for ensuring any paths in the XML are absolute or resolvable.
		const std::filesystem::path base_dir = ".";

		xml_parser_utils::ParserContext ctx;
		ctx.world = world;
		ctx.master_seeder = &masterSeeder;
		ctx.base_dir = base_dir;
		ctx.loaders = xml_parser_utils::createDefaultAssetLoaders();

		xml_parser_utils::processParsedDocument(doc, ctx);

		// Push the isolated context parameters into global application parameters
		params::params = ctx.parameters;
	}
}
