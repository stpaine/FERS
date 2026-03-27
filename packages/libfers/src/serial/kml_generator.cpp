// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2023-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file kml_generator.cpp
 * @brief Source file for KML file generation from FERS simulation scenarios.
 */

#include "serial/kml_generator.h"

#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/LocalCartesian.hpp>
#include <GeographicLib/UTMUPS.hpp>
#include <fstream>
#include <memory>

#include "core/logging.h"
#include "core/parameters.h"
#include "serial/kml_generator_utils.h"

namespace serial
{
	bool KmlGenerator::generateKml(const core::World& world, const std::string& outputKmlPath)
	{
		try
		{
			kml_generator_utils::KmlContext ctx;
			ctx.parameters = params::params;

			switch (ctx.parameters.coordinate_frame)
			{
			case params::CoordinateFrame::ENU:
				{
					auto proj = std::make_shared<GeographicLib::LocalCartesian>(ctx.parameters.origin_latitude,
																				ctx.parameters.origin_longitude,
																				ctx.parameters.origin_altitude);
					ctx.converter = [proj](const math::Vec3& pos, double& lat, double& lon, double& alt)
					{ proj->Reverse(pos.x, pos.y, pos.z, lat, lon, alt); };
					break;
				}
			case params::CoordinateFrame::UTM:
				{
					const int zone = ctx.parameters.utm_zone;
					const bool northp = ctx.parameters.utm_north_hemisphere;
					ctx.converter = [zone, northp](const math::Vec3& pos, double& lat, double& lon, double& alt)
					{
						double gamma, k;
						GeographicLib::UTMUPS::Reverse(zone, northp, pos.x, pos.y, lat, lon, gamma, k);
						alt = pos.z;
					};
					break;
				}
			case params::CoordinateFrame::ECEF:
				{
					const auto& earth = GeographicLib::Geocentric::WGS84();
					ctx.converter = [&earth](const math::Vec3& pos, double& lat, double& lon, double& alt)
					{ earth.Reverse(pos.x, pos.y, pos.z, lat, lon, alt); };
					break;
				}
			}

			std::ofstream kml_file(outputKmlPath.c_str());
			if (!kml_file.is_open())
			{
				LOG(logging::Level::ERROR, "Error opening output KML file {}", outputKmlPath);
				return false;
			}

			kml_generator_utils::generateKmlToStream(kml_file, world, ctx);
			kml_file.close();
		}
		catch (const std::exception& e)
		{
			LOG(logging::Level::ERROR, "Error generating KML file: {}", e.what());
			return false;
		}
		catch (...)
		{
			LOG(logging::Level::ERROR, "Unknown error occurred while generating KML file.");
			return false;
		}
		return true;
	}
}
