// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file kml_generator.h
 * @brief KML file generator for geographical visualization of FERS scenarios.
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
	 * @class KmlGenerator
	 * @brief Generates KML files from FERS simulation scenarios for geographical visualization.
	 *
	 * This class generates KML files for geographical visualization of FERS scenarios. It
	 * interprets the simulation coordinates based on the user-specified coordinate system
	 * in the XML file, which can be one of:
	 *
	 * - ENU (East-North-Up): Default. Local Cartesian coordinates (x, y, z) are
	 *   treated as meters in an ENU tangent plane centered at a geodetic `<origin>`.
	 *
	 * - UTM (Universal Transverse Mercator): Coordinates (x, y, z) are treated as
	 *   easting (m), northing (m), and altitude (m) within a specified UTM zone and
	 *   hemisphere.
	 *
	 * - ECEF (Earth-Centered, Earth-Fixed): Coordinates (x, y, z) are treated as
	 *   geocentric X, Y, Z values in meters.
	 *
	 * All input coordinates are converted to WGS84 geodetic coordinates (latitude,
	 * longitude, altitude) for the final KML output. The KML is written with
	 * `<altitudeMode>absolute</altitudeMode>`, where altitude is relative to Mean Sea
	 * Level (MSL).
	 */
	class KmlGenerator
	{
	public:
		/**
		 * @brief Generates a KML file from a pre-built simulation world.
		 *
		 * @param world The simulation world containing all objects and paths.
		 * @param outputKmlPath The path for the output KML file.
		 * @return True on success, false on failure.
		 */
		static bool generateKml(const core::World& world, const std::string& outputKmlPath);
	};
}
