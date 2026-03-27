// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file kml_generator_utils.h
 * @brief Utility definitions and functions for generating KML files from simulation scenarios.
 *
 * This file provides the internal mechanisms required to export FERS scenarios
 * to KML format for visualization in Google Earth or similar tools. It defines
 * functions for formatting coordinates, rendering styling, creating models,
 * and assembling the top-level document structure.
 */

#pragma once

#include <functional>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include "core/parameters.h"
#include "math/geometry_ops.h"

namespace core
{
	class World;
}

namespace radar
{
	class Platform;
	class Radar;
	class Object;
}

namespace antenna
{
	class Gaussian;
	class Parabolic;
	class SquareHorn;
}

namespace serial::kml_generator_utils
{
	/**
	 * @brief Callback signature used for converting math coordinates to lat/lon/altitude.
	 */
	using ConverterFunc = std::function<void(const math::Vec3&, double&, double&, double&)>;

	/**
	 * @struct KmlContext
	 * @brief Context data required during KML generation.
	 */
	struct KmlContext
	{
		params::Parameters parameters; ///< A copy of the global simulation parameters.
		ConverterFunc converter; ///< Function used to translate simulation Cartesian space into geographic coords.
	};

	/**
	 * @brief Calculates a normalized sinc-based antenna gain mapping.
	 * @param theta The angle to evaluate.
	 * @param alpha Primary dimensional scale factor.
	 * @param beta Secondary dimensional scale factor.
	 * @param gamma Falloff rate.
	 * @return The calculated gain.
	 */
	double sincAntennaGain(double theta, double alpha, double beta, double gamma);

	/**
	 * @brief Numerically determines the 3dB drop angle for a parameterized generic antenna.
	 * @param alpha Primary dimensional scale factor.
	 * @param beta Secondary dimensional scale factor.
	 * @param gamma Falloff rate.
	 * @return The calculated 3dB drop angle in radians.
	 */
	double find3DbDropAngle(double alpha, double beta, double gamma);

	/**
	 * @brief Calculates the 3dB drop angle for a Gaussian antenna.
	 * @param gaussianAnt Pointer to the Gaussian antenna.
	 * @return The 3dB drop angle in degrees.
	 */
	double findGaussian3DbDropAngle(const antenna::Gaussian* gaussianAnt);

	/**
	 * @brief Calculates the 3dB drop angle for a Parabolic antenna.
	 * @param parabolicAnt Pointer to the Parabolic antenna.
	 * @param wavelength The operating wavelength.
	 * @return The 3dB drop angle in degrees.
	 */
	double findParabolic3DbDropAngle(const antenna::Parabolic* parabolicAnt, double wavelength);

	/**
	 * @brief Calculates the 3dB drop angle for a Square Horn antenna.
	 * @param squarehornAnt Pointer to the Square Horn antenna.
	 * @param wavelength The operating wavelength.
	 * @return The 3dB drop angle in degrees.
	 */
	double findSquareHorn3DbDropAngle(const antenna::SquareHorn* squarehornAnt, double wavelength);

	/**
	 * @brief Formats coordinates into a comma-separated string suitable for KML `<coordinates>`.
	 * @param lon Longitude in degrees.
	 * @param lat Latitude in degrees.
	 * @param alt Altitude in meters.
	 * @return Formatted coordinate string containing "lon,lat,alt".
	 */
	std::string formatCoordinates(double lon, double lat, double alt);

	/**
	 * @brief Calculates a destination coordinate given a starting position, bearing, and distance.
	 * @param startLatitude Starting latitude in degrees.
	 * @param startLongitude Starting longitude in degrees.
	 * @param angle Bearing angle in degrees.
	 * @param distance Distance to travel in meters.
	 * @param destLatitude Output destination latitude in degrees.
	 * @param destLongitude Output destination longitude in degrees.
	 */
	void calculateDestinationCoordinate(double startLatitude, double startLongitude, double angle, double distance,
										double& destLatitude, double& destLongitude);

	/**
	 * @brief Generates a collection of points tracing a circle around a center coordinate.
	 * @param lat Center latitude in degrees.
	 * @param lon Center longitude in degrees.
	 * @param radius_km Circle radius in kilometers.
	 * @return A vector of latitude-longitude pairs defining the circle.
	 */
	std::vector<std::pair<double, double>> generateCircleCoordinates(double lat, double lon, double radius_km);

	/**
	 * @brief Writes the standard KML preamble and style definitions to the output stream.
	 * @param out The stream to write to.
	 * @param ctx The current KML generation context.
	 */
	void writeKmlHeaderAndStyles(std::ostream& out, const KmlContext& ctx);

	/**
	 * @brief Writes a KML `<Point>` placemark to the output stream.
	 * @param out The stream to write to.
	 * @param indent The indentation string.
	 * @param name The name of the point.
	 * @param styleUrl The ID of the style to apply.
	 * @param coordinates The pre-formatted coordinate string.
	 * @param objectAltitude The absolute altitude of the point.
	 * @param referenceAltitude The reference altitude (e.g. ground level) for relative extrusions.
	 */
	void writePoint(std::ostream& out, const std::string& indent, const std::string& name, const std::string& styleUrl,
					const std::string& coordinates, double objectAltitude, double referenceAltitude);

	/**
	 * @brief Writes a visual cone/beam line representing the antenna look direction.
	 * @param out The stream to write to.
	 * @param indent The indentation string.
	 * @param name The name of the beam projection line.
	 * @param style The ID of the line style to apply.
	 * @param startCoords The pre-formatted origin coordinate string.
	 * @param endCoords The pre-formatted destination coordinate string.
	 */
	void writeAntennaBeamLine(std::ostream& out, const std::string& indent, const std::string& name,
							  const std::string& style, const std::string& startCoords, const std::string& endCoords);

	/**
	 * @brief Determines the proper KML style definition ID to use based on the platform's attached objects.
	 * @param objects A collection of objects attached to the platform.
	 * @return The style URL string (e.g., "#RadarPlatformStyle").
	 */
	std::string getPlacemarkStyleForPlatform(const std::vector<const radar::Object*>& objects);

	/**
	 * @brief Identifies the primary radar object within a platform for styling and direction tracking.
	 * @param objects A collection of objects attached to the platform.
	 * @return Pointer to the primary radar, or nullptr if none exists.
	 */
	const radar::Radar* getPrimaryRadar(const std::vector<const radar::Object*>& objects);

	/**
	 * @brief Renders the visual representation for an isotropic antenna into the output stream.
	 * @param out The stream to write to.
	 * @param position The position of the antenna in simulation space.
	 * @param ctx The current KML context.
	 * @param indent The indentation string.
	 */
	void generateIsotropicAntennaKml(std::ostream& out, const math::Vec3& position, const KmlContext& ctx,
									 const std::string& indent);

	/**
	 * @brief Renders the visual pointing representation for a directional (beam) antenna.
	 * @param out The stream to write to.
	 * @param platform The platform the antenna is mounted on.
	 * @param ctx The current KML context.
	 * @param angle3DbDropDeg An optional computed angular drop-off to dictate the beam width rendering.
	 * @param indent The indentation string.
	 */
	void generateDirectionalAntennaKml(std::ostream& out, const radar::Platform* platform, const KmlContext& ctx,
									   const std::optional<double>& angle3DbDropDeg, const std::string& indent);

	/**
	 * @brief Dispatch function that selects and generates the appropriate KML for a given radar's antenna.
	 * @param out The stream to write to.
	 * @param platform The platform containing the radar.
	 * @param radar The radar to visualize.
	 * @param ctx The current KML context.
	 * @param indent The indentation string.
	 */
	void generateAntennaKml(std::ostream& out, const radar::Platform* platform, const radar::Radar* radar,
							const KmlContext& ctx, const std::string& indent);

	/**
	 * @brief Generates KML for a continuously moving dynamic platform path.
	 * @param out The stream to write to.
	 * @param platform The platform undergoing dynamic motion.
	 * @param styleUrl The style ID for the path line.
	 * @param refAlt The baseline reference altitude.
	 * @param ctx The current KML context.
	 * @param indent The indentation string.
	 */
	void generateDynamicPathKml(std::ostream& out, const radar::Platform* platform, const std::string& styleUrl,
								double refAlt, const KmlContext& ctx, const std::string& indent);

	/**
	 * @brief Generates KML rendering start and end pushpins for a moving platform's track.
	 * @param out The stream to write to.
	 * @param platform The platform containing the path boundaries.
	 * @param refAlt The baseline reference altitude.
	 * @param ctx The current KML context.
	 * @param indent The indentation string.
	 */
	void generateTrackEndpointsKml(std::ostream& out, const radar::Platform* platform, double refAlt,
								   const KmlContext& ctx, const std::string& indent);

	/**
	 * @brief Generates a simple static placemark KML for a non-moving platform.
	 * @param out The stream to write to.
	 * @param platform The stationary platform.
	 * @param styleUrl The style ID to apply to the placemark.
	 * @param refAlt The baseline reference altitude.
	 * @param ctx The current KML context.
	 * @param indent The indentation string.
	 */
	void generateStaticPlacemarkKml(std::ostream& out, const radar::Platform* platform, const std::string& styleUrl,
									double refAlt, const KmlContext& ctx, const std::string& indent);

	/**
	 * @brief Dispatches the generation of a platform's path representation (static vs dynamic).
	 * @param out The stream to write to.
	 * @param platform The platform whose path needs to be rendered.
	 * @param style The path tracing styling.
	 * @param refAlt The baseline reference altitude.
	 * @param ctx The current KML context.
	 * @param indent The indentation string.
	 */
	void generatePlatformPathKml(std::ostream& out, const radar::Platform* platform, const std::string& style,
								 double refAlt, const KmlContext& ctx, const std::string& indent);

	/**
	 * @brief Orchestrates full processing and rendering of an individual platform into the KML stream.
	 * @param out The stream to write to.
	 * @param platform The platform to process.
	 * @param objects Attached radar or target objects corresponding to the platform.
	 * @param ctx The current KML context.
	 * @param referenceAltitude Baseline reference altitude for the entire platform.
	 * @param indent The indentation string.
	 */
	void processPlatform(std::ostream& out, const radar::Platform* platform,
						 const std::vector<const radar::Object*>& objects, const KmlContext& ctx,
						 double referenceAltitude, const std::string& indent);

	/**
	 * @brief Master entry point designed to convert the comprehensive simulation world state into a valid KML document.
	 * @param out The output stream where the resulting KML markup goes.
	 * @param world The aggregated global simulation state.
	 * @param ctx The configured KML context containing simulation parameters and a coordinate converter.
	 */
	void generateKmlToStream(std::ostream& out, const core::World& world, const KmlContext& ctx);
}
