// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

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
	using ConverterFunc = std::function<void(const math::Vec3&, double&, double&, double&)>;

	struct KmlContext
	{
		params::Parameters parameters;
		ConverterFunc converter;
	};

	double sincAntennaGain(double theta, double alpha, double beta, double gamma);
	double find3DbDropAngle(double alpha, double beta, double gamma);
	double findGaussian3DbDropAngle(const antenna::Gaussian* gaussianAnt);
	double findParabolic3DbDropAngle(const antenna::Parabolic* parabolicAnt, double wavelength);
	double findSquareHorn3DbDropAngle(const antenna::SquareHorn* squarehornAnt, double wavelength);

	std::string formatCoordinates(double lon, double lat, double alt);
	void calculateDestinationCoordinate(double startLatitude, double startLongitude, double angle, double distance,
										double& destLatitude, double& destLongitude);
	std::vector<std::pair<double, double>> generateCircleCoordinates(double lat, double lon, double radius_km);

	void writeKmlHeaderAndStyles(std::ostream& out, const KmlContext& ctx);
	void writePoint(std::ostream& out, const std::string& indent, const std::string& name, const std::string& styleUrl,
					const std::string& coordinates, double objectAltitude, double referenceAltitude);
	void writeAntennaBeamLine(std::ostream& out, const std::string& indent, const std::string& name,
							  const std::string& style, const std::string& startCoords, const std::string& endCoords);

	std::string getPlacemarkStyleForPlatform(const std::vector<const radar::Object*>& objects);
	const radar::Radar* getPrimaryRadar(const std::vector<const radar::Object*>& objects);

	void generateIsotropicAntennaKml(std::ostream& out, const math::Vec3& position, const KmlContext& ctx,
									 const std::string& indent);
	void generateDirectionalAntennaKml(std::ostream& out, const radar::Platform* platform, const KmlContext& ctx,
									   const std::optional<double>& angle3DbDropDeg, const std::string& indent);
	void generateAntennaKml(std::ostream& out, const radar::Platform* platform, const radar::Radar* radar,
							const KmlContext& ctx, const std::string& indent);

	void generateDynamicPathKml(std::ostream& out, const radar::Platform* platform, const std::string& styleUrl,
								double refAlt, const KmlContext& ctx, const std::string& indent);
	void generateTrackEndpointsKml(std::ostream& out, const radar::Platform* platform, double refAlt,
								   const KmlContext& ctx, const std::string& indent);
	void generateStaticPlacemarkKml(std::ostream& out, const radar::Platform* platform, const std::string& styleUrl,
									double refAlt, const KmlContext& ctx, const std::string& indent);
	void generatePlatformPathKml(std::ostream& out, const radar::Platform* platform, const std::string& style,
								 double refAlt, const KmlContext& ctx, const std::string& indent);

	void processPlatform(std::ostream& out, const radar::Platform* platform,
						 const std::vector<const radar::Object*>& objects, const KmlContext& ctx,
						 double referenceAltitude, const std::string& indent);

	void generateKmlToStream(std::ostream& out, const core::World& world, const KmlContext& ctx);
}
