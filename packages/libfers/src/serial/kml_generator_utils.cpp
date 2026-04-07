// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "serial/kml_generator_utils.h"

#include <GeographicLib/Geodesic.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <map>
#include <ranges>
#include <sstream>

#include "antenna/antenna_factory.h"
#include "core/logging.h"
#include "core/world.h"
#include "math/coord.h"
#include "math/path.h"
#include "radar/platform.h"
#include "radar/radar_obj.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/rotation_angle_utils.h"
#include "signal/radar_signal.h"

namespace serial::kml_generator_utils
{
	constexpr int TRACK_NUM_DIVISIONS = 100;
	constexpr int ISOTROPIC_PATTERN_POINTS = 100;
	constexpr double ISOTROPIC_PATTERN_RADIUS_KM = 20.0;
	constexpr double DIRECTIONAL_ANTENNA_ARROW_LENGTH_M = 20000.0;

	double sincAntennaGain(const double theta, const double alpha, const double beta, const double gamma)
	{
		if (theta == 0.0)
		{
			return alpha;
		}
		const double gain = alpha * std::pow(std::sin(beta * theta) / (beta * theta), gamma);
		return gain;
	}

	double find3DbDropAngle(const double alpha, const double beta, const double gamma)
	{
		constexpr std::size_t num_points = 1000;
		const auto midpoint = static_cast<std::ptrdiff_t>(num_points / 2);
		std::vector<double> theta(num_points);
		std::vector<double> gain(num_points);
		for (std::size_t i = 0; i < num_points; ++i)
		{
			theta[i] = -PI + 2.0 * PI * static_cast<double>(i) / static_cast<double>(num_points - 1);
			gain[i] = sincAntennaGain(theta[i], alpha, beta, gamma);
		}
		const auto search_begin = gain.begin() + midpoint;
		const double max_gain = *std::max_element(search_begin, gain.end());
		const double max_gain_db = 10.0 * std::log10(max_gain);
		const double target_gain_db = max_gain_db - 3.0;
		const double target_gain = std::pow(10.0, target_gain_db / 10.0);
		const auto min_gain = std::min_element(search_begin, gain.end(), [target_gain](const double a, const double b)
											   { return std::abs(a - target_gain) < std::abs(b - target_gain); });
		const auto idx = static_cast<std::size_t>(std::distance(search_begin, min_gain));
		const double angle_3db_drop = theta[static_cast<std::size_t>(midpoint) + idx];
		return angle_3db_drop * 180.0 / PI;
	}

	double findGaussian3DbDropAngle(const antenna::Gaussian* gaussianAnt)
	{
		if (gaussianAnt->getAzimuthScale() <= 0.0)
		{
			LOG(logging::Level::WARNING,
				"Gaussian antenna '{}' has a non-positive azimuth scale ({}). 3dB beamwidth is undefined. KML will "
				"only show boresight.",
				gaussianAnt->getName(), gaussianAnt->getAzimuthScale());
			return 0.0;
		}
		const double half_angle_rad = std::sqrt(std::log(2.0) / gaussianAnt->getAzimuthScale());
		return half_angle_rad * 180.0 / PI;
	}

	double findParabolic3DbDropAngle(const antenna::Parabolic* parabolicAnt, const double wavelength)
	{
		if (parabolicAnt->getDiameter() <= 0.0)
		{
			LOG(logging::Level::WARNING,
				"Parabolic antenna '{}' has a non-positive diameter ({}). This is physically impossible. KML will only "
				"show boresight.",
				parabolicAnt->getName(), parabolicAnt->getDiameter());
			return 0.0;
		}
		const double arg = 1.6 * wavelength / (PI * parabolicAnt->getDiameter());
		if (arg > 1.0)
		{
			LOG(logging::Level::INFO,
				"Parabolic antenna '{}': The operating wavelength ({:.4f}m) is very large compared to its diameter "
				"({:.4f}m), resulting in a nearly omnidirectional pattern. KML visualization will cap the 3dB "
				"half-angle at 90 degrees.",
				parabolicAnt->getName(), wavelength, parabolicAnt->getDiameter());
			return 90.0;
		}
		const double half_angle_rad = std::asin(arg);
		return half_angle_rad * 180.0 / PI;
	}

	double findSquareHorn3DbDropAngle(const antenna::SquareHorn* squarehornAnt, const double wavelength)
	{
		if (squarehornAnt->getDimension() <= 0.0)
		{
			LOG(logging::Level::WARNING,
				"SquareHorn antenna '{}' has a non-positive dimension ({}). This is physically impossible. KML will "
				"only show boresight.",
				squarehornAnt->getName(), squarehornAnt->getDimension());
			return 0.0;
		}
		const double arg = 1.39155 * wavelength / (PI * squarehornAnt->getDimension());
		if (arg > 1.0)
		{
			LOG(logging::Level::INFO,
				"SquareHorn antenna '{}': The operating wavelength ({:.4f}m) is very large compared to its dimension "
				"({:.4f}m), resulting in a nearly omnidirectional pattern. KML visualization will cap the 3dB "
				"half-angle at 90 degrees.",
				squarehornAnt->getName(), wavelength, squarehornAnt->getDimension());
			return 90.0;
		}
		const double half_angle_rad = std::asin(arg);
		return half_angle_rad * 180.0 / PI;
	}

	std::string formatCoordinates(const double lon, const double lat, const double alt)
	{
		std::stringstream ss;
		ss << std::fixed << std::setprecision(6) << lon << "," << lat << "," << alt;
		return ss.str();
	}

	void calculateDestinationCoordinate(const double startLatitude, const double startLongitude, const double angle,
										const double distance, double& destLatitude, double& destLongitude)
	{
		const GeographicLib::Geodesic& geod = GeographicLib::Geodesic::WGS84();
		geod.Direct(startLatitude, startLongitude, angle, distance, destLatitude, destLongitude);
	}

	std::vector<std::pair<double, double>> generateCircleCoordinates(const double lat, const double lon,
																	 const double radius_km)
	{
		std::vector<std::pair<double, double>> circle_coordinates;
		for (int i = 0; i < ISOTROPIC_PATTERN_POINTS; i++)
		{
			const double bearing = i * 360.0 / ISOTROPIC_PATTERN_POINTS;
			double new_lat, new_lon;
			calculateDestinationCoordinate(lat, lon, bearing, radius_km * 1000.0, new_lat, new_lon);
			circle_coordinates.emplace_back(new_lat, new_lon);
		}
		return circle_coordinates;
	}

	void writeKmlHeaderAndStyles(std::ostream& out, const KmlContext& ctx)
	{
		out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
		out << "<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\">\n";
		out << "<Document>\n";
		out << "  <name>";
		if (ctx.parameters.simulation_name.empty())
		{
			out << "FERS Simulation Visualization";
		}
		else
		{
			out << ctx.parameters.simulation_name;
		}
		out << "</name>\n";
		out << "  <Style "
			   "id=\"receiver\"><IconStyle><Icon><href>https://cdn-icons-png.flaticon.com/512/645/645436.png</href></"
			   "Icon></IconStyle></Style>\n";
		out << "  <Style "
			   "id=\"transmitter\"><IconStyle><Icon><href>https://cdn-icons-png.flaticon.com/128/224/224666.png</"
			   "href></Icon></IconStyle></Style>\n";
		out << "  <Style "
			   "id=\"target\"><IconStyle><Icon><href>https://upload.wikimedia.org/wikipedia/commons/thumb/a/ad/"
			   "Target_red_dot1.svg/1200px-Target_red_dot1.svg.png</href></Icon></IconStyle><LineStyle><width>2</"
			   "width></LineStyle></Style>\n";
		out << "  <Style "
			   "id=\"translucentPolygon\"><LineStyle><color>ff0000ff</color><width>2</width></"
			   "LineStyle><PolyStyle><color>00ffffff</color></PolyStyle></Style>\n";
		out << "  <Style "
			   "id=\"arrowStyle\"><IconStyle><Icon><href>http://maps.google.com/mapfiles/kml/shapes/arrow.png</href></"
			   "Icon><scale>0.5</scale></IconStyle></Style>\n";
		out << "  <Style id=\"lineStyle\"><LineStyle><color>ff0000ff</color><width>2</width></LineStyle></Style>\n";
		out << "  <Style id=\"lineStyleBlue\"><LineStyle><color>ffff0000</color><width>2</width></LineStyle></Style>\n";
	}

	void writePoint(std::ostream& out, const std::string& indent, const std::string& name, const std::string& styleUrl,
					const std::string& coordinates, const double objectAltitude, const double referenceAltitude)
	{
		out << indent << "<Placemark>\n";
		out << indent << "  <name>" << name << "</name>\n";
		out << indent << "  <styleUrl>" << styleUrl << "</styleUrl>\n";
		out << indent << "  <Point>\n";
		out << indent << "    <coordinates>" << coordinates << "</coordinates>\n";
		out << indent << "    <altitudeMode>absolute</altitudeMode>\n";
		if (objectAltitude > referenceAltitude)
		{
			out << indent << "    <extrude>1</extrude>\n";
		}
		out << indent << "  </Point>\n";
		out << indent << "</Placemark>\n";
	}

	void writeAntennaBeamLine(std::ostream& out, const std::string& indent, const std::string& name,
							  const std::string& style, const std::string& startCoords, const std::string& endCoords)
	{
		out << indent << "<Placemark>\n";
		out << indent << "  <name>" << name << "</name>\n";
		out << indent << "  <styleUrl>" << style << "</styleUrl>\n";
		out << indent << "  <LineString>\n";
		out << indent << "    <altitudeMode>absolute</altitudeMode>\n";
		out << indent << "    <tessellate>1</tessellate>\n";
		out << indent << "    <coordinates>" << startCoords << " " << endCoords << "</coordinates>\n";
		out << indent << "  </LineString>\n";
		out << indent << "</Placemark>\n";
	}

	std::string getPlacemarkStyleForPlatform(const std::vector<const radar::Object*>& objects)
	{
		bool has_receiver = false;
		bool has_transmitter = false;
		for (const auto* obj : objects)
		{
			if (dynamic_cast<const radar::Receiver*>(obj) != nullptr)
			{
				has_receiver = true;
			}
			if (dynamic_cast<const radar::Transmitter*>(obj) != nullptr)
			{
				has_transmitter = true;
			}
		}

		if (has_receiver)
		{
			return "#receiver";
		}
		if (has_transmitter)
		{
			return "#transmitter";
		}
		return "#target";
	}

	const radar::Radar* getPrimaryRadar(const std::vector<const radar::Object*>& objects)
	{
		for (const auto* obj : objects)
		{
			if (const auto* const r = dynamic_cast<const radar::Radar*>(obj))
			{
				return r;
			}
		}
		return nullptr;
	}

	void generateIsotropicAntennaKml(std::ostream& out, const math::Vec3& position, const KmlContext& ctx,
									 const std::string& indent)
	{
		double lat, lon, alt_abs;
		ctx.converter(position, lat, lon, alt_abs);

		const auto circle_coordinates = generateCircleCoordinates(lat, lon, ISOTROPIC_PATTERN_RADIUS_KM);

		out << indent << "<Placemark>\n";
		out << indent << "  <name>Isotropic pattern range</name>\n";
		out << indent << "  <styleUrl>#translucentPolygon</styleUrl>\n";
		out << indent << "  <Polygon>\n";
		out << indent << "    <extrude>1</extrude>\n";
		out << indent << "    <altitudeMode>absolute</altitudeMode>\n";
		out << indent << "    <outerBoundaryIs><LinearRing><coordinates>\n";
		for (const auto& [pt_lat, pt_lon] : circle_coordinates)
		{
			out << indent << "      " << formatCoordinates(pt_lon, pt_lat, alt_abs) << "\n";
		}
		out << indent << "      "
			<< formatCoordinates(circle_coordinates[0].second, circle_coordinates[0].first, alt_abs) << "\n";
		out << indent << "    </coordinates></LinearRing></outerBoundaryIs>\n";
		out << indent << "  </Polygon>\n";
		out << indent << "</Placemark>\n";
	}

	void generateDirectionalAntennaKml(std::ostream& out, const radar::Platform* platform, const KmlContext& ctx,
									   const std::optional<double>& angle3DbDropDeg, const std::string& indent)
	{
		const auto& first_wp_pos = platform->getMotionPath()->getCoords().front().pos;
		double start_lat, start_lon, start_alt;
		ctx.converter(first_wp_pos, start_lat, start_lon, start_alt);
		const std::string start_coords_str = formatCoordinates(start_lon, start_lat, start_alt);

		const math::SVec3 initial_rotation = platform->getRotationPath()->getPosition(ctx.parameters.start);
		const double display_azimuth = rotation_angle_utils::internal_azimuth_to_external(
			initial_rotation.azimuth, ctx.parameters.rotation_angle_unit);
		const double display_elevation = rotation_angle_utils::internal_elevation_to_external(
			initial_rotation.elevation, ctx.parameters.rotation_angle_unit);

		const double fers_azimuth_deg = initial_rotation.azimuth * 180.0 / PI;
		double start_azimuth_deg_kml = 90.0 - fers_azimuth_deg;
		start_azimuth_deg_kml = std::fmod(start_azimuth_deg_kml, 360.0);
		if (start_azimuth_deg_kml < 0.0)
		{
			start_azimuth_deg_kml += 360.0;
		}

		const double horizontal_distance = DIRECTIONAL_ANTENNA_ARROW_LENGTH_M * std::cos(initial_rotation.elevation);
		const double delta_altitude = DIRECTIONAL_ANTENNA_ARROW_LENGTH_M * std::sin(initial_rotation.elevation);
		const double end_alt = start_alt + delta_altitude;

		double dest_lat, dest_lon;
		calculateDestinationCoordinate(start_lat, start_lon, start_azimuth_deg_kml, horizontal_distance, dest_lat,
									   dest_lon);
		const std::string end_coords_str = formatCoordinates(dest_lon, dest_lat, end_alt);
		out << indent << "<Placemark>\n";
		out << indent << "  <name>Antenna Boresight</name>\n";
		out << indent << "  <ExtendedData>\n";
		out << indent << "    <Data name=\"rotationangleunit\"><value>"
			<< params::rotationAngleUnitToken(ctx.parameters.rotation_angle_unit) << "</value></Data>\n";
		out << indent << "    <Data name=\"azimuth\"><value>" << display_azimuth << "</value></Data>\n";
		out << indent << "    <Data name=\"elevation\"><value>" << display_elevation << "</value></Data>\n";
		out << indent << "  </ExtendedData>\n";
		out << indent << "  <styleUrl>#lineStyle</styleUrl>\n";
		out << indent << "  <LineString>\n";
		out << indent << "    <altitudeMode>absolute</altitudeMode>\n";
		out << indent << "    <tessellate>1</tessellate>\n";
		out << indent << "    <coordinates>" << start_coords_str << " " << end_coords_str << "</coordinates>\n";
		out << indent << "  </LineString>\n";
		out << indent << "</Placemark>\n";

		if (angle3DbDropDeg.has_value() && *angle3DbDropDeg > EPSILON)
		{
			double side1_lat, side1_lon;
			calculateDestinationCoordinate(start_lat, start_lon, start_azimuth_deg_kml - *angle3DbDropDeg,
										   horizontal_distance, side1_lat, side1_lon);
			const std::string side1_coords_str = formatCoordinates(side1_lon, side1_lat, end_alt);
			writeAntennaBeamLine(out, indent, "Antenna 3dB Beamwidth", "#lineStyleBlue", start_coords_str,
								 side1_coords_str);

			double side2_lat, side2_lon;
			calculateDestinationCoordinate(start_lat, start_lon, start_azimuth_deg_kml + *angle3DbDropDeg,
										   horizontal_distance, side2_lat, side2_lon);
			const std::string side2_coords_str = formatCoordinates(side2_lon, side2_lat, end_alt);
			writeAntennaBeamLine(out, indent, "Antenna 3dB Beamwidth", "#lineStyleBlue", start_coords_str,
								 side2_coords_str);
		}

		const double arrow_heading = std::fmod(start_azimuth_deg_kml + 180.0, 360.0);
		out << indent << "<Placemark>\n";
		out << indent << "  <name>Antenna Arrow</name>\n";
		out << indent << "  <styleUrl>#arrowStyle</styleUrl>\n";
		out << indent << "  <Point><coordinates>" << end_coords_str
			<< "</coordinates><altitudeMode>absolute</altitudeMode></Point>\n";
		out << indent << "  <Style>\n";
		out << indent << "    <IconStyle><heading>" << arrow_heading << "</heading></IconStyle>\n";
		out << indent << "  </Style>\n";
		out << indent << "</Placemark>\n";
	}

	void generateAntennaKml(std::ostream& out, const radar::Platform* platform, const radar::Radar* radar,
							const KmlContext& ctx, const std::string& indent)
	{
		const antenna::Antenna* ant = radar->getAntenna();
		if ((ant == nullptr) || platform->getMotionPath()->getCoords().empty())
		{
			return;
		}

		if (dynamic_cast<const antenna::Isotropic*>(ant) != nullptr)
		{
			const math::Vec3 initial_pos = platform->getMotionPath()->getCoords().front().pos;
			generateIsotropicAntennaKml(out, initial_pos, ctx, indent);
		}
		else
		{
			std::optional<double> angle_3db_drop_deg;

			std::optional<double> wavelength;
			if (const auto* tx = dynamic_cast<const radar::Transmitter*>(radar))
			{
				if (tx->getSignal() != nullptr)
				{
					wavelength = ctx.parameters.c / tx->getSignal()->getCarrier();
				}
			}
			else if (const auto* rx = dynamic_cast<const radar::Receiver*>(radar))
			{
				if (const auto* attached_tx = dynamic_cast<const radar::Transmitter*>(rx->getAttached()))
				{
					if (attached_tx->getSignal() != nullptr)
					{
						wavelength = ctx.parameters.c / attached_tx->getSignal()->getCarrier();
					}
				}
			}

			if (const auto* sinc_ant = dynamic_cast<const antenna::Sinc*>(ant))
			{
				angle_3db_drop_deg = find3DbDropAngle(sinc_ant->getAlpha(), sinc_ant->getBeta(), sinc_ant->getGamma());
			}
			else if (const auto* gaussian_ant = dynamic_cast<const antenna::Gaussian*>(ant))
			{
				angle_3db_drop_deg = findGaussian3DbDropAngle(gaussian_ant);
			}
			else if (const auto* parabolic_ant = dynamic_cast<const antenna::Parabolic*>(ant))
			{
				if (wavelength)
				{
					angle_3db_drop_deg = findParabolic3DbDropAngle(parabolic_ant, *wavelength);
				}
			}
			else if (const auto* squarehorn_ant = dynamic_cast<const antenna::SquareHorn*>(ant))
			{
				if (wavelength)
				{
					angle_3db_drop_deg = findSquareHorn3DbDropAngle(squarehorn_ant, *wavelength);
				}
			}
			else if ((dynamic_cast<const antenna::XmlAntenna*>(ant) != nullptr) ||
					 (dynamic_cast<const antenna::H5Antenna*>(ant) != nullptr))
			{
				LOG(logging::Level::INFO,
					"KML visualization for antenna '{}' ('{}') is symbolic. "
					"Only the boresight direction is shown, as a 3dB beamwidth is not calculated from file-based "
					"patterns.",
					ant->getName(), dynamic_cast<const antenna::XmlAntenna*>(ant) ? "xml" : "file");
			}

			generateDirectionalAntennaKml(out, platform, ctx, angle_3db_drop_deg, indent);
		}
	}

	void generateDynamicPathKml(std::ostream& out, const radar::Platform* platform, const std::string& styleUrl,
								const double refAlt, const KmlContext& ctx, const std::string& indent)
	{
		const math::Path* path = platform->getMotionPath();
		const auto& waypoints = path->getCoords();

		double first_alt_abs;
		{
			double lat, lon;
			ctx.converter(waypoints.front().pos, lat, lon, first_alt_abs);
		}

		out << indent << "<Placemark>\n";
		out << indent << "  <name>" << platform->getName() << " Path</name>\n";
		out << indent << "  <styleUrl>" << styleUrl << "</styleUrl>\n";
		out << indent << "  <gx:Track>\n";
		out << indent << "    <altitudeMode>absolute</altitudeMode>\n";
		if (first_alt_abs > refAlt)
		{
			out << indent << "    <extrude>1</extrude>\n";
		}

		const double start_time = waypoints.front().t;
		const double end_time = waypoints.back().t;

		if (const double time_diff = end_time - start_time; time_diff <= 0.0)
		{
			const math::Vec3 p_pos = path->getPosition(start_time);
			double p_lon, p_lat, p_alt_abs;
			ctx.converter(p_pos, p_lat, p_lon, p_alt_abs);
			out << indent << "    <when>" << start_time << "</when>\n";
			out << indent << "    <gx:coord>" << p_lon << " " << p_lat << " " << p_alt_abs << "</gx:coord>\n";
		}
		else
		{
			const double time_step = time_diff / TRACK_NUM_DIVISIONS;
			for (int i = 0; i <= TRACK_NUM_DIVISIONS; ++i)
			{
				const double current_time = start_time + i * time_step;
				const math::Vec3 p_pos = path->getPosition(current_time);
				double p_lon, p_lat, p_alt_abs;
				ctx.converter(p_pos, p_lat, p_lon, p_alt_abs);
				out << indent << "    <when>" << current_time << "</when>\n";
				out << indent << "    <gx:coord>" << p_lon << " " << p_lat << " " << p_alt_abs << "</gx:coord>\n";
			}
		}

		out << indent << "  </gx:Track>\n";
		out << indent << "</Placemark>\n";
	}

	void generateTrackEndpointsKml(std::ostream& out, const radar::Platform* platform, const double refAlt,
								   const KmlContext& ctx, const std::string& indent)
	{
		const math::Path* path = platform->getMotionPath();
		if (path->getCoords().size() <= 1)
		{
			return;
		}

		const auto& [start_wp_pos, start_wp_t] = path->getCoords().front();
		const auto& [end_wp_pos, end_wp_t] = path->getCoords().back();

		double start_lat, start_lon, start_alt_abs;
		ctx.converter(start_wp_pos, start_lat, start_lon, start_alt_abs);
		const std::string start_coordinates = formatCoordinates(start_lon, start_lat, start_alt_abs);

		double end_lat, end_lon, end_alt_abs;
		ctx.converter(end_wp_pos, end_lat, end_lon, end_alt_abs);
		const std::string end_coordinates = formatCoordinates(end_lon, end_lat, end_alt_abs);

		writePoint(out, indent, "Start: " + platform->getName(), "#target", start_coordinates, start_alt_abs, refAlt);
		writePoint(out, indent, "End: " + platform->getName(), "#target", end_coordinates, end_alt_abs, refAlt);
	}

	void generateStaticPlacemarkKml(std::ostream& out, const radar::Platform* platform, const std::string& styleUrl,
									const double refAlt, const KmlContext& ctx, const std::string& indent)
	{
		const auto& [first_wp_pos, first_wp_t] = platform->getMotionPath()->getCoords().front();
		double lat, lon, alt_abs;
		ctx.converter(first_wp_pos, lat, lon, alt_abs);
		const std::string coordinates = formatCoordinates(lon, lat, alt_abs);

		out << indent << "<Placemark>\n";
		out << indent << "  <name>" << platform->getName() << "</name>\n";
		out << indent << "  <styleUrl>" << styleUrl << "</styleUrl>\n";
		out << indent << "  <LookAt>\n";
		out << indent << "    <longitude>" << lon << "</longitude>\n";
		out << indent << "    <latitude>" << lat << "</latitude>\n";
		out << indent << "    <altitude>" << alt_abs << "</altitude>\n";
		out << indent << "    <heading>-148.41</heading><tilt>40.55</tilt><range>500.65</range>\n";
		out << indent << "  </LookAt>\n";
		out << indent << "  <Point>\n";
		out << indent << "    <coordinates>" << coordinates << "</coordinates>\n";
		out << indent << "    <altitudeMode>absolute</altitudeMode>\n";
		if (alt_abs > refAlt)
		{
			out << indent << "    <extrude>1</extrude>\n";
		}
		out << indent << "  </Point>\n";
		out << indent << "</Placemark>\n";
	}

	void generatePlatformPathKml(std::ostream& out, const radar::Platform* platform, const std::string& style,
								 const double refAlt, const KmlContext& ctx, const std::string& indent)
	{
		const auto path_type = platform->getMotionPath()->getType();
		const bool is_dynamic =
			path_type == math::Path::InterpType::INTERP_LINEAR || path_type == math::Path::InterpType::INTERP_CUBIC;

		if (is_dynamic)
		{
			generateDynamicPathKml(out, platform, style, refAlt, ctx, indent);
			generateTrackEndpointsKml(out, platform, refAlt, ctx, indent);
		}
		else
		{
			generateStaticPlacemarkKml(out, platform, style, refAlt, ctx, indent);
		}
	}

	void processPlatform(std::ostream& out, const radar::Platform* platform,
						 const std::vector<const radar::Object*>& objects, const KmlContext& ctx,
						 const double referenceAltitude, const std::string& indent)
	{
		if (platform->getMotionPath()->getCoords().empty())
		{
			return;
		}

		out << indent << "<Folder>\n";
		out << indent << "  <name>" << platform->getName() << "</name>\n";

		const std::string inner_indent = indent + "  ";
		const auto placemark_style = getPlacemarkStyleForPlatform(objects);

		if (const auto* radar_obj = getPrimaryRadar(objects))
		{
			generateAntennaKml(out, platform, radar_obj, ctx, inner_indent);
		}

		generatePlatformPathKml(out, platform, placemark_style, referenceAltitude, ctx, inner_indent);

		out << indent << "</Folder>\n";
	}

	void generateKmlToStream(std::ostream& out, const core::World& world, const KmlContext& ctx)
	{
		std::map<const radar::Platform*, std::vector<const radar::Object*>> platform_to_objects;
		const auto group_objects = [&](const auto& objectCollection)
		{
			for (const auto& obj_ptr : objectCollection)
			{
				platform_to_objects[obj_ptr->getPlatform()].push_back(obj_ptr.get());
			}
		};

		group_objects(world.getReceivers());
		group_objects(world.getTransmitters());
		group_objects(world.getTargets());

		double reference_latitude = ctx.parameters.origin_latitude;
		double reference_longitude = ctx.parameters.origin_longitude;
		double reference_altitude = ctx.parameters.origin_altitude;

		if (ctx.parameters.coordinate_frame != params::CoordinateFrame::ENU)
		{
			bool ref_set = false;
			for (const auto& platform : platform_to_objects | std::views::keys)
			{
				if (!platform->getMotionPath()->getCoords().empty())
				{
					const math::Vec3& first_pos = platform->getMotionPath()->getCoords().front().pos;
					ctx.converter(first_pos, reference_latitude, reference_longitude, reference_altitude);
					ref_set = true;
					break;
				}
			}
			if (!ref_set)
			{
				reference_latitude = ctx.parameters.origin_latitude;
				reference_longitude = ctx.parameters.origin_longitude;
				reference_altitude = ctx.parameters.origin_altitude;
			}
		}

		writeKmlHeaderAndStyles(out, ctx);

		out << "  <Folder>\n";
		out << "    <name>Reference Coordinate</name>\n";
		out << "    <description>Placemarks for various elements in the FERSXML file. All Placemarks are "
			   "situated relative to this reference point.</description>\n";
		out << "    <LookAt>\n";
		out << "      <longitude>" << reference_longitude << "</longitude>\n";
		out << "      <latitude>" << reference_latitude << "</latitude>\n";
		out << "      <altitude>" << reference_altitude << "</altitude>\n";
		out << "      <heading>-148.41</heading><tilt>40.55</tilt><range>10000</range>\n";
		out << "    </LookAt>\n";

		const std::string platform_indent = "    ";
		for (const auto& [platform, objects] : platform_to_objects)
		{
			processPlatform(out, platform, objects, ctx, reference_altitude, platform_indent);
		}

		out << "  </Folder>\n";
		out << "</Document>\n";
		out << "</kml>\n";
	}
}
