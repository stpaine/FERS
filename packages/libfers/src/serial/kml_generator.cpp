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
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/LocalCartesian.hpp>
#include <GeographicLib/UTMUPS.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <memory>
#include <optional>
#include <ranges>
#include <sstream>
#include <string>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/logging.h"
#include "core/parameters.h"
#include "core/world.h"
#include "math/coord.h"
#include "math/path.h"
#include "radar/platform.h"
#include "radar/radar_obj.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "signal/radar_signal.h"


namespace kmlGen
{
	using namespace std;
	using ConverterFunc = std::function<void(const math::Vec3&, double&, double&, double&)>;

	// --- Constants ---
	constexpr int TRACK_NUM_DIVISIONS = 100;

	constexpr int ISOTROPIC_PATTERN_POINTS = 100;

	constexpr double ISOTROPIC_PATTERN_RADIUS_KM = 20.0;

	constexpr double DIRECTIONAL_ANTENNA_ARROW_LENGTH_M = 20000.0;

	// --- Geodetic and Coordinate Helpers ---
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
		constexpr int num_points = 1000;
		std::vector<double> theta(num_points);
		std::vector<double> gain(num_points);
		for (int i = 0; i < num_points; ++i)
		{
			theta[i] = -PI + 2.0 * PI * i / (num_points - 1);
			gain[i] = sincAntennaGain(theta[i], alpha, beta, gamma);
		}
		const double max_gain = *std::max_element(gain.begin() + num_points / 2, gain.end());
		const double max_gain_db = 10.0 * std::log10(max_gain);
		const double target_gain_db = max_gain_db - 3.0;
		double target_gain = std::pow(10.0, target_gain_db / 10.0);
		const int idx = std::distance(
			gain.begin() + num_points / 2,
			std::min_element(gain.begin() + num_points / 2, gain.end(), [target_gain](const double a, const double b)
							 { return std::abs(a - target_gain) < std::abs(b - target_gain); }));
		const double angle_3db_drop = theta[idx + num_points / 2];
		return angle_3db_drop * 180.0 / PI;
	}

	/**
	 * @brief Calculates the half-power (-3dB) beamwidth angle for a Gaussian antenna.
	 * @param gaussianAnt Pointer to the Gaussian antenna object.
	 * @return The 3dB drop half-angle in degrees.
	 */
	double findGaussian3DbDropAngle(const antenna::Gaussian* gaussianAnt)
	{
		// 3dB drop is when gain is 0.5. For gaussian, G = exp(-theta^2 * scale)
		// 0.5 = exp(-theta^2 * scale) => ln(0.5) = -theta^2 * scale
		// theta = sqrt(-ln(0.5) / scale) = sqrt(ln(2) / scale)
		// We use the azimuth scale for the visualization in the horizontal plane.
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

	/**
	 * @brief Calculates the half-power (-3dB) beamwidth angle for a Parabolic antenna.
	 * @param parabolicAnt Pointer to the Parabolic antenna object.
	 * @param wavelength The operating wavelength in meters.
	 * @return The 3dB drop half-angle in degrees.
	 */
	double findParabolic3DbDropAngle(const antenna::Parabolic* parabolicAnt, const double wavelength)
	{
		// For a parabolic antenna, the gain pattern is related to (2*J1(x)/x)^2,
		// where J1 is the Bessel function of the first kind of order one.
		// The 3dB point occurs approximately when x = 1.6.
		// x = PI * diameter * sin(theta) / wavelength
		// sin(theta) = 1.6 * wavelength / (PI * diameter)
		if (parabolicAnt->getDiameter() <= 0.0)
		{
			LOG(logging::Level::WARNING,
				"Parabolic antenna '{}' has a non-positive diameter ({}). This is physically impossible. KML will only "
				"show boresight.",
				parabolicAnt->getName(), parabolicAnt->getDiameter());
			return 0.0;
		}
		// TODO: magic numbers
		const double arg = 1.6 * wavelength / (PI * parabolicAnt->getDiameter());
		// For physically realizable antennas, arg should be <= 1.
		if (arg > 1.0)
		{
			LOG(logging::Level::INFO,
				"Parabolic antenna '{}': The operating wavelength ({:.4f}m) is very large compared to its diameter "
				"({:.4f}m), resulting in a nearly omnidirectional pattern. KML visualization will cap the 3dB "
				"half-angle at 90 degrees.",
				parabolicAnt->getName(), wavelength, parabolicAnt->getDiameter());
			return 90.0; // Extremely wide beam, cap at 90 degrees.
		}
		const double half_angle_rad = std::asin(arg);
		return half_angle_rad * 180.0 / PI;
	}

	/**
	 * @brief Calculates the half-power (-3dB) beamwidth angle for a SquareHorn antenna.
	 * @param squarehornAnt Pointer to the SquareHorn antenna object.
	 * @param wavelength The operating wavelength in meters.
	 * @return The 3dB drop half-angle in degrees.
	 */
	double findSquareHorn3DbDropAngle(const antenna::SquareHorn* squarehornAnt, const double wavelength)
	{
		// The gain pattern for a square horn is related to sinc(x)^2.
		// The 3dB point occurs when sinc(x) = sqrt(0.5), which is approx. x = 1.39155.
		// x = PI * dimension * sin(theta) / wavelength
		// sin(theta) = 1.39155 * wavelength / (PI * dimension)
		if (squarehornAnt->getDimension() <= 0.0)
		{
			LOG(logging::Level::WARNING,
				"SquareHorn antenna '{}' has a non-positive dimension ({}). This is physically impossible. KML will "
				"only show boresight.",
				squarehornAnt->getName(), squarehornAnt->getDimension());
			return 0.0;
		}
		// TODO: magic numbers
		const double arg = 1.39155 * wavelength / (PI * squarehornAnt->getDimension());
		if (arg > 1.0)
		{
			LOG(logging::Level::INFO,
				"SquareHorn antenna '{}': The operating wavelength ({:.4f}m) is very large compared to its dimension "
				"({:.4f}m), resulting in a nearly omnidirectional pattern. KML visualization will cap the 3dB "
				"half-angle at 90 degrees.",
				squarehornAnt->getName(), wavelength, squarehornAnt->getDimension());
			return 90.0; // Extremely wide beam, cap at 90 degrees.
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

	// --- KML Generation Helpers ---
	void writeKmlHeaderAndStyles(std::ofstream& kmlFile)
	{
		kmlFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
		kmlFile << "<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\">\n";
		kmlFile << "<Document>\n";
		kmlFile << "  <name>";
		if (params::params.simulation_name.empty())
		{
			kmlFile << "FERS Simulation Visualization";
		}
		else
		{
			kmlFile << params::params.simulation_name;
		}
		kmlFile << "</name>\n";
		kmlFile << "  <Style "
				   "id=\"receiver\"><IconStyle><Icon><href>https://cdn-icons-png.flaticon.com/512/645/645436.png</"
				   "href></Icon></IconStyle></Style>\n";
		kmlFile << "  <Style "
				   "id=\"transmitter\"><IconStyle><Icon><href>https://cdn-icons-png.flaticon.com/128/224/224666.png</"
				   "href></Icon></IconStyle></Style>\n";
		kmlFile << "  <Style "
				   "id=\"target\"><IconStyle><Icon><href>https://upload.wikimedia.org/wikipedia/commons/thumb/a/ad/"
				   "Target_red_dot1.svg/1200px-Target_red_dot1.svg.png</href></Icon></IconStyle><LineStyle><width>2</"
				   "width></LineStyle></Style>\n";
		kmlFile << "  <Style "
				   "id=\"translucentPolygon\"><LineStyle><color>ff0000ff</color><width>2</width></"
				   "LineStyle><PolyStyle><color>00ffffff</color></PolyStyle></Style>\n";
		kmlFile << "  <Style "
				   "id=\"arrowStyle\"><IconStyle><Icon><href>http://maps.google.com/mapfiles/kml/shapes/arrow.png</"
				   "href></Icon><scale>0.5</scale></IconStyle></Style>\n";
		kmlFile << "  <Style id=\"lineStyle\"><LineStyle><color>ff0000ff</color><width>2</width></LineStyle></Style>\n";
		kmlFile
			<< "  <Style id=\"lineStyleBlue\"><LineStyle><color>ffff0000</color><width>2</width></LineStyle></Style>\n";
	}

	void writePoint(std::ofstream& kmlFile, const std::string& indent, const std::string& name,
					const std::string& styleUrl, const std::string& coordinates, const double objectAltitude,
					const double referenceAltitude)
	{
		kmlFile << indent << "<Placemark>\n";
		kmlFile << indent << "  <name>" << name << "</name>\n";
		kmlFile << indent << "  <styleUrl>" << styleUrl << "</styleUrl>\n";
		kmlFile << indent << "  <Point>\n";
		kmlFile << indent << "    <coordinates>" << coordinates << "</coordinates>\n";
		kmlFile << indent << "    <altitudeMode>absolute</altitudeMode>\n";
		if (objectAltitude > referenceAltitude)
		{
			kmlFile << indent << "    <extrude>1</extrude>\n";
		}
		kmlFile << indent << "  </Point>\n";
		kmlFile << indent << "</Placemark>\n";
	}

	void writeAntennaBeamLine(std::ofstream& kmlFile, const std::string& indent, const std::string& name,
							  const std::string& style, const std::string& startCoords, const std::string& endCoords)
	{
		kmlFile << indent << "<Placemark>\n";
		kmlFile << indent << "  <name>" << name << "</name>\n";
		kmlFile << indent << "  <styleUrl>" << style << "</styleUrl>\n";
		kmlFile << indent << "  <LineString>\n";
		kmlFile << indent << "    <altitudeMode>absolute</altitudeMode>\n";
		kmlFile << indent << "    <tessellate>1</tessellate>\n";
		kmlFile << indent << "    <coordinates>" << startCoords << " " << endCoords << "</coordinates>\n";
		kmlFile << indent << "  </LineString>\n";
		kmlFile << indent << "</Placemark>\n";
	}

	// --- Platform Processing Logic ---
	std::string getPlacemarkStyleForPlatform(const std::vector<const radar::Object*>& objects)
	{
		bool has_receiver = false;
		bool has_transmitter = false;
		for (const auto* obj : objects)
		{
			if (dynamic_cast<const radar::Receiver*>(obj))
			{
				has_receiver = true;
			}
			if (dynamic_cast<const radar::Transmitter*>(obj))
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
			if (const auto r = dynamic_cast<const radar::Radar*>(obj))
			{
				return r;
			}
		}
		return nullptr;
	}

	void generateIsotropicAntennaKml(std::ofstream& kmlFile, const math::Vec3& position, const ConverterFunc& converter,
									 const std::string& indent)
	{
		double lat, lon, alt_abs;
		converter(position, lat, lon, alt_abs);

		const auto circle_coordinates = generateCircleCoordinates(lat, lon, ISOTROPIC_PATTERN_RADIUS_KM);

		kmlFile << indent << "<Placemark>\n";
		kmlFile << indent << "  <name>Isotropic pattern range</name>\n";
		kmlFile << indent << "  <styleUrl>#translucentPolygon</styleUrl>\n";
		kmlFile << indent << "  <Polygon>\n";
		kmlFile << indent << "    <extrude>1</extrude>\n";
		kmlFile << indent << "    <altitudeMode>absolute</altitudeMode>\n";
		kmlFile << indent << "    <outerBoundaryIs><LinearRing><coordinates>\n";
		for (const auto& [pt_lat, pt_lon] : circle_coordinates)
		{
			kmlFile << indent << "      " << formatCoordinates(pt_lon, pt_lat, alt_abs) << "\n";
		}
		kmlFile << indent << "      "
				<< formatCoordinates(circle_coordinates[0].second, circle_coordinates[0].first, alt_abs) << "\n";
		kmlFile << indent << "    </coordinates></LinearRing></outerBoundaryIs>\n";
		kmlFile << indent << "  </Polygon>\n";
		kmlFile << indent << "</Placemark>\n";
	}

	void generateDirectionalAntennaKml(std::ofstream& kmlFile, const radar::Platform* platform,
									   const ConverterFunc& converter, const std::optional<double>& angle3DbDropDeg,
									   const std::string& indent)
	{
		const auto& first_wp_pos = platform->getMotionPath()->getCoords().front().pos;
		double start_lat, start_lon, start_alt;
		converter(first_wp_pos, start_lat, start_lon, start_alt);
		const std::string start_coords_str = formatCoordinates(start_lon, start_lat, start_alt);

		const math::SVec3 initial_rotation = platform->getRotationPath()->getPosition(params::startTime());

		// The parser now handles the conversion from compass heading to the internal
		// FERS format (radians, CCW from East). The KML generator needs to convert
		// this back to a standard KML heading (degrees, CW from North).
		const double fers_azimuth_deg = initial_rotation.azimuth * 180.0 / PI;
		double start_azimuth_deg_kml = 90.0 - fers_azimuth_deg;
		// Normalize to [0, 360)
		start_azimuth_deg_kml = std::fmod(start_azimuth_deg_kml, 360.0);
		if (start_azimuth_deg_kml < 0.0)
		{
			start_azimuth_deg_kml += 360.0;
		}

		// Project the arrow length onto the horizontal plane for the geodetic calculation
		// and calculate the change in altitude separately.
		const double horizontal_distance = DIRECTIONAL_ANTENNA_ARROW_LENGTH_M * std::cos(initial_rotation.elevation);
		const double delta_altitude = DIRECTIONAL_ANTENNA_ARROW_LENGTH_M * std::sin(initial_rotation.elevation);
		const double end_alt = start_alt + delta_altitude;

		// TODO: Antenna beam visualization is static, showing only the orientation at the simulation's start time.
		//       This does not represent dynamic scanning defined by a platform's <rotationpath>. The KML should
		//       ideally visualize the scan volume or animate the beam's movement using a <gx:Track> to match FERS's
		//       capabilities.
		// Main beam
		double dest_lat, dest_lon;
		calculateDestinationCoordinate(start_lat, start_lon, start_azimuth_deg_kml, horizontal_distance, dest_lat,
									   dest_lon);
		const std::string end_coords_str = formatCoordinates(dest_lon, dest_lat, end_alt);
		writeAntennaBeamLine(kmlFile, indent, "Antenna Boresight", "#lineStyle", start_coords_str, end_coords_str);

		// 3dB beamwidth lines, if angle is provided and is greater than a small epsilon
		if (angle3DbDropDeg.has_value() && *angle3DbDropDeg > EPSILON)
		{
			double side1_lat, side1_lon;
			calculateDestinationCoordinate(start_lat, start_lon, start_azimuth_deg_kml - *angle3DbDropDeg,
										   horizontal_distance, side1_lat, side1_lon);
			const std::string side1_coords_str = formatCoordinates(side1_lon, side1_lat, end_alt);
			writeAntennaBeamLine(kmlFile, indent, "Antenna 3dB Beamwidth", "#lineStyleBlue", start_coords_str,
								 side1_coords_str);

			double side2_lat, side2_lon;
			calculateDestinationCoordinate(start_lat, start_lon, start_azimuth_deg_kml + *angle3DbDropDeg,
										   horizontal_distance, side2_lat, side2_lon);
			const std::string side2_coords_str = formatCoordinates(side2_lon, side2_lat, end_alt);
			writeAntennaBeamLine(kmlFile, indent, "Antenna 3dB Beamwidth", "#lineStyleBlue", start_coords_str,
								 side2_coords_str);
		}

		// Arrow placemark
		const double arrow_heading = std::fmod(start_azimuth_deg_kml + 180.0, 360.0);
		kmlFile << indent << "<Placemark>\n";
		kmlFile << indent << "  <name>Antenna Arrow</name>\n";
		kmlFile << indent << "  <styleUrl>#arrowStyle</styleUrl>\n";
		kmlFile << indent << "  <Point><coordinates>" << end_coords_str
				<< "</coordinates><altitudeMode>absolute</altitudeMode></Point>\n";
		kmlFile << indent << "  <Style>\n";
		kmlFile << indent << "    <IconStyle><heading>" << arrow_heading << "</heading></IconStyle>\n";
		kmlFile << indent << "  </Style>\n";
		kmlFile << indent << "</Placemark>\n";
	}

	void generateAntennaKml(std::ofstream& kmlFile, const radar::Platform* platform, const radar::Radar* radar,
							const ConverterFunc& converter, const std::string& indent)
	{
		const antenna::Antenna* ant = radar->getAntenna();
		if (!ant || platform->getMotionPath()->getCoords().empty())
		{
			return;
		}

		if (dynamic_cast<const antenna::Isotropic*>(ant))
		{
			const math::Vec3 initial_pos = platform->getMotionPath()->getCoords().front().pos;
			generateIsotropicAntennaKml(kmlFile, initial_pos, converter, indent);
		}
		else // Handle all directional antennas
		{
			std::optional<double> angle_3db_drop_deg;

			// Attempt to find wavelength for wavelength-dependent patterns
			std::optional<double> wavelength;
			if (const auto* tx = dynamic_cast<const radar::Transmitter*>(radar))
			{
				if (tx->getSignal())
				{
					wavelength = params::c() / tx->getSignal()->getCarrier();
				}
			}
			else if (const auto* rx = dynamic_cast<const radar::Receiver*>(radar))
			{
				if (const auto* attached_tx = dynamic_cast<const radar::Transmitter*>(rx->getAttached()))
				{
					if (attached_tx->getSignal())
					{
						wavelength = params::c() / attached_tx->getSignal()->getCarrier();
					}
				}
			}

			// Calculate 3dB drop angle based on antenna type
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
			else if (dynamic_cast<const antenna::XmlAntenna*>(ant) || dynamic_cast<const antenna::H5Antenna*>(ant))
			{
				// For XmlAntenna and H5Antenna, angle_3db_drop_deg remains nullopt,
				// resulting in only the boresight arrow being drawn. This is an intentional
				// symbolic representation. Alert the user about this.
				LOG(logging::Level::INFO,
					"KML visualization for antenna '{}' ('{}') is symbolic. "
					"Only the boresight direction is shown, as a 3dB beamwidth is not calculated from file-based "
					"patterns.",
					ant->getName(), dynamic_cast<const antenna::XmlAntenna*>(ant) ? "xml" : "file");
			}

			generateDirectionalAntennaKml(kmlFile, platform, converter, angle_3db_drop_deg, indent);
		}
	}

	void generateDynamicPathKml(std::ofstream& kmlFile, const radar::Platform* platform, const std::string& styleUrl,
								const double refAlt, const ConverterFunc& converter, const std::string& indent)
	{
		const math::Path* path = platform->getMotionPath();
		const auto& waypoints = path->getCoords();

		double first_alt_abs;
		{
			double lat, lon;
			converter(waypoints.front().pos, lat, lon, first_alt_abs);
		}

		kmlFile << indent << "<Placemark>\n";
		kmlFile << indent << "  <name>" << platform->getName() << " Path</name>\n";
		kmlFile << indent << "  <styleUrl>" << styleUrl << "</styleUrl>\n";
		kmlFile << indent << "  <gx:Track>\n";
		kmlFile << indent << "    <altitudeMode>absolute</altitudeMode>\n";
		if (first_alt_abs > refAlt)
		{
			kmlFile << indent << "    <extrude>1</extrude>\n";
		}

		// The sampling time range is now based on the platform's specific motion path duration,
		// ensuring accurate track resolution for objects with short lifespans.
		const double start_time = waypoints.front().t;
		const double end_time = waypoints.back().t;

		// Handle single-point paths or paths with zero duration by emitting a single coordinate.
		if (const double time_diff = end_time - start_time; time_diff <= 0.0)
		{
			const math::Vec3 p_pos = path->getPosition(start_time);
			double p_lon, p_lat, p_alt_abs;
			converter(p_pos, p_lat, p_lon, p_alt_abs);
			kmlFile << indent << "    <when>" << start_time << "</when>\n";
			kmlFile << indent << "    <gx:coord>" << p_lon << " " << p_lat << " " << p_alt_abs << "</gx:coord>\n";
		}
		else
		{
			const double time_step = time_diff / TRACK_NUM_DIVISIONS;
			for (int i = 0; i <= TRACK_NUM_DIVISIONS; ++i)
			{
				const double current_time = start_time + i * time_step;
				const math::Vec3 p_pos = path->getPosition(current_time);
				double p_lon, p_lat, p_alt_abs;
				converter(p_pos, p_lat, p_lon, p_alt_abs);
				kmlFile << indent << "    <when>" << current_time << "</when>\n";
				kmlFile << indent << "    <gx:coord>" << p_lon << " " << p_lat << " " << p_alt_abs << "</gx:coord>\n";
			}
		}

		kmlFile << indent << "  </gx:Track>\n";
		kmlFile << indent << "</Placemark>\n";
	}

	void generateTrackEndpointsKml(std::ofstream& kmlFile, const radar::Platform* platform, const double refAlt,
								   const ConverterFunc& converter, const std::string& indent)
	{
		const math::Path* path = platform->getMotionPath();
		if (path->getCoords().size() <= 1)
		{
			return;
		}

		const auto& [start_wp_pos, start_wp_t] = path->getCoords().front();
		const auto& [end_wp_pos, end_wp_t] = path->getCoords().back();

		double start_lat, start_lon, start_alt_abs;
		converter(start_wp_pos, start_lat, start_lon, start_alt_abs);
		const std::string start_coordinates = formatCoordinates(start_lon, start_lat, start_alt_abs);

		double end_lat, end_lon, end_alt_abs;
		converter(end_wp_pos, end_lat, end_lon, end_alt_abs);
		const std::string end_coordinates = formatCoordinates(end_lon, end_lat, end_alt_abs);

		writePoint(kmlFile, indent, "Start: " + platform->getName(), "#target", start_coordinates, start_alt_abs,
				   refAlt);
		writePoint(kmlFile, indent, "End: " + platform->getName(), "#target", end_coordinates, end_alt_abs, refAlt);
	}

	void generateStaticPlacemarkKml(std::ofstream& kmlFile, const radar::Platform* platform,
									const std::string& styleUrl, const double refAlt, const ConverterFunc& converter,
									const std::string& indent)
	{
		const auto& [first_wp_pos, first_wp_t] = platform->getMotionPath()->getCoords().front();
		double lat, lon, alt_abs;
		converter(first_wp_pos, lat, lon, alt_abs);
		const std::string coordinates = formatCoordinates(lon, lat, alt_abs);

		kmlFile << indent << "<Placemark>\n";
		kmlFile << indent << "  <name>" << platform->getName() << "</name>\n";
		kmlFile << indent << "  <styleUrl>" << styleUrl << "</styleUrl>\n";
		kmlFile << indent << "  <LookAt>\n";
		kmlFile << indent << "    <longitude>" << lon << "</longitude>\n";
		kmlFile << indent << "    <latitude>" << lat << "</latitude>\n";
		kmlFile << indent << "    <altitude>" << alt_abs << "</altitude>\n";
		kmlFile << indent << "    <heading>-148.41</heading><tilt>40.55</tilt><range>500.65</range>\n";
		kmlFile << indent << "  </LookAt>\n";
		kmlFile << indent << "  <Point>\n";
		kmlFile << indent << "    <coordinates>" << coordinates << "</coordinates>\n";
		kmlFile << indent << "    <altitudeMode>absolute</altitudeMode>\n";
		if (alt_abs > refAlt)
		{
			kmlFile << indent << "    <extrude>1</extrude>\n";
		}
		kmlFile << indent << "  </Point>\n";
		kmlFile << indent << "</Placemark>\n";
	}

	void generatePlatformPathKml(std::ofstream& kmlFile, const radar::Platform* platform, const std::string& style,
								 const double refAlt, const ConverterFunc& converter, const std::string& indent)
	{
		const auto path_type = platform->getMotionPath()->getType();
		const bool is_dynamic =
			path_type == math::Path::InterpType::INTERP_LINEAR || path_type == math::Path::InterpType::INTERP_CUBIC;

		if (is_dynamic)
		{
			generateDynamicPathKml(kmlFile, platform, style, refAlt, converter, indent);
			generateTrackEndpointsKml(kmlFile, platform, refAlt, converter, indent);
		}
		else
		{
			generateStaticPlacemarkKml(kmlFile, platform, style, refAlt, converter, indent);
		}
	}

	void processPlatform(const radar::Platform* platform, const std::vector<const radar::Object*>& objects,
						 std::ofstream& kmlFile, const ConverterFunc& converter, const double referenceAltitude,
						 const std::string& indent)
	{
		if (platform->getMotionPath()->getCoords().empty())
		{
			return;
		}

		kmlFile << indent << "<Folder>\n";
		kmlFile << indent << "  <name>" << platform->getName() << "</name>\n";

		const std::string inner_indent = indent + "  ";
		const auto placemark_style = getPlacemarkStyleForPlatform(objects);

		if (const auto* radar_obj = getPrimaryRadar(objects))
		{
			generateAntennaKml(kmlFile, platform, radar_obj, converter, inner_indent);
		}

		generatePlatformPathKml(kmlFile, platform, placemark_style, referenceAltitude, converter, inner_indent);

		kmlFile << indent << "</Folder>\n";
	}

}

namespace serial
{
	bool KmlGenerator::generateKml(const core::World& world, const std::string& outputKmlPath)
	{
		try
		{
			// Setup coordinate conversion based on global parameters
			kmlGen::ConverterFunc converter;
			double reference_latitude, reference_longitude, reference_altitude;

			switch (params::coordinateFrame())
			{
			case params::CoordinateFrame::ENU:
				{
					reference_latitude = params::originLatitude();
					reference_longitude = params::originLongitude();
					reference_altitude = params::originAltitude();
					auto proj = std::make_shared<GeographicLib::LocalCartesian>(reference_latitude, reference_longitude,
																				reference_altitude);
					converter = [proj](const math::Vec3& pos, double& lat, double& lon, double& alt)
					{ proj->Reverse(pos.x, pos.y, pos.z, lat, lon, alt); };
					break;
				}
			case params::CoordinateFrame::UTM:
				{
					const int zone = params::utmZone();
					const bool northp = params::utmNorthHemisphere();
					converter = [zone, northp](const math::Vec3& pos, double& lat, double& lon, double& alt)
					{
						double gamma, k;
						GeographicLib::UTMUPS::Reverse(zone, northp, pos.x, pos.y, lat, lon, gamma, k);
						alt = pos.z; // Altitude is given directly in the z-coordinate
					};
					break;
				}
			case params::CoordinateFrame::ECEF:
				{
					const auto& earth = GeographicLib::Geocentric::WGS84();
					converter = [&earth](const math::Vec3& pos, double& lat, double& lon, double& alt)
					{ earth.Reverse(pos.x, pos.y, pos.z, lat, lon, alt); };
					break;
				}
			}

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

			if (params::coordinateFrame() != params::CoordinateFrame::ENU)
			{
				bool ref_set = false;
				for (const auto& platform : platform_to_objects | std::views::keys)
				{
					if (!platform->getMotionPath()->getCoords().empty())
					{
						const math::Vec3& first_pos = platform->getMotionPath()->getCoords().front().pos;
						converter(first_pos, reference_latitude, reference_longitude, reference_altitude);
						ref_set = true;
						break;
					}
				}
				if (!ref_set) // Fallback if no platforms or no waypoints
				{
					reference_latitude = params::originLatitude(); // UCT default
					reference_longitude = params::originLongitude();
					reference_altitude = params::originAltitude();
				}
			}

			std::ofstream kml_file(outputKmlPath.c_str());
			if (!kml_file.is_open())
			{
				LOG(logging::Level::ERROR, "Error opening output KML file {}", outputKmlPath);
				return false;
			}

			kmlGen::writeKmlHeaderAndStyles(kml_file);

			kml_file << "  <Folder>\n";
			kml_file << "    <name>Reference Coordinate</name>\n";
			kml_file << "    <description>Placemarks for various elements in the FERSXML file. All Placemarks are "
						"situated relative to this reference point.</description>\n";
			kml_file << "    <LookAt>\n";
			kml_file << "      <longitude>" << reference_longitude << "</longitude>\n";
			kml_file << "      <latitude>" << reference_latitude << "</latitude>\n";
			kml_file << "      <altitude>" << reference_altitude << "</altitude>\n";
			kml_file << "      <heading>-148.41</heading><tilt>40.55</tilt><range>10000</range>\n";
			kml_file << "    </LookAt>\n";

			const std::string platform_indent = "    ";
			for (const auto& [platform, objects] : platform_to_objects)
			{
				kmlGen::processPlatform(platform, objects, kml_file, converter, reference_altitude, platform_indent);
			}

			kml_file << "  </Folder>\n";
			kml_file << "</Document>\n";
			kml_file << "</kml>\n";
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
