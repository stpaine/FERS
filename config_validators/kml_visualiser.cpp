// FERS input Validator Function sub-system
// Outputs KML GIS readable file.
// Script written by Michael Altshuler
// University of Cape Town: ALTMIC003

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>
#include <algorithm>
#include <vector>
#include <iomanip> // Include the iomanip header for setprecision

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/XMLString.hpp>

using namespace std;
using namespace xercesc;

// The following two functions calculate the 3dB drop from the max gain of the antenna gain pattern
/*Start*/
double sinc_antenna_gain(double theta, double alpha, double beta, double gamma)
{
    double gain = alpha * std::pow(std::sin(beta * theta) / (beta * theta), gamma);
    return gain;
}

double find_3db_drop_angle(double alpha, double beta, double gamma)
{
    const int num_points = 1000;
    const double pi = 3.14159265358979323846;
    std::vector<double> theta(num_points);
    std::vector<double> gain(num_points);

    // Calculate gain values for each angle
    for (int i = 0; i < num_points; ++i)
    {
        theta[i] = -pi + 2.0 * pi * i / (num_points - 1);
        gain[i] = sinc_antenna_gain(theta[i], alpha, beta, gamma);
    }

    // Find the maximum gain
    double max_gain = *std::max_element(gain.begin() + num_points / 2, gain.end());

    // Convert the maximum gain to dB
    double max_gain_dB = 10.0 * std::log10(max_gain);

    // Find the target gain (3dB drop)
    double target_gain_dB = max_gain_dB - 3.0;
    double target_gain = std::pow(10.0, target_gain_dB / 10.0);

    // Find the angle index where the gain is closest to the target gain
    // Considering only positive angles (from 0 to pi)
    int idx = std::distance(gain.begin() + num_points / 2, std::min_element(gain.begin() + num_points / 2, gain.end(),
                                                                            [target_gain](double a, double b)
                                                                            { return std::abs(a - target_gain) < std::abs(b - target_gain); }));

    // Get the angle at which the 3dB drop occurs
    double angle_3dB_drop = theta[idx + num_points / 2];

    // Convert the angle to degrees
    double angle_3dB_drop_deg = angle_3dB_drop * 180.0 / pi;

    return angle_3dB_drop_deg;
}
/*End*/

// Function returns coordinates from positionWayPoint values
std::string getCoordinatesFromPositionWaypoint(const DOMElement *positionWaypointElement, double referenceLatitude, double referenceLongitude, double referenceAltitude)
{
    const XMLCh *xTag = XMLString::transcode("x");
    const XMLCh *yTag = XMLString::transcode("y");
    const XMLCh *altitudeTag = XMLString::transcode("altitude");

    double x = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(xTag)->item(0)->getTextContent()));
    double y = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(yTag)->item(0)->getTextContent()));
    double altitude = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(altitudeTag)->item(0)->getTextContent()));

    double longitude = referenceLongitude + x / (cos(referenceLatitude * M_PI / 180) * 111319.9);
    double latitude = referenceLatitude + y / 111319.9;
    double altitudeAboveGround = altitude - referenceAltitude;

    std::stringstream coordinates;
    coordinates << std::fixed << std::setprecision(6) << longitude << "," << latitude << "," << altitudeAboveGround;

    return coordinates.str();
}

// Function to calculate the destination coordinate given starting coordinate, angle (in degrees), and distance (in meters)
void calculateDestinationCoordinate(double startLatitude, double startLongitude, double angle, double distance, double &destLatitude, double &destLongitude)
{
    const double R = 6371000; // Earth's radius in meters
    double d = distance / R;  // Angular distance in radians

    // Convert degrees to radians
    double startLatRad = startLatitude * M_PI / 180;
    double startLonRad = startLongitude * M_PI / 180;
    double angleRad = angle * M_PI / 180;

    // Calculate destination latitude and longitude in radians
    double destLatRad = asin(sin(startLatRad) * cos(d) + cos(startLatRad) * sin(d) * cos(angleRad));
    double destLonRad = startLonRad + atan2(sin(angleRad) * sin(d) * cos(startLatRad), cos(d) - sin(startLatRad) * sin(destLatRad));

    // Convert radians to degrees
    destLatitude = destLatRad * 180 / M_PI;
    destLongitude = destLonRad * 180 / M_PI;
}

// Function calculates the hyperbolic path and updates the longitude and latitude values. *Not a valid interpolation path at this time (04/05/2023)
void updateLongitudeLatitudeHyperbolic(double &longitude, double &latitude, double t, double a, double b)
{
    double x_hyperbolic = a * std::cosh(t);
    double y_hyperbolic = b * std::sinh(t);

    // Assuming the start point of the hyperbolic path is at the origin
    longitude += x_hyperbolic / (cos(latitude * M_PI / 180) * 111319.9);
    latitude += y_hyperbolic / 111319.9;
}

// Older implementation (deprecated)
void updateLongitudeLatitudeCubic2(double &x, double &y, double t, double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
    double t2 = t * t;
    double t3 = t2 * t;

    double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    double h10 = -2.0 * t3 + 3.0 * t2;
    double h01 = t3 - 2.0 * t2 + t;
    double h11 = t3 - t2;

    x = h00 * x1 + h10 * x2 + h01 * x3 + h11 * x4;
    y = h00 * y1 + h10 * y2 + h01 * y3 + h11 * y4;
}

// Function calculates the cubic path and updates the longitude and latitude values.
void updateLongitudeLatitudeCubic(double &newLongitude, double &newLatitude, double t, double longitude1, double latitude1, double longitude4, double latitude4)
{
    double t2 = t * t;
    double t3 = t2 * t;

    double controlPointAngle = 45.0 * M_PI / 180.0; // 45 degrees in radians
    double controlPointDistance = 111319.9;         // Fixed distance for control points (e.g. 111319.9 meters or 1 degree)

    // Calculate control points based on fixed angle and distance
    double x2 = longitude1 + controlPointDistance * cos(controlPointAngle) / (cos(latitude1 * M_PI / 180) * 111319.9);
    double y2 = latitude1 + controlPointDistance * sin(controlPointAngle) / 111319.9;

    double x3 = longitude4 - controlPointDistance * cos(controlPointAngle) / (cos(latitude4 * M_PI / 180) * 111319.9);
    double y3 = latitude4 - controlPointDistance * sin(controlPointAngle) / 111319.9;

    // Perform cubic interpolation
    double one_minus_t = 1 - t;
    double one_minus_t2 = one_minus_t * one_minus_t;
    double one_minus_t3 = one_minus_t2 * one_minus_t;

    newLongitude = one_minus_t3 * longitude1 + 3 * one_minus_t2 * t * x2 + 3 * one_minus_t * t2 * x3 + t3 * longitude4;
    newLatitude = one_minus_t3 * latitude1 + 3 * one_minus_t2 * t * y2 + 3 * one_minus_t * t2 * y3 + t3 * latitude4;
}

// Function to populate antenna maps
void populateAntennaMaps(const DOMElement *element, std::map<std::string, const DOMElement *> &isotropic_antennas, std::map<std::string, const DOMElement *> &patterned_antennas)
{
    DOMNodeList *antennaElements = element->getElementsByTagName(XMLString::transcode("antenna"));
    for (XMLSize_t i = 0; i < antennaElements->getLength(); i++)
    {
        const DOMElement *antennaElement = dynamic_cast<const DOMElement *>(antennaElements->item(i));
        const XMLCh *nameAttr = XMLString::transcode("name");
        const XMLCh *patternAttr = XMLString::transcode("pattern");
        const XMLCh *nameValue = antennaElement->getAttribute(nameAttr);
        const XMLCh *patternValue = antennaElement->getAttribute(patternAttr);
        std::string nameStr = XMLString::transcode(nameValue);

        if (XMLString::equals(patternValue, XMLString::transcode("isotropic")))
        {
            isotropic_antennas[nameStr] = antennaElement;
        }
        else
        {
            patterned_antennas[nameStr] = antennaElement;
        }
    }
}

// Function to check if an antenna name is isotropic
bool isAntennaIsotropic(const std::string &antennaName, const std::map<std::string, const DOMElement *> &isotropic_antennas)
{
    return isotropic_antennas.find(antennaName) != isotropic_antennas.end();
}

// Function that converts degrees to radians
double deg2rad(double degrees)
{
    return degrees * M_PI / 180.0;
}

// Function to generate coordinates for a circle around a given latitude, longitude and radius
std::vector<std::pair<double, double>> generate_circle_coordinates(double lat, double lon, double radius_km, int num_points = 100)
{
    std::vector<std::pair<double, double>> circle_coordinates;
    double radius_earth = 6371.0; // Earth's radius in km

    for (int i = 0; i < num_points; i++)
    {
        double bearing = deg2rad(i * 360.0 / num_points);
        double lat_rad = deg2rad(lat);
        double lon_rad = deg2rad(lon);
        double angular_distance = radius_km / radius_earth;

        double new_lat_rad = asin(sin(lat_rad) * cos(angular_distance) +
                                  cos(lat_rad) * sin(angular_distance) * cos(bearing));
        double new_lon_rad = lon_rad + atan2(sin(bearing) * sin(angular_distance) * cos(lat_rad),
                                             cos(angular_distance) - sin(lat_rad) * sin(new_lat_rad));

        double new_lat = new_lat_rad * 180.0 / M_PI;
        double new_lon = new_lon_rad * 180.0 / M_PI;

        circle_coordinates.push_back(std::make_pair(new_lat, new_lon));
    }

    return circle_coordinates;
}

// Function to return antenna element with 'sinc' pattern
const DOMElement *getAntennaElementWithSincPattern(const DOMElement *rootElement)
{
    const XMLCh *antennaTag = XMLString::transcode("antenna");
    const XMLCh *patternAttribute = XMLString::transcode("pattern");
    DOMNodeList *antennaList = rootElement->getElementsByTagName(antennaTag);

    for (XMLSize_t i = 0; i < antennaList->getLength(); ++i)
    {
        const DOMElement *currentAntennaElement = dynamic_cast<const DOMElement *>(antennaList->item(i));
        if (XMLString::equals(currentAntennaElement->getAttribute(patternAttribute), XMLString::transcode("sinc")))
        {
            return currentAntennaElement;
        }
    }
    return nullptr;
}

// Function to process the DOMElement, extract necessary data from FERSXML file, and output accordingly to KML file
void processElement(const DOMElement *element, std::ofstream &kmlFile, double referenceLatitude, double referenceLongitude, double referenceAltitude, DOMDocument *document)
{

    // Defining constants
    const XMLCh *platformTag = XMLString::transcode("platform");
    const XMLCh *receiverTag = XMLString::transcode("receiver");
    const XMLCh *transmitterTag = XMLString::transcode("transmitter");
    const XMLCh *targetTag = XMLString::transcode("target");
    const XMLCh *positionWaypointTag = XMLString::transcode("positionwaypoint");
    const XMLCh *xTag = XMLString::transcode("x");
    const XMLCh *yTag = XMLString::transcode("y");
    const XMLCh *altitudeTag = XMLString::transcode("altitude");
    const XMLCh *motionPathTag = XMLString::transcode("motionpath");
    const XMLCh *interpolationAttr = XMLString::transcode("interpolation");
    const XMLCh *alphaTag = XMLString::transcode("alpha");
    const XMLCh *betaTag = XMLString::transcode("beta");
    const XMLCh *gammaTag = XMLString::transcode("gamma");

    // Define maps to store isotropic and patterned antennas
    std::map<std::string, const DOMElement *> isotropic_antennas;
    std::map<std::string, const DOMElement *> patterned_antennas;
    populateAntennaMaps(document->getDocumentElement(), isotropic_antennas, patterned_antennas);

    // Check if the element is a platform
    if (XMLString::equals(element->getTagName(), platformTag))
    {
        // Get the positionwaypoint element
        // Check if the getElementsByTagName() function returns a valid result before using it
        DOMNodeList *positionWaypointList = element->getElementsByTagName(positionWaypointTag);
        if (positionWaypointList->getLength() == 0)
        {
            return;
        }

        const DOMElement *positionWaypointElement = dynamic_cast<const DOMElement *>(element->getElementsByTagName(positionWaypointTag)->item(0));

        // Extract the position coordinates
        double x = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(xTag)->item(0)->getTextContent()));
        double y = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(yTag)->item(0)->getTextContent()));
        double altitude = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(altitudeTag)->item(0)->getTextContent()));

        // Rough estimation equirectangular projection method.
        double longitude = referenceLongitude + x / (cos(referenceLatitude * M_PI / 180) * 111319.9);
        double latitude = referenceLatitude + y / 111319.9;
        double altitudeAboveGround = altitude - referenceAltitude;

        // Get the motionpath element
        const DOMElement *motionPathElement = dynamic_cast<const DOMElement *>(element->getElementsByTagName(motionPathTag)->item(0));

        // Extract the interpolation attribute
        const XMLCh *interpolation = motionPathElement->getAttribute(interpolationAttr);

        // Determine if the interpolation is linear, hyperbolic or cubic
        bool isLinear = (XMLString::equals(interpolation, XMLString::transcode("linear")));
        bool isHyperbolic = (XMLString::equals(interpolation, XMLString::transcode("hyperbolic")));
        bool isCubic = (XMLString::equals(interpolation, XMLString::transcode("cubic")));

        // Determine the type of placemark to use
        std::string placemarkStyle;
        if (element->getElementsByTagName(receiverTag)->getLength() > 0)
        {
            placemarkStyle = "receiver";
        }
        else if (element->getElementsByTagName(transmitterTag)->getLength() > 0)
        {
            placemarkStyle = "transmitter";
        }
        else if (element->getElementsByTagName(targetTag)->getLength() > 0)
        {
            placemarkStyle = "target";
        }

        // Determine if antenna 'pattern' is isotropic or patterned
        bool isIsotropic = false;

        if (element->getElementsByTagName(receiverTag)->getLength() > 0)
        {
            const DOMElement *receiverElement = dynamic_cast<const DOMElement *>(element->getElementsByTagName(receiverTag)->item(0));
            const XMLCh *antennaAttr = XMLString::transcode("antenna");
            const XMLCh *antennaValue = receiverElement->getAttribute(antennaAttr);
            std::string antennaName = XMLString::transcode(antennaValue);

            isIsotropic = isAntennaIsotropic(antennaName, isotropic_antennas);
        }
        else if (element->getElementsByTagName(transmitterTag)->getLength() > 0)
        {
            const DOMElement *transmitterElement = dynamic_cast<const DOMElement *>(element->getElementsByTagName(transmitterTag)->item(0));
            const XMLCh *antennaAttr = XMLString::transcode("antenna");
            const XMLCh *antennaValue = transmitterElement->getAttribute(antennaAttr);
            std::string antennaName = XMLString::transcode(antennaValue);

            isIsotropic = isAntennaIsotropic(antennaName, isotropic_antennas);
        }

        // Testing if 'isIsotropic' set to True if pattern = isotropic
        // std::cout << isIsotropic << std::endl;

        // If the associated pattern is isotropic, add a circular ring of radius 20 km
        if (isIsotropic)
        {
            double circle_radius = 20; // Radius in km
            int num_points = 100;      // Number of points to form the circle

            kmlFile << "<Placemark>\n";
            kmlFile << "    <name>Isotropic pattern range</name>\n";
            kmlFile << "    <styleUrl>#translucentPolygon</styleUrl>\n";
            std::vector<std::pair<double, double>> circle_coordinates = generate_circle_coordinates(latitude, longitude, circle_radius, num_points);
            kmlFile << "    <Polygon>\n";
            kmlFile << "        <extrude>1</extrude>\n";
            kmlFile << "        <altitudeMode>relativeToGround</altitudeMode>\n";
            kmlFile << "        <outerBoundaryIs>\n";
            kmlFile << "            <LinearRing>\n";
            kmlFile << "                <coordinates>\n";

            for (const auto &coord : circle_coordinates)
            {
                kmlFile << "                    " << coord.second << "," << coord.first << "," << altitudeAboveGround << "\n";
            }
            // Close the circle by repeating the first point
            kmlFile << "                    " << circle_coordinates[0].second << "," << circle_coordinates[0].first << "," << altitudeAboveGround << "\n";

            kmlFile << "                </coordinates>\n";
            kmlFile << "            </LinearRing>\n";
            kmlFile << "        </outerBoundaryIs>\n";
            kmlFile << "    </Polygon>\n";
            kmlFile << "</Placemark>\n";
        }
        else if (!isIsotropic && (element->getElementsByTagName(transmitterTag)->getLength() > 0 || element->getElementsByTagName(receiverTag)->getLength() > 0))
        {
            // Extract necessary values
            const XMLCh *startAzimuthTag = XMLString::transcode("startazimuth");
            const XMLCh *positionWaypointTag = XMLString::transcode("positionwaypoint");

            double startAzimuth = std::stod(XMLString::transcode(element->getElementsByTagName(startAzimuthTag)->item(0)->getTextContent()));
            const DOMElement *positionWaypointElement = dynamic_cast<DOMElement *>(element->getElementsByTagName(positionWaypointTag)->item(0));

            // Convert positionWaypoint to coordinates
            std::string coordinates = getCoordinatesFromPositionWaypoint(positionWaypointElement, referenceLatitude, referenceLongitude, referenceAltitude);

            // Calculate end coordinates
            double arrowLength = 20000; // Adjust this value according to the desired length of the arrow

            // Parse coordinates from the positionWaypoint
            double startLatitude, startLongitude, startAltitude;
            std::istringstream coordinatesStream(coordinates);
            coordinatesStream >> startLongitude;
            coordinatesStream.ignore(1); // skip the comma
            coordinatesStream >> startLatitude;
            coordinatesStream.ignore(1); // skip the comma
            coordinatesStream >> startAltitude;

            // Adjust startAzimuth to make 0 degrees point East
            startAzimuth = startAzimuth + 180;

            // Calculate end coordinates
            double destLatitude, destLongitude;
            calculateDestinationCoordinate(startLatitude, startLongitude, startAzimuth, arrowLength, destLatitude, destLongitude);

            std::stringstream endCoordinatesStream;
            endCoordinatesStream << std::fixed << std::setprecision(6) << destLongitude << "," << destLatitude << "," << startAltitude;
            std::string endCoordinates = endCoordinatesStream.str();

            // Define values for testing
            // double alpha = 1;
            // double beta = 2;
            // double gamma = 3.6;

            // Extract the antenna element with pattern="sinc"
            const DOMElement *sincAntennaElement = getAntennaElementWithSincPattern(document->getDocumentElement());

            // Extract alpha, beta, and gamma values if the antenna element was found
            if (sincAntennaElement != nullptr)
            {
                double alpha = std::stod(XMLString::transcode(sincAntennaElement->getElementsByTagName(alphaTag)->item(0)->getTextContent()));
                double beta = std::stod(XMLString::transcode(sincAntennaElement->getElementsByTagName(betaTag)->item(0)->getTextContent()));
                double gamma = std::stod(XMLString::transcode(sincAntennaElement->getElementsByTagName(gammaTag)->item(0)->getTextContent()));

                double angle_3dB_drop_deg = find_3db_drop_angle(alpha, beta, gamma);

                // Calculate end coordinates for both side lines
                double sideLine1Azimuth = startAzimuth - angle_3dB_drop_deg;
                double sideLine2Azimuth = startAzimuth + angle_3dB_drop_deg;
                double sideLine1DestLatitude, sideLine1DestLongitude;
                double sideLine2DestLatitude, sideLine2DestLongitude;

                calculateDestinationCoordinate(startLatitude, startLongitude, sideLine1Azimuth, arrowLength, sideLine1DestLatitude, sideLine1DestLongitude);
                calculateDestinationCoordinate(startLatitude, startLongitude, sideLine2Azimuth, arrowLength, sideLine2DestLatitude, sideLine2DestLongitude);

                std::stringstream sideLine1EndCoordinatesStream, sideLine2EndCoordinatesStream;
                sideLine1EndCoordinatesStream << std::fixed << std::setprecision(6) << sideLine1DestLongitude << "," << sideLine1DestLatitude << "," << startAltitude;
                sideLine2EndCoordinatesStream << std::fixed << std::setprecision(6) << sideLine2DestLongitude << "," << sideLine2DestLatitude << "," << startAltitude;
                std::string sideLine1EndCoordinates = sideLine1EndCoordinatesStream.str();
                std::string sideLine2EndCoordinates = sideLine2EndCoordinatesStream.str();

                // Add placemarks for side lines
                for (int i = 1; i <= 2; ++i)
                {
                    std::string sideLineName = "Antenna Side Line " + std::to_string(i);
                    std::string sideLineEndCoordinates = (i == 1) ? sideLine1EndCoordinates : sideLine2EndCoordinates;

                    kmlFile << "<Placemark>\n";
                    kmlFile << "      <name>" << sideLineName << "</name>\n";
                    kmlFile << "      <styleUrl>#lineStyleBlue</styleUrl>\n";
                    kmlFile << "      <LineString>\n";
                    kmlFile << "            <tessellate>1</tessellate>\n";
                    kmlFile << "            <coordinates>\n";
                    kmlFile << "            " + coordinates + " " + sideLineEndCoordinates + "\n";
                    kmlFile << "            </coordinates>\n";
                    kmlFile << "      </LineString>\n";
                    kmlFile << "</Placemark>\n";
                }
            }
            else
            {
                std::cerr << "Error: Antenna element with pattern='sinc' not found in the XML file" << std::endl;
            }

            kmlFile << "<Placemark>\n";
            kmlFile << "      <name>Antenna Direction</name>\n";
            kmlFile << "      <styleUrl>#lineStyle</styleUrl>\n";
            kmlFile << "      <LineString>\n";
            kmlFile << "            <tessellate>1</tessellate>\n";
            kmlFile << "            <coordinates>\n";
            kmlFile << "            " + coordinates + " " + endCoordinates + "\n";
            kmlFile << "            </coordinates>\n";
            kmlFile << "      </LineString>\n";
            kmlFile << "</Placemark>\n";

            kmlFile << "<Placemark>\n";
            kmlFile << "      <name>Antenna Arrow</name>\n";
            kmlFile << "      <styleUrl>#arrowStyle</styleUrl>\n";
            kmlFile << "      <Point>\n";
            kmlFile << "          <coordinates>" + endCoordinates + "</coordinates>\n";
            kmlFile << "      </Point>\n";
            kmlFile << "      <IconStyle>\n";
            kmlFile << "          <heading>" << startAzimuth << "</heading>\n";
            kmlFile << "      </IconStyle>\n";
            kmlFile << "</Placemark>\n";
        }

        // Write the placemark data to the KML file
        kmlFile << "<Placemark>\n";
        kmlFile << "    <name>" << XMLString::transcode(element->getAttribute(XMLString::transcode("name"))) << "</name>\n";
        kmlFile << "    <description>" << XMLString::transcode(element->getAttribute(XMLString::transcode("description"))) << "</description>\n";

        if (element->getElementsByTagName(receiverTag)->getLength() > 0)
        {
            kmlFile << "    <styleUrl>#receiver</styleUrl>\n";
        }
        else if (element->getElementsByTagName(transmitterTag)->getLength() > 0)
        {
            kmlFile << "    <styleUrl>#transmitter</styleUrl>\n";
        }
        else if (element->getElementsByTagName(targetTag)->getLength() > 0)
        {
            kmlFile << "    <styleUrl>#target</styleUrl>\n";
        }

        // If the interpolation is linear, hyperbolic or exponential, use the gx:Track element
        if (isLinear || isHyperbolic || isCubic)
        {
            kmlFile << "    <gx:Track>\n";
            if (altitudeAboveGround > 0)
            {
                kmlFile << "        <altitudeMode>relativeToGround</altitudeMode>\n";
                kmlFile << "        <extrude>1</extrude>\n";
            }
            else
            {
                kmlFile << "        <altitudeMode>clampToGround</altitudeMode>\n";
            }

            // Iterate through the position waypoints
            for (XMLSize_t i = 0; i < positionWaypointList->getLength(); ++i)
            {
                const DOMElement *positionWaypointElement = dynamic_cast<const DOMElement *>(positionWaypointList->item(i));

                // Extract the position coordinates
                double x = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(xTag)->item(0)->getTextContent()));
                double y = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(yTag)->item(0)->getTextContent()));
                double altitude = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(altitudeTag)->item(0)->getTextContent()));

                // Convert the position coordinates to geographic coordinates
                double longitude = referenceLongitude + x / (cos(referenceLatitude * M_PI / 180) * 111319.9);
                double latitude = referenceLatitude + y / 111319.9;
                double altitudeAboveGround = altitude - referenceAltitude;

                // Extract the time value
                const XMLCh *timeTag = XMLString::transcode("time");
                double time = std::stod(XMLString::transcode(positionWaypointElement->getElementsByTagName(timeTag)->item(0)->getTextContent()));

                // Check if interpolation is hyperbolic
                if (isHyperbolic)
                {
                    // Calculate the hyperbolic path and update longitude and latitude values accordingly
                    double a = 0.5;                                                              // Set the desired value for 'a' based on the shape of the hyperbola
                    double b = 0.5;                                                              // Set the desired value for 'b' based on the shape of the hyperbola
                    double t = (double)i / (positionWaypointList->getLength() - 1) * 2.0 * M_PI; // Parameter 't' varies from 0 to 2 * PI
                    updateLongitudeLatitudeHyperbolic(longitude, latitude, t, a, b);
                }

                // Check if interpolation is cubic
                if (isCubic && i + 1 < positionWaypointList->getLength())
                {
                    // Calculate time difference between two consecutive position waypoints
                    const DOMElement *nextPositionWaypointElement = dynamic_cast<const DOMElement *>(positionWaypointList->item(i + 1));
                    double nextTime = std::stod(XMLString::transcode(nextPositionWaypointElement->getElementsByTagName(timeTag)->item(0)->getTextContent()));
                    double time_diff = nextTime - time;

                    // Extract the position coordinates for the next waypoint
                    double nextX = std::stod(XMLString::transcode(nextPositionWaypointElement->getElementsByTagName(xTag)->item(0)->getTextContent()));
                    double nextY = std::stod(XMLString::transcode(nextPositionWaypointElement->getElementsByTagName(yTag)->item(0)->getTextContent()));
                    double nextAltitude = std::stod(XMLString::transcode(nextPositionWaypointElement->getElementsByTagName(altitudeTag)->item(0)->getTextContent()));

                    // Convert the position coordinates to geographic coordinates
                    double nextLongitude = referenceLongitude + nextX / (cos(referenceLatitude * M_PI / 180) * 111319.9);
                    double nextLatitude = referenceLatitude + nextY / 111319.9;
                    double nextAltitudeAboveGround = nextAltitude - referenceAltitude;

                    // Calculate control points for cubic interpolation
                    double x1 = longitude;
                    double y1 = latitude;
                    double x4 = nextLongitude;
                    double y4 = nextLatitude;
                    double newX, newY;
                    updateLongitudeLatitudeCubic(newX, newY, 0.0, x1, y1, x4, y4); // Calculate first point on cubic curve

                    int num_divisions = 100;
                    for (int j = 0; j <= num_divisions; ++j)
                    {
                        double t = (double)j / num_divisions;

                        double newLongitude, newLatitude;
                        updateLongitudeLatitudeCubic(newLongitude, newLatitude, t, x1, y1, x4, y4);

                        double newAltitudeAboveGround = altitudeAboveGround + t * (nextAltitudeAboveGround - altitudeAboveGround);

                        kmlFile << "        <when>" << time + (double)(j * time_diff) / num_divisions << "</when>\n";
                        kmlFile << "        <gx:coord>" << newLongitude << " " << newLatitude << " " << newAltitudeAboveGround << "</gx:coord>\n";
                    }
                }

                else
                {
                    // Write the time and coordinates to the gx:Track element
                    kmlFile << "        <when>" << time << "</when>\n";
                    kmlFile << "        <gx:coord>" << longitude << " " << latitude << " " << altitudeAboveGround << "</gx:coord>\n";
                }
            }

            kmlFile << "    </gx:Track>\n";
        }

        else
        {
            kmlFile << "    <LookAt>\n";
            kmlFile << "        <longitude>" << longitude << "</longitude>\n";
            kmlFile << "        <latitude>" << latitude << "</latitude>\n";
            kmlFile << "        <altitude>" << altitudeAboveGround << "</altitude>\n";
            kmlFile << "        <heading>-148.4122922628044</heading>\n";
            kmlFile << "        <tilt>40.5575073395506</tilt>\n";
            kmlFile << "        <range>500.6566641072245</range>\n";
            kmlFile << "    </LookAt>\n";

            kmlFile << "    <Point>\n";
            kmlFile << "        <coordinates>" << longitude << "," << latitude << "," << altitudeAboveGround << "</coordinates>\n";

            if (altitudeAboveGround > 0)
            {
                kmlFile << "        <altitudeMode>relativeToGround</altitudeMode>\n";
                kmlFile << "        <extrude>1</extrude>\n";
            }
            else
            {
                kmlFile << "        <altitudeMode>clampToGround</altitudeMode>\n";
            }

            kmlFile << "    </Point>\n";
        }

        kmlFile << "</Placemark>\n";

        if (isLinear || isHyperbolic || isCubic)
        {
            // Get the first and last position waypoints
            const DOMElement *firstPositionWaypointElement = dynamic_cast<const DOMElement *>(positionWaypointList->item(0));
            const DOMElement *lastPositionWaypointElement = dynamic_cast<const DOMElement *>(positionWaypointList->item(positionWaypointList->getLength() - 1));

            // Extract the start and end coordinates
            std::string startCoordinates = getCoordinatesFromPositionWaypoint(firstPositionWaypointElement, referenceLatitude, referenceLongitude, referenceAltitude);
            std::string endCoordinates = getCoordinatesFromPositionWaypoint(lastPositionWaypointElement, referenceLatitude, referenceLongitude, referenceAltitude);

            // Start point placemark
            kmlFile << "<Placemark>\n";
            kmlFile << "    <name>Start: " << XMLString::transcode(element->getAttribute(XMLString::transcode("name"))) << "</name>\n";
            kmlFile << "    <styleUrl>#target</styleUrl>\n"; // Replace with your desired style URL for the start icon
            kmlFile << "    <Point>\n";
            kmlFile << "        <coordinates>" << startCoordinates << "</coordinates>\n";
            if (altitudeAboveGround > 0)
            {
                kmlFile << "        <altitudeMode>relativeToGround</altitudeMode>\n";
                kmlFile << "        <extrude>1</extrude>\n";
            }
            else
            {
                kmlFile << "        <altitudeMode>clampToGround</altitudeMode>\n";
            }
            kmlFile << "    </Point>\n";
            kmlFile << "</Placemark>\n";

            // End point placemark
            kmlFile << "<Placemark>\n";
            kmlFile << "    <name>End: " << XMLString::transcode(element->getAttribute(XMLString::transcode("name"))) << "</name>\n";
            kmlFile << "    <styleUrl>#target</styleUrl>\n"; // Replace with your desired style URL for the end icon
            kmlFile << "    <Point>\n";
            kmlFile << "        <coordinates>" << endCoordinates << "</coordinates>\n";
            if (altitudeAboveGround > 0)
            {
                kmlFile << "        <altitudeMode>relativeToGround</altitudeMode>\n";
                kmlFile << "        <extrude>1</extrude>\n";
            }
            else
            {
                kmlFile << "        <altitudeMode>clampToGround</altitudeMode>\n";
            }
            kmlFile << "    </Point>\n";
            kmlFile << "</Placemark>\n";
        }
    }
}

// Function to traverse the DOMNode by recursively calling itself and processElement()
void traverseDOMNode(const DOMNode *node, std::ofstream &kmlFile, double referenceLatitude, double referenceLongitude, double referenceAltitude, DOMDocument *document)
{
    if (node->getNodeType() == DOMNode::ELEMENT_NODE)
    {
        const DOMElement *element = dynamic_cast<const DOMElement *>(node);
        processElement(element, kmlFile, referenceLatitude, referenceLongitude, referenceAltitude, document);
    }

    for (DOMNode *child = node->getFirstChild(); child != nullptr; child = child->getNextSibling())
    {
        traverseDOMNode(child, kmlFile, referenceLatitude, referenceLongitude, referenceAltitude, document);
    }
}

// Main function
int main(int argc, char *argv[])
{

    // double alpha = 1;
    // double beta = 2;
    // double gamma = 3.6;

    // double angle_3dB_drop_deg = find_3db_drop_angle(alpha, beta, gamma);
    // std::cout << "3dB drop occurs at: " << angle_3dB_drop_deg << " degrees" << std::endl;

    if (argc > 3 && argc < 6)
    {
        std::cerr << "Usage: " << argv[0] << " <input XML file> <output KML file> [<referenceLatitude> <referenceLongitude> <referenceAltitude>]" << std::endl;
        return 1;
    }

    // Setting default geographical and altitude coordinates
    double referenceLatitude = -33.9545;
    double referenceLongitude = 18.4563;
    double referenceAltitude = 0;

    // Update file_path with command line argument
    string file_path = argv[1];
    // Setting mode to evironment variable from command line
    string output_file = argv[2];
    // Setting georgraphical coordinates to command line input
    if (argc == 6)
    {
        try
        {
            referenceLatitude = std::stod(argv[3]);
            referenceLongitude = std::stod(argv[4]);
            referenceAltitude = std::stod(argv[5]);
        }
        catch (const std::invalid_argument &e)
        {
            std::cerr << "Error: Invalid argument. Please provide valid numbers for referenceLatitude, referenceLongitude, and referenceAltitude.\n";
            return 1;
        }
        catch (const std::out_of_range &e)
        {
            std::cerr << "Error: Out of range. Please provide valid numbers for referenceLatitude, referenceLongitude, and referenceAltitude.\n";
            return 1;
        }
    }

    try
    {

        // Initializing Xerces-C++ library:
        XMLPlatformUtils::Initialize();

        // A XercesDOMParser object is set along with its features:
        XercesDOMParser parser;

        // Error handler configuration
        ErrorHandler *errorHandler = (ErrorHandler *)new HandlerBase();
        parser.setErrorHandler(errorHandler);

        // Disables validation during parsing.
        parser.setValidationScheme(XercesDOMParser::Val_Never);

        // Namespace set to false
        parser.setDoNamespaces(false);

        // Validation against schema set to false
        parser.setDoSchema(false);

        parser.setLoadExternalDTD(false);

        // Use file_path from command line argument
        parser.parse(file_path.c_str());

        // Creating DOMDocument and checking if document pointer is valid
        DOMDocument *document = parser.getDocument();
        if (!document)
        {
            std::cerr << "Error: document not found" << std::endl;
            XMLPlatformUtils::Terminate();
            return 1;
        }

        // Creating rootElement and checking if rootElement pointer is valid
        DOMElement *rootElement = document->getDocumentElement();
        if (!rootElement)
        {
            std::cerr << "Error: root element not found" << std::endl;
            XMLPlatformUtils::Terminate();
            return 1;
        }

        std::ofstream kmlFile(output_file.c_str());
        if (!kmlFile.is_open())
        {
            std::cerr << "Error opening output KML file" << std::endl;
            XMLPlatformUtils::Terminate();
            return 1;
        }

        // Write the KML header
        kmlFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        kmlFile << "<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\">\n";
        kmlFile << "<Document>\n";
        kmlFile << "<name>" << file_path << "</name>\n";

        // KML styles appended to document
        kmlFile << "<Style id=\"receiver\">\n";
        kmlFile << "  <IconStyle>\n";
        kmlFile << "    <Icon>\n";
        kmlFile << "      <href>https://cdn-icons-png.flaticon.com/512/645/645436.png</href>\n";
        kmlFile << "    </Icon>\n";
        kmlFile << "  </IconStyle>\n";
        kmlFile << "</Style>\n";
        kmlFile << "<Style id=\"transmitter\">\n";
        kmlFile << "  <IconStyle>\n";
        kmlFile << "    <Icon>\n";
        kmlFile << "      <href>https://cdn-icons-png.flaticon.com/128/224/224666.png</href>\n";
        kmlFile << "    </Icon>\n";
        kmlFile << "  </IconStyle>\n";
        kmlFile << "</Style>\n";
        kmlFile << "<Style id=\"target\">\n";
        kmlFile << "  <IconStyle>\n";
        kmlFile << "    <Icon>\n";
        kmlFile << "      <href>https://upload.wikimedia.org/wikipedia/commons/thumb/a/ad/Target_red_dot1.svg/1200px-Target_red_dot1.svg.png</href>\n";
        kmlFile << "    </Icon>\n";
        kmlFile << "  </IconStyle>\n";
        kmlFile << "  <LineStyle>\n";
        kmlFile << "    <width>2</width>\n";
        kmlFile << "  </LineStyle>\n";
        kmlFile << "</Style>\n";
        kmlFile << "<Style id=\"translucentPolygon\">\n";
        kmlFile << "    <LineStyle>\n";
        kmlFile << "        <color>ff0000ff</color>\n";
        kmlFile << "        <width>2</width>\n";
        kmlFile << "    </LineStyle>\n";
        kmlFile << "    <PolyStyle>\n";
        kmlFile << "        <color>00ffffff</color> <!-- RGBA: 50% transparent white --> \n";
        kmlFile << "     </PolyStyle>\n";
        kmlFile << "</Style>\n";
        kmlFile << "<Style id=\"arrowStyle\">\n";
        kmlFile << "    <IconStyle>\n";
        kmlFile << "        <Icon>\n";
        kmlFile << "            <href>http://maps.google.com/mapfiles/kml/shapes/arrow.png</href>\n";
        kmlFile << "        </Icon>\n";
        kmlFile << "        <scale>0.5</scale>\n";
        kmlFile << "    </IconStyle>\n";
        kmlFile << "</Style>\n";
        kmlFile << "<Style id=\"lineStyle\">\n";
        kmlFile << "    <LineStyle>\n";
        kmlFile << "        <color>ff0000ff</color>\n";
        kmlFile << "        <width>2</width>\n";
        kmlFile << "     </LineStyle>\n";
        kmlFile << "</Style>\n";
        kmlFile << "<Style id=\"lineStyleBlue\">\n";
        kmlFile << "    <LineStyle>\n";
        kmlFile << "        <color>ffff0000</color>\n";
        kmlFile << "        <width>2</width>\n";
        kmlFile << "     </LineStyle>\n";
        kmlFile << "</Style>\n";

        // Folder element appended
        kmlFile << "<Folder>\n";
        kmlFile << "  <name>Reference Coordinate</name>\n";
        kmlFile << "  <description>Placemarks for various elements in the FERSXML file. All Placemarks are situated relative to this reference point.</description>\n";

        // Add the LookAt element with given values
        kmlFile << "  <LookAt>\n";
        kmlFile << "    <longitude>" << referenceLongitude << "</longitude>\n";
        kmlFile << "    <latitude>" << referenceLatitude << "</latitude>\n";
        kmlFile << "    <altitude>" << referenceAltitude << "</altitude>\n";
        kmlFile << "    <heading>-148.4122922628044</heading>\n";
        kmlFile << "    <tilt>40.5575073395506</tilt>\n";
        kmlFile << "    <range>10000</range>\n";
        kmlFile << "  </LookAt>\n";

        // Traverse DOMNode and output extracted FERSXML data:
        traverseDOMNode(rootElement, kmlFile, referenceLatitude, referenceLongitude, referenceAltitude, document);

        // Close the Folder and Document elements
        kmlFile << "</Folder>\n";
        kmlFile << "</Document>\n";
        kmlFile << "</kml>\n";

        kmlFile.close();

        delete errorHandler; // Clean up the error handler
    }
    catch (const XMLException &e)
    {
        cerr << "Error initializing Xerces-C++: " << XMLString::transcode(e.getMessage()) << endl;
        return 1;
    }
    catch (const DOMException &e)
    {
        cerr << "Error parsing XML: " << XMLString::transcode(e.getMessage()) << endl;
        return 1;
    }
    catch (const SAXException &e)
    {
        cerr << "Error parsing XML: " << XMLString::transcode(e.getMessage()) << endl;
        return 1;
    }
    catch (...)
    {
        cerr << "Unknown error occurred while parsing XML." << endl;
        return 1;
    }

    XMLPlatformUtils::Terminate();
    return 0;
}