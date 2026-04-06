#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <highfive/highfive.hpp>
#include <iomanip>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/sim_id.h"
#include "math/geometry_ops.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace
{
	using Sample = std::pair<RealType, RealType>;

	std::string uniqueFileName(const std::string& prefix, const std::string& extension)
	{
		return prefix + "_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + extension;
	}

	std::filesystem::path tempFilePath(const std::string& filename)
	{
		return std::filesystem::temp_directory_path() / filename;
	}

	void removeIfExists(const std::filesystem::path& path)
	{
		std::error_code ec;
		std::filesystem::remove(path, ec);
	}

	void writeXmlAntennaFile(const std::filesystem::path& path, const std::string& body)
	{
		std::ofstream out(path, std::ios::binary);
		REQUIRE(out.is_open());
		out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << body;
	}

	void writeHdf5PatternFile(const std::filesystem::path& path, const std::vector<std::vector<RealType>>& pattern)
	{
		HighFive::File file(path.string(), HighFive::File::Overwrite);
		file.createDataSet<RealType>("antenna", HighFive::DataSpace::From(pattern)).write(pattern);
	}

	math::SVec3 unitDirection(const RealType azimuth, const RealType elevation) { return {1.0, azimuth, elevation}; }

	std::string formatReal(const RealType value)
	{
		std::ostringstream out;
		out << std::setprecision(17) << value;
		return out.str();
	}

	RealType radToDeg(const RealType value) { return value * 180.0 / PI; }

	RealType linearToDbi(const RealType value) { return 10.0 * std::log10(value); }

	std::string makeAxisXml(const std::string_view axis_name, const std::vector<Sample>& samples,
							const std::string_view attributes = {})
	{
		std::ostringstream out;
		out << "<" << axis_name;
		if (!attributes.empty())
		{
			out << " " << attributes;
		}
		out << ">";
		for (const auto& [angle, gain] : samples)
		{
			out << "<gainsample><angle>" << formatReal(angle) << "</angle><gain>" << formatReal(gain)
				<< "</gain></gainsample>";
		}
		out << "</" << axis_name << ">";
		return out.str();
	}

	std::string xmlPatternFixture()
	{
		return R"(<antenna>
		<azimuth>
			<gainsample><angle>0.0</angle><gain>8.0</gain></gainsample>
			<gainsample><angle>1.0</angle><gain>4.0</gain></gainsample>
			<gainsample><angle>2.0</angle><gain>0.0</gain></gainsample>
		</azimuth>
		<elevation>
			<gainsample><angle>0.0</angle><gain>8.0</gain></gainsample>
			<gainsample><angle>1.0</angle><gain>0.0</gain></gainsample>
			<gainsample><angle>2.0</angle><gain>0.0</gain></gainsample>
		</elevation>
	</antenna>)";
	}

	std::string xmlPositiveLinearPatternFixture()
	{
		const std::vector<Sample> azimuth = {{0.0, 8.0}, {1.0, 4.0}, {2.0, 2.0}};
		const std::vector<Sample> elevation = {{0.0, 8.0}, {1.0, 2.0}, {2.0, 1.0}};
		return "<antenna>" + makeAxisXml("azimuth", azimuth) + makeAxisXml("elevation", elevation) + "</antenna>";
	}

	std::string xmlPatternFixtureDegreesDbi()
	{
		const std::vector<Sample> azimuth = {
			{radToDeg(0.0), linearToDbi(8.0)},
			{radToDeg(1.0), linearToDbi(4.0)},
			{radToDeg(2.0), linearToDbi(2.0)},
		};
		const std::vector<Sample> elevation = {
			{radToDeg(0.0), linearToDbi(8.0)},
			{radToDeg(1.0), linearToDbi(2.0)},
			{radToDeg(2.0), linearToDbi(1.0)},
		};
		return "<antenna>" + makeAxisXml("azimuth", azimuth, "unit=\"deg\" format=\"dBi\" symmetry=\"mirrored\"") +
			makeAxisXml("elevation", elevation, "unit=\"deg\" format=\"dBi\" symmetry=\"mirrored\"") + "</antenna>";
	}

	std::string xmlAsymmetricPatternFixture(const bool explicit_symmetry)
	{
		const std::string symmetry = explicit_symmetry ? " symmetry=\"full\"" : "";
		const std::vector<Sample> azimuth = {{-90.0, 10.0}, {0.0, 0.0}, {90.0, -10.0}};
		const std::vector<Sample> elevation = {{0.0, 0.0}, {90.0, 0.0}};
		return "<antenna>" + makeAxisXml("azimuth", azimuth, std::string("unit=\"deg\" format=\"dBi\"") + symmetry) +
			makeAxisXml("elevation", elevation, "unit=\"deg\" format=\"dBi\" symmetry=\"mirrored\"") + "</antenna>";
	}
}

TEST_CASE("XmlAntenna loads gain axes and metadata", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("xml_antenna", ".xml"));
	removeIfExists(path);
	writeXmlAntennaFile(path, xmlPatternFixture());

	constexpr SimId explicitId = 314;
	antenna::XmlAntenna antenna("xml", path.string(), explicitId);
	antenna.setEfficiencyFactor(0.5);

	REQUIRE(antenna.getId() == explicitId);
	REQUIRE(antenna.getName() == "xml");
	REQUIRE(antenna.getFilename() == path.string());
	REQUIRE(antenna.getAzimuthSamples() != nullptr);
	REQUIRE(antenna.getElevationSamples() != nullptr);
	REQUIRE_THAT(antenna.getMaxGain(), WithinAbs(8.0, 1e-12));
	REQUIRE_THAT(antenna.getGain(unitDirection(0.0, 0.0), unitDirection(0.0, 0.0), 1.0), WithinAbs(4.0, 1e-12));

	removeIfExists(path);
}

TEST_CASE("XmlAntenna interpolates normalized axis gains", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("xml_antenna_interp", ".xml"));
	removeIfExists(path);
	writeXmlAntennaFile(path, xmlPatternFixture());

	antenna::XmlAntenna antenna("xml", path.string(), 99);
	antenna.setEfficiencyFactor(0.5);

	const auto azimuth_value = antenna.getAzimuthSamples()->getValueAt(0.5);
	const auto elevation_value = antenna.getElevationSamples()->getValueAt(0.25);
	REQUIRE(azimuth_value.has_value());
	REQUIRE(elevation_value.has_value());
	REQUIRE_THAT(*azimuth_value, WithinAbs(0.75, 1e-12));
	REQUIRE_THAT(*elevation_value, WithinAbs(0.75, 1e-12));

	const RealType expected = 0.75 * 0.75 * 8.0 * 0.5;
	REQUIRE_THAT(antenna.getGain(unitDirection(-0.5, 0.25), unitDirection(0.0, 0.0), 1.0), WithinAbs(expected, 1e-12));
	REQUIRE_THAT(antenna.getGain(unitDirection(-1.0, 0.0), unitDirection(0.0, 0.0), 1.0),
				 WithinAbs(antenna.getGain(unitDirection(1.0, 0.0), unitDirection(0.0, 0.0), 1.0), 1e-12));

	removeIfExists(path);
}

TEST_CASE("XmlAntenna converts degree and dBi axes to internal units", "[antenna][io]")
{
	const std::filesystem::path linear_path = tempFilePath(uniqueFileName("xml_antenna_linear_ref", ".xml"));
	const std::filesystem::path converted_path = tempFilePath(uniqueFileName("xml_antenna_dbi_deg", ".xml"));
	removeIfExists(linear_path);
	removeIfExists(converted_path);
	writeXmlAntennaFile(linear_path, xmlPositiveLinearPatternFixture());
	writeXmlAntennaFile(converted_path, xmlPatternFixtureDegreesDbi());

	antenna::XmlAntenna linear("linear", linear_path.string(), 101);
	antenna::XmlAntenna converted("converted", converted_path.string(), 102);
	linear.setEfficiencyFactor(0.5);
	converted.setEfficiencyFactor(0.5);

	REQUIRE_THAT(converted.getMaxGain(), WithinAbs(linear.getMaxGain(), 1e-10));
	REQUIRE_THAT(converted.getGain(unitDirection(0.0, 0.0), unitDirection(0.0, 0.0), 1.0),
				 WithinAbs(linear.getGain(unitDirection(0.0, 0.0), unitDirection(0.0, 0.0), 1.0), 1e-10));
	REQUIRE_THAT(converted.getGain(unitDirection(-0.5, 0.25), unitDirection(0.0, 0.0), 1.0),
				 WithinAbs(linear.getGain(unitDirection(-0.5, 0.25), unitDirection(0.0, 0.0), 1.0), 1e-10));
	REQUIRE_THAT(converted.getGain(unitDirection(1.0, 1.0), unitDirection(0.0, 0.0), 1.0),
				 WithinAbs(linear.getGain(unitDirection(1.0, 1.0), unitDirection(0.0, 0.0), 1.0), 1e-10));

	removeIfExists(linear_path);
	removeIfExists(converted_path);
}

TEST_CASE("XmlAntenna supports explicit full-range asymmetric lookup", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("xml_antenna_asym_explicit", ".xml"));
	removeIfExists(path);
	writeXmlAntennaFile(path, xmlAsymmetricPatternFixture(true));

	antenna::XmlAntenna antenna("xml", path.string(), 15);

	const RealType negative_gain = antenna.getGain(unitDirection(-PI / 2.0, 0.0), unitDirection(0.0, 0.0), 1.0);
	const RealType positive_gain = antenna.getGain(unitDirection(PI / 2.0, 0.0), unitDirection(0.0, 0.0), 1.0);

	REQUIRE(negative_gain > positive_gain);

	removeIfExists(path);
}

TEST_CASE("XmlAntenna auto-detects full-range lookup from negative sample angles", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("xml_antenna_asym_auto", ".xml"));
	removeIfExists(path);
	writeXmlAntennaFile(path, xmlAsymmetricPatternFixture(false));

	antenna::XmlAntenna antenna("xml", path.string(), 16);

	const RealType negative_gain = antenna.getGain(unitDirection(-PI / 2.0, 0.0), unitDirection(0.0, 0.0), 1.0);
	const RealType positive_gain = antenna.getGain(unitDirection(PI / 2.0, 0.0), unitDirection(0.0, 0.0), 1.0);

	REQUIRE(negative_gain > positive_gain);

	removeIfExists(path);
}

TEST_CASE("XmlAntenna throws for missing antenna description file", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("xml_antenna_missing", ".xml"));
	removeIfExists(path);

	REQUIRE_THROWS_AS(antenna::XmlAntenna("xml", path.string(), 1), std::runtime_error);
}

TEST_CASE("XmlAntenna throws for malformed antenna description file", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("xml_antenna_malformed", ".xml"));
	removeIfExists(path);

	{
		std::ofstream out(path, std::ios::binary);
		REQUIRE(out.is_open());
		out << "<antenna><azimuth><gainsample></antenna>";
	}

	REQUIRE_THROWS_AS(antenna::XmlAntenna("xml", path.string(), 1), std::runtime_error);

	removeIfExists(path);
}

TEST_CASE("XmlAntenna throws when gain samples are unavailable", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("xml_antenna_empty_axis", ".xml"));
	removeIfExists(path);
	writeXmlAntennaFile(path,
						R"(<antenna>
		<azimuth>
			<gainsample><angle>0.0</angle><gain>6.0</gain></gainsample>
			<gainsample><angle>1.0</angle><gain>3.0</gain></gainsample>
		</azimuth>
		<elevation></elevation>
	</antenna>)");

	REQUIRE_THROWS_AS(antenna::XmlAntenna("xml", path.string(), 7), std::runtime_error);

	removeIfExists(path);
}

TEST_CASE("XmlAntenna rejects non-positive peak gain during load", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("xml_antenna_zero_peak", ".xml"));
	removeIfExists(path);
	writeXmlAntennaFile(path,
						R"(<antenna>
		<azimuth>
			<gainsample><angle>0.0</angle><gain>0.0</gain></gainsample>
			<gainsample><angle>1.0</angle><gain>0.0</gain></gainsample>
		</azimuth>
		<elevation>
			<gainsample><angle>0.0</angle><gain>0.0</gain></gainsample>
			<gainsample><angle>1.0</angle><gain>0.0</gain></gainsample>
		</elevation>
	</antenna>)");

	REQUIRE_THROWS_WITH(antenna::XmlAntenna("xml", path.string(), 8), ContainsSubstring("peak linear gain"));

	removeIfExists(path);
}

TEST_CASE("H5Antenna loads pattern metadata", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("h5_antenna", ".h5"));
	removeIfExists(path);

	const std::vector<std::vector<RealType>> pattern = {
		{1.0, 2.0, 3.0, 4.0},
		{5.0, 6.0, 7.0, 8.0},
		{9.0, 10.0, 11.0, 12.0},
		{13.0, 14.0, 15.0, 16.0},
	};
	writeHdf5PatternFile(path, pattern);

	antenna::H5Antenna antenna("h5", path.string(), 808);
	antenna.setEfficiencyFactor(0.5);

	REQUIRE(antenna.getId() == 808);
	REQUIRE(antenna.getName() == "h5");
	REQUIRE(antenna.getFilename() == path.string());
	REQUIRE(antenna.getPattern().size() == 4u);
	REQUIRE(antenna.getPattern()[0].size() == 4u);
	REQUIRE_THAT(antenna.getPattern()[0][0], WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(antenna.getPattern()[3][3], WithinAbs(16.0, 1e-12));
	const RealType exact_grid_azimuth = -PI / 3.0;
	const RealType exact_grid_elevation = -PI / 3.0;
	REQUIRE_THAT(antenna.getGain(unitDirection(exact_grid_azimuth, exact_grid_elevation), unitDirection(0.0, 0.0), 1.0),
				 WithinAbs(pattern[1][1] * 0.5, 1e-12));

	removeIfExists(path);
}

TEST_CASE("H5Antenna interpolates bilinearly across the pattern grid", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("h5_antenna_interp", ".h5"));
	removeIfExists(path);

	const std::vector<std::vector<RealType>> pattern = {
		{1.0, 2.0, 3.0, 4.0},
		{5.0, 6.0, 7.0, 8.0},
		{9.0, 10.0, 11.0, 12.0},
		{13.0, 14.0, 15.0, 16.0},
	};
	writeHdf5PatternFile(path, pattern);

	antenna::H5Antenna antenna("h5", path.string(), 9);
	antenna.setEfficiencyFactor(0.5);

	const RealType p00 = pattern[0][0];
	const RealType p10 = pattern[1][0];
	const RealType p11 = pattern[1][1];
	const RealType p01 = pattern[0][1];
	const RealType expected = 0.25 * (p00 + p10 + p11 + p01) * 0.5;

	REQUIRE_THAT(antenna.getGain(unitDirection(-3.0 * PI / 4.0, -3.0 * PI / 4.0), unitDirection(0.0, 0.0), 1.0),
				 WithinAbs(expected, 1e-12));

	removeIfExists(path);
}

TEST_CASE("H5Antenna preserves explicit antenna id", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("h5_antenna_id", ".h5"));
	removeIfExists(path);
	writeHdf5PatternFile(path, {{2.0, 4.0}, {6.0, 8.0}});

	constexpr SimId explicitId = (static_cast<SimId>(0xABCD) << 32) | 0x1234;
	antenna::H5Antenna antenna("h5", path.string(), explicitId);

	REQUIRE(antenna.getId() == explicitId);

	removeIfExists(path);
}

TEST_CASE("H5Antenna throws for missing pattern file", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("h5_antenna_missing", ".h5"));
	removeIfExists(path);

	REQUIRE_THROWS_AS(antenna::H5Antenna("h5", path.string(), 1), std::runtime_error);
}

TEST_CASE("H5Antenna throws when antenna dataset is missing", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("h5_antenna_missing_dataset", ".h5"));
	removeIfExists(path);

	{
		HighFive::File file(path.string(), HighFive::File::Overwrite);
		const std::vector<std::vector<RealType>> pattern = {{1.0, 2.0}, {3.0, 4.0}};
		file.createDataSet<RealType>("other", HighFive::DataSpace::From(pattern)).write(pattern);
	}

	REQUIRE_THROWS_AS(antenna::H5Antenna("h5", path.string(), 1), std::runtime_error);

	removeIfExists(path);
}

TEST_CASE("H5Antenna throws when antenna dataset is not two-dimensional", "[antenna][io]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("h5_antenna_bad_dims", ".h5"));
	removeIfExists(path);

	{
		HighFive::File file(path.string(), HighFive::File::Overwrite);
		const std::vector<RealType> values = {1.0, 2.0, 3.0};
		file.createDataSet<RealType>("antenna", HighFive::DataSpace::From(values)).write(values);
	}

	REQUIRE_THROWS_AS(antenna::H5Antenna("h5", path.string(), 1), std::runtime_error);

	removeIfExists(path);
}

// TODO: Add a stable wraparound-boundary interpolation test for H5Antenna once
// the upper-edge interpolation path avoids the current x2 == x1 / y2 == y1 case
// at exact +pi boundaries, which makes boundary assertions non-deterministic.
