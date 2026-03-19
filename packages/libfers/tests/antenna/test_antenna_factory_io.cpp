#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <highfive/highfive.hpp>
#include <string>
#include <system_error>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/sim_id.h"
#include "math/geometry_ops.h"

using Catch::Matchers::WithinAbs;

namespace
{
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

	antenna::XmlAntenna antenna("xml", path.string(), 7);
	REQUIRE_THROWS_AS(antenna.getGain(unitDirection(0.0, 0.0), unitDirection(0.0, 0.0), 1.0), std::runtime_error);

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
