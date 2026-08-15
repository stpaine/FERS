#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <highfive/highfive.hpp>
#include <string>
#include <vector>

#include "core/parameters.h"
#include "core/world.h"
#include "serial/xml_parser.h"

using Catch::Matchers::ContainsSubstring;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		ParamGuard(const ParamGuard&) = delete;
		ParamGuard& operator=(const ParamGuard&) = delete;
		ParamGuard(ParamGuard&&) = delete;
		ParamGuard& operator=(ParamGuard&&) = delete;
		~ParamGuard() { params::params = saved; }
	};

	std::string getMinimalValidXml()
	{
		return R"(<?xml version="1.0" encoding="UTF-8"?>
		<simulation name="TestSim">
		  <parameters>
		    <starttime>0</starttime>
		    <endtime>1</endtime>
		    <rate>1000</rate>
		  </parameters>
		  <waveform name="w1">
		    <power>1000</power>
		    <carrier_frequency>1e9</carrier_frequency>
		    <cw/>
		  </waveform>
		  <timing name="t1">
		    <frequency>1e6</frequency>
		  </timing>
		  <antenna name="a1" pattern="isotropic"/>
		  <platform name="p1">
		    <motionpath interpolation="static">
		      <positionwaypoint><x>0</x><y>0</y><altitude>0</altitude><time>0</time></positionwaypoint>
		    </motionpath>
		    <rotationpath interpolation="static">
		      <rotationwaypoint><azimuth>0</azimuth><elevation>0</elevation><time>0</time></rotationwaypoint>
		    </rotationpath>
		    <transmitter name="tx1" waveform="w1" antenna="a1" timing="t1">
		      <cw_mode/>
		    </transmitter>
		  </platform>
		</simulation>)";
	}

	std::filesystem::path writeTempXml(const std::string& filename, const std::string& content)
	{
		auto path = std::filesystem::temp_directory_path() / filename;
		std::ofstream out(path);
		out << content;
		return path;
	}

	void writeHdf5Waveform(const std::filesystem::path& path)
	{
		HighFive::File file(path.string(), HighFive::File::Overwrite);
		const std::vector<RealType> i_values{1.0, 0.0, -1.0, 0.0};
		const std::vector<RealType> q_values{0.0, 1.0, 0.0, -1.0};
		file.createGroup("/I").createDataSet<RealType>("value", HighFive::DataSpace::From(i_values)).write(i_values);
		file.createGroup("/Q").createDataSet<RealType>("value", HighFive::DataSpace::From(q_values)).write(q_values);
	}
}

TEST_CASE("parseSimulationFromString successfully parses a valid scenario", "[serial][xml_parser]")
{
	ParamGuard const guard;
	core::World world;
	std::mt19937 seeder(42);

	REQUIRE_NOTHROW(serial::parseSimulationFromString(getMinimalValidXml(), &world, true, seeder));

	REQUIRE(params::params.simulation_name == "TestSim");
	REQUIRE(world.getPlatforms().size() == 1);
	REQUIRE(world.getTransmitters().size() == 1);
	REQUIRE(world.getWaveforms().size() == 1);
	REQUIRE(world.getAntennas().size() == 1);
	REQUIRE(world.getTimings().size() == 1);
}

TEST_CASE("parseSimulationFromString throws on malformed XML", "[serial][xml_parser]")
{
	ParamGuard const guard;
	core::World world;
	std::mt19937 seeder(42);

	std::string const bad_xml = "<simulation><parameters><starttime>0</starttime></parameters>"; // Missing closing tags

	REQUIRE_THROWS_WITH(serial::parseSimulationFromString(bad_xml, &world, false, seeder),
						ContainsSubstring("Failed to parse XML from memory string"));
}

TEST_CASE("parseSimulationFromString throws on schema validation failure", "[serial][xml_parser]")
{
	ParamGuard const guard;
	core::World world;
	std::mt19937 seeder(42);

	// Missing required <parameters> block
	std::string const invalid_schema_xml = R"(<?xml version="1.0" encoding="UTF-8"?>
		<simulation name="TestSim">
		  <waveform name="w1"><power>1000</power><carrier_frequency>1e9</carrier_frequency><cw/></waveform>
		</simulation>)";

	REQUIRE_THROWS_WITH(serial::parseSimulationFromString(invalid_schema_xml, &world, true, seeder),
						ContainsSubstring("XML failed DTD validation"));
}

TEST_CASE("parseSimulation handles file loading and includes", "[serial][xml_parser]")
{
	ParamGuard const guard;
	core::World world;
	std::mt19937 seeder(42);

	std::string const include_xml = R"(<?xml version="1.0" encoding="UTF-8"?>
		<simulation name="IncludeSim">
		  <antenna name="a1" pattern="isotropic"/>
		</simulation>)";

	std::string const main_xml = R"(<?xml version="1.0" encoding="UTF-8"?>
		<simulation name="MainSim">
		  <parameters>
		    <starttime>0</starttime>
		    <endtime>1</endtime>
		    <rate>1000</rate>
		  </parameters>
		  <include>test_include.xml</include>
		</simulation>)";

	auto include_path = writeTempXml("test_include.xml", include_xml);
	auto main_path = writeTempXml("test_main.xml", main_xml);

	// We disable validation here because the split files individually might not pass strict XSD/DTD
	// until they are merged, but the test ensures the merge logic works.
	REQUIRE_NOTHROW(serial::parseSimulation(main_path.string(), &world, false, seeder));

	REQUIRE(world.getAntennas().size() == 1); // Proves the include was merged

	std::filesystem::remove(include_path);
	std::filesystem::remove(main_path);
}

TEST_CASE("parseSimulation resolves relative CW and FMCW HDF5 waveform files",
		  "[serial][xml_parser][file-waveform][hdf5]")
{
	ParamGuard const guard;
	const auto directory = std::filesystem::temp_directory_path() /
		("fers_relative_file_waveforms_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
	std::filesystem::create_directories(directory);
	const auto hdf5_path = directory / "sampled_iq.h5";
	const auto scenario_path = directory / "scenario.fersxml";
	writeHdf5Waveform(hdf5_path);
	{
		std::ofstream scenario(scenario_path);
		scenario << R"(<?xml version="1.0" encoding="UTF-8"?>
<simulation name="RelativeFileWaveforms">
  <parameters><starttime>0</starttime><endtime>1</endtime><rate>8</rate></parameters>
  <waveform name="FileCw"><power>1</power><carrier_frequency>100</carrier_frequency><cw_from_file filename="sampled_iq.h5"/></waveform>
  <waveform name="FileFmcw"><power>1</power><carrier_frequency>100</carrier_frequency><fmcw_from_file filename="sampled_iq.h5"/></waveform>
</simulation>)";
	}

	core::World world;
	std::mt19937 seeder(42);
	REQUIRE_NOTHROW(serial::parseSimulation(scenario_path.string(), &world, true, seeder));
	REQUIRE(world.getWaveforms().size() == 2u);
	REQUIRE(world.findWaveformByName("FileCw")->isCw());
	REQUIRE(world.findWaveformByName("FileFmcw")->isFmcwFamily());
	REQUIRE(world.findWaveformByName("FileCw")->getFilename() == hdf5_path.string());
	REQUIRE(world.findWaveformByName("FileFmcw")->getFilename() == hdf5_path.string());
	std::error_code ec;
	std::filesystem::remove_all(directory, ec);
}

TEST_CASE("parseSimulation throws on missing file", "[serial][xml_parser]")
{
	ParamGuard const guard;
	core::World world;
	std::mt19937 seeder(42);

	REQUIRE_THROWS_WITH(serial::parseSimulation("non_existent_file_12345.xml", &world, false, seeder),
						ContainsSubstring("Failed to load main XML file"));
}
