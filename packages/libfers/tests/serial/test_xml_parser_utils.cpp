#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>

#include "antenna/antenna_factory.h"
#include "core/logging.h"
#include "core/parameters.h"
#include "core/world.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/libxml_wrapper.h"
#include "serial/xml_parser_utils.h"
#include "signal/radar_signal.h"
#include "timing/prototype_timing.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
	};

	struct CerrCapture
	{
		std::ostringstream buffer;
		std::streambuf* old{nullptr};
		CerrCapture() { old = std::cerr.rdbuf(buffer.rdbuf()); }
		~CerrCapture() { std::cerr.rdbuf(old); }
		[[nodiscard]] std::string str() const { return buffer.str(); }
	};

	struct LogLevelGuard
	{
		explicit LogLevelGuard(const logging::Level level) { logging::logger.setLevel(level); }
		~LogLevelGuard() { logging::logger.setLevel(logging::Level::INFO); }
	};

	XmlDocument loadXml(const std::string& xmlString)
	{
		XmlDocument doc;
		REQUIRE(doc.loadString("<?xml version=\"1.0\" encoding=\"UTF-8\"?>" + xmlString));
		return doc;
	}

	serial::xml_parser_utils::AssetLoaders createMockLoaders()
	{
		serial::xml_parser_utils::AssetLoaders loaders;
		loaders.loadWaveform =
			[](const std::string& name, const std::filesystem::path&, RealType power, RealType carrierFreq, SimId id)
		{
			auto sig = std::make_unique<fers_signal::CwSignal>();
			return std::make_unique<fers_signal::RadarSignal>(name, power, carrierFreq, 1.0, std::move(sig), id);
		};
		loaders.loadXmlAntenna = [](const std::string& name, const std::string&, SimId id)
		{ return std::make_unique<antenna::Isotropic>(name, id); };
		loaders.loadH5Antenna = [](const std::string& name, const std::string&, SimId id)
		{ return std::make_unique<antenna::Isotropic>(name, id); };
		loaders.loadFileTarget =
			[](radar::Platform* platform, const std::string& name, const std::string&, unsigned seed, SimId id)
		{ return radar::createIsoTarget(platform, name, 1.0, seed, id); };
		return loaders;
	}
}

TEST_CASE("get_child_real_type extracts floating point values", "[serial][xml_parser_utils]")
{
	auto doc = loadXml("<root><val>3.14159</val></root>");
	REQUIRE_THAT(serial::xml_parser_utils::get_child_real_type(doc.getRootElement(), "val"), WithinAbs(3.14159, 1e-5));

	auto empty_doc = loadXml("<root><val></val></root>");
	REQUIRE_THROWS_AS(serial::xml_parser_utils::get_child_real_type(empty_doc.getRootElement(), "val"), XmlException);

	auto missing_doc = loadXml("<root></root>");
	REQUIRE_THROWS_AS(serial::xml_parser_utils::get_child_real_type(missing_doc.getRootElement(), "val"), XmlException);
}

TEST_CASE("get_attribute_bool extracts boolean values safely", "[serial][xml_parser_utils]")
{
	LogLevelGuard log_level(logging::Level::WARNING);
	CerrCapture capture;
	auto doc = loadXml("<root flag_true=\"true\" flag_false=\"false\" flag_invalid=\"yes\"/>");
	auto root = doc.getRootElement();

	REQUIRE(serial::xml_parser_utils::get_attribute_bool(root, "flag_true", false) == true);
	REQUIRE(serial::xml_parser_utils::get_attribute_bool(root, "flag_false", true) == false);
	REQUIRE(serial::xml_parser_utils::get_attribute_bool(root, "flag_invalid", true) == true);
	REQUIRE(serial::xml_parser_utils::get_attribute_bool(root, "missing", true) == true); // Default fallback
	REQUIRE_THAT(capture.str(), ContainsSubstring("Invalid boolean value"));
}

TEST_CASE("resolve_reference_id maps string names to SimIds", "[serial][xml_parser_utils]")
{
	auto doc = loadXml("<root ref=\"target_a\"/>");
	std::unordered_map<std::string, SimId> map = {{"target_a", 42}, {"target_b", 99}};

	REQUIRE(serial::xml_parser_utils::resolve_reference_id(doc.getRootElement(), "ref", "owner", map) == 42);
	REQUIRE_THROWS_AS(serial::xml_parser_utils::resolve_reference_id(doc.getRootElement(), "missing", "owner", map),
					  XmlException);

	auto bad_doc = loadXml("<root ref=\"unknown\"/>");
	REQUIRE_THROWS_AS(serial::xml_parser_utils::resolve_reference_id(bad_doc.getRootElement(), "ref", "owner", map),
					  XmlException);
}

TEST_CASE("parseSchedule handles valid and invalid periods", "[serial][xml_parser_utils]")
{
	ParamGuard guard;
	params::setTime(0.0, 10.0);

	auto doc = loadXml("<parent>"
					   "  <schedule>"
					   "    <period start=\"1.0\" end=\"2.0\"/>"
					   "    <period start=\"bad\" end=\"3.0\"/>"
					   "  </schedule>"
					   "</parent>");

	// The bad period should be caught and logged, but the valid one should be parsed
	auto periods = serial::xml_parser_utils::parseSchedule(doc.getRootElement(), "test_owner", false);
	REQUIRE(periods.size() == 1);
	REQUIRE_THAT(periods[0].start, WithinAbs(1.0, 1e-5));
	REQUIRE_THAT(periods[0].end, WithinAbs(2.0, 1e-5));
}

TEST_CASE("parseParameters extracts simulation parameters", "[serial][xml_parser_utils]")
{
	ParamGuard guard;

	SECTION("Full parameters with UTM South")
	{
		auto doc = loadXml("<parameters>"
						   "  <starttime>1.5</starttime>"
						   "  <endtime>10.0</endtime>"
						   "  <rate>2000</rate>"
						   "  <c>3e8</c>"
						   "  <adc_bits>12</adc_bits>"
						   "  <oversample>4</oversample>"
						   "  <origin latitude=\"-33.0\" longitude=\"18.0\" altitude=\"100.0\"/>"
						   "  <coordinatesystem frame=\"UTM\" zone=\"34\" hemisphere=\"S\"/>"
						   "</parameters>");

		params::Parameters p;
		serial::xml_parser_utils::parseParameters(doc.getRootElement(), p);

		REQUIRE_THAT(p.start, WithinAbs(1.5, 1e-5));
		REQUIRE_THAT(p.end, WithinAbs(10.0, 1e-5));
		REQUIRE_THAT(p.rate, WithinAbs(2000.0, 1e-5));
		REQUIRE_THAT(p.c, WithinAbs(3e8, 1e-5));
		REQUIRE(p.adc_bits == 12);
		REQUIRE(p.oversample_ratio == 4);
		REQUIRE_THAT(p.origin_latitude, WithinAbs(-33.0, 1e-5));
		REQUIRE_THAT(p.origin_longitude, WithinAbs(18.0, 1e-5));
		REQUIRE_THAT(p.origin_altitude, WithinAbs(100.0, 1e-5));
		REQUIRE(p.coordinate_frame == params::CoordinateFrame::UTM);
		REQUIRE(p.utm_zone == 34);
		REQUIRE(p.utm_north_hemisphere == false);
	}

	SECTION("Optional parameters (simSamplingRate, randomseed) and ECEF")
	{
		auto doc = loadXml("<parameters>"
						   "  <starttime>0</starttime><endtime>1</endtime><rate>1000</rate>"
						   "  <simSamplingRate>500.0</simSamplingRate>"
						   "  <randomseed>42</randomseed>"
						   "  <rotationangleunit>rad</rotationangleunit>"
						   "  <coordinatesystem frame=\"ECEF\"/>"
						   "</parameters>");

		params::Parameters p;
		serial::xml_parser_utils::parseParameters(doc.getRootElement(), p);

		REQUIRE_THAT(p.sim_sampling_rate, WithinAbs(500.0, 1e-5));
		REQUIRE(p.random_seed.has_value());
		REQUIRE(p.random_seed.value() == 42);
		REQUIRE(p.rotation_angle_unit == params::RotationAngleUnit::Radians);
		REQUIRE(p.coordinate_frame == params::CoordinateFrame::ECEF);
	}

	SECTION("Unsigned optional parameters floor positive fractional values")
	{
		auto doc = loadXml("<parameters>"
						   "  <starttime>0</starttime><endtime>1</endtime><rate>1000</rate>"
						   "  <randomseed>42.9</randomseed>"
						   "  <adc_bits>12.8</adc_bits>"
						   "  <oversample>4.2</oversample>"
						   "</parameters>");

		params::Parameters p;
		serial::xml_parser_utils::parseParameters(doc.getRootElement(), p);

		REQUIRE(p.random_seed.has_value());
		REQUIRE(p.random_seed.value() == 42);
		REQUIRE(p.adc_bits == 12);
		REQUIRE(p.oversample_ratio == 4);
	}

	SECTION("Missing optional parameters use defaults without warnings")
	{
		LogLevelGuard log_level(logging::Level::WARNING);
		CerrCapture capture;
		auto doc = loadXml("<parameters>"
						   "  <starttime>0</starttime><endtime>1</endtime><rate>1000</rate>"
						   "</parameters>");

		params::Parameters p;
		serial::xml_parser_utils::parseParameters(doc.getRootElement(), p);

		REQUIRE_THAT(p.c, WithinAbs(params::Parameters::DEFAULT_C, 1e-5));
		REQUIRE_THAT(p.sim_sampling_rate, WithinAbs(1000.0, 1e-5));
		REQUIRE(p.adc_bits == 0);
		REQUIRE(p.oversample_ratio == 1);
		REQUIRE(capture.str().empty());
	}

	SECTION("Origin altitude defaults to zero when omitted")
	{
		LogLevelGuard log_level(logging::Level::WARNING);
		CerrCapture capture;
		auto doc = loadXml("<parameters>"
						   "  <starttime>0</starttime><endtime>1</endtime><rate>1000</rate>"
						   "  <origin latitude=\"-33.0\" longitude=\"18.0\"/>"
						   "  <coordinatesystem frame=\"ENU\"/>"
						   "</parameters>");

		params::Parameters p;
		serial::xml_parser_utils::parseParameters(doc.getRootElement(), p);

		REQUIRE_THAT(p.origin_latitude, WithinAbs(-33.0, 1e-5));
		REQUIRE_THAT(p.origin_longitude, WithinAbs(18.0, 1e-5));
		REQUIRE_THAT(p.origin_altitude, WithinAbs(0.0, 1e-5));
		REQUIRE(capture.str().empty());
	}

	SECTION("Unsigned optional parameters reject invalid values")
	{
		const auto parse_invalid = [](const std::string& xml)
		{
			params::Parameters p;
			auto doc = loadXml(xml);
			serial::xml_parser_utils::parseParameters(doc.getRootElement(), p);
		};
		const auto too_large_unsigned =
			std::to_string(static_cast<unsigned long long>(std::numeric_limits<unsigned>::max()) + 1ULL);

		REQUIRE_THROWS_AS(parse_invalid("<parameters>"
										"  <starttime>0</starttime><endtime>1</endtime><rate>1000</rate>"
										"  <randomseed>-1</randomseed>"
										"</parameters>"),
						  XmlException);

		REQUIRE_THROWS_AS(parse_invalid("<parameters>"
										"  <starttime>0</starttime><endtime>1</endtime><rate>1000</rate>"
										"  <adc_bits>nan</adc_bits>"
										"</parameters>"),
						  XmlException);

		REQUIRE_THROWS_AS(parse_invalid("<parameters>"
										"  <starttime>0</starttime><endtime>1</endtime><rate>1000</rate>"
										"  <oversample>0</oversample>"
										"</parameters>"),
						  std::runtime_error);

		REQUIRE_THROWS_WITH(parse_invalid("<parameters>"
										  "  <starttime>0</starttime><endtime>1</endtime><rate>1000</rate>"
										  "  <oversample>9</oversample>"
										  "</parameters>"),
							ContainsSubstring("Oversampling ratios > 8 are not supported"));

		REQUIRE_THROWS_AS(parse_invalid(std::string("<parameters>") +
										"  <starttime>0</starttime><endtime>1</endtime><rate>1000</rate>" +
										"  <adc_bits>" + too_large_unsigned + "</adc_bits>" + "</parameters>"),
						  XmlException);
	}

	SECTION("UTM North Hemisphere")
	{
		auto doc = loadXml("<parameters>"
						   "  <starttime>0</starttime><endtime>1</endtime><rate>1000</rate>"
						   "  <coordinatesystem frame=\"UTM\" zone=\"34\" hemisphere=\"N\"/>"
						   "</parameters>");

		params::Parameters p;
		serial::xml_parser_utils::parseParameters(doc.getRootElement(), p);
		REQUIRE(p.coordinate_frame == params::CoordinateFrame::UTM);
		REQUIRE(p.utm_north_hemisphere == true);
	}

	SECTION("ENU without origin logs warning but succeeds")
	{
		auto doc = loadXml("<parameters>"
						   "  <starttime>0</starttime><endtime>1</endtime><rate>1000</rate>"
						   "  <coordinatesystem frame=\"ENU\"/>"
						   "</parameters>");

		params::Parameters p;
		serial::xml_parser_utils::parseParameters(doc.getRootElement(), p);
		REQUIRE(p.coordinate_frame == params::CoordinateFrame::ENU);
	}
}

TEST_CASE("parseParameters throws on invalid UTM zones", "[serial][xml_parser_utils]")
{
	ParamGuard guard;
	auto doc = loadXml("<parameters>"
					   "  <starttime>0</starttime>"
					   "  <endtime>1</endtime>"
					   "  <rate>1000</rate>"
					   "  <coordinatesystem frame=\"UTM\" zone=\"99\" hemisphere=\"N\"/>"
					   "</parameters>");
	params::Parameters p;

	// Should catch the exception internally and fallback to ENU
	serial::xml_parser_utils::parseParameters(doc.getRootElement(), p);
	REQUIRE(p.coordinate_frame == params::CoordinateFrame::ENU);
}

TEST_CASE("parseWaveform handles CW and delegates file loading", "[serial][xml_parser_utils]")
{
	core::World world;
	std::mt19937 seeder(42);
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;
	ctx.master_seeder = &seeder;
	ctx.parameters.start = 0;
	ctx.parameters.end = 1;
	ctx.loaders = createMockLoaders();

	SECTION("CW Waveform")
	{
		auto doc = loadXml(
			"<waveform name=\"cw1\"><power>10</power><carrier_frequency>1e9</carrier_frequency><cw/></waveform>");
		serial::xml_parser_utils::parseWaveform(doc.getRootElement(), ctx);
		REQUIRE(world.getWaveforms().size() == 1);
		REQUIRE(world.getWaveforms().begin()->second->getName() == "cw1");
	}

	SECTION("Pulsed from file")
	{
		auto doc =
			loadXml("<waveform name=\"p1\"><power>5</power><carrier_frequency>2e9</carrier_frequency><pulsed_from_file "
					"filename=\"dummy.csv\"/></waveform>");
		serial::xml_parser_utils::parseWaveform(doc.getRootElement(), ctx);
		REQUIRE(world.getWaveforms().size() == 1);
		REQUIRE(world.getWaveforms().begin()->second->getName() == "p1");
	}
}

TEST_CASE("parseWaveform warns for large FMCW streaming allocation", "[serial][xml_parser_utils]")
{
	ParamGuard guard;
	LogLevelGuard log_level(logging::Level::WARNING);
	CerrCapture capture;

	params::setRate(1.0e9);
	params::setOversampleRatio(1);
	params::setTime(0.0, 5.0);

	core::World world;
	std::mt19937 seeder(42);
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;
	ctx.master_seeder = &seeder;
	ctx.parameters.start = 0.0;
	ctx.parameters.end = 5.0;
	ctx.loaders = createMockLoaders();

	auto doc = loadXml("<waveform name=\"huge_fmcw\">"
					   "  <power>10</power>"
					   "  <carrier_frequency>2.4e9</carrier_frequency>"
					   "  <fmcw_up_chirp>"
					   "    <chirp_bandwidth>1e6</chirp_bandwidth>"
					   "    <chirp_duration>1e-3</chirp_duration>"
					   "    <chirp_period>1e-3</chirp_period>"
					   "  </fmcw_up_chirp>"
					   "</waveform>");

	serial::xml_parser_utils::parseWaveform(doc.getRootElement(), ctx);
	REQUIRE(world.getWaveforms().size() == 1);
	REQUIRE_THAT(capture.str(), ContainsSubstring("GiB of FMCW streaming IQ data"));
}

TEST_CASE("parseWaveform validates FMCW chirp schema constraints", "[serial][xml_parser_utils][fmcw]")
{
	ParamGuard guard;
	params::setRate(2.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 1.0);

	core::World world;
	std::mt19937 seeder(42);
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;
	ctx.master_seeder = &seeder;
	ctx.parameters.start = 0.0;
	ctx.parameters.end = 1.0;
	ctx.loaders = createMockLoaders();

	SECTION("well-formed FMCW waveform")
	{
		auto doc = loadXml("<waveform name=\"fmcw_good\">"
						   "  <power>10</power>"
						   "  <carrier_frequency>2.4e9</carrier_frequency>"
						   "  <fmcw_up_chirp>"
						   "    <chirp_bandwidth>1e6</chirp_bandwidth>"
						   "    <chirp_duration>1e-3</chirp_duration>"
						   "    <chirp_period>2e-3</chirp_period>"
						   "    <start_frequency_offset>-2.5e5</start_frequency_offset>"
						   "    <chirp_count>4</chirp_count>"
						   "  </fmcw_up_chirp>"
						   "</waveform>");

		serial::xml_parser_utils::parseWaveform(doc.getRootElement(), ctx);
		REQUIRE(world.getWaveforms().size() == 1);
		const auto* wave = world.getWaveforms().begin()->second.get();
		REQUIRE(wave->getFmcwChirpSignal() != nullptr);
		REQUIRE_THAT(wave->getLength(), WithinAbs(1.0e-3, 1.0e-12));
	}

	SECTION("chirp period shorter than duration is rejected")
	{
		auto doc = loadXml("<waveform name=\"fmcw_bad_period\">"
						   "  <power>10</power>"
						   "  <carrier_frequency>2.4e9</carrier_frequency>"
						   "  <fmcw_up_chirp>"
						   "    <chirp_bandwidth>1e6</chirp_bandwidth>"
						   "    <chirp_duration>2e-3</chirp_duration>"
						   "    <chirp_period>1e-3</chirp_period>"
						   "  </fmcw_up_chirp>"
						   "</waveform>");

		REQUIRE_THROWS_WITH(serial::xml_parser_utils::parseWaveform(doc.getRootElement(), ctx),
							ContainsSubstring("FMCW requires T_rep >= T_c"));
	}

	SECTION("complex-baseband aliasing constraint uses both sweep edges")
	{
		auto doc = loadXml("<waveform name=\"fmcw_bad_alias\">"
						   "  <power>10</power>"
						   "  <carrier_frequency>2.4e9</carrier_frequency>"
						   "  <fmcw_up_chirp>"
						   "    <chirp_bandwidth>5e5</chirp_bandwidth>"
						   "    <chirp_duration>1e-3</chirp_duration>"
						   "    <chirp_period>1e-3</chirp_period>"
						   "    <start_frequency_offset>-3e6</start_frequency_offset>"
						   "  </fmcw_up_chirp>"
						   "</waveform>");

		REQUIRE_THROWS_AS(serial::xml_parser_utils::parseWaveform(doc.getRootElement(), ctx), XmlException);
	}
}

TEST_CASE("parseTiming extracts clock parameters and noise entries", "[serial][xml_parser_utils]")
{
	core::World world;
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;

	auto doc = loadXml("<timing name=\"clk\" synconpulse=\"true\">"
					   "  <frequency>1e6</frequency>"
					   "  <freq_offset>10.0</freq_offset>"
					   "  <phase_offset>3.14</phase_offset>"
					   "  <noise_entry><alpha>1.0</alpha><weight>0.5</weight></noise_entry>"
					   "  <noise_entry><alpha>2.0</alpha><weight>0.25</weight></noise_entry>"
					   "</timing>");

	serial::xml_parser_utils::parseTiming(doc.getRootElement(), ctx);
	REQUIRE(world.getTimings().size() == 1);

	auto* timing = world.getTimings().begin()->second.get();
	REQUIRE(timing->getName() == "clk");
	REQUIRE_THAT(timing->getFrequency(), WithinAbs(1e6, 1e-5));

	REQUIRE(timing->getFreqOffset().has_value());
	REQUIRE_THAT(timing->getFreqOffset().value(), WithinAbs(10.0, 1e-5));

	REQUIRE(timing->getPhaseOffset().has_value());
	REQUIRE_THAT(timing->getPhaseOffset().value(), WithinAbs(3.14, 1e-5));

	REQUIRE(timing->getSyncOnPulse() == true);
	// Noise entries are added internally, we just verify it didn't crash
}

TEST_CASE("parseTiming accepts omitted defaults without warnings", "[serial][xml_parser_utils]")
{
	LogLevelGuard log_level(logging::Level::WARNING);
	CerrCapture capture;
	core::World world;
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;

	auto doc = loadXml("<timing name=\"quiet\"><frequency>1e6</frequency></timing>");

	serial::xml_parser_utils::parseTiming(doc.getRootElement(), ctx);
	REQUIRE(world.getTimings().size() == 1);

	auto* timing = world.getTimings().begin()->second.get();
	REQUIRE(timing->getName() == "quiet");
	REQUIRE_FALSE(timing->getFreqOffset().has_value());
	REQUIRE_FALSE(timing->getRandomFreqOffsetStdev().has_value());
	REQUIRE_FALSE(timing->getPhaseOffset().has_value());
	REQUIRE_FALSE(timing->getRandomPhaseOffsetStdev().has_value());
	REQUIRE_FALSE(timing->getSyncOnPulse());
	REQUIRE(capture.str().empty());
}

TEST_CASE("parseAntenna instantiates correct antenna types", "[serial][xml_parser_utils]")
{
	core::World world;
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;
	ctx.loaders = createMockLoaders();

	SECTION("Isotropic")
	{
		auto doc = loadXml("<antenna name=\"iso\" pattern=\"isotropic\"><efficiency>0.9</efficiency></antenna>");
		serial::xml_parser_utils::parseAntenna(doc.getRootElement(), ctx);
		REQUIRE(world.getAntennas().begin()->second->getEfficiencyFactor() == 0.9);
	}

	SECTION("Sinc")
	{
		auto doc = loadXml(
			"<antenna name=\"sinc1\" pattern=\"sinc\"><alpha>1</alpha><beta>2</beta><gamma>3</gamma></antenna>");
		serial::xml_parser_utils::parseAntenna(doc.getRootElement(), ctx);
		auto* ant = dynamic_cast<antenna::Sinc*>(world.getAntennas().begin()->second.get());
		REQUIRE(ant != nullptr);
		REQUIRE(ant->getAlpha() == 1.0);
	}

	SECTION("Gaussian")
	{
		auto doc = loadXml(
			"<antenna name=\"gauss\" pattern=\"gaussian\"><azscale>1.5</azscale><elscale>2.5</elscale></antenna>");
		serial::xml_parser_utils::parseAntenna(doc.getRootElement(), ctx);
		auto* ant = dynamic_cast<antenna::Gaussian*>(world.getAntennas().begin()->second.get());
		REQUIRE(ant != nullptr);
		REQUIRE_THAT(ant->getAzimuthScale(), WithinAbs(1.5, 1e-5));
		REQUIRE_THAT(ant->getElevationScale(), WithinAbs(2.5, 1e-5));
	}

	SECTION("SquareHorn")
	{
		auto doc = loadXml("<antenna name=\"horn\" pattern=\"squarehorn\"><diameter>0.5</diameter></antenna>");
		serial::xml_parser_utils::parseAntenna(doc.getRootElement(), ctx);
		auto* ant = dynamic_cast<antenna::SquareHorn*>(world.getAntennas().begin()->second.get());
		REQUIRE(ant != nullptr);
		REQUIRE_THAT(ant->getDimension(), WithinAbs(0.5, 1e-5));
	}

	SECTION("Parabolic")
	{
		auto doc = loadXml("<antenna name=\"dish\" pattern=\"parabolic\"><diameter>2.0</diameter></antenna>");
		serial::xml_parser_utils::parseAntenna(doc.getRootElement(), ctx);
		auto* ant = dynamic_cast<antenna::Parabolic*>(world.getAntennas().begin()->second.get());
		REQUIRE(ant != nullptr);
		REQUIRE_THAT(ant->getDiameter(), WithinAbs(2.0, 1e-5));
	}

	SECTION("XML File")
	{
		auto doc = loadXml("<antenna name=\"xml_ant\" pattern=\"xml\" filename=\"dummy.xml\"></antenna>");
		serial::xml_parser_utils::parseAntenna(doc.getRootElement(), ctx);
		REQUIRE(world.getAntennas().begin()->second->getName() == "xml_ant");
	}

	SECTION("HDF5 File")
	{
		auto doc = loadXml("<antenna name=\"h5_ant\" pattern=\"file\" filename=\"dummy.h5\"></antenna>");
		serial::xml_parser_utils::parseAntenna(doc.getRootElement(), ctx);
		REQUIRE(world.getAntennas().begin()->second->getName() == "h5_ant");
	}

	SECTION("Unsupported pattern throws")
	{
		auto doc = loadXml("<antenna name=\"bad\" pattern=\"magic\"></antenna>");
		REQUIRE_THROWS_AS(serial::xml_parser_utils::parseAntenna(doc.getRootElement(), ctx), XmlException);
	}
}

TEST_CASE("parseAntenna accepts omitted efficiency without warnings", "[serial][xml_parser_utils]")
{
	LogLevelGuard log_level(logging::Level::WARNING);
	CerrCapture capture;
	core::World world;
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;
	ctx.loaders = createMockLoaders();

	auto doc = loadXml("<antenna name=\"iso\" pattern=\"isotropic\"></antenna>");

	serial::xml_parser_utils::parseAntenna(doc.getRootElement(), ctx);
	REQUIRE(world.getAntennas().size() == 1);
	REQUIRE_THAT(world.getAntennas().begin()->second->getEfficiencyFactor(), WithinAbs(1.0, 1e-12));
	REQUIRE(capture.str().empty());
}

TEST_CASE("parseMotionPath handles all interpolation types", "[serial][xml_parser_utils]")
{
	radar::Platform platform("plat");

	SECTION("Linear")
	{
		auto doc =
			loadXml("<motionpath interpolation=\"linear\">"
					"  <positionwaypoint><x>1</x><y>2</y><altitude>3</altitude><time>0</time></positionwaypoint>"
					"  <positionwaypoint><x>4</x><y>5</y><altitude>6</altitude><time>10</time></positionwaypoint>"
					"</motionpath>");
		serial::xml_parser_utils::parseMotionPath(doc.getRootElement(), &platform);
		REQUIRE(platform.getMotionPath()->getType() == math::Path::InterpType::INTERP_LINEAR);
		REQUIRE_THAT(platform.getPosition(5.0).x, WithinAbs(2.5, 1e-5));
	}

	SECTION("Cubic")
	{
		auto doc = loadXml("<motionpath interpolation=\"cubic\">"
						   "  <positionwaypoint><x>0</x><y>0</y><altitude>0</altitude><time>0</time></positionwaypoint>"
						   "  <positionwaypoint><x>1</x><y>1</y><altitude>1</altitude><time>1</time></positionwaypoint>"
						   "  <positionwaypoint><x>2</x><y>2</y><altitude>2</altitude><time>2</time></positionwaypoint>"
						   "</motionpath>");
		serial::xml_parser_utils::parseMotionPath(doc.getRootElement(), &platform);
		REQUIRE(platform.getMotionPath()->getType() == math::Path::InterpType::INTERP_CUBIC);
	}

	SECTION("Unsupported defaults to static")
	{
		auto doc = loadXml("<motionpath interpolation=\"magic\"></motionpath>");
		serial::xml_parser_utils::parseMotionPath(doc.getRootElement(), &platform);
		REQUIRE(platform.getMotionPath()->getType() == math::Path::InterpType::INTERP_STATIC);
	}

	SECTION("Missing interpolation defaults to static without warnings")
	{
		LogLevelGuard log_level(logging::Level::WARNING);
		CerrCapture capture;
		auto doc = loadXml("<motionpath>"
						   "  <positionwaypoint><x>1</x><y>2</y><altitude>3</altitude><time>0</time></positionwaypoint>"
						   "</motionpath>");
		serial::xml_parser_utils::parseMotionPath(doc.getRootElement(), &platform);
		REQUIRE(platform.getMotionPath()->getType() == math::Path::InterpType::INTERP_STATIC);
		REQUIRE(capture.str().empty());
	}
}

TEST_CASE("parseRotationPath handles all interpolation types", "[serial][xml_parser_utils]")
{
	radar::Platform platform("plat");

	SECTION("Linear")
	{
		auto doc = loadXml(
			"<rotationpath interpolation=\"linear\">"
			"  <rotationwaypoint><azimuth>0</azimuth><elevation>0</elevation><time>0</time></rotationwaypoint>"
			"  <rotationwaypoint><azimuth>10</azimuth><elevation>10</elevation><time>1</time></rotationwaypoint>"
			"</rotationpath>");
		serial::xml_parser_utils::parseRotationPath(doc.getRootElement(), &platform,
													params::RotationAngleUnit::Degrees);
		REQUIRE(platform.getRotationPath()->getType() == math::RotationPath::InterpType::INTERP_LINEAR);
	}

	SECTION("Cubic")
	{
		auto doc = loadXml(
			"<rotationpath interpolation=\"cubic\">"
			"  <rotationwaypoint><azimuth>0</azimuth><elevation>0</elevation><time>0</time></rotationwaypoint>"
			"  <rotationwaypoint><azimuth>10</azimuth><elevation>10</elevation><time>1</time></rotationwaypoint>"
			"  <rotationwaypoint><azimuth>20</azimuth><elevation>20</elevation><time>2</time></rotationwaypoint>"
			"</rotationpath>");
		serial::xml_parser_utils::parseRotationPath(doc.getRootElement(), &platform,
													params::RotationAngleUnit::Degrees);
		REQUIRE(platform.getRotationPath()->getType() == math::RotationPath::InterpType::INTERP_CUBIC);
	}
}

TEST_CASE("parseFixedRotation sets constant rate rotation", "[serial][xml_parser_utils]")
{
	radar::Platform platform("plat");

	SECTION("Valid")
	{
		auto doc = loadXml("<fixedrotation>"
						   "  <startazimuth>90</startazimuth>"
						   "  <startelevation>0</startelevation>"
						   "  <azimuthrate>10</azimuthrate>"
						   "  <elevationrate>5</elevationrate>"
						   "</fixedrotation>");
		serial::xml_parser_utils::parseFixedRotation(doc.getRootElement(), &platform,
													 params::RotationAngleUnit::Degrees);

		// 90 deg azimuth -> 0 rad (since az_rad = (90 - az_deg) * PI/180)
		// rate 10 deg/s -> -10 * PI/180 rad/s
		REQUIRE_THAT(platform.getRotation(1.0).azimuth, WithinAbs(-10.0 * PI / 180.0, 1e-5));
		REQUIRE_THAT(platform.getRotation(1.0).elevation, WithinAbs(5.0 * PI / 180.0, 1e-5));
	}

	SECTION("Invalid throws")
	{
		auto doc = loadXml("<fixedrotation><startazimuth>90</startazimuth></fixedrotation>");
		REQUIRE_THROWS_AS(serial::xml_parser_utils::parseFixedRotation(doc.getRootElement(), &platform,
																	   params::RotationAngleUnit::Degrees),
						  XmlException);
	}
}

TEST_CASE("parseRotationPath warns when values look like the opposite unit", "[serial][xml_parser_utils]")
{
	LogLevelGuard log_guard(logging::Level::WARNING);
	CerrCapture capture;
	radar::Platform platform("warning-platform", 77);
	auto doc =
		loadXml("<rotationpath interpolation=\"static\">"
				"  <rotationwaypoint><azimuth>90</azimuth><elevation>0</elevation><time>0</time></rotationwaypoint>"
				"</rotationpath>");

	serial::xml_parser_utils::parseRotationPath(doc.getRootElement(), &platform, params::RotationAngleUnit::Radians);

	REQUIRE_THAT(capture.str(), ContainsSubstring("platform 'warning-platform' rotation waypoint 0"));
	REQUIRE_THAT(capture.str(), ContainsSubstring("'azimuth'"));
	REQUIRE_THAT(capture.str(), ContainsSubstring("declared"));
}

TEST_CASE("parseTransmitter resolves references and builds object with schedule", "[serial][xml_parser_utils]")
{
	ParamGuard guard;
	params::setRate(10000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	core::World world;
	std::mt19937 seeder(42);
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;
	ctx.master_seeder = &seeder;

	// Populate dependencies
	auto wave =
		std::make_unique<fers_signal::RadarSignal>("w1", 1.0, 1e9, 1.0, std::make_unique<fers_signal::CwSignal>(), 10);
	auto ant = std::make_unique<antenna::Isotropic>("a1", 20);
	auto tim = std::make_unique<timing::PrototypeTiming>("t1", 30);
	tim->setFrequency(1e6);

	world.add(std::move(wave));
	world.add(std::move(ant));
	world.add(std::move(tim));

	std::unordered_map<std::string, SimId> w_refs = {{"w1", 10}};
	std::unordered_map<std::string, SimId> a_refs = {{"a1", 20}};
	std::unordered_map<std::string, SimId> t_refs = {{"t1", 30}};
	serial::xml_parser_utils::ReferenceLookup refs{&w_refs, &a_refs, &t_refs};

	radar::Platform platform("plat");

	auto doc = loadXml("<transmitter name=\"tx1\" waveform=\"w1\" antenna=\"a1\" timing=\"t1\">"
					   "  <cw_mode/>"
					   "  <schedule><period start=\"0.1\" end=\"0.5\"/></schedule>"
					   "</transmitter>");

	auto* tx = serial::xml_parser_utils::parseTransmitter(doc.getRootElement(), &platform, ctx, refs);

	REQUIRE(tx->getName() == "tx1");
	REQUIRE(tx->getMode() == radar::OperationMode::CW_MODE);
	REQUIRE(tx->getSignal()->getName() == "w1");
	REQUIRE(tx->getAntenna()->getName() == "a1");
	REQUIRE(tx->getSchedule().size() == 1);
}

TEST_CASE("parseTransmitter rejects FMCW waveform and mode mismatches", "[serial][xml_parser_utils][fmcw]")
{
	ParamGuard guard;
	params::setRate(2.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	core::World world;
	std::mt19937 seeder(42);
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;
	ctx.master_seeder = &seeder;

	auto fmcw_signal = std::make_unique<fers_signal::FmcwChirpSignal>(1.0e6, 1.0e-3, 1.0e-3);
	world.add(std::make_unique<fers_signal::RadarSignal>("fmcw_wave", 1.0, 1e9, 1.0e-3, std::move(fmcw_signal), 10));
	world.add(std::make_unique<antenna::Isotropic>("a1", 20));
	auto timing_proto = std::make_unique<timing::PrototypeTiming>("t1", 30);
	timing_proto->setFrequency(1e6);
	world.add(std::move(timing_proto));

	std::unordered_map<std::string, SimId> w_refs = {{"fmcw_wave", 10}};
	std::unordered_map<std::string, SimId> a_refs = {{"a1", 20}};
	std::unordered_map<std::string, SimId> t_refs = {{"t1", 30}};
	serial::xml_parser_utils::ReferenceLookup refs{&w_refs, &a_refs, &t_refs};

	radar::Platform platform("plat");

	auto mismatch_doc = loadXml("<transmitter name=\"tx1\" waveform=\"fmcw_wave\" antenna=\"a1\" timing=\"t1\">"
								"  <cw_mode/>"
								"</transmitter>");
	REQUIRE_THROWS_AS(serial::xml_parser_utils::parseTransmitter(mismatch_doc.getRootElement(), &platform, ctx, refs),
					  XmlException);

	auto ambiguous_doc = loadXml("<transmitter name=\"tx2\" waveform=\"fmcw_wave\" antenna=\"a1\" timing=\"t1\">"
								 "  <cw_mode/>"
								 "  <fmcw_mode/>"
								 "</transmitter>");
	REQUIRE_THROWS_AS(serial::xml_parser_utils::parseTransmitter(ambiguous_doc.getRootElement(), &platform, ctx, refs),
					  XmlException);
}

TEST_CASE("parseTransmitter validates FMCW schedule duration against chirp timing", "[serial][xml_parser_utils][fmcw]")
{
	ParamGuard guard;
	params::setRate(2.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	core::World world;
	std::mt19937 seeder(42);
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;
	ctx.master_seeder = &seeder;

	auto fmcw_signal = std::make_unique<fers_signal::FmcwChirpSignal>(1.0e6, 1.0e-3, 2.0e-3);
	world.add(std::make_unique<fers_signal::RadarSignal>("fmcw_wave", 1.0, 1e9, 1.0e-3, std::move(fmcw_signal), 10));
	world.add(std::make_unique<antenna::Isotropic>("a1", 20));
	auto timing_proto = std::make_unique<timing::PrototypeTiming>("t1", 30);
	timing_proto->setFrequency(1e6);
	world.add(std::move(timing_proto));

	std::unordered_map<std::string, SimId> w_refs = {{"fmcw_wave", 10}};
	std::unordered_map<std::string, SimId> a_refs = {{"a1", 20}};
	std::unordered_map<std::string, SimId> t_refs = {{"t1", 30}};
	serial::xml_parser_utils::ReferenceLookup refs{&w_refs, &a_refs, &t_refs};

	radar::Platform platform("plat");

	SECTION("period shorter than T_c is rejected")
	{
		auto doc = loadXml("<transmitter name=\"tx_short_tc\" waveform=\"fmcw_wave\" antenna=\"a1\" timing=\"t1\">"
						   "  <fmcw_mode/>"
						   "  <schedule><period start=\"0.1\" end=\"0.1005\"/></schedule>"
						   "</transmitter>");

		REQUIRE_THROWS_WITH(serial::xml_parser_utils::parseTransmitter(doc.getRootElement(), &platform, ctx, refs),
							ContainsSubstring("shorter than FMCW chirp_duration T_c"));
	}

	SECTION("period shorter than T_rep but at least T_c only warns")
	{
		LogLevelGuard log_level(logging::Level::WARNING);
		CerrCapture capture;
		auto doc = loadXml("<transmitter name=\"tx_short_trep\" waveform=\"fmcw_wave\" antenna=\"a1\" timing=\"t1\">"
						   "  <fmcw_mode/>"
						   "  <schedule><period start=\"0.1\" end=\"0.1015\"/></schedule>"
						   "</transmitter>");

		auto* tx = serial::xml_parser_utils::parseTransmitter(doc.getRootElement(), &platform, ctx, refs);
		REQUIRE(tx->getSchedule().size() == 1);
		REQUIRE_THAT(capture.str(), ContainsSubstring("shorter than FMCW chirp_period"));
	}
}

TEST_CASE("parseReceiver resolves references and builds object with flags and schedule", "[serial][xml_parser_utils]")
{
	ParamGuard guard;
	params::setRate(10000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	core::World world;
	std::mt19937 seeder(42);
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;
	ctx.master_seeder = &seeder;

	auto ant = std::make_unique<antenna::Isotropic>("a1", 20);
	auto tim = std::make_unique<timing::PrototypeTiming>("t1", 30);
	tim->setFrequency(1e6);
	world.add(std::move(ant));
	world.add(std::move(tim));

	std::unordered_map<std::string, SimId> w_refs;
	std::unordered_map<std::string, SimId> a_refs = {{"a1", 20}};
	std::unordered_map<std::string, SimId> t_refs = {{"t1", 30}};
	serial::xml_parser_utils::ReferenceLookup refs{&w_refs, &a_refs, &t_refs};

	radar::Platform platform("plat");

	SECTION("Valid pulsed receiver with flags and schedule")
	{
		// Use a window_skip that is an exact multiple of the sample period (1/10000 = 1e-4)
		// to avoid quantization down to 0.0
		auto doc =
			loadXml("<receiver name=\"rx1\" antenna=\"a1\" timing=\"t1\" nodirect=\"true\" nopropagationloss=\"true\">"
					"  "
					"<pulsed_mode><prf>1000</prf><window_length>1e-4</window_length><window_skip>2e-4</window_skip></"
					"pulsed_mode>"
					"  <schedule><period start=\"0.1\" end=\"0.5\"/></schedule>"
					"</receiver>");

		auto* rx = serial::xml_parser_utils::parseReceiver(doc.getRootElement(), &platform, ctx, refs);

		REQUIRE(rx->getName() == "rx1");
		REQUIRE(rx->getMode() == radar::OperationMode::PULSED_MODE);
		REQUIRE_THAT(rx->getWindowPrf(), WithinAbs(1000.0, 1e-5));
		REQUIRE_THAT(rx->getWindowLength(), WithinAbs(1e-4, 1e-9));
		REQUIRE_THAT(rx->getWindowSkip(), WithinAbs(2e-4, 1e-9));
		REQUIRE(rx->checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT));
		REQUIRE(rx->checkFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS));
		REQUIRE(rx->getSchedule().size() == 1);
	}

	SECTION("Missing receiver defaults do not warn")
	{
		LogLevelGuard log_level(logging::Level::WARNING);
		CerrCapture capture;
		auto doc = loadXml("<receiver name=\"rx_defaults\" antenna=\"a1\" timing=\"t1\">"
						   "  <cw_mode/>"
						   "</receiver>");

		auto* rx = serial::xml_parser_utils::parseReceiver(doc.getRootElement(), &platform, ctx, refs);

		REQUIRE(rx->getName() == "rx_defaults");
		REQUIRE(rx->getMode() == radar::OperationMode::CW_MODE);
		REQUIRE_FALSE(rx->checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT));
		REQUIRE_FALSE(rx->checkFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS));
		REQUIRE(capture.str().empty());
	}

	SECTION("Invalid pulsed parameters throw")
	{
		auto doc1 = loadXml(
			"<receiver name=\"rx1\" antenna=\"a1\" "
			"timing=\"t1\"><pulsed_mode><prf>1000</prf><window_length>-1</window_length></pulsed_mode></receiver>");
		REQUIRE_THROWS_AS(serial::xml_parser_utils::parseReceiver(doc1.getRootElement(), &platform, ctx, refs),
						  XmlException);

		auto doc2 = loadXml(
			"<receiver name=\"rx1\" antenna=\"a1\" "
			"timing=\"t1\"><pulsed_mode><prf>-1000</prf><window_length>1e-4</window_length></pulsed_mode></receiver>");
		REQUIRE_THROWS_AS(serial::xml_parser_utils::parseReceiver(doc2.getRootElement(), &platform, ctx, refs),
						  XmlException);

		auto doc3 = loadXml("<receiver name=\"rx1\" antenna=\"a1\" "
							"timing=\"t1\"><pulsed_mode><prf>1000</prf><window_length>1e-4</"
							"window_length><window_skip>-1</window_skip></pulsed_mode></receiver>");
		REQUIRE_THROWS_AS(serial::xml_parser_utils::parseReceiver(doc3.getRootElement(), &platform, ctx, refs),
						  XmlException);
	}
}

TEST_CASE("parseMonostatic reuses one shared timing instance for a common timing id", "[serial][xml_parser_utils]")
{
	ParamGuard guard;
	params::setRate(10000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	core::World world;
	std::mt19937 seeder(42);
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;
	ctx.master_seeder = &seeder;

	world.add(
		std::make_unique<fers_signal::RadarSignal>("w1", 1.0, 1e9, 1.0, std::make_unique<fers_signal::CwSignal>(), 10));
	world.add(std::make_unique<antenna::Isotropic>("a1", 20));
	auto timing_proto = std::make_unique<timing::PrototypeTiming>("t1", 30);
	timing_proto->setFrequency(1e6);
	world.add(std::move(timing_proto));

	std::unordered_map<std::string, SimId> w_refs = {{"w1", 10}};
	std::unordered_map<std::string, SimId> a_refs = {{"a1", 20}};
	std::unordered_map<std::string, SimId> t_refs = {{"t1", 30}};
	serial::xml_parser_utils::ReferenceLookup refs{&w_refs, &a_refs, &t_refs};

	radar::Platform platform("plat");

	auto doc = loadXml("<monostatic name=\"mono\" waveform=\"w1\" antenna=\"a1\" timing=\"t1\">"
					   "  <cw_mode/>"
					   "</monostatic>");

	serial::xml_parser_utils::parseMonostatic(doc.getRootElement(), &platform, ctx, refs);

	REQUIRE(world.getTransmitters().size() == 1);
	REQUIRE(world.getReceivers().size() == 1);
	REQUIRE(world.getTransmitters().front()->getTiming().get() == world.getReceivers().front()->getTiming().get());
}

TEST_CASE("parseMonostatic derives FMCW mode from the monostatic block without receiver override",
		  "[serial][xml_parser_utils][fmcw]")
{
	ParamGuard guard;
	params::setRate(2.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	core::World world;
	std::mt19937 seeder(42);
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;
	ctx.master_seeder = &seeder;

	auto fmcw_signal = std::make_unique<fers_signal::FmcwChirpSignal>(1.0e6, 1.0e-3, 1.0e-3);
	world.add(std::make_unique<fers_signal::RadarSignal>("w1", 1.0, 1e9, 1.0e-3, std::move(fmcw_signal), 10));
	world.add(std::make_unique<antenna::Isotropic>("a1", 20));
	auto timing_proto = std::make_unique<timing::PrototypeTiming>("t1", 30);
	timing_proto->setFrequency(1e6);
	world.add(std::move(timing_proto));

	std::unordered_map<std::string, SimId> w_refs = {{"w1", 10}};
	std::unordered_map<std::string, SimId> a_refs = {{"a1", 20}};
	std::unordered_map<std::string, SimId> t_refs = {{"t1", 30}};
	serial::xml_parser_utils::ReferenceLookup refs{&w_refs, &a_refs, &t_refs};

	radar::Platform platform("plat");

	auto doc = loadXml("<monostatic name=\"mono\" waveform=\"w1\" antenna=\"a1\" timing=\"t1\">"
					   "  <fmcw_mode/>"
					   "</monostatic>");

	serial::xml_parser_utils::parseMonostatic(doc.getRootElement(), &platform, ctx, refs);

	REQUIRE(world.getTransmitters().size() == 1);
	REQUIRE(world.getReceivers().size() == 1);
	REQUIRE(world.getTransmitters().front()->getMode() == radar::OperationMode::FMCW_MODE);
	REQUIRE(world.getReceivers().front()->getMode() == radar::OperationMode::FMCW_MODE);
	REQUIRE(world.getTransmitters().front()->getAttached() == world.getReceivers().front().get());
}

TEST_CASE("parseTarget handles chisquare model", "[serial][xml_parser_utils]")
{
	core::World world;
	std::mt19937 seeder(42);
	const unsigned expected_seed = static_cast<unsigned>(seeder());
	seeder.seed(42);
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;
	ctx.master_seeder = &seeder;
	radar::Platform platform("plat");

	auto doc = loadXml("<target name=\"tgt1\">"
					   "  <rcs type=\"isotropic\"><value>10.0</value></rcs>"
					   "  <model type=\"chisquare\"><k>2.0</k></model>"
					   "</target>");

	serial::xml_parser_utils::parseTarget(doc.getRootElement(), &platform, ctx);
	REQUIRE(world.getTargets().size() == 1);

	auto* tgt = world.getTargets().front().get();
	auto* model = dynamic_cast<const radar::RcsChiSquare*>(tgt->getFluctuationModel());
	REQUIRE(model != nullptr);
	REQUIRE(tgt->getSeed() == expected_seed);
	REQUIRE_THAT(model->getK(), WithinAbs(2.0, 1e-5));
}

TEST_CASE("parsePlatform prefers rotationpath over fixedrotation and supports fixed-only", "[serial][xml_parser_utils]")
{
	core::World world;
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;

	std::unordered_map<std::string, SimId> w_refs;
	std::unordered_map<std::string, SimId> a_refs;
	std::unordered_map<std::string, SimId> t_refs;
	serial::xml_parser_utils::ReferenceLookup refs{&w_refs, &a_refs, &t_refs};

	auto register_name = [](const XmlElement&, std::string_view) {};

	SECTION("Both rotationpath and fixedrotation: rotationpath is used")
	{
		auto doc = loadXml(
			"<platform name=\"plat_both\">"
			"  <rotationpath interpolation=\"linear\">"
			"    <rotationwaypoint><azimuth>0</azimuth><elevation>0</elevation><time>0</time></rotationwaypoint>"
			"    <rotationwaypoint><azimuth>90</azimuth><elevation>0</elevation><time>1</time></rotationwaypoint>"
			"  </rotationpath>"
			"  <fixedrotation>"
			"    <startazimuth>90</startazimuth><startelevation>0</startelevation>"
			"    <azimuthrate>10</azimuthrate><elevationrate>0</elevationrate>"
			"  </fixedrotation>"
			"</platform>");

		serial::xml_parser_utils::parsePlatform(doc.getRootElement(), ctx, register_name, refs);

		REQUIRE(world.getPlatforms().size() == 1);
		auto* plat = world.getPlatforms().front().get();
		REQUIRE(plat->getRotationPath()->getType() == math::RotationPath::InterpType::INTERP_LINEAR);
		REQUIRE_THAT(plat->getRotation(1.0).azimuth, WithinAbs(0.0, 1e-5));
	}

	SECTION("Only fixedrotation: constant-rate rotation is used")
	{
		world.clear();

		auto doc = loadXml("<platform name=\"plat_fixed\">"
						   "  <fixedrotation>"
						   "    <startazimuth>90</startazimuth><startelevation>0</startelevation>"
						   "    <azimuthrate>10</azimuthrate><elevationrate>5</elevationrate>"
						   "  </fixedrotation>"
						   "</platform>");

		serial::xml_parser_utils::parsePlatform(doc.getRootElement(), ctx, register_name, refs);

		REQUIRE(world.getPlatforms().size() == 1);
		auto* plat = world.getPlatforms().front().get();
		REQUIRE(plat->getRotationPath()->getType() == math::RotationPath::InterpType::INTERP_CONSTANT);
		REQUIRE_THAT(plat->getRotation(1.0).azimuth, WithinAbs(-10.0 * PI / 180.0, 1e-5));
		REQUIRE_THAT(plat->getRotation(1.0).elevation, WithinAbs(5.0 * PI / 180.0, 1e-5));
	}
}

TEST_CASE("collectIncludeElements skips empty and unreadable includes", "[serial][xml_parser_utils]")
{
	SECTION("Empty include text is ignored")
	{
		auto doc = loadXml("<simulation><include></include></simulation>");
		std::vector<std::filesystem::path> include_paths;

		serial::xml_parser_utils::collectIncludeElements(doc, std::filesystem::temp_directory_path(), include_paths);

		REQUIRE(include_paths.empty());
	}

	SECTION("Unreadable include file is collected once and skipped")
	{
		const std::string missing_name = "collect_include_missing_file_7f9a8115.xml";
		auto doc = loadXml("<simulation><include>" + missing_name + "</include></simulation>");
		const auto base_dir = std::filesystem::temp_directory_path();
		std::vector<std::filesystem::path> include_paths;

		serial::xml_parser_utils::collectIncludeElements(doc, base_dir, include_paths);

		REQUIRE(include_paths.size() == 1);
		REQUIRE(include_paths[0] == (base_dir / missing_name));
	}
}
