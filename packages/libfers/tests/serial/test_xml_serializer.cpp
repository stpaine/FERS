// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <fstream>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/parameters.h"
#include "core/world.h"
#include "math/path.h"
#include "math/rotation_path.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/libxml_wrapper.h"
#include "serial/xml_parser_utils.h"
#include "serial/xml_serializer.h"
#include "serial/xml_serializer_utils.h"
#include "signal/radar_signal.h"
#include "timing/prototype_timing.h"
#include "timing/timing.h"

using Catch::Matchers::ContainsSubstring;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
	};

	std::string dumpElement(const XmlElement& elem)
	{
		xmlBufferPtr buf = xmlBufferCreate();
		xmlNodeDump(buf, nullptr, elem.getNode(), 0, 0);
		std::string result(reinterpret_cast<const char*>(xmlBufferContent(buf)));
		xmlBufferFree(buf);
		return result;
	}
}

TEST_CASE("serializeParameters creates correct tags across frames", "[serial][xml_serializer]")
{
	XmlDocument doc;
	XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);

	params::Parameters p;
	p.start = 1.0;
	p.end = 2.0;
	p.rate = 1000.0;
	p.c = 3e8;
	p.sim_sampling_rate = 500.0;
	p.random_seed = 42;
	p.adc_bits = 12;
	p.oversample_ratio = 4;
	p.origin_latitude = -33.0;
	p.origin_longitude = 18.0;
	p.origin_altitude = 100.0;

	SECTION("UTM North")
	{
		p.coordinate_frame = params::CoordinateFrame::UTM;
		p.utm_zone = 34;
		p.utm_north_hemisphere = true;

		serial::xml_serializer_utils::serializeParameters(root, p);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("<starttime>1</starttime>"));
		REQUIRE_THAT(s, ContainsSubstring("<c>3e+08</c>"));
		REQUIRE_THAT(s, ContainsSubstring("<simSamplingRate>500</simSamplingRate>"));
		REQUIRE_THAT(s, ContainsSubstring("<randomseed>42</randomseed>"));
		REQUIRE_THAT(s, ContainsSubstring("<adc_bits>12</adc_bits>"));
		REQUIRE_THAT(s, ContainsSubstring("<oversample>4</oversample>"));
		REQUIRE_THAT(s, ContainsSubstring("<coordinatesystem frame=\"UTM\" zone=\"34\" hemisphere=\"N\"/>"));
	}

	SECTION("UTM South")
	{
		p.coordinate_frame = params::CoordinateFrame::UTM;
		p.utm_zone = 35;
		p.utm_north_hemisphere = false;

		serial::xml_serializer_utils::serializeParameters(root, p);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("<coordinatesystem frame=\"UTM\" zone=\"35\" hemisphere=\"S\"/>"));
	}

	SECTION("ECEF")
	{
		p.coordinate_frame = params::CoordinateFrame::ECEF;
		p.rotation_angle_unit = params::RotationAngleUnit::Radians;
		serial::xml_serializer_utils::serializeParameters(root, p);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("<coordinatesystem frame=\"ECEF\"/>"));
		REQUIRE_THAT(s, ContainsSubstring("<rotationangleunit>rad</rotationangleunit>"));
	}

	SECTION("ENU")
	{
		p.coordinate_frame = params::CoordinateFrame::ENU;
		serial::xml_serializer_utils::serializeParameters(root, p);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("<coordinatesystem frame=\"ENU\"/>"));
	}

	SECTION("Defaults parameters omitted gracefully")
	{
		params::Parameters p2; // Default constructor
		serial::xml_serializer_utils::serializeParameters(root, p2);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, !ContainsSubstring("<c>"));
		REQUIRE_THAT(s, !ContainsSubstring("<simSamplingRate>"));
		REQUIRE_THAT(s, !ContainsSubstring("<randomseed>"));
		REQUIRE_THAT(s, !ContainsSubstring("<adc_bits>"));
		REQUIRE_THAT(s, !ContainsSubstring("<oversample>"));
		REQUIRE_THAT(s, !ContainsSubstring("<rotationangleunit>"));
	}
}

TEST_CASE("serializeWaveform processes CW and Pulsed correctly", "[serial][xml_serializer]")
{
	XmlDocument doc;
	XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);

	SECTION("CW Mode")
	{
		auto sig = std::make_unique<fers_signal::CwSignal>();
		fers_signal::RadarSignal wave("w1", 10.0, 1e9, 1.0, std::move(sig));
		serial::xml_serializer_utils::serializeWaveform(wave, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("name=\"w1\""));
		REQUIRE_THAT(s, ContainsSubstring("<cw/>"));
	}

	SECTION("Pulsed Mode File")
	{
		auto sig = std::make_unique<fers_signal::Signal>(); // Implies NOT CwSignal
		fers_signal::RadarSignal wave("w2", 20.0, 2e9, 1.0, std::move(sig));
		wave.setFilename("pulse.csv");
		serial::xml_serializer_utils::serializeWaveform(wave, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("name=\"w2\""));
		REQUIRE_THAT(s, ContainsSubstring("<pulsed_from_file filename=\"pulse.csv\"/>"));
	}

	SECTION("Pulsed Missing Filename safely defaults")
	{
		auto sig = std::make_unique<fers_signal::Signal>();
		fers_signal::RadarSignal wave("w3", 20.0, 2e9, 1.0, std::move(sig));
		serial::xml_serializer_utils::serializeWaveform(wave, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("<pulsed_from_file filename=\"\"/>"));
	}
}

TEST_CASE("serializeWaveform round trips FMCW linear chirp direction", "[serial][xml_serializer][fmcw]")
{
	ParamGuard guard;
	params::setRate(2.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 1.0);

	XmlDocument doc;
	XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("waveform")));
	doc.setRootElement(root);

	auto sig = std::make_unique<fers_signal::FmcwChirpSignal>(1.0e6, 1.0e-3, 1.0e-3, 0.0, std::nullopt,
															  fers_signal::FmcwChirpDirection::Down);
	fers_signal::RadarSignal wave("down", 20.0, 2e9, 1.0e-3, std::move(sig));

	serial::xml_serializer_utils::serializeWaveform(wave, root);
	const std::string serialized = dumpElement(root);
	REQUIRE_THAT(serialized, ContainsSubstring("<fmcw_linear_chirp direction=\"down\">"));
	REQUIRE_THAT(serialized, !ContainsSubstring("fmcw_up_chirp"));

	core::World world;
	std::mt19937 seeder(42);
	serial::xml_parser_utils::ParserContext ctx;
	ctx.world = &world;
	ctx.master_seeder = &seeder;
	ctx.parameters.start = 0.0;
	ctx.parameters.end = 1.0;
	serial::xml_parser_utils::parseWaveform(root, ctx);

	REQUIRE(world.getWaveforms().size() == 1);
	const auto* parsed_wave = world.getWaveforms().begin()->second.get();
	REQUIRE(parsed_wave->getFmcwChirpSignal() != nullptr);
	REQUIRE(parsed_wave->getFmcwChirpSignal()->isDownChirp());
}

TEST_CASE("serializeTiming preserves clock phase and jitter characteristics", "[serial][xml_serializer]")
{
	XmlDocument doc;
	XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);

	timing::PrototypeTiming t("t1");
	t.setFrequency(123456.0);
	t.setFreqOffset(10.0);
	t.setPhaseOffset(3.14);
	t.setRandomFreqOffsetStdev(5.0);
	t.setRandomPhaseOffsetStdev(1.0);
	t.setSyncOnPulse();
	t.setAlpha(1.0, 0.5);

	serial::xml_serializer_utils::serializeTiming(t, root);
	std::string s = dumpElement(root);
	REQUIRE_THAT(s, ContainsSubstring("name=\"t1\""));
	REQUIRE_THAT(s, ContainsSubstring("synconpulse=\"true\""));
	REQUIRE_THAT(s, ContainsSubstring("<frequency>123456</frequency>"));
	REQUIRE_THAT(s, ContainsSubstring("<freq_offset>10</freq_offset>"));
	REQUIRE_THAT(s, ContainsSubstring("<phase_offset>3.14</phase_offset>"));
	REQUIRE_THAT(s, ContainsSubstring("<random_freq_offset_stdev>5</random_freq_offset_stdev>"));
	REQUIRE_THAT(s, ContainsSubstring("<random_phase_offset_stdev>1</random_phase_offset_stdev>"));
	REQUIRE_THAT(s, ContainsSubstring("<noise_entry><alpha>1</alpha><weight>0.5</weight></noise_entry>"));
}

TEST_CASE("serializeAntenna dispatches to correct pattern forms", "[serial][xml_serializer]")
{
	XmlDocument doc;
	XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);

	SECTION("Isotropic")
	{
		antenna::Isotropic ant("a1");
		ant.setEfficiencyFactor(0.9);
		serial::xml_serializer_utils::serializeAntenna(ant, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("pattern=\"isotropic\""));
		REQUIRE_THAT(s, ContainsSubstring("<efficiency>0.9</efficiency>"));
	}

	SECTION("Sinc")
	{
		antenna::Sinc ant("a2", 1.0, 2.0, 3.0);
		serial::xml_serializer_utils::serializeAntenna(ant, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("pattern=\"sinc\""));
		REQUIRE_THAT(s, ContainsSubstring("<alpha>1</alpha>"));
	}

	SECTION("Gaussian")
	{
		antenna::Gaussian ant("a3", 1.5, 2.5);
		serial::xml_serializer_utils::serializeAntenna(ant, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("pattern=\"gaussian\""));
		REQUIRE_THAT(s, ContainsSubstring("<azscale>1.5</azscale>"));
	}

	SECTION("SquareHorn")
	{
		antenna::SquareHorn ant("a4", 0.5);
		serial::xml_serializer_utils::serializeAntenna(ant, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("pattern=\"squarehorn\""));
		REQUIRE_THAT(s, ContainsSubstring("<diameter>0.5</diameter>"));
	}

	SECTION("Parabolic")
	{
		antenna::Parabolic ant("a5", 2.0);
		serial::xml_serializer_utils::serializeAntenna(ant, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("pattern=\"parabolic\""));
		REQUIRE_THAT(s, ContainsSubstring("<diameter>2</diameter>"));
	}

	// TODO: XmlAntenna and H5Antenna coverage skipped as they require structurally-valid backing files
}

TEST_CASE("serializeMotionPath sets right interpolation models", "[serial][xml_serializer]")
{
	XmlDocument doc;
	XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);

	math::Path path;
	path.addCoord({math::Vec3(1, 2, 3), 0.0});
	path.addCoord({math::Vec3(4, 5, 6), 10.0});

	SECTION("Static")
	{
		path.setInterp(math::Path::InterpType::INTERP_STATIC);
		serial::xml_serializer_utils::serializeMotionPath(path, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("interpolation=\"static\""));
	}
	SECTION("Linear")
	{
		path.setInterp(math::Path::InterpType::INTERP_LINEAR);
		serial::xml_serializer_utils::serializeMotionPath(path, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("interpolation=\"linear\""));
		REQUIRE_THAT(s, ContainsSubstring("<x>1</x>"));
	}
	SECTION("Cubic")
	{
		path.setInterp(math::Path::InterpType::INTERP_CUBIC);
		serial::xml_serializer_utils::serializeMotionPath(path, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("interpolation=\"cubic\""));
	}
}

TEST_CASE("serializeRotation translates internal math to compass correctly", "[serial][xml_serializer]")
{
	XmlDocument doc;
	XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);

	math::RotationPath rot;

	SECTION("Fixed Constant")
	{
		// Internal math equivalent to Compass 90 (East), Elevation 0
		math::RotationCoord start{(90.0 - 90.0) * PI / 180.0, 0, 0};
		// Rotate +10 deg/sec in compass, up 5 deg/sec in elevation
		math::RotationCoord rate{-10.0 * PI / 180.0, 5.0 * PI / 180.0, 0};
		rot.setConstantRate(start, rate);

		serial::xml_serializer_utils::serializeRotation(rot, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("<fixedrotation>"));
		REQUIRE_THAT(s, ContainsSubstring("<startazimuth>90</startazimuth>"));
		REQUIRE_THAT(s, ContainsSubstring("<startelevation>0</startelevation>"));
		REQUIRE_THAT(s, ContainsSubstring("<azimuthrate>10</azimuthrate>"));
		REQUIRE_THAT(s, ContainsSubstring("<elevationrate>5</elevationrate>"));
	}

	SECTION("Linear")
	{
		rot.setInterp(math::RotationPath::InterpType::INTERP_LINEAR);
		rot.addCoord({(90.0 - 45.0) * PI / 180.0, 10.0 * PI / 180.0, 0});
		serial::xml_serializer_utils::serializeRotation(rot, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("<rotationpath interpolation=\"linear\">"));
		REQUIRE_THAT(s, ContainsSubstring("<azimuth>45</azimuth>"));
		REQUIRE_THAT(s, ContainsSubstring("<elevation>10</elevation>"));
	}

	SECTION("Static")
	{
		rot.setInterp(math::RotationPath::InterpType::INTERP_STATIC);
		serial::xml_serializer_utils::serializeRotation(rot, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("<rotationpath interpolation=\"static\""));
	}

	SECTION("Cubic")
	{
		rot.setInterp(math::RotationPath::InterpType::INTERP_CUBIC);
		rot.addCoord({(90.0 - 45.0) * PI / 180.0, 10.0 * PI / 180.0, 0});
		serial::xml_serializer_utils::serializeRotation(rot, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("<rotationpath interpolation=\"cubic\">"));
		REQUIRE_THAT(s, ContainsSubstring("<azimuth>45</azimuth>"));
	}

	SECTION("Default/Unknown")
	{
		rot.setInterp(static_cast<math::RotationPath::InterpType>(999));
		serial::xml_serializer_utils::serializeRotation(rot, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("<rotationpath/>"));
	}
}

TEST_CASE("serializeTransmitter", "[serial][xml_serializer]")
{
	ParamGuard guard;
	params::setRate(10000.0);

	XmlDocument doc;
	XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);

	radar::Platform plat("p1");
	auto timing = std::make_shared<timing::Timing>("t1", 42);

	SECTION("CW w/ Schedule")
	{
		radar::Transmitter tx(&plat, "tx1", radar::OperationMode::CW_MODE);
		tx.setTiming(timing);
		tx.setSchedule({{1.0, 2.0}});
		serial::xml_serializer_utils::serializeTransmitter(tx, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("name=\"tx1\""));
		REQUIRE_THAT(s, ContainsSubstring("<cw_mode/>"));
		REQUIRE_THAT(s, ContainsSubstring("<period start=\"1.000000\" end=\"2.000000\"/>"));
	}

	SECTION("Pulsed")
	{
		radar::Transmitter tx(&plat, "tx2", radar::OperationMode::PULSED_MODE);
		tx.setTiming(timing);
		tx.setPrf(1000);
		serial::xml_serializer_utils::serializeTransmitter(tx, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("<pulsed_mode><prf>1000</prf></pulsed_mode>"));
	}
}

TEST_CASE("serializeReceiver", "[serial][xml_serializer]")
{
	XmlDocument doc;
	XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);

	radar::Platform plat("p1");
	auto timing = std::make_shared<timing::Timing>("t1", 42);

	SECTION("CW Noise + Flags")
	{
		radar::Receiver rx(&plat, "rx1", 42, radar::OperationMode::CW_MODE);
		rx.setTiming(timing);
		rx.setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
		rx.setNoiseTemperature(290.0);
		serial::xml_serializer_utils::serializeReceiver(rx, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("name=\"rx1\""));
		REQUIRE_THAT(s, ContainsSubstring("nodirect=\"true\""));
		REQUIRE_THAT(s, ContainsSubstring("nopropagationloss=\"false\""));
		REQUIRE_THAT(s, ContainsSubstring("<cw_mode/>"));
		REQUIRE_THAT(s, ContainsSubstring("<noise_temp>290</noise_temp>"));
	}

	SECTION("Pulsed Options")
	{
		radar::Receiver rx(&plat, "rx2", 42, radar::OperationMode::PULSED_MODE);
		rx.setTiming(timing);
		params::setRate(10000.0);
		params::setOversampleRatio(1);
		rx.setWindowProperties(1e-4, 1000.0, 0.0);
		serial::xml_serializer_utils::serializeReceiver(rx, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("<pulsed_mode>"));
		REQUIRE_THAT(s, ContainsSubstring("<prf>1000</prf>"));
		REQUIRE_THAT(s, ContainsSubstring("<window_skip>0</window_skip>"));
		REQUIRE_THAT(s, ContainsSubstring("<window_length>1e-04</window_length>"));
	}

	SECTION("FMCW IF-chain options")
	{
		radar::Receiver rx(&plat, "rx3", 42, radar::OperationMode::FMCW_MODE);
		rx.setTiming(timing);
		rx.setDechirpMode(radar::Receiver::DechirpMode::Physical);
		rx.setDechirpReference({.source = radar::Receiver::DechirpReferenceSource::Attached,
								.name = "",
								.transmitter_name = "",
								.waveform_name = ""});
		rx.setFmcwIfChainRequest(
			{.sample_rate_hz = 1.0e6, .filter_bandwidth_hz = 4.0e5, .filter_transition_width_hz = 1.0e5});
		serial::xml_serializer_utils::serializeReceiver(rx, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("dechirp_mode=\"physical\""));
		REQUIRE_THAT(s, ContainsSubstring("<dechirp_reference source=\"attached\"/>"));
		REQUIRE_THAT(s, ContainsSubstring("<if_sample_rate>1e+06</if_sample_rate>"));
		REQUIRE_THAT(s, ContainsSubstring("<if_filter_bandwidth>4e+05</if_filter_bandwidth>"));
		REQUIRE_THAT(s, ContainsSubstring("<if_filter_transition_width>1e+05</if_filter_transition_width>"));
		REQUIRE(s.find("<dechirp_reference") < s.find("<if_sample_rate>"));
	}
}

TEST_CASE("serializeMonostatic pairs matching attached Transmitter/Receiver", "[serial][xml_serializer]")
{
	XmlDocument doc;
	XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);

	radar::Platform plat("p1");
	auto timing = std::make_shared<timing::Timing>("t1", 42);
	radar::Transmitter tx(&plat, "mono1", radar::OperationMode::PULSED_MODE);
	tx.setTiming(timing);
	tx.setPrf(1000.0);
	radar::Receiver rx(&plat, "rx1", 42, radar::OperationMode::PULSED_MODE);
	rx.setTiming(timing);
	params::setRate(10000.0);
	params::setOversampleRatio(1);
	rx.setWindowProperties(1e-4, 1000.0, 0.0);
	rx.setNoiseTemperature(290.0);

	serial::xml_serializer_utils::serializeMonostatic(tx, rx, root);
	std::string s = dumpElement(root);
	REQUIRE_THAT(s, ContainsSubstring("<monostatic"));
	REQUIRE_THAT(s, ContainsSubstring("name=\"mono1\""));
	REQUIRE_THAT(s, ContainsSubstring("<noise_temp>290</noise_temp>"));
}

TEST_CASE("serializeTarget models", "[serial][xml_serializer]")
{
	XmlDocument doc;
	XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);

	radar::Platform plat("p1");

	SECTION("Isotropic")
	{
		auto t = radar::createIsoTarget(&plat, "tgt1", 10.0, 42);
		serial::xml_serializer_utils::serializeTarget(*t, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("name=\"tgt1\""));
		REQUIRE_THAT(s, ContainsSubstring("<rcs type=\"isotropic\"><value>10</value></rcs>"));
	}

	SECTION("Fluctuation Model ChiSquare")
	{
		auto t = radar::createIsoTarget(&plat, "tgt2", 10.0, 42);
		t->setFluctuationModel(std::make_unique<radar::RcsChiSquare>(t->getRngEngine(), 2.0));
		serial::xml_serializer_utils::serializeTarget(*t, root);
		std::string s = dumpElement(root);
		REQUIRE_THAT(s, ContainsSubstring("<model type=\"chisquare\"><k>2</k></model>"));
	}

	// TODO: Missing Coverage for `FileTarget` due to inability to mock generic external valid file loads cleanly
}

TEST_CASE("serializePlatform iterates correctly linked nodes", "[serial][xml_serializer]")
{
	XmlDocument doc;
	XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);

	core::World world;
	auto plat = std::make_unique<radar::Platform>("p1");
	auto timing = std::make_shared<timing::Timing>("t1", 42);

	auto tx = std::make_unique<radar::Transmitter>(plat.get(), "tx1", radar::OperationMode::CW_MODE);
	tx->setTiming(timing);
	auto rx = std::make_unique<radar::Receiver>(plat.get(), "rx1", 42, radar::OperationMode::CW_MODE);
	rx->setTiming(timing);
	auto tgt = radar::createIsoTarget(plat.get(), "tgt1", 10.0, 42);

	world.add(std::make_unique<radar::Transmitter>(plat.get(), "mono_tx", radar::OperationMode::CW_MODE));
	world.getTransmitters().back()->setTiming(timing);
	world.add(std::make_unique<radar::Receiver>(plat.get(), "mono_rx", 42, radar::OperationMode::CW_MODE));
	world.getReceivers().back()->setTiming(timing);
	world.getTransmitters().back()->setAttached(world.getReceivers().back().get());
	world.getReceivers().back()->setAttached(world.getTransmitters().back().get());

	world.add(std::move(tx));
	world.add(std::move(rx));
	world.add(std::move(tgt));

	serial::xml_serializer_utils::serializePlatform(*plat, world, root);
	std::string s = dumpElement(root);

	REQUIRE_THAT(s, ContainsSubstring("name=\"p1\""));
	REQUIRE_THAT(s, ContainsSubstring("<transmitter name=\"tx1\""));
	REQUIRE_THAT(s, ContainsSubstring("<receiver name=\"rx1\""));
	REQUIRE_THAT(s, ContainsSubstring("<monostatic name=\"mono_tx\""));
	REQUIRE_THAT(s, ContainsSubstring("<target name=\"tgt1\""));
}

TEST_CASE("world_to_xml_string outputs a full well-formed simulation", "[serial][xml_serializer]")
{
	ParamGuard guard;
	core::World world;
	SECTION("Custom Name")
	{
		params::params.simulation_name = "Custom Name";
		std::string xml = serial::world_to_xml_string(world);

		REQUIRE_THAT(xml, ContainsSubstring("<simulation name=\"Custom Name\">"));
		REQUIRE_THAT(xml, ContainsSubstring("<parameters>"));
		REQUIRE_THAT(xml, ContainsSubstring("<origin latitude"));
	}

	SECTION("Default Name")
	{
		params::params.simulation_name = "";
		std::string xml = serial::world_to_xml_string(world);

		REQUIRE_THAT(xml, ContainsSubstring("<simulation name=\"FERS Scenario\">"));
		REQUIRE_THAT(xml, ContainsSubstring("<parameters>"));
	}
}

TEST_CASE("addChildWithNumber handles various floating point values", "[serial][xml_serializer]")
{
	XmlDocument doc;
	XmlElement root(xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("root")));
	doc.setRootElement(root);

	serial::xml_serializer_utils::addChildWithNumber(root, "nan_val", std::numeric_limits<double>::quiet_NaN());
	serial::xml_serializer_utils::addChildWithNumber(root, "inf_val", std::numeric_limits<double>::infinity());
	serial::xml_serializer_utils::addChildWithNumber(root, "float_val", 3.14f);
	serial::xml_serializer_utils::addChildWithNumber(root, "double_val", 3.14159);
	serial::xml_serializer_utils::addChildWithNumber(root, "ldouble_val", 3.1415926535L);

	std::string s = dumpElement(root);
	REQUIRE_THAT(s, ContainsSubstring("<nan_val>nan</nan_val>"));
	REQUIRE_THAT(s, ContainsSubstring("<inf_val>inf</inf_val>"));
	REQUIRE_THAT(s, ContainsSubstring("<float_val>3.14"));
	REQUIRE_THAT(s, ContainsSubstring("<double_val>3.14159"));
	REQUIRE_THAT(s, ContainsSubstring("<ldouble_val>3.14159"));
}
