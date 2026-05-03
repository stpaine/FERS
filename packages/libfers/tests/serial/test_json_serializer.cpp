// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <iostream>
#include <nlohmann/json.hpp>
#include <random>
#include <sstream>
#include <string_view>

#include "antenna/antenna_factory.h"
#include "core/logging.h"
#include "core/parameters.h"
#include "core/world.h"
#include "math/coord.h"
#include "math/path.h"
#include "math/rotation_path.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/json_serializer.h"
#include "signal/radar_signal.h"
#include "timing/prototype_timing.h"
#include "timing/timing.h"

using json = nlohmann::json;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace math
{
	// Forward-declare ADL functions to make them visible to nlohmann::json in this TU.
	// Their definitions are in json_serializer.cpp.
	void to_json(nlohmann::json& j, const Vec3& v);
	void from_json(const nlohmann::json& j, Vec3& v);
}

namespace fers_signal
{
	void to_json(nlohmann::json& j, const RadarSignal& rs);
}

namespace
{
	// Protects the global `params::params` from bleeding state across tests
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
}

TEST_CASE("JSON: Granular parsing of Antenna and Waveform", "[serial][json]")
{
	ParamGuard guard;

	SECTION("Parse Antenna")
	{
		json ant_json = {{"id", 123}, {"name", "TestAntenna"}, {"pattern", "isotropic"}, {"efficiency", 0.85}};

		auto ant = serial::parse_antenna_from_json(ant_json);
		REQUIRE(ant != nullptr);
		REQUIRE(ant->getId() == 123);
		REQUIRE(ant->getName() == "TestAntenna");
		REQUIRE_THAT(ant->getEfficiencyFactor(), WithinAbs(0.85, 1e-9));
	}

	SECTION("Parse Waveform")
	{
		json wf_json = {{"id", 456},
						{"name", "TestWaveform"},
						{"power", 500.0},
						{"carrier_frequency", 2e9},
						{"cw", json::object()}};

		auto wf = serial::parse_waveform_from_json(wf_json);
		REQUIRE(wf != nullptr);
		REQUIRE(wf->getId() == 456);
		REQUIRE(wf->getName() == "TestWaveform");
		REQUIRE_THAT(wf->getPower(), WithinAbs(500.0, 1e-9));
		REQUIRE_THAT(wf->getCarrier(), WithinAbs(2e9, 1e-9));
	}
}

TEST_CASE("JSON: FMCW waveform emits large-buffer warning", "[serial][json]")
{
	ParamGuard guard;
	LogLevelGuard log_level(logging::Level::WARNING);
	CerrCapture capture;

	params::setRate(1.0e9);
	params::setOversampleRatio(1);
	params::setTime(0.0, 5.0);

	json wf_json = {
		{"id", 456},
		{"name", "HugeFmcw"},
		{"power", 500.0},
		{"carrier_frequency", 2.4e9},
		{"fmcw_linear_chirp",
		 {{"direction", "up"}, {"chirp_bandwidth", 1.0e6}, {"chirp_duration", 1.0e-3}, {"chirp_period", 1.0e-3}}}};

	auto wf = serial::parse_waveform_from_json(wf_json);
	REQUIRE(wf != nullptr);
	REQUIRE_THAT(capture.str(), ContainsSubstring("GiB of FMCW streaming IQ data"));
}

TEST_CASE("JSON: FMCW linear chirp direction round trips", "[serial][json][fmcw]")
{
	ParamGuard guard;
	params::setRate(2.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 1.0);

	for (const auto direction : {"up", "down"})
	{
		json wf_json = {{"id", 457},
						{"name", std::string("Fmcw") + direction},
						{"power", 500.0},
						{"carrier_frequency", 2.4e9},
						{"fmcw_linear_chirp",
						 {{"direction", direction},
						  {"chirp_bandwidth", 1.0e6},
						  {"chirp_duration", 1.0e-3},
						  {"chirp_period", 1.0e-3}}}};

		auto wf = serial::parse_waveform_from_json(wf_json);
		REQUIRE(wf != nullptr);
		REQUIRE(wf->getFmcwChirpSignal() != nullptr);
		REQUIRE(wf->getFmcwChirpSignal()->isDownChirp() == (std::string_view(direction) == "down"));

		json serialized;
		fers_signal::to_json(serialized, *wf);
		REQUIRE(serialized.contains("fmcw_linear_chirp"));
		REQUIRE_FALSE(serialized.contains("fmcw_up_chirp"));
		REQUIRE(serialized.at("fmcw_linear_chirp").at("direction") == direction);

		auto reparsed = serial::parse_waveform_from_json(serialized);
		REQUIRE(reparsed != nullptr);
		REQUIRE(reparsed->getFmcwChirpSignal() != nullptr);
		REQUIRE(reparsed->getFmcwChirpSignal()->isDownChirp() == (std::string_view(direction) == "down"));
	}
}

TEST_CASE("JSON: FMCW triangle round trips", "[serial][json][fmcw]")
{
	ParamGuard guard;
	params::setRate(2.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 1.0);

	json wf_json = {{"id", 458},
					{"name", "Tri"},
					{"power", 500.0},
					{"carrier_frequency", 2.4e9},
					{"fmcw_triangle",
					 {{"chirp_bandwidth", 1.0e6},
					  {"chirp_duration", 1.0e-3},
					  {"start_frequency_offset", 1.0e3},
					  {"triangle_count", 3}}}};

	auto wf = serial::parse_waveform_from_json(wf_json);
	REQUIRE(wf != nullptr);
	REQUIRE(wf->isFmcwFamily());
	REQUIRE(wf->isFmcwTriangle());
	REQUIRE(wf->getFmcwTriangleSignal() != nullptr);
	REQUIRE(wf->getFmcwTriangleSignal()->getTriangleCount().value() == 3);

	json serialized;
	fers_signal::to_json(serialized, *wf);
	REQUIRE(serialized.contains("fmcw_triangle"));
	REQUIRE_FALSE(serialized.contains("fmcw_linear_chirp"));
	REQUIRE(serialized.at("fmcw_triangle").at("triangle_count") == 3);

	auto reparsed = serial::parse_waveform_from_json(serialized);
	REQUIRE(reparsed != nullptr);
	REQUIRE(reparsed->getFmcwTriangleSignal() != nullptr);
	REQUIRE(reparsed->getFmcwTriangleSignal()->getTriangleCount().value() == 3);
}

TEST_CASE("JSON: FMCW triangle rejects fractional triangle count", "[serial][json][fmcw]")
{
	ParamGuard guard;
	params::setRate(2.0e6);
	params::setOversampleRatio(1);
	params::setTime(0.0, 1.0);

	json wf_json = {
		{"id", 459},
		{"name", "TriFractional"},
		{"power", 500.0},
		{"carrier_frequency", 2.4e9},
		{"fmcw_triangle", {{"chirp_bandwidth", 1.0e6}, {"chirp_duration", 1.0e-3}, {"triangle_count", 1.5}}}};

	REQUIRE_THROWS_WITH(serial::parse_waveform_from_json(wf_json), ContainsSubstring("invalid triangle_count"));
}

TEST_CASE("JSON: Serialization of Math and Timing Structures", "[serial][json]")
{
	ParamGuard guard;
	core::World w;

	// 1. Setup Platform with Math Paths
	auto p = std::make_unique<radar::Platform>("p1", 1);
	p->getMotionPath()->setInterp(math::Path::InterpType::INTERP_CUBIC);
	p->getMotionPath()->addCoord({math::Vec3(1.1, 2.2, 3.3), 4.5});
	p->getMotionPath()->addCoord({math::Vec3(2.1, 3.2, 4.3), 5.5});
	p->getMotionPath()->addCoord({math::Vec3(3.1, 4.2, 5.3), 6.5});
	p->getMotionPath()->addCoord({math::Vec3(4.1, 5.2, 6.3), 7.5});
	p->getMotionPath()->finalize();

	math::RotationCoord start(0.0, 0.0, 0.0);
	math::RotationCoord rate(-PI / 2.0, PI / 4.0, 0.0);
	p->getRotationPath()->setInterp(math::RotationPath::InterpType::INTERP_CONSTANT);
	p->getRotationPath()->setConstantRate(start, rate);
	p->getRotationPath()->finalize();
	w.add(std::move(p));

	// 2. Setup Timing
	auto pt = std::make_unique<timing::PrototypeTiming>("TestTiming", 77);
	pt->setFrequency(10e6);
	pt->setSyncOnPulse();
	pt->setFreqOffset(5.0);
	pt->setRandomFreqOffsetStdev(1.2);
	pt->setPhaseOffset(0.5);
	pt->setRandomPhaseOffsetStdev(0.1);
	pt->setAlpha(1.0, 0.8);
	w.add(std::move(pt));

	// Serialize
	json j = serial::world_to_json(w);

	SECTION("Math Paths Serialize Correctly")
	{
		auto& plat_json = j["simulation"]["platforms"][0];
		REQUIRE(j["simulation"]["parameters"]["rotationangleunit"] == "deg");

		REQUIRE(plat_json["motionpath"]["interpolation"] == "cubic");
		REQUIRE_THAT(plat_json["motionpath"]["positionwaypoints"][0]["x"].get<double>(), WithinAbs(1.1, 1e-9));
		REQUIRE_THAT(plat_json["motionpath"]["positionwaypoints"][0]["time"].get<double>(), WithinAbs(4.5, 1e-9));

		REQUIRE(plat_json["fixedrotation"]["interpolation"] == "constant");
		REQUIRE_THAT(plat_json["fixedrotation"]["startazimuth"].get<double>(), WithinAbs(90.0, 1e-9));
		REQUIRE_THAT(plat_json["fixedrotation"]["azimuthrate"].get<double>(), WithinAbs(90.0, 1e-9));
	}

	SECTION("Timing Structures Serialize Correctly")
	{
		auto& timing_json = j["simulation"]["timings"][0];

		REQUIRE(timing_json["id"] == "77");
		REQUIRE(timing_json["name"] == "TestTiming");
		REQUIRE_THAT(timing_json["frequency"].get<double>(), WithinAbs(10e6, 1e-9));
		REQUIRE(timing_json["synconpulse"] == true);
		REQUIRE_THAT(timing_json["freq_offset"].get<double>(), WithinAbs(5.0, 1e-9));
		REQUIRE_THAT(timing_json["random_freq_offset_stdev"].get<double>(), WithinAbs(1.2, 1e-9));
		REQUIRE_THAT(timing_json["noise_entries"][0]["alpha"].get<double>(), WithinAbs(1.0, 1e-9));
	}
}

TEST_CASE("JSON: Serialization of Assets and Radar Components", "[serial][json]")
{
	ParamGuard guard;

	// MUST initialize global simulation parameters to avoid division-by-zero
	// Without this, params::rate() is 0.0, causing the PRF math to evaluate to NaN.
	params::setRate(1000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	core::World w;

	// Assets
	auto cw = std::make_unique<fers_signal::CwSignal>();
	w.add(std::make_unique<fers_signal::RadarSignal>("CwWave", 10.0, 1e9, 1.0, std::move(cw), 10));
	w.add(std::make_unique<antenna::Gaussian>("gauss", 1.5, 2.5, 20));

	// PrototypeTiming AND Timing must both be properly initialized
	auto proto_tim = std::make_unique<timing::PrototypeTiming>("dummy_proto", 104);
	proto_tim->setFrequency(10e6);
	auto tim = std::make_shared<timing::Timing>("dummy_time", 42, 103);
	tim->initializeModel(proto_tim.get());
	w.add(std::move(proto_tim));

	// Platform & Components
	auto p = std::make_unique<radar::Platform>("p1", 100);
	auto tx = std::make_unique<radar::Transmitter>(p.get(), "tx", radar::OperationMode::PULSED_MODE, 101);

	// Set timing FIRST, then PRF, then Schedule, to guarantee internal math has all variables
	tx->setTiming(tim);
	tx->setPrf(1000.0);
	tx->setSchedule({{0.1, 0.5}});

	auto tgt = radar::createIsoTarget(p.get(), "tgt", 10.0, 42, 102);
	tgt->setFluctuationModel(std::make_unique<radar::RcsChiSquare>(tgt->getRngEngine(), 2.5));

	auto rx = std::make_unique<radar::Receiver>(p.get(), "rx", 88, radar::OperationMode::PULSED_MODE, 105);
	rx->setTiming(tim);
	rx->setWindowProperties(1e-4, 1000.0, 1e-5);
	rx->setNoiseTemperature(300.0);
	rx->setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
	rx->setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);

	w.add(std::move(tx));
	w.add(std::move(tgt));
	w.add(std::move(rx));
	w.add(std::move(p));

	json j = serial::world_to_json(w);

	SECTION("Assets Serialize Correctly")
	{
		REQUIRE(j["simulation"]["waveforms"][0]["name"] == "CwWave");
		REQUIRE_THAT(j["simulation"]["waveforms"][0]["power"].get<double>(), WithinAbs(10.0, 1e-9));
		REQUIRE(j["simulation"]["waveforms"][0].contains("cw"));

		REQUIRE(j["simulation"]["antennas"][0]["pattern"] == "gaussian");
		REQUIRE_THAT(j["simulation"]["antennas"][0]["azscale"].get<double>(), WithinAbs(1.5, 1e-9));
	}

	SECTION("Radar Components Serialize Correctly")
	{
		auto& comps = j["simulation"]["platforms"][0]["components"];
		REQUIRE(comps.size() == 3);

		bool found_tx = false, found_tgt = false, found_rx = false;
		for (const auto& c : comps)
		{
			if (c.contains("transmitter"))
			{
				found_tx = true;
				REQUIRE_THAT(c["transmitter"]["pulsed_mode"]["prf"].get<double>(), WithinAbs(1000.0, 1e-9));
				REQUIRE_THAT(c["transmitter"]["schedule"][0]["start"].get<double>(), WithinAbs(0.1, 1e-9));
			}
			if (c.contains("target"))
			{
				found_tgt = true;
				REQUIRE(c["target"]["rcs"]["type"] == "isotropic");
				REQUIRE_THAT(c["target"]["rcs"]["value"].get<double>(), WithinAbs(10.0, 1e-9));
				REQUIRE(c["target"]["model"]["type"] == "chisquare");
				REQUIRE_THAT(c["target"]["model"]["k"].get<double>(), WithinAbs(2.5, 1e-9));
			}
			if (c.contains("receiver"))
			{
				found_rx = true;
				auto& rx_json = c["receiver"];
				REQUIRE_THAT(rx_json["noise_temp"].get<double>(), WithinAbs(300.0, 1e-9));
				REQUIRE(rx_json["nodirect"] == true);
				REQUIRE(rx_json["nopropagationloss"] == true);
				REQUIRE_THAT(rx_json["pulsed_mode"]["window_length"].get<double>(), WithinAbs(1e-4, 1e-9));
			}
		}
		REQUIRE(found_tx);
		REQUIRE(found_tgt);
		REQUIRE(found_rx);
	}
}

TEST_CASE("JSON: Full World Scenario Deserialization", "[serial][json]")
{
	ParamGuard guard;
	core::World world;
	std::mt19937 seeder(42);

	// Safely construct JSON using initializer lists
	json scenario = {
		{"simulation",
		 {{"parameters",
		   {{"starttime", 0.0},
			{"endtime", 1.0},
			{"rate", 1000.0},
			{"rotationangleunit", "rad"},
			{"origin", {{"latitude", 0.0}, {"longitude", 0.0}, {"altitude", 0.0}}},
			{"coordinatesystem", {{"frame", "ENU"}}},
			{"randomseed", 777}}},
		  {"waveforms",
		   json::array(
			   {{{"id", 10}, {"name", "w1"}, {"power", 10.0}, {"carrier_frequency", 1e9}, {"cw", json::object()}}})},
		  {"antennas", json::array({{{"id", 20}, {"name", "a1"}, {"pattern", "isotropic"}}})},
		  {"timings", json::array({{{"id", 30}, {"name", "t1"}, {"frequency", 1e6}}})},
		  {"platforms",
		   json::array({{{"id", 100},
						 {"name", "p1"},
						 {"components",
						  json::array({{{"monostatic",
										 {{"name", "mono1"},
										  {"tx_id", 101},
										  {"rx_id", 102},
										  {"waveform", 10},
										  {"antenna", 20},
										  {"timing", 30},
										  {"cw_mode", json::object()},
										  {"noise_temp", 290.0},
										  {"nodirect", true},
										  {"nopropagationloss", true},
										  {"schedule", json::array({{{"start", 0.1}, {"end", 0.5}}})}}}},
									   {{"target",
										 {{"id", 103},
										  {"name", "tgt1"},
										  {"rcs", {{"type", "isotropic"}, {"value", 1.0}}},
										  {"model", {{"type", "chisquare"}, {"k", 2.0}}}}}}})}}})}}}};

	REQUIRE_NOTHROW(serial::json_to_world(scenario, world, seeder));

	REQUIRE(world.getPlatforms().size() == 1);
	REQUIRE(world.getTransmitters().size() == 1);
	REQUIRE(world.getReceivers().size() == 1);
	REQUIRE(world.getTargets().size() == 1);
	REQUIRE(params::params.random_seed == 777);
	REQUIRE(params::params.rotation_angle_unit == params::RotationAngleUnit::Radians);

	auto tx = world.getTransmitters()[0].get();
	auto rx = world.getReceivers()[0].get();
	REQUIRE(tx->getAttached() == rx);
	REQUIRE(rx->getAttached() == tx);
	REQUIRE(tx->getTiming().get() == rx->getTiming().get());
	REQUIRE(rx->checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT));
	REQUIRE(rx->checkFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS));
	REQUIRE(tx->getSchedule().size() == 1);
}

TEST_CASE("JSON: Full world skips incomplete numeric-placeholder radar components", "[serial][json]")
{
	ParamGuard guard;
	core::World world;
	std::mt19937 seeder(42);

	json scenario = {
		{"simulation",
		 {{"parameters",
		   {{"starttime", 0.0},
			{"endtime", 1.0},
			{"rate", 1000.0},
			{"origin", {{"latitude", 0.0}, {"longitude", 0.0}, {"altitude", 0.0}}},
			{"coordinatesystem", {{"frame", "ENU"}}}}},
		  {"platforms",
		   json::array(
			   {{{"id", 100},
				 {"name", "draft-platform"},
				 {"components",
				  json::array({{{"transmitter",
								 {{"id", 101},
								  {"name", "draft-tx"},
								  {"waveform", 0},
								  {"antenna", 0},
								  {"timing", 0},
								  {"pulsed_mode", {{"prf", 1000.0}}}}}},
							   {{"receiver",
								 {{"id", 102},
								  {"name", "draft-rx"},
								  {"antenna", 0},
								  {"timing", 0},
								  {"pulsed_mode", {{"prf", 1000.0}, {"window_length", 1.0e-5}, {"window_skip", 0.0}}},
								  {"noise_temp", 290.0}}}},
							   {{"monostatic",
								 {{"name", "draft-mono"},
								  {"tx_id", 103},
								  {"rx_id", 104},
								  {"waveform", 0},
								  {"antenna", 0},
								  {"timing", 0},
								  {"pulsed_mode", {{"prf", 1000.0}, {"window_length", 1.0e-5}, {"window_skip", 0.0}}},
								  {"noise_temp", 290.0}}}}})}}})}}}};

	REQUIRE_NOTHROW(serial::json_to_world(scenario, world, seeder));
	REQUIRE(world.getPlatforms().size() == 1);
	REQUIRE(world.findPlatform(100) != nullptr);
	REQUIRE(world.getTransmitters().empty());
	REQUIRE(world.getReceivers().empty());
}

TEST_CASE("JSON: Full world rejects duplicate scenario object names", "[serial][json]")
{
	ParamGuard guard;
	core::World world;
	std::mt19937 seeder(42);

	json scenario = {{"simulation",
					  {{"parameters",
						{{"starttime", 0.0},
						 {"endtime", 1.0},
						 {"rate", 1000.0},
						 {"origin", {{"latitude", 0.0}, {"longitude", 0.0}, {"altitude", 0.0}}},
						 {"coordinatesystem", {{"frame", "ENU"}}}}},
					   {"platforms",
						json::array({{{"id", 100},
									  {"name", "MonoRadar"},
									  {"components",
									   json::array({{{"monostatic",
													  {{"name", "MonoRadar Monostatic"},
													   {"tx_id", 101},
													   {"rx_id", 102},
													   {"waveform", 10},
													   {"antenna", 20},
													   {"timing", 30},
													   {"cw_mode", json::object()}}}}})}},
									 {{"id", 103},
									  {"name", "MonoRadar Copy"},
									  {"components",
									   json::array({{{"monostatic",
													  {{"name", "MonoRadar Monostatic"},
													   {"tx_id", 104},
													   {"rx_id", 105},
													   {"waveform", 11},
													   {"antenna", 21},
													   {"timing", 31},
													   {"cw_mode", json::object()}}}}})}}})}}}};

	REQUIRE_THROWS_WITH(
		serial::json_to_world(scenario, world, seeder),
		ContainsSubstring(
			"Duplicate name 'MonoRadar Monostatic' found for monostatic; previously used by monostatic."));
}

TEST_CASE("JSON: Rotation parsing warns when values look like the opposite unit", "[serial][json]")
{
	ParamGuard guard;
	LogLevelGuard log_guard(logging::Level::WARNING);
	CerrCapture capture;
	core::World world;
	std::mt19937 seeder(42);

	json scenario = {
		{"simulation",
		 {{"name", "Warning Scenario"},
		  {"parameters",
		   {{"starttime", 0.0},
			{"endtime", 1.0},
			{"rate", 1000.0},
			{"rotationangleunit", "rad"},
			{"origin", {{"latitude", 0.0}, {"longitude", 0.0}, {"altitude", 0.0}}},
			{"coordinatesystem", {{"frame", "ENU"}}}}},
		  {"platforms",
		   json::array({{{"id", 100},
						 {"name", "warn-platform"},
						 {"rotationpath",
						  {{"interpolation", "static"},
						   {"rotationwaypoints",
							json::array({{{"time", 0.0}, {"azimuth", 90.0}, {"elevation", 0.0}}})}}}}})}}}};

	REQUIRE_NOTHROW(serial::json_to_world(scenario, world, seeder));
	REQUIRE_THAT(capture.str(), ContainsSubstring("platform 'warn-platform' rotation waypoint 0"));
	REQUIRE_THAT(capture.str(), ContainsSubstring("'azimuth'"));
	REQUIRE_THAT(capture.str(), ContainsSubstring("declared"));
}

TEST_CASE("JSON: Deserialization Error Paths", "[serial][json]")
{
	ParamGuard guard;
	core::World world;
	std::mt19937 seeder(42);

	auto run_bad_scenario = [&](const json& test_comps)
	{
		json scenario = {
			{"simulation",
			 {{"parameters",
			   {{"starttime", 0.0},
				{"endtime", 1.0},
				{"rate", 1000.0},
				{"origin", {{"latitude", 0.0}, {"longitude", 0.0}, {"altitude", 0.0}}},
				{"coordinatesystem", {{"frame", "ENU"}}}}},
			  {"waveforms",
			   json::array(
				   {{{"id", 10}, {"name", "w1"}, {"power", 1.0}, {"carrier_frequency", 1.0}, {"cw", json::object()}}})},
			  {"antennas", json::array({{{"id", 20}, {"name", "a1"}, {"pattern", "isotropic"}}})},
			  {"timings", json::array({{{"id", 30}, {"name", "t1"}, {"frequency", 1e6}}})},
			  {"platforms", json::array({{{"id", 100}, {"name", "p1"}, {"components", test_comps}}})}}}};
		serial::json_to_world(scenario, world, seeder);
	};

	SECTION("Missing mode throws")
	{
		json test_comps = json::array(
			{{{"transmitter", {{"id", 1}, {"name", "tx1"}, {"waveform", 10}, {"antenna", 20}, {"timing", 30}}}}});
		REQUIRE_THROWS_WITH(run_bad_scenario(test_comps),
							ContainsSubstring("must have a 'pulsed_mode', 'cw_mode', or 'fmcw_mode' block"));
	}

	SECTION("Unsupported RCS type throws")
	{
		json test_comps =
			json::array({{{"target", {{"id", 1}, {"name", "bad-target"}, {"rcs", {{"type", "magic"}}}}}}});
		REQUIRE_THROWS_WITH(run_bad_scenario(test_comps), ContainsSubstring("Unsupported target RCS type: magic"));
	}

	SECTION("Unsupported Fluctuation model type throws")
	{
		json test_comps = json::array({{{"target",
										 {{"id", 1},
										  {"name", "bad-target"},
										  {"rcs", {{"type", "isotropic"}, {"value", 1.0}}},
										  {"model", {{"type", "magic"}}}}}}});
		REQUIRE_THROWS_WITH(run_bad_scenario(test_comps),
							ContainsSubstring("Unsupported fluctuation model type: magic"));
	}

	SECTION("Negative ID throws")
	{
		json test_comps = json::array(
			{{{"target", {{"id", -5}, {"name", "bad-target"}, {"rcs", {{"type", "isotropic"}, {"value", 1.0}}}}}}});
		REQUIRE_THROWS_WITH(run_bad_scenario(test_comps), ContainsSubstring("negative id"));
	}
}

TEST_CASE("JSON: FMCW schedule validation matches chirp timing", "[serial][json][fmcw]")
{
	ParamGuard guard;
	std::mt19937 seeder(42);

	const auto make_scenario = [](const json& schedule)
	{
		json scenario;
		scenario["simulation"]["parameters"] = {{"starttime", 0.0},
												{"endtime", 10.0},
												{"rate", 2.0e6},
												{"origin", {{"latitude", 0.0}, {"longitude", 0.0}, {"altitude", 0.0}}},
												{"coordinatesystem", {{"frame", "ENU"}}}};
		scenario["simulation"]["waveforms"] = json::array({{{"id", 10},
															{"name", "fmcw_wave"},
															{"power", 1.0},
															{"carrier_frequency", 1.0e9},
															{"fmcw_linear_chirp",
															 {{"direction", "up"},
															  {"chirp_bandwidth", 1.0e6},
															  {"chirp_duration", 1.0e-3},
															  {"chirp_period", 2.0e-3}}}}});
		scenario["simulation"]["antennas"] = json::array({{{"id", 20}, {"name", "a1"}, {"pattern", "isotropic"}}});
		scenario["simulation"]["timings"] = json::array({{{"id", 30}, {"name", "t1"}, {"frequency", 1.0e6}}});
		scenario["simulation"]["platforms"] = json::array({{{"id", 100},
															{"name", "p1"},
															{"components",
															 json::array({{{"transmitter",
																			{{"id", 101},
																			 {"name", "tx1"},
																			 {"waveform", 10},
																			 {"antenna", 20},
																			 {"timing", 30},
																			 {"fmcw_mode", json::object()},
																			 {"schedule", schedule}}}}})}}});
		return scenario;
	};

	SECTION("period shorter than T_c is rejected")
	{
		core::World world;
		const auto scenario = make_scenario(json::array({{{"start", 0.1}, {"end", 0.1005}}}));
		REQUIRE_THROWS_WITH(serial::json_to_world(scenario, world, seeder),
							ContainsSubstring("shorter than FMCW chirp_duration T_c"));
	}

	SECTION("period shorter than T_rep but at least T_c only warns")
	{
		core::World world;
		LogLevelGuard log_level(logging::Level::WARNING);
		CerrCapture capture;
		const auto scenario = make_scenario(json::array({{{"start", 0.1}, {"end", 0.1015}}}));
		REQUIRE_NOTHROW(serial::json_to_world(scenario, world, seeder));
		REQUIRE(world.getTransmitters().size() == 1);
		REQUIRE_THAT(capture.str(), ContainsSubstring("shorter than FMCW chirp_period"));
	}
}

TEST_CASE("JSON: FMCW dechirp configuration validates and round-trips", "[serial][json][fmcw][dechirp]")
{
	ParamGuard guard;
	std::mt19937 seeder(42);

	const auto make_scenario = [](json fmcw_mode)
	{
		json scenario;
		scenario["simulation"]["parameters"] = {{"starttime", 0.0},
												{"endtime", 1.0e-3},
												{"rate", 4.0e6},
												{"origin", {{"latitude", 0.0}, {"longitude", 0.0}, {"altitude", 0.0}}},
												{"coordinatesystem", {{"frame", "ENU"}}}};
		scenario["simulation"]["waveforms"] = json::array({{{"id", 10},
															{"name", "fmcw_wave"},
															{"power", 1.0},
															{"carrier_frequency", 1.0e9},
															{"fmcw_linear_chirp",
															 {{"direction", "up"},
															  {"chirp_bandwidth", 1.0e6},
															  {"chirp_duration", 1.0e-4},
															  {"chirp_period", 1.0e-4}}}}});
		scenario["simulation"]["antennas"] = json::array({{{"id", 20}, {"name", "a1"}, {"pattern", "isotropic"}}});
		scenario["simulation"]["timings"] = json::array({{{"id", 30}, {"name", "t1"}, {"frequency", 1.0e6}}});
		scenario["simulation"]["platforms"] = json::array({{{"id", 100},
															{"name", "p1"},
															{"components",
															 json::array({{{"monostatic",
																			{{"name", "mono1"},
																			 {"tx_id", 101},
																			 {"rx_id", 102},
																			 {"waveform", 10},
																			 {"antenna", 20},
																			 {"timing", 30},
																			 {"fmcw_mode", fmcw_mode}}}}})}}});
		return scenario;
	};

	SECTION("attached reference is resolved and serialized")
	{
		core::World world;
		const auto scenario =
			make_scenario({{"dechirp_mode", "physical"}, {"dechirp_reference", {{"source", "attached"}}}});

		REQUIRE_NOTHROW(serial::json_to_world(scenario, world, seeder));
		REQUIRE(world.getReceivers().size() == 1);
		const auto* rx = world.getReceivers().front().get();
		REQUIRE(rx->getDechirpMode() == radar::Receiver::DechirpMode::Physical);
		REQUIRE(rx->getDechirpReference().source == radar::Receiver::DechirpReferenceSource::Attached);
		REQUIRE_FALSE(rx->getDechirpSources().empty());

		const json serialized = serial::world_to_json(world);
		const auto& mode_json =
			serialized.at("simulation").at("platforms").at(0).at("components").at(0).at("monostatic").at("fmcw_mode");
		REQUIRE(mode_json.at("dechirp_mode") == "physical");
		REQUIRE(mode_json.at("dechirp_reference").at("source") == "attached");
	}

	SECTION("IF-chain fields are parsed and serialized")
	{
		core::World world;
		const auto scenario = make_scenario({{"dechirp_mode", "physical"},
											 {"dechirp_reference", {{"source", "attached"}}},
											 {"if_sample_rate", 1.0e6},
											 {"if_filter_bandwidth", 4.0e5},
											 {"if_filter_transition_width", 1.0e5}});

		REQUIRE_NOTHROW(serial::json_to_world(scenario, world, seeder));
		const auto* rx = world.getReceivers().front().get();
		const auto& if_chain = rx->getFmcwIfChainRequest();
		REQUIRE(if_chain.sample_rate_hz.has_value());
		REQUIRE(if_chain.filter_bandwidth_hz.has_value());
		REQUIRE(if_chain.filter_transition_width_hz.has_value());
		REQUIRE_THAT(*if_chain.sample_rate_hz, WithinAbs(1.0e6, 1.0e-9));
		REQUIRE_THAT(*if_chain.filter_bandwidth_hz, WithinAbs(4.0e5, 1.0e-9));
		REQUIRE_THAT(*if_chain.filter_transition_width_hz, WithinAbs(1.0e5, 1.0e-9));

		const json serialized = serial::world_to_json(world);
		const auto& mode_json =
			serialized.at("simulation").at("platforms").at(0).at("components").at(0).at("monostatic").at("fmcw_mode");
		REQUIRE(mode_json.at("if_sample_rate") == 1.0e6);
		REQUIRE(mode_json.at("if_filter_bandwidth") == 4.0e5);
		REQUIRE(mode_json.at("if_filter_transition_width") == 1.0e5);
	}

	SECTION("orphan dechirp_reference is rejected")
	{
		core::World world;
		const auto scenario = make_scenario({{"dechirp_reference", {{"source", "attached"}}}});
		REQUIRE_THROWS_WITH(serial::json_to_world(scenario, world, seeder), ContainsSubstring("dechirp_reference"));
	}

	SECTION("IF sample rate requires dechirp")
	{
		core::World world;
		const auto scenario = make_scenario({{"if_sample_rate", 1.0e6}});
		REQUIRE_THROWS_WITH(serial::json_to_world(scenario, world, seeder), ContainsSubstring("IF-chain fields"));
	}

	SECTION("IF-chain values must be positive")
	{
		core::World world;
		const auto scenario = make_scenario(
			{{"dechirp_mode", "physical"}, {"dechirp_reference", {{"source", "attached"}}}, {"if_sample_rate", -1.0}});
		REQUIRE_THROWS_WITH(serial::json_to_world(scenario, world, seeder), ContainsSubstring("finite positive"));
	}

	SECTION("IF filter bandwidth must be below IF Nyquist")
	{
		core::World world;
		const auto scenario = make_scenario({{"dechirp_mode", "physical"},
											 {"dechirp_reference", {{"source", "attached"}}},
											 {"if_sample_rate", 1.0e6},
											 {"if_filter_bandwidth", 5.0e5}});
		REQUIRE_THROWS_WITH(serial::json_to_world(scenario, world, seeder),
							ContainsSubstring("less than half if_sample_rate"));
	}

	SECTION("conflicting mode blocks cannot hide dechirp configuration")
	{
		core::World world;
		auto scenario = make_scenario({{"dechirp_mode", "physical"}, {"dechirp_reference", {{"source", "attached"}}}});
		auto& monostatic = scenario["simulation"]["platforms"][0]["components"][0]["monostatic"];
		monostatic["cw_mode"] = json::object();

		REQUIRE_THROWS_WITH(serial::json_to_world(scenario, world, seeder), ContainsSubstring("at most one"));
	}

	SECTION("unknown dechirp reference keys are rejected")
	{
		core::World world;
		const auto scenario = make_scenario(
			{{"dechirp_mode", "physical"}, {"dechirp_reference", {{"source", "attached"}, {"bogus", true}}}});
		REQUIRE_THROWS_WITH(serial::json_to_world(scenario, world, seeder), ContainsSubstring("unsupported key"));
	}
}

TEST_CASE("JSON: Vec3 Serialization and Deserialization", "[serial][json]")
{
	math::Vec3 v_orig(1.2, 3.4, 5.6);
	json j = v_orig;

	REQUIRE_THAT(j["x"].get<double>(), WithinAbs(1.2, 1e-9));
	REQUIRE_THAT(j["y"].get<double>(), WithinAbs(3.4, 1e-9));
	REQUIRE_THAT(j["z"].get<double>(), WithinAbs(5.6, 1e-9));

	math::Vec3 v_new = j.get<math::Vec3>();
	REQUIRE_THAT(v_new.x, WithinAbs(v_orig.x, 1e-9));
	REQUIRE_THAT(v_new.y, WithinAbs(v_orig.y, 1e-9));
	REQUIRE_THAT(v_new.z, WithinAbs(v_orig.z, 1e-9));
}

TEST_CASE("JSON: Focused PrototypeTiming Serialization", "[serial][json]")
{
	ParamGuard guard;
	core::World w;

	auto pt = std::make_unique<timing::PrototypeTiming>("FocusTiming", 99);
	pt->setFrequency(50e6);
	w.add(std::move(pt));

	json j = serial::world_to_json(w);

	auto& timing_json = j["simulation"]["timings"][0];
	REQUIRE(timing_json["name"] == "FocusTiming");
	REQUIRE_THAT(timing_json["frequency"].get<double>(), WithinAbs(50e6, 1e-9));
}

TEST_CASE("JSON: Monostatic Radar Serialization", "[serial][json]")
{
	ParamGuard guard;
	params::setRate(10000.0); // 10 kHz -> 1e-4s sample period
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	core::World w;

	// Assets
	auto cw = std::make_unique<fers_signal::CwSignal>();
	w.add(std::make_unique<fers_signal::RadarSignal>("MonoWave", 1.0, 1e9, 1.0, std::move(cw), 210));
	w.add(std::make_unique<antenna::Isotropic>("MonoAnt", 220));
	auto proto_tim = std::make_unique<timing::PrototypeTiming>("MonoProto", 230);
	proto_tim->setFrequency(10e6);
	w.add(std::move(proto_tim));

	// Platform & Monostatic Components
	auto p = std::make_unique<radar::Platform>("p_mono", 200);
	auto tx = std::make_unique<radar::Transmitter>(p.get(), "mono", radar::OperationMode::PULSED_MODE, 201);
	auto rx = std::make_unique<radar::Receiver>(p.get(), "mono", 99, radar::OperationMode::PULSED_MODE, 202);

	auto tim_inst = std::make_shared<timing::Timing>("mono_timing_inst", 42, 230);
	tim_inst->initializeModel(w.findTiming(230));

	tx->setTiming(tim_inst);
	tx->setPrf(2000.0);
	tx->setWave(w.findWaveform(210));
	tx->setAntenna(w.findAntenna(220));
	tx->setSchedule({{0.2, 0.8}});

	rx->setTiming(tim_inst);
	rx->setAntenna(w.findAntenna(220));
	// Use window_length and window_skip that are exact multiples of the 1e-4s sample period
	// to prevent the engine from quantizing them down to 0.0.
	rx->setWindowProperties(1e-4, 2000.0, 2e-4); // length, prf, skip
	rx->setNoiseTemperature(295.0);
	rx->setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
	rx->setSchedule({{0.2, 0.8}});

	// Link them
	tx->setAttached(rx.get());
	rx->setAttached(tx.get());

	w.add(std::move(tx));
	w.add(std::move(rx));
	w.add(std::move(p));

	json j = serial::world_to_json(w);

	SECTION("Monostatic component serializes correctly")
	{
		auto& comps = j["simulation"]["platforms"][0]["components"];
		REQUIRE(comps.size() == 1);
		REQUIRE(comps[0].contains("monostatic"));

		auto& mono_json = comps[0]["monostatic"];
		REQUIRE(mono_json["name"] == "mono");
		REQUIRE(mono_json["tx_id"] == "201");
		REQUIRE(mono_json["rx_id"] == "202");
		REQUIRE_THAT(mono_json["noise_temp"].get<double>(), WithinAbs(295.0, 1e-9));
		REQUIRE(mono_json["nodirect"] == true);
		REQUIRE(mono_json["nopropagationloss"] == false); // default

		auto& pulsed_json = mono_json["pulsed_mode"];
		REQUIRE_THAT(pulsed_json["prf"].get<double>(), WithinAbs(2000.0, 1e-9));
		REQUIRE_THAT(pulsed_json["window_length"].get<double>(), WithinAbs(1e-4, 1e-9));
		REQUIRE_THAT(pulsed_json["window_skip"].get<double>(), WithinAbs(2e-4, 1e-9));
	}
}

TEST_CASE("JSON: Granular updates of Radar Components and Timing", "[serial][json]")
{
	ParamGuard guard;
	params::setRate(10000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	core::World w;
	std::mt19937 seeder(42);

	// Setup basic assets
	auto wf = std::make_unique<fers_signal::RadarSignal>("wf1", 10.0, 1e9, 1.0,
														 std::make_unique<fers_signal::CwSignal>(), 10);
	w.add(std::move(wf));

	auto ant = std::make_unique<antenna::Isotropic>("ant1", 20);
	w.add(std::move(ant));

	auto pt = std::make_unique<timing::PrototypeTiming>("tim1", 30);
	pt->setFrequency(1e6);
	w.add(std::move(pt));

	auto p = std::make_unique<radar::Platform>("p1", 100);

	auto tx = std::make_unique<radar::Transmitter>(p.get(), "tx1", radar::OperationMode::CW_MODE, 101);
	auto rx = std::make_unique<radar::Receiver>(p.get(), "rx1", 42, radar::OperationMode::CW_MODE, 102);
	auto tgt = radar::createIsoTarget(p.get(), "tgt1", 1.0, 54321, 103); // Explicit seed 54321

	// Give tx and rx an initial timing object to verify seed preservation
	auto initial_timing = std::make_shared<timing::Timing>("tim1", 12345, 30);
	initial_timing->initializeModel(w.findTiming(30));
	tx->setTiming(initial_timing);
	rx->setTiming(initial_timing);

	auto* tx_ptr = tx.get();
	auto* rx_ptr = rx.get();
	auto* tgt_ptr = tgt.get();

	w.add(std::move(tx));
	w.add(std::move(rx));
	w.add(std::move(tgt));
	w.add(std::move(p));

	SECTION("Update Parameters")
	{
		json j = {{"starttime", 1.0},
				  {"endtime", 20.0},
				  {"rate", 2000.0},
				  {"origin", {{"latitude", 10.0}, {"longitude", 20.0}, {"altitude", 30.0}}},
				  {"coordinatesystem", {{"frame", "ECEF"}}}};
		serial::update_parameters_from_json(j, seeder);

		REQUIRE_THAT(params::startTime(), WithinAbs(1.0, 1e-9));
		REQUIRE_THAT(params::endTime(), WithinAbs(20.0, 1e-9));
		REQUIRE_THAT(params::rate(), WithinAbs(2000.0, 1e-9));
		REQUIRE(params::coordinateFrame() == params::CoordinateFrame::ECEF);
	}

	SECTION("Update Parameters rejects unsupported oversample ratios")
	{
		json j = {{"starttime", 1.0},
				  {"endtime", 20.0},
				  {"rate", 2000.0},
				  {"oversample", 9},
				  {"origin", {{"latitude", 10.0}, {"longitude", 20.0}, {"altitude", 30.0}}},
				  {"coordinatesystem", {{"frame", "ECEF"}}}};

		REQUIRE_THROWS_WITH(serial::update_parameters_from_json(j, seeder),
							ContainsSubstring("Oversampling ratios > 8 are not supported"));
	}

	SECTION("Update Antenna In-Place")
	{
		auto* ant_ptr = w.findAntenna(20);
		json j = {{"id", 20}, {"name", "ant_updated"}, {"pattern", "isotropic"}, {"efficiency", 0.5}};
		serial::update_antenna_from_json(j, ant_ptr, w);

		REQUIRE(ant_ptr->getName() == "ant_updated");
		REQUIRE_THAT(ant_ptr->getEfficiencyFactor(), WithinAbs(0.5, 1e-9));
	}

	SECTION("Update Antenna Replacement")
	{
		auto* ant_ptr = w.findAntenna(20);
		json j = {{"id", 20},	 {"name", "ant_replaced"}, {"pattern", "sinc"}, {"alpha", 1.0}, {"beta", 2.0},
				  {"gamma", 3.0}};
		serial::update_antenna_from_json(j, ant_ptr, w);

		auto* new_ant = w.findAntenna(20);
		REQUIRE(new_ant != nullptr);
		REQUIRE(new_ant != ant_ptr); // Pointer changed
		REQUIRE(new_ant->getName() == "ant_replaced");
		REQUIRE(dynamic_cast<antenna::Sinc*>(new_ant) != nullptr);
	}

	SECTION("Update Antenna Rejects Empty File Replacement")
	{
		auto* ant_ptr = w.findAntenna(20);
		json j = {{"id", 20}, {"name", "draft_h5"}, {"pattern", "file"}, {"filename", ""}};

		REQUIRE_THROWS_WITH(serial::update_antenna_from_json(j, ant_ptr, w), ContainsSubstring("without a filename"));
		REQUIRE(w.findAntenna(20) == ant_ptr);
		REQUIRE(ant_ptr->getName() == "ant1");
		REQUIRE(dynamic_cast<antenna::Isotropic*>(ant_ptr) != nullptr);
	}

	SECTION("Update Transmitter")
	{
		json j = {{"name", "tx_updated"}, {"cw_mode", json::object()},
				  {"waveform", 10},		  {"antenna", 20},
				  {"timing", 30},		  {"schedule", json::array({{{"start", 0.1}, {"end", 0.5}}})}};
		serial::update_transmitter_from_json(j, tx_ptr, w, seeder);

		REQUIRE(tx_ptr->getName() == "tx_updated");
		REQUIRE(tx_ptr->getMode() == radar::OperationMode::CW_MODE);
		REQUIRE(tx_ptr->getSignal() != nullptr);
		REQUIRE(tx_ptr->getAntenna() != nullptr);
		REQUIRE(tx_ptr->getTiming() != nullptr);
		REQUIRE(tx_ptr->getTiming()->getSeed() == 12345); // Seed preserved
		REQUIRE(tx_ptr->getSchedule().size() == 1);
	}

	SECTION("Update Transmitter PRF")
	{
		json j = {{"pulsed_mode", {{"prf", 1250.0}}}};
		serial::update_transmitter_from_json(j, tx_ptr, w, seeder);

		REQUIRE(tx_ptr->getMode() == radar::OperationMode::PULSED_MODE);
		REQUIRE_THAT(tx_ptr->getPrf(), WithinAbs(1250.0, 1e-9));
	}

	SECTION("Update Receiver")
	{
		json j = {{"name", "rx_updated"},
				  {"pulsed_mode", {{"prf", 1250.0}, {"window_length", 1e-4}, {"window_skip", 1e-5}}},
				  {"noise_temp", 400.0},
				  {"nodirect", true},
				  {"nopropagationloss", true},
				  {"antenna", 20},
				  {"timing", 30}};
		serial::update_receiver_from_json(j, rx_ptr, w, seeder);

		REQUIRE(rx_ptr->getName() == "rx_updated");
		REQUIRE(rx_ptr->getMode() == radar::OperationMode::PULSED_MODE);
		REQUIRE_THAT(rx_ptr->getWindowPrf(), WithinAbs(1250.0, 1e-9));
		REQUIRE_THAT(rx_ptr->getNoiseTemperature(), WithinAbs(400.0, 1e-9));
		REQUIRE(rx_ptr->checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT));
		REQUIRE(rx_ptr->checkFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS));
		REQUIRE(rx_ptr->getTiming()->getSeed() == 12345); // Seed preserved
	}

	SECTION("Update Target")
	{
		json j = {{"name", "tgt_updated"},
				  {"rcs", {{"type", "isotropic"}, {"value", 50.0}}},
				  {"model", {{"type", "chisquare"}, {"k", 3.0}}}};
		serial::update_target_from_json(j, tgt_ptr, w, seeder);

		// Note: update_target_from_json replaces the target in the world.
		// We need to find it again.
		auto* new_tgt = w.findTarget(103);
		REQUIRE(new_tgt != nullptr);
		REQUIRE(new_tgt->getName() == "tgt_updated");
		REQUIRE(new_tgt->getSeed() == 54321); // Seed preserved
		auto* iso_tgt = dynamic_cast<radar::IsoTarget*>(new_tgt);
		REQUIRE(iso_tgt != nullptr);
		REQUIRE_THAT(iso_tgt->getConstRcs(), WithinAbs(50.0, 1e-9));
		REQUIRE(new_tgt->getFluctuationModel() != nullptr);
	}

	SECTION("Update Timing")
	{
		auto* pt_ptr = w.findTiming(30);
		json j = {{"name", "tim_updated"},
				  {"frequency", 2e6},
				  {"synconpulse", true},
				  {"freq_offset", 10.0},
				  {"random_freq_offset_stdev", 2.0},
				  {"phase_offset", 0.5},
				  {"random_phase_offset_stdev", 0.1},
				  {"noise_entries", json::array({{{"alpha", 1.0}, {"weight", 0.5}}})}};
		serial::update_timing_from_json(j, w, 30);

		auto* updated_timing = w.findTiming(30);
		REQUIRE(updated_timing != nullptr);
		REQUIRE(updated_timing != pt_ptr);

		REQUIRE(updated_timing->getName() == "tim_updated");
		REQUIRE_THAT(updated_timing->getFrequency(), WithinAbs(2e6, 1e-9));
		REQUIRE(updated_timing->getSyncOnPulse() == true);
		REQUIRE_THAT(updated_timing->getFreqOffset().value(), WithinAbs(10.0, 1e-9));
		REQUIRE_THAT(updated_timing->getRandomFreqOffsetStdev().value(), WithinAbs(2.0, 1e-9));
		REQUIRE_THAT(updated_timing->getPhaseOffset().value(), WithinAbs(0.5, 1e-9));
		REQUIRE_THAT(updated_timing->getRandomPhaseOffsetStdev().value(), WithinAbs(0.1, 1e-9));

		std::vector<RealType> alphas, weights;
		updated_timing->copyAlphas(alphas, weights);
		REQUIRE(alphas.size() == 1);
		REQUIRE_THAT(alphas[0], WithinAbs(1.0, 1e-9));
		REQUIRE(tx_ptr->getTiming().get() != initial_timing.get());
		REQUIRE(rx_ptr->getTiming().get() != initial_timing.get());
		REQUIRE(tx_ptr->getTiming()->getSeed() == 12345);
		REQUIRE(rx_ptr->getTiming()->getSeed() == 12345);
		REQUIRE(tx_ptr->getTiming()->getName() == "tim_updated");
		REQUIRE(rx_ptr->getTiming()->getName() == "tim_updated");
		REQUIRE_THAT(tx_ptr->getTiming()->getFrequency(), WithinAbs(2e6, 1e-9));
		REQUIRE_THAT(rx_ptr->getTiming()->getFrequency(), WithinAbs(2e6, 1e-9));
		REQUIRE(tx_ptr->getTiming()->getSyncOnPulse() == true);
		REQUIRE(rx_ptr->getTiming()->getSyncOnPulse() == true);
		REQUIRE_THAT(tx_ptr->getTiming()->getFreqOffset(), WithinAbs(rx_ptr->getTiming()->getFreqOffset(), 1e-9));
		REQUIRE_THAT(tx_ptr->getTiming()->getPhaseOffset(), WithinAbs(rx_ptr->getTiming()->getPhaseOffset(), 1e-9));
	}
}

TEST_CASE("JSON: Granular updates of Monostatic Radar", "[serial][json]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	core::World w;
	std::mt19937 seeder(42);

	auto pt = std::make_unique<timing::PrototypeTiming>("tim1", 30);
	pt->setFrequency(1e6);
	w.add(std::move(pt));

	auto ant = std::make_unique<antenna::Isotropic>("ant1", 20);
	w.add(std::move(ant));

	auto wf = std::make_unique<fers_signal::RadarSignal>("wf1", 10.0, 1e9, 1.0,
														 std::make_unique<fers_signal::CwSignal>(), 10);
	w.add(std::move(wf));
	auto fmcw_wf = std::make_unique<fers_signal::RadarSignal>(
		"wf_fmcw", 10.0, 1e9, 0.1, std::make_unique<fers_signal::FmcwChirpSignal>(100.0, 0.1, 0.2), 11);
	w.add(std::move(fmcw_wf));

	auto p = std::make_unique<radar::Platform>("p1", 100);
	auto tx = std::make_unique<radar::Transmitter>(p.get(), "mono_tx", radar::OperationMode::CW_MODE, 101);
	auto rx = std::make_unique<radar::Receiver>(p.get(), "mono_rx", 42, radar::OperationMode::CW_MODE, 102);

	tx->setAttached(rx.get());
	rx->setAttached(tx.get());

	// Initialize with existing timing
	auto initial_timing = std::make_shared<timing::Timing>("tim1", 12345, 30);
	initial_timing->initializeModel(w.findTiming(30));
	tx->setTiming(initial_timing);
	rx->setTiming(initial_timing);

	auto* tx_ptr = tx.get();
	auto* rx_ptr = rx.get();

	w.add(std::move(tx));
	w.add(std::move(rx));
	w.add(std::move(p));

	json j = {{"name", "mono_updated"},	   {"tx_id", 101},		  {"rx_id", 102},	  {"waveform", 10}, {"antenna", 20},
			  {"cw_mode", json::object()}, {"noise_temp", 300.0}, {"nodirect", true}, {"timing", 30}};

	serial::update_monostatic_from_json(j, tx_ptr, rx_ptr, w, seeder);

	REQUIRE(tx_ptr->getName() == "mono_updated");
	REQUIRE(rx_ptr->getName() == "mono_updated");

	REQUIRE(tx_ptr->getMode() == radar::OperationMode::CW_MODE);
	REQUIRE(rx_ptr->getMode() == radar::OperationMode::CW_MODE);
	REQUIRE_THAT(rx_ptr->getNoiseTemperature(), WithinAbs(300.0, 1e-9));
	REQUIRE(rx_ptr->checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT));
	REQUIRE(tx_ptr->getTiming()->getSeed() == 12345);
	REQUIRE(rx_ptr->getTiming()->getSeed() == 12345);
	REQUIRE(tx_ptr->getTiming().get() == rx_ptr->getTiming().get());

	json fmcw_update = {{"name", "mono_fmcw"},
						{"tx_id", 101},
						{"rx_id", 102},
						{"waveform", 11},
						{"antenna", 20},
						{"fmcw_mode",
						 {{"dechirp_mode", "physical"},
						  {"dechirp_reference", {{"source", "attached"}}},
						  {"if_sample_rate", 100.0},
						  {"if_filter_bandwidth", 40.0},
						  {"if_filter_transition_width", 10.0}}},
						{"timing", 30}};

	REQUIRE_NOTHROW(serial::update_monostatic_from_json(fmcw_update, tx_ptr, rx_ptr, w, seeder));
	REQUIRE(tx_ptr->getMode() == radar::OperationMode::FMCW_MODE);
	REQUIRE(rx_ptr->getMode() == radar::OperationMode::FMCW_MODE);
	REQUIRE(rx_ptr->getDechirpMode() == radar::Receiver::DechirpMode::Physical);
	REQUIRE_THAT(*rx_ptr->getIfSampleRate(), WithinAbs(100.0, 1e-12));
	REQUIRE_THAT(*rx_ptr->getIfFilterBandwidth(), WithinAbs(40.0, 1e-12));
	REQUIRE_THAT(*rx_ptr->getIfFilterTransitionWidth(), WithinAbs(10.0, 1e-12));
}
