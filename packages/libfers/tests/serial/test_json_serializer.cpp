// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <nlohmann/json.hpp>
#include <random>

#include "antenna/antenna_factory.h"
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

namespace
{
	// Protects the global `params::params` from bleeding state across tests
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
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
			{"origin", {{"latitude", 0.0}, {"longitude", 0.0}, {"altitude", 0.0}}},
			{"coordinatesystem", {{"frame", "ENU"}}},
			{"randomseed", 777}}},
		  {"waveforms",
		   json::array(
			   {{{"id", 10}, {"name", "w1"}, {"power", 10.0}, {"carrier_frequency", 1e9}, {"cw", json::object()}}})},
		  {"antennas", json::array({{{"id", 20}, {"name", "a1"}, {"pattern", "isotropic"}}})},
		  {"timings", json::array({{{"id", 30}, {"name", "t1"}, {"frequency", 1e6}}})},
		  {"platforms",
		   json::array(
			   {{{"id", 100},
				 {"name", "p1"},
				 {"components",
				  json::array({{{"monostatic",
								 {{"name", "mono1"},
								  {"tx_id", 101},
								  {"rx_id", 102},
								  {"waveform", 10},
								  {"antenna", 20},
								  {"timing", 30},
								  {"pulsed_mode", {{"prf", 1000.0}, {"window_length", 1e-4}, {"window_skip", 0.0}}},
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

	auto tx = world.getTransmitters()[0].get();
	auto rx = world.getReceivers()[0].get();
	REQUIRE(tx->getAttached() == rx);
	REQUIRE(rx->getAttached() == tx);
	REQUIRE(rx->checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT));
	REQUIRE(rx->checkFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS));
	REQUIRE(tx->getSchedule().size() == 1);
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
							ContainsSubstring("must have a 'pulsed_mode' or 'cw_mode' block"));
	}

	SECTION("Unsupported RCS type throws")
	{
		json test_comps = json::array({{{"target", {{"id", 1}, {"name", "t1"}, {"rcs", {{"type", "magic"}}}}}}});
		REQUIRE_THROWS_WITH(run_bad_scenario(test_comps), ContainsSubstring("Unsupported target RCS type: magic"));
	}

	SECTION("Unsupported Fluctuation model type throws")
	{
		json test_comps = json::array({{{"target",
										 {{"id", 1},
										  {"name", "t1"},
										  {"rcs", {{"type", "isotropic"}, {"value", 1.0}}},
										  {"model", {{"type", "magic"}}}}}}});
		REQUIRE_THROWS_WITH(run_bad_scenario(test_comps),
							ContainsSubstring("Unsupported fluctuation model type: magic"));
	}

	SECTION("Negative ID throws")
	{
		json test_comps =
			json::array({{{"target", {{"id", -5}, {"name", "t1"}, {"rcs", {{"type", "isotropic"}, {"value", 1.0}}}}}}});
		REQUIRE_THROWS_WITH(run_bad_scenario(test_comps), ContainsSubstring("negative id"));
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
