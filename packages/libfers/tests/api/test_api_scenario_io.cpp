// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <filesystem>
#include <libfers/api.h>

#include "api_test_helpers.h"

using api_test::json;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

TEST_CASE("API loading scenario from XML string rejects null arguments", "[api][scenario]")
{
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);
	const std::string xml = api_test::minimalScenarioXml();

	SECTION("null context")
	{
		REQUIRE(fers_load_scenario_from_xml_string(nullptr, xml.c_str(), 0) == -1);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("context or xml_content is NULL"));
	}

	SECTION("null xml")
	{
		REQUIRE(fers_load_scenario_from_xml_string(context.get(), nullptr, 0) == -1);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("context or xml_content is NULL"));
	}
}

TEST_CASE("API set output directory rejects null arguments", "[api][scenario]")
{
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	SECTION("null context")
	{
		REQUIRE(fers_set_output_directory(nullptr, "/tmp") == -1);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("context or out_dir is NULL"));
	}

	SECTION("null out_dir")
	{
		REQUIRE(fers_set_output_directory(context.get(), nullptr) == -1);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("context or out_dir is NULL"));
	}
}

TEST_CASE("API set output directory succeeds with valid arguments", "[api][scenario]")
{
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	REQUIRE(fers_set_output_directory(context.get(), "/tmp/custom_dir") == 0);
	api_test::ApiString error = api_test::lastError();
	REQUIRE(error.get() == nullptr);
}

TEST_CASE("API loading scenario from XML file rejects null arguments", "[api][scenario]")
{
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);
	const auto xml_path = api_test::uniqueTempPath("api_minimal_scenario", ".xml");
	api_test::ScopedPath xml_guard(xml_path);
	api_test::writeTextFile(xml_path, api_test::minimalScenarioXml());
	const std::string xml_path_string = api_test::pathString(xml_path);

	SECTION("null context")
	{
		REQUIRE(fers_load_scenario_from_xml_file(nullptr, xml_path_string.c_str(), 0) == -1);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("context or xml_filepath is NULL"));
	}

	SECTION("null path")
	{
		REQUIRE(fers_load_scenario_from_xml_file(context.get(), nullptr, 0) == -1);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("context or xml_filepath is NULL"));
	}
}

TEST_CASE("API loads a minimal scenario from XML string and exports JSON", "[api][scenario]")
{
	api_test::ParamGuard guard;
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::minimalScenarioXml("XML String Scenario");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	api_test::ApiString json_text = api_test::scenarioAsJson(context.get());
	REQUIRE(json_text.get() != nullptr);

	const json scenario = json::parse(json_text.str());
	REQUIRE(scenario.contains("simulation"));
	REQUIRE(scenario.at("simulation").is_object());
	REQUIRE(scenario.at("simulation").contains("name"));
}

TEST_CASE("API loads a minimal scenario from XML file and exports XML", "[api][scenario]")
{
	api_test::ParamGuard guard;
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const auto xml_path = api_test::uniqueTempPath("api_minimal_scenario", ".xml");
	api_test::ScopedPath xml_guard(xml_path);
	api_test::writeTextFile(xml_path, api_test::minimalScenarioXml("XML File Scenario"));
	const std::string xml_path_string = api_test::pathString(xml_path);

	REQUIRE(fers_load_scenario_from_xml_file(context.get(), xml_path_string.c_str(), 0) == 0);

	api_test::ApiString xml_text = api_test::scenarioAsXml(context.get());
	REQUIRE(xml_text.get() != nullptr);
	REQUIRE_FALSE(xml_text.view().empty());
	REQUIRE_THAT(xml_text.str(), ContainsSubstring("<simulation"));
}

TEST_CASE("API exposes memory projection JSON", "[api][scenario][memory_projection]")
{
	api_test::ParamGuard guard;
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::previewScenarioXml("Memory Projection Scenario");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	api_test::ApiString projection_text(fers_get_memory_projection_json(context.get()));
	REQUIRE(projection_text.get() != nullptr);

	const json projection = json::parse(projection_text.str());
	REQUIRE(projection["phase_noise_lookups"].contains("bytes"));
	REQUIRE(projection["streaming_iq_buffers"]["bytes"].get<std::uint64_t>() > 0);
	REQUIRE(projection["rendered_hdf5_dataset_payload"]["bytes"].get<std::uint64_t>() > 0);
	REQUIRE(projection["resident_baseline"].contains("bytes"));
	REQUIRE(projection["projected_total_footprint"].contains("bytes"));

	REQUIRE(fers_get_memory_projection_json(nullptr) == nullptr);
	api_test::ApiString error = api_test::lastError();
	REQUIRE(error.get() != nullptr);
	REQUIRE_THAT(error.str(), ContainsSubstring("Invalid context provided to fers_get_memory_projection_json"));
}

TEST_CASE("API scenario serialization and update helpers reject null arguments", "[api][scenario]")
{
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	SECTION("json export rejects null context")
	{
		REQUIRE(fers_get_scenario_as_json(nullptr) == nullptr);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("Invalid context provided to fers_get_scenario_as_json"));
	}

	SECTION("xml export rejects null context")
	{
		REQUIRE(fers_get_scenario_as_xml(nullptr) == nullptr);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("Invalid context provided to fers_get_scenario_as_xml"));
	}

	SECTION("json update rejects null context")
	{
		REQUIRE(fers_update_scenario_from_json(nullptr, "{}") == -1);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("context or scenario_json is NULL"));
	}

	SECTION("json update rejects null content")
	{
		REQUIRE(fers_update_scenario_from_json(context.get(), nullptr) == -1);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("context or scenario_json is NULL"));
	}

	SECTION("kml generation rejects null context")
	{
		REQUIRE(fers_generate_kml(nullptr, "out.kml") == -1);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("context or output_kml_filepath is NULL"));
	}

	SECTION("kml generation rejects null path")
	{
		REQUIRE(fers_generate_kml(context.get(), nullptr) == -1);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("context or output_kml_filepath is NULL"));
	}
}

TEST_CASE("API update scenario from JSON modifies serialized state", "[api][scenario]")
{
	api_test::ParamGuard guard;
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::minimalScenarioXml("Before Update");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	auto scenario = api_test::parseScenarioJson(context.get());
	scenario["simulation"]["name"] = "After Update";

	REQUIRE(fers_update_scenario_from_json(context.get(), scenario.dump().c_str()) == 0);

	const auto updated = api_test::parseScenarioJson(context.get());
	REQUIRE(updated.at("simulation").at("name") == "After Update");
}

TEST_CASE("API update scenario from malformed JSON returns JSON-specific error code", "[api][scenario]")
{
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	REQUIRE(fers_update_scenario_from_json(context.get(), api_test::malformedJson().c_str()) == 2);

	api_test::ApiString error = api_test::lastError();
	REQUIRE(error.get() != nullptr);
	REQUIRE_THAT(error.str(), ContainsSubstring("JSON parsing/deserialization error:"));
}

TEST_CASE("API generate KML writes a file for a loaded scenario", "[api][scenario]")
{
	api_test::ParamGuard guard;
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::minimalScenarioXml("KML Scenario");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	const auto kml_path = api_test::uniqueTempPath("api_scenario", ".kml");
	api_test::ScopedPath kml_guard(kml_path);
	const std::string kml_path_string = api_test::pathString(kml_path);
	REQUIRE(fers_generate_kml(context.get(), kml_path_string.c_str()) == 0);
	REQUIRE(std::filesystem::exists(kml_path));
	REQUIRE(std::filesystem::file_size(kml_path) > 0);
}

TEST_CASE("API missing XML file returns a populated last error", "[api][scenario]")
{
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const auto missing_path = api_test::uniqueTempPath("api_missing", ".xml");
	api_test::ScopedPath missing_guard(missing_path);
	const std::string missing_path_string = api_test::pathString(missing_path);

	REQUIRE(fers_load_scenario_from_xml_file(context.get(), missing_path_string.c_str(), 0) == 1);

	api_test::ApiString error = api_test::lastError();
	REQUIRE(error.get() != nullptr);
	REQUIRE_THAT(error.str(), ContainsSubstring("Failed to load main XML file"));
}

TEST_CASE("API malformed XML string returns a populated last error", "[api][scenario]")
{
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	REQUIRE(fers_load_scenario_from_xml_string(context.get(), api_test::malformedXml().c_str(), 0) == 1);

	api_test::ApiString error = api_test::lastError();
	REQUIRE(error.get() != nullptr);
	REQUIRE_THAT(error.str(), ContainsSubstring("Failed to parse XML from memory string"));
}

TEST_CASE("API KML generation reports non-creatable output paths", "[api][scenario]")
{
	api_test::ParamGuard guard;
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::minimalScenarioXml("Bad KML Path");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	const auto missing_parent = api_test::uniqueTempPath("api_missing_kml_dir");
	api_test::ScopedPath missing_guard(missing_parent);
	const auto kml_path = missing_parent / "out.kml";
	const std::string kml_path_string = api_test::pathString(kml_path);

	REQUIRE(fers_generate_kml(context.get(), kml_path_string.c_str()) == 2);

	api_test::ApiString error = api_test::lastError();
	REQUIRE(error.get() != nullptr);
	REQUIRE_THAT(error.str(), ContainsSubstring("KML generation failed"));
}

TEST_CASE("API granular updates modify specific objects", "[api][scenario]")
{
	api_test::ParamGuard guard;
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::minimalScenarioXml("Granular Update Test");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	auto scenario = api_test::parseScenarioJson(context.get());

	SECTION("Parameters update")
	{
		json params_json = {{"starttime", 2.0},
							{"endtime", 50.0},
							{"rate", 5000.0},
							{"origin", {{"latitude", 0.0}, {"longitude", 0.0}, {"altitude", 0.0}}},
							{"coordinatesystem", {{"frame", "ENU"}}}};
		REQUIRE(fers_update_parameters_from_json(context.get(), params_json.dump().c_str()) == 0);

		auto updated = api_test::parseScenarioJson(context.get());
		REQUIRE_THAT(updated["simulation"]["parameters"]["starttime"].get<double>(), WithinAbs(2.0, 1e-9));
		REQUIRE_THAT(updated["simulation"]["parameters"]["rate"].get<double>(), WithinAbs(5000.0, 1e-9));
	}

	SECTION("Platform update")
	{
		auto plat = scenario["simulation"]["platforms"][0];
		uint64_t plat_id = api_test::parseId(plat["id"]);
		plat["name"] = "UpdatedPlatform";

		REQUIRE(fers_update_platform_from_json(context.get(), plat_id, plat.dump().c_str()) == 0);

		auto updated = api_test::parseScenarioJson(context.get());
		REQUIRE(updated["simulation"]["platforms"][0]["name"] == "UpdatedPlatform");
	}

	SECTION("Antenna update")
	{
		auto ant = scenario["simulation"]["antennas"][0];
		ant["name"] = "UpdatedAntenna";
		ant["efficiency"] = 0.5;

		REQUIRE(fers_update_antenna_from_json(context.get(), ant.dump().c_str()) == 0);

		auto updated = api_test::parseScenarioJson(context.get());
		REQUIRE(updated["simulation"]["antennas"][0]["name"] == "UpdatedAntenna");
		REQUIRE_THAT(updated["simulation"]["antennas"][0]["efficiency"].get<double>(), WithinAbs(0.5, 1e-9));
	}

	SECTION("Waveform update")
	{
		auto wf = scenario["simulation"]["waveforms"][0];
		wf["name"] = "UpdatedWaveform";
		wf["power"] = 999.0;

		REQUIRE(fers_update_waveform_from_json(context.get(), wf.dump().c_str()) == 0);

		auto updated = api_test::parseScenarioJson(context.get());
		REQUIRE(updated["simulation"]["waveforms"][0]["name"] == "UpdatedWaveform");
		REQUIRE_THAT(updated["simulation"]["waveforms"][0]["power"].get<double>(), WithinAbs(999.0, 1e-9));
	}
}

TEST_CASE("API granular updates modify specific objects (Preview Scenario)", "[api][scenario]")
{
	api_test::ParamGuard guard;
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::previewScenarioXml("Granular Update Preview");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	auto scenario = api_test::parseScenarioJson(context.get());

	SECTION("Transmitter update")
	{
		auto tx = scenario["simulation"]["platforms"][0]["components"][0]["transmitter"];
		uint64_t tx_id = api_test::parseId(tx["id"]);
		tx["name"] = "UpdatedTx";
		tx["cw_mode"] = json::object();
		tx.erase("pulsed_mode");

		REQUIRE(fers_update_transmitter_from_json(context.get(), tx_id, tx.dump().c_str()) == 0);

		auto updated = api_test::parseScenarioJson(context.get());
		auto updated_tx = updated["simulation"]["platforms"][0]["components"][0]["transmitter"];
		REQUIRE(updated_tx["name"] == "UpdatedTx");
		REQUIRE(updated_tx["cw_mode"].is_object());
	}

	SECTION("Receiver update")
	{
		auto rx = scenario["simulation"]["platforms"][1]["components"][0]["receiver"];
		uint64_t rx_id = api_test::parseId(rx["id"]);
		rx["name"] = "UpdatedRx";
		rx["noise_temp"] = 400.0;

		REQUIRE(fers_update_receiver_from_json(context.get(), rx_id, rx.dump().c_str()) == 0);

		auto updated = api_test::parseScenarioJson(context.get());
		auto updated_rx = updated["simulation"]["platforms"][1]["components"][0]["receiver"];
		REQUIRE(updated_rx["name"] == "UpdatedRx");
		REQUIRE_THAT(updated_rx["noise_temp"].get<double>(), WithinAbs(400.0, 1e-9));
	}

	SECTION("Target update")
	{
		auto tgt = scenario["simulation"]["platforms"][2]["components"][0]["target"];
		uint64_t tgt_id = api_test::parseId(tgt["id"]);
		tgt["name"] = "UpdatedTgt";
		tgt["rcs"]["value"] = 99.0;

		REQUIRE(fers_update_target_from_json(context.get(), tgt_id, tgt.dump().c_str()) == 0);

		auto updated = api_test::parseScenarioJson(context.get());
		auto updated_tgt = updated["simulation"]["platforms"][2]["components"][0]["target"];
		REQUIRE(updated_tgt["name"] == "UpdatedTgt");
		REQUIRE_THAT(updated_tgt["rcs"]["value"].get<double>(), WithinAbs(99.0, 1e-9));
	}

	SECTION("Timing update")
	{
		auto timing = scenario["simulation"]["timings"][0];
		uint64_t timing_id = api_test::parseId(timing["id"]);
		timing["name"] = "UpdatedTiming";
		timing["frequency"] = 2e6;

		REQUIRE(fers_update_timing_from_json(context.get(), timing_id, timing.dump().c_str()) == 0);

		auto updated = api_test::parseScenarioJson(context.get());
		auto updated_timing = updated["simulation"]["timings"][0];
		REQUIRE(updated_timing["name"] == "UpdatedTiming");
		REQUIRE_THAT(updated_timing["frequency"].get<double>(), WithinAbs(2e6, 1e-9));
	}
}

TEST_CASE("API granular updates modify specific objects (Monostatic Scenario)", "[api][scenario]")
{
	api_test::ParamGuard guard;
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::monostaticPreviewScenarioXml("Granular Update Monostatic");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	auto scenario = api_test::parseScenarioJson(context.get());

	SECTION("Monostatic update")
	{
		auto mono = scenario["simulation"]["platforms"][0]["components"][0]["monostatic"];
		mono["name"] = "UpdatedMono";
		mono["noise_temp"] = 500.0;

		REQUIRE(fers_update_monostatic_from_json(context.get(), mono.dump().c_str()) == 0);

		auto updated = api_test::parseScenarioJson(context.get());
		auto updated_mono = updated["simulation"]["platforms"][0]["components"][0]["monostatic"];
		REQUIRE(updated_mono["name"] == "UpdatedMono");
		REQUIRE_THAT(updated_mono["noise_temp"].get<double>(), WithinAbs(500.0, 1e-9));
	}
}

TEST_CASE("API granular updates handle errors gracefully", "[api][scenario]")
{
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);
	const std::string xml = api_test::minimalScenarioXml();
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	SECTION("Null context or JSON")
	{
		REQUIRE(fers_update_parameters_from_json(nullptr, "{}") == -1);
		REQUIRE(fers_update_parameters_from_json(context.get(), nullptr) == -1);
		REQUIRE(fers_update_platform_from_json(nullptr, 1, "{}") == -1);
		REQUIRE(fers_update_platform_from_json(context.get(), 1, nullptr) == -1);
		REQUIRE(fers_update_antenna_from_json(nullptr, "{}") == -1);
		REQUIRE(fers_update_antenna_from_json(context.get(), nullptr) == -1);
		REQUIRE(fers_update_waveform_from_json(nullptr, "{}") == -1);
		REQUIRE(fers_update_waveform_from_json(context.get(), nullptr) == -1);
		REQUIRE(fers_update_transmitter_from_json(nullptr, 1, "{}") == -1);
		REQUIRE(fers_update_transmitter_from_json(context.get(), 1, nullptr) == -1);
		REQUIRE(fers_update_receiver_from_json(nullptr, 1, "{}") == -1);
		REQUIRE(fers_update_receiver_from_json(context.get(), 1, nullptr) == -1);
		REQUIRE(fers_update_target_from_json(nullptr, 1, "{}") == -1);
		REQUIRE(fers_update_target_from_json(context.get(), 1, nullptr) == -1);
		REQUIRE(fers_update_monostatic_from_json(nullptr, "{}") == -1);
		REQUIRE(fers_update_monostatic_from_json(context.get(), nullptr) == -1);
		REQUIRE(fers_update_timing_from_json(nullptr, 1, "{}") == -1);
		REQUIRE(fers_update_timing_from_json(context.get(), 1, nullptr) == -1);
	}

	SECTION("Invalid JSON")
	{
		REQUIRE(fers_update_parameters_from_json(context.get(), "{bad") == 1);
		REQUIRE(fers_update_platform_from_json(context.get(), 1, "{bad") == 1);
		REQUIRE(fers_update_antenna_from_json(context.get(), "{bad") == 1);
		REQUIRE(fers_update_waveform_from_json(context.get(), "{bad") == 1);
		REQUIRE(fers_update_transmitter_from_json(context.get(), 1, "{bad") == 1);
		REQUIRE(fers_update_receiver_from_json(context.get(), 1, "{bad") == 1);
		REQUIRE(fers_update_target_from_json(context.get(), 1, "{bad") == 1);
		REQUIRE(fers_update_monostatic_from_json(context.get(), "{bad") == 1);
		REQUIRE(fers_update_timing_from_json(context.get(), 1, "{bad") == 1);
	}

	SECTION("Object not found")
	{
		REQUIRE(fers_update_platform_from_json(context.get(), 999999, "{}") == 1);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("Platform not found"));

		REQUIRE(fers_update_transmitter_from_json(context.get(), 999999, "{}") == 1);
		REQUIRE(fers_update_receiver_from_json(context.get(), 999999, "{}") == 1);
		REQUIRE(fers_update_target_from_json(context.get(), 999999, "{}") == 1);
		REQUIRE(fers_update_monostatic_from_json(context.get(), "{\"tx_id\": 999999, \"rx_id\": 888888}") == 1);
		REQUIRE(fers_update_timing_from_json(context.get(), 999999, "{}") == 1);
	}
}
