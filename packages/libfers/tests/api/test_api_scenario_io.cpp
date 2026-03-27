#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <libfers/api.h>

#include "api_test_helpers.h"

using Catch::Matchers::ContainsSubstring;

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

TEST_CASE("API loading scenario from XML file rejects null arguments", "[api][scenario]")
{
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);
	const auto xml_path = api_test::uniqueTempPath("api_minimal_scenario", ".xml");
	api_test::ScopedPath xml_guard(xml_path);
	api_test::writeTextFile(xml_path, api_test::minimalScenarioXml());

	SECTION("null context")
	{
		REQUIRE(fers_load_scenario_from_xml_file(nullptr, xml_path.c_str(), 0) == -1);
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

	const api_test::json scenario = api_test::json::parse(json_text.str());
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

	REQUIRE(fers_load_scenario_from_xml_file(context.get(), xml_path.c_str(), 0) == 0);

	api_test::ApiString xml_text = api_test::scenarioAsXml(context.get());
	REQUIRE(xml_text.get() != nullptr);
	REQUIRE_FALSE(xml_text.view().empty());
	REQUIRE_THAT(xml_text.str(), ContainsSubstring("<simulation"));
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
	REQUIRE(fers_generate_kml(context.get(), kml_path.c_str()) == 0);
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

	REQUIRE(fers_load_scenario_from_xml_file(context.get(), missing_path.c_str(), 0) == 1);

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

	REQUIRE(fers_generate_kml(context.get(), kml_path.c_str()) == 2);

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
		REQUIRE_THAT(updated["simulation"]["antennas"][0]["efficiency"].get<double>(),
					 Catch::Matchers::WithinAbs(0.5, 1e-9));
	}

	SECTION("Waveform update")
	{
		auto wf = scenario["simulation"]["waveforms"][0];
		wf["name"] = "UpdatedWaveform";
		wf["power"] = 999.0;

		REQUIRE(fers_update_waveform_from_json(context.get(), wf.dump().c_str()) == 0);

		auto updated = api_test::parseScenarioJson(context.get());
		REQUIRE(updated["simulation"]["waveforms"][0]["name"] == "UpdatedWaveform");
		REQUIRE_THAT(updated["simulation"]["waveforms"][0]["power"].get<double>(),
					 Catch::Matchers::WithinAbs(999.0, 1e-9));
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
		REQUIRE(fers_update_platform_from_json(nullptr, 1, "{}") == -1);
		REQUIRE(fers_update_platform_from_json(context.get(), 1, nullptr) == -1);
		REQUIRE(fers_update_antenna_from_json(nullptr, "{}") == -1);
		REQUIRE(fers_update_antenna_from_json(context.get(), nullptr) == -1);
		REQUIRE(fers_update_waveform_from_json(nullptr, "{}") == -1);
		REQUIRE(fers_update_waveform_from_json(context.get(), nullptr) == -1);
	}

	SECTION("Invalid JSON")
	{
		REQUIRE(fers_update_platform_from_json(context.get(), 1, "{bad") == 1);
		REQUIRE(fers_update_antenna_from_json(context.get(), "{bad") == 1);
		REQUIRE(fers_update_waveform_from_json(context.get(), "{bad") == 1);
	}

	SECTION("Platform not found")
	{
		REQUIRE(fers_update_platform_from_json(context.get(), 999999, "{}") == 1);
		api_test::ApiString error = api_test::lastError();
		REQUIRE_THAT(error.str(), ContainsSubstring("Platform not found"));
	}
}

TEST_CASE("API XML JSON round-trip keeps the context usable", "[api][scenario]")
{
	api_test::ParamGuard guard;
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::minimalScenarioXml("Round Trip Original");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	auto scenario = api_test::parseScenarioJson(context.get());
	scenario["simulation"]["name"] = "Round Trip Updated";
	REQUIRE(fers_update_scenario_from_json(context.get(), scenario.dump().c_str()) == 0);

	api_test::ApiString xml_text = api_test::scenarioAsXml(context.get());
	REQUIRE(xml_text.get() != nullptr);
	REQUIRE_THAT(xml_text.str(), ContainsSubstring("Round Trip Updated"));

	const auto after_update = api_test::parseScenarioJson(context.get());
	REQUIRE(after_update.at("simulation").at("name") == "Round Trip Updated");
}
