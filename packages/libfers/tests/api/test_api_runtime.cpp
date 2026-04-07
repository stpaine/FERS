#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <libfers/api.h>
#include <string>
#include <vector>

#include "api_test_helpers.h"

using Catch::Matchers::ContainsSubstring;

namespace
{
	std::string uniqueLogMessage(const std::string& prefix)
	{
		return prefix + "_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
	}

	struct CallbackState
	{
		int calls = 0;
		void* seen_user_data = nullptr;
		std::vector<std::string> messages;
	};

	void recordProgress(const char* message, int /*current*/, int /*total*/, void* user_data)
	{
		auto* state = static_cast<CallbackState*>(user_data);
		++state->calls;
		state->seen_user_data = user_data;
		state->messages.emplace_back(message ? message : "");
	}
}

TEST_CASE("API log level mapping covers all exported enum values", "[api][runtime]")
{
	const auto log_path = api_test::uniqueTempPath("api_log_levels", ".log");
	api_test::ScopedPath log_guard(log_path);
	const auto rollover_path = api_test::uniqueTempPath("api_log_levels_rollover", ".log");
	api_test::ScopedPath rollover_guard(rollover_path);

	const struct
	{
		fers_log_level_t level;
		const char* label;
	} cases[] = {
		{FERS_LOG_TRACE, "trace"},
		{FERS_LOG_DEBUG, "debug"},
		{FERS_LOG_WARNING, "warning"},
		{FERS_LOG_ERROR, "error"},
		{FERS_LOG_FATAL, "fatal"},
		// Exercise the fallback branch with a value outside the exported enum set but still representable.
		{static_cast<fers_log_level_t>(FERS_LOG_FATAL + 1), "default"},
	};

	for (const auto& entry : cases)
	{
		const std::string message = uniqueLogMessage(std::string("api_level_") + entry.label);
		REQUIRE(fers_configure_logging(entry.level, log_path.c_str()) == 0);
		fers_log(entry.level, message.c_str());
		REQUIRE(fers_configure_logging(FERS_LOG_INFO, rollover_path.c_str()) == 0);

		const std::string log_text = api_test::readTextFile(log_path);
		REQUIRE_THAT(log_text, ContainsSubstring(message));
	}
}

TEST_CASE("API log ignores null messages", "[api][runtime]")
{
	const auto log_path = api_test::uniqueTempPath("api_log_null_message", ".log");
	api_test::ScopedPath log_guard(log_path);
	const auto rollover_path = api_test::uniqueTempPath("api_log_null_message_rollover", ".log");
	api_test::ScopedPath rollover_guard(rollover_path);

	REQUIRE(fers_configure_logging(FERS_LOG_INFO, log_path.c_str()) == 0);
	fers_log(FERS_LOG_INFO, nullptr);
	REQUIRE(fers_configure_logging(FERS_LOG_INFO, rollover_path.c_str()) == 0);

	REQUIRE(api_test::readTextFile(log_path).empty());
}

TEST_CASE("API configure logging accepts null and writable file destinations", "[api][runtime]")
{
	api_test::clearLastError();

	SECTION("null log path")
	{
		REQUIRE(fers_configure_logging(FERS_LOG_INFO, nullptr) == 0);
		api_test::ApiString error = api_test::lastError();
		REQUIRE(error.get() == nullptr);
	}

	SECTION("temp log path")
	{
		const auto log_path = api_test::uniqueTempPath("api_runtime", ".log");
		api_test::ScopedPath log_guard(log_path);

		REQUIRE(fers_configure_logging(FERS_LOG_INFO, log_path.c_str()) == 0);
		REQUIRE(std::filesystem::exists(log_path));

		api_test::ApiString error = api_test::lastError();
		REQUIRE(error.get() == nullptr);
	}
}

TEST_CASE("API configure logging reports file open failures", "[api][runtime]")
{
	api_test::clearLastError();
	const auto missing_parent = api_test::uniqueTempPath("api_log_dir");
	api_test::ScopedPath missing_guard(missing_parent);
	const auto log_path = missing_parent / "runtime.log";

	REQUIRE(fers_configure_logging(FERS_LOG_INFO, log_path.c_str()) == 1);

	api_test::ApiString error = api_test::lastError();
	REQUIRE(error.get() != nullptr);
	REQUIRE_THAT(error.str(), ContainsSubstring("Unable to open log file:"));
}

TEST_CASE("API log writes to configured file", "[api][runtime]")
{
	api_test::clearLastError();
	const auto log_path = api_test::uniqueTempPath("api_log_smoke", ".log");
	api_test::ScopedPath log_guard(log_path);
	const auto rollover_path = api_test::uniqueTempPath("api_log_rollover", ".log");
	api_test::ScopedPath rollover_guard(rollover_path);

	REQUIRE(fers_configure_logging(FERS_LOG_INFO, log_path.c_str()) == 0);
	const std::string message = uniqueLogMessage("api runtime logging smoke message");
	fers_log(FERS_LOG_INFO, message.c_str());
	REQUIRE(fers_configure_logging(FERS_LOG_INFO, rollover_path.c_str()) == 0);

	const std::string log_text = api_test::readTextFile(log_path);
	REQUIRE_THAT(log_text, ContainsSubstring(message));
	REQUIRE_THAT(log_text, ContainsSubstring("INFO"));
}

TEST_CASE("API warning getter returns deduplicated rotation-unit warnings", "[api][runtime]")
{
	api_test::ParamGuard guard;
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::minimalScenarioXml("API Warning Runtime");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	auto scenario = api_test::parseScenarioJson(context.get());
	scenario["simulation"]["parameters"]["rotationangleunit"] = "rad";
	scenario["simulation"]["platforms"][0]["rotationpath"]["rotationwaypoints"][0]["azimuth"] = 90.0;

	REQUIRE(fers_update_scenario_from_json(context.get(), scenario.dump().c_str()) == 0);

	api_test::ApiString warnings_json(fers_get_last_warning_messages_json());
	REQUIRE(warnings_json.get() != nullptr);

	const auto warnings = api_test::json::parse(warnings_json.str());
	REQUIRE(warnings.is_array());
	REQUIRE(warnings.size() == 1);
	REQUIRE_THAT(warnings[0].get<std::string>(), ContainsSubstring("platform 'api_sensor' rotation waypoint 0"));
	REQUIRE_THAT(warnings[0].get<std::string>(), ContainsSubstring("'azimuth'"));

	api_test::ApiString cleared(fers_get_last_warning_messages_json());
	REQUIRE(cleared.get() == nullptr);
}

TEST_CASE("API run simulation rejects null context", "[api][runtime]")
{
	api_test::clearLastError();

	REQUIRE(fers_run_simulation(nullptr, nullptr, nullptr) == -1);

	api_test::ApiString error = api_test::lastError();
	REQUIRE(error.get() != nullptr);
	REQUIRE_THAT(error.str(), ContainsSubstring("Invalid context provided to fers_run_simulation"));
}

TEST_CASE("API run simulation accepts a minimal valid scenario", "[api][runtime]")
{
	api_test::ParamGuard guard;
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const auto out_dir = api_test::uniqueTempPath("api_out_dir");
	std::filesystem::create_directories(out_dir);
	api_test::ScopedPath dir_guard(out_dir);

	const std::string unique_rx_name =
		"api_preview_rx_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
	const auto output_path = out_dir / (unique_rx_name + "_results.h5");

	REQUIRE(fers_set_output_directory(context.get(), out_dir.string().c_str()) == 0);

	std::string xml = api_test::previewScenarioXml("Runtime Scenario");
	size_t pos = xml.find("api_preview_rx");
	while (pos != std::string::npos)
	{
		xml.replace(pos, 14, unique_rx_name);
		pos = xml.find("api_preview_rx", pos + unique_rx_name.length());
	}
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);
	REQUIRE(fers_run_simulation(context.get(), nullptr, nullptr) == 0);

	api_test::ApiString error = api_test::lastError();
	REQUIRE(error.get() == nullptr);
	REQUIRE(std::filesystem::exists(output_path));
}

TEST_CASE("API run simulation invokes progress callbacks with caller user data", "[api][runtime]")
{
	api_test::ParamGuard guard;
	api_test::clearLastError();
	api_test::Context context;
	REQUIRE(context.get() != nullptr);

	const auto out_dir = api_test::uniqueTempPath("api_out_dir_cb");
	std::filesystem::create_directories(out_dir);
	api_test::ScopedPath dir_guard(out_dir);

	const std::string unique_rx_name =
		"api_preview_rx_cb_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
	const auto output_path = out_dir / (unique_rx_name + "_results.h5");

	REQUIRE(fers_set_output_directory(context.get(), out_dir.string().c_str()) == 0);

	std::string xml = api_test::previewScenarioXml("Runtime Callback Scenario");
	size_t pos = xml.find("api_preview_rx");
	while (pos != std::string::npos)
	{
		xml.replace(pos, 14, unique_rx_name);
		pos = xml.find("api_preview_rx", pos + unique_rx_name.length());
	}
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	CallbackState state;
	REQUIRE(fers_run_simulation(context.get(), recordProgress, &state) == 0);
	REQUIRE(state.calls > 0);
	REQUIRE(state.seen_user_data == &state);

	bool saw_expected_message = false;
	for (const auto& message : state.messages)
	{
		if (message.find("Initializing event-driven simulation") != std::string::npos ||
			message.find("Simulation complete") != std::string::npos)
		{
			saw_expected_message = true;
			break;
		}
	}
	REQUIRE(saw_expected_message);
	REQUIRE(std::filesystem::exists(output_path));
}
