#include <atomic>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <chrono>
#include <filesystem>
#include <libfers/api.h>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include "api_test_helpers.h"
#include "core/output_metadata.h"

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
		std::vector<std::chrono::steady_clock::time_point> message_times;
	};

	struct LogCallbackState
	{
		int calls = 0;
		void* seen_user_data = nullptr;
		std::vector<fers_log_level_t> levels;
		std::vector<std::string> lines;
	};

	void recordProgress(const char* message, int /*current*/, int /*total*/, void* user_data)
	{
		auto* state = static_cast<CallbackState*>(user_data);
		++state->calls;
		state->seen_user_data = user_data;
		state->messages.emplace_back((message != nullptr) ? message : "");
		state->message_times.push_back(std::chrono::steady_clock::now());
	}

	void recordLog(fers_log_level_t level, const char* line, void* user_data)
	{
		auto* state = static_cast<LogCallbackState*>(user_data);
		++state->calls;
		state->seen_user_data = user_data;
		state->levels.push_back(level);
		state->lines.emplace_back((line != nullptr) ? line : "");
	}

	int requestCancel(void* user_data)
	{
		auto* cancel = static_cast<std::atomic_bool*>(user_data);
		return cancel->load() ? 1 : 0;
	}

	void replaceOnce(std::string& value, const std::string& from, const std::string& to)
	{
		const auto pos = value.find(from);
		REQUIRE(pos != std::string::npos);
		value.replace(pos, from.size(), to);
	}

	[[nodiscard]] std::string vita49DrainTimingScenarioXml()
	{
		std::string xml = api_test::previewScenarioXml("API VITA Drain Timing");
		replaceOnce(xml, "<endtime>0.1</endtime>", "<endtime>0.2</endtime>");
		replaceOnce(xml, "<rate>1000000</rate>", "<rate>100000</rate>");
		return xml;
	}

	[[nodiscard]] std::chrono::steady_clock::time_point messageTime(const CallbackState& state,
																	const std::string_view message)
	{
		for (std::size_t i = 0; i < state.messages.size(); ++i)
		{
			if (state.messages.at(i) == message)
			{
				return state.message_times.at(i);
			}
		}
		return {};
	}

	void requireVita49DrainCompletionTiming(const CallbackState& state,
											const std::chrono::steady_clock::time_point start)
	{
		const auto drain_time = messageTime(state, "Waiting for VITA output stream drain...");
		const auto completion_time = messageTime(state, "Simulation complete");

		REQUIRE(drain_time != std::chrono::steady_clock::time_point{});
		REQUIRE(completion_time != std::chrono::steady_clock::time_point{});
		CHECK(drain_time <= completion_time);
		CHECK(completion_time >= start + std::chrono::milliseconds(180));
	}
}

TEST_CASE("API log level mapping writes emitted enum values", "[api][runtime]")
{
	const auto log_path = api_test::uniqueTempPath("api_log_levels", ".log");
	api_test::ScopedPath const log_guard(log_path);
	const auto rollover_path = api_test::uniqueTempPath("api_log_levels_rollover", ".log");
	api_test::ScopedPath const rollover_guard(rollover_path);

	const struct
	{
		fers_log_level_t level;
		const char* label;
	} cases[] = {
		{FERS_LOG_TRACE, "trace"},
		{FERS_LOG_DEBUG, "debug"},
		{FERS_LOG_INFO, "info"},
		{FERS_LOG_WARNING, "warning"},
		{FERS_LOG_ERROR, "error"},
		{FERS_LOG_FATAL, "fatal"},
		// Exercise the fallback branch with a value outside the exported enum set but still representable.
		{static_cast<fers_log_level_t>(FERS_LOG_OFF + 1), "default"},
	};

	for (const auto& entry : cases)
	{
		const std::string message = uniqueLogMessage(std::string("api_level_") + entry.label);
		const std::string log_path_string = api_test::pathString(log_path);
		const std::string rollover_path_string = api_test::pathString(rollover_path);
		REQUIRE(fers_configure_logging(entry.level, log_path_string.c_str()) == 0);
		fers_log(entry.level, message.c_str());
		REQUIRE(fers_configure_logging(FERS_LOG_INFO, rollover_path_string.c_str()) == 0);

		const std::string log_text = api_test::readTextFile(log_path);
		REQUIRE_THAT(log_text, ContainsSubstring(message));
	}
}

TEST_CASE("API log level getter round-trips configured levels", "[api][runtime]")
{
	const fers_log_level_t cases[] = {
		FERS_LOG_TRACE, FERS_LOG_DEBUG, FERS_LOG_INFO, FERS_LOG_WARNING, FERS_LOG_ERROR, FERS_LOG_FATAL, FERS_LOG_OFF,
	};

	for (const auto level : cases)
	{
		REQUIRE(fers_configure_logging(level, nullptr) == 0);
		REQUIRE(fers_get_log_level() == level);
	}

	REQUIRE(fers_configure_logging(FERS_LOG_INFO, nullptr) == 0);
}

TEST_CASE("API log ignores null messages", "[api][runtime]")
{
	const auto log_path = api_test::uniqueTempPath("api_log_null_message", ".log");
	api_test::ScopedPath const log_guard(log_path);
	const auto rollover_path = api_test::uniqueTempPath("api_log_null_message_rollover", ".log");
	api_test::ScopedPath const rollover_guard(rollover_path);
	const std::string log_path_string = api_test::pathString(log_path);
	const std::string rollover_path_string = api_test::pathString(rollover_path);

	REQUIRE(fers_configure_logging(FERS_LOG_INFO, log_path_string.c_str()) == 0);
	fers_log(FERS_LOG_INFO, nullptr);
	REQUIRE(fers_configure_logging(FERS_LOG_INFO, rollover_path_string.c_str()) == 0);

	REQUIRE(api_test::readTextFile(log_path).empty());
}

TEST_CASE("API configure logging accepts null and writable file destinations", "[api][runtime]")
{
	api_test::clearLastError();

	SECTION("null log path")
	{
		REQUIRE(fers_configure_logging(FERS_LOG_INFO, nullptr) == 0);
		api_test::ApiString const error = api_test::lastError();
		REQUIRE(error.get() == nullptr);
	}

	SECTION("temp log path")
	{
		const auto log_path = api_test::uniqueTempPath("api_runtime", ".log");
		api_test::ScopedPath const log_guard(log_path);
		const std::string log_path_string = api_test::pathString(log_path);

		REQUIRE(fers_configure_logging(FERS_LOG_INFO, log_path_string.c_str()) == 0);
		REQUIRE(std::filesystem::exists(log_path));

		api_test::ApiString const error = api_test::lastError();
		REQUIRE(error.get() == nullptr);
	}
}

TEST_CASE("API configure logging reports file open failures", "[api][runtime]")
{
	api_test::clearLastError();
	const auto missing_parent = api_test::uniqueTempPath("api_log_dir");
	api_test::ScopedPath const missing_guard(missing_parent);
	const auto log_path = missing_parent / "runtime.log";
	const std::string log_path_string = api_test::pathString(log_path);

	REQUIRE(fers_configure_logging(FERS_LOG_INFO, log_path_string.c_str()) == 1);

	api_test::ApiString const error = api_test::lastError();
	REQUIRE(error.get() != nullptr);
	REQUIRE_THAT(error.str(), ContainsSubstring("Unable to open log file:"));
}

TEST_CASE("API log writes to configured file", "[api][runtime]")
{
	api_test::clearLastError();
	const auto log_path = api_test::uniqueTempPath("api_log_smoke", ".log");
	api_test::ScopedPath const log_guard(log_path);
	const auto rollover_path = api_test::uniqueTempPath("api_log_rollover", ".log");
	api_test::ScopedPath const rollover_guard(rollover_path);
	const std::string log_path_string = api_test::pathString(log_path);
	const std::string rollover_path_string = api_test::pathString(rollover_path);

	REQUIRE(fers_configure_logging(FERS_LOG_INFO, log_path_string.c_str()) == 0);
	const std::string message = uniqueLogMessage("api runtime logging smoke message");
	fers_log(FERS_LOG_INFO, message.c_str());
	REQUIRE(fers_configure_logging(FERS_LOG_INFO, rollover_path_string.c_str()) == 0);

	const std::string log_text = api_test::readTextFile(log_path);
	REQUIRE_THAT(log_text, ContainsSubstring(message));
	REQUIRE_THAT(log_text, ContainsSubstring("INFO"));
}

TEST_CASE("API OFF log level suppresses file and callback output", "[api][runtime]")
{
	api_test::clearLastError();
	const auto log_path = api_test::uniqueTempPath("api_log_off", ".log");
	api_test::ScopedPath const log_guard(log_path);
	const auto rollover_path = api_test::uniqueTempPath("api_log_off_rollover", ".log");
	api_test::ScopedPath const rollover_guard(rollover_path);
	const std::string log_path_string = api_test::pathString(log_path);
	const std::string rollover_path_string = api_test::pathString(rollover_path);

	LogCallbackState state;
	fers_set_log_callback(recordLog, &state);

	REQUIRE(fers_configure_logging(FERS_LOG_OFF, log_path_string.c_str()) == 0);
	fers_log(FERS_LOG_FATAL, "suppressed fatal message");
	fers_log(FERS_LOG_OFF, "suppressed off message");
	REQUIRE(fers_configure_logging(FERS_LOG_INFO, rollover_path_string.c_str()) == 0);
	fers_set_log_callback(nullptr, nullptr);

	REQUIRE(state.calls == 0);
	REQUIRE(api_test::readTextFile(log_path).empty());
}

TEST_CASE("API log callback receives formatted accepted lines", "[api][runtime]")
{
	api_test::clearLastError();
	REQUIRE(fers_configure_logging(FERS_LOG_WARNING, nullptr) == 0);

	LogCallbackState state;
	fers_set_log_callback(recordLog, &state);

	fers_log(FERS_LOG_INFO, "ignored by callback level");
	REQUIRE(state.calls == 0);

	const std::string message = uniqueLogMessage("api callback message");
	fers_log(FERS_LOG_ERROR, message.c_str());

	REQUIRE(state.calls == 1);
	REQUIRE(state.seen_user_data == &state);
	REQUIRE(state.levels.front() == FERS_LOG_ERROR);
	REQUIRE_THAT(state.lines.front(), ContainsSubstring(message));
	REQUIRE_THAT(state.lines.front(), ContainsSubstring("ERROR"));
	REQUIRE(state.lines.front().back() != '\n');

	fers_set_log_callback(nullptr, nullptr);
	fers_log(FERS_LOG_ERROR, "callback disabled");
	REQUIRE(state.calls == 1);

	REQUIRE(fers_configure_logging(FERS_LOG_INFO, nullptr) == 0);
}

TEST_CASE("API warning getter returns deduplicated rotation-unit warnings", "[api][runtime]")
{
	api_test::ParamGuard const guard;
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::minimalScenarioXml("API Warning Runtime");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	auto scenario = api_test::parseScenarioJson(context.get());
	scenario["simulation"]["parameters"]["rotationangleunit"] = "rad";
	scenario["simulation"]["platforms"][0]["rotationpath"]["rotationwaypoints"][0]["azimuth"] = 90.0;

	REQUIRE(fers_update_scenario_from_json(context.get(), scenario.dump().c_str()) == 0);

	api_test::ApiString const warnings_json(fers_get_last_warning_messages_json());
	REQUIRE(warnings_json.get() != nullptr);

	const auto warnings = api_test::json::parse(warnings_json.str());
	REQUIRE(warnings.is_array());
	REQUIRE(warnings.size() == 1);
	REQUIRE_THAT(warnings[0].get<std::string>(), ContainsSubstring("platform 'api_sensor' rotation waypoint 0"));
	REQUIRE_THAT(warnings[0].get<std::string>(), ContainsSubstring("'azimuth'"));

	api_test::ApiString const cleared(fers_get_last_warning_messages_json());
	REQUIRE(cleared.get() == nullptr);
}

TEST_CASE("API run simulation rejects null context", "[api][runtime]")
{
	api_test::clearLastError();

	REQUIRE(fers_run_simulation(nullptr, nullptr, nullptr) == -1);

	api_test::ApiString const error = api_test::lastError();
	REQUIRE(error.get() != nullptr);
	REQUIRE_THAT(error.str(), ContainsSubstring("Invalid context provided to fers_run_simulation"));
}

TEST_CASE("API VITA49 setters validate control-plane inputs", "[api][runtime][vita49]")
{
	api_test::clearLastError();
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	REQUIRE(fers_enable_vita49_udp_output(nullptr, "127.0.0.1", 4991) == -1);
	api_test::ApiString const null_context_error = api_test::lastError();
	REQUIRE(null_context_error.get() != nullptr);
	REQUIRE_THAT(null_context_error.str(), ContainsSubstring("context is NULL"));

	REQUIRE(fers_enable_vita49_udp_output(context.get(), nullptr, 4991) == -1);
	api_test::ApiString const null_host_error = api_test::lastError();
	REQUIRE(null_host_error.get() != nullptr);
	REQUIRE_THAT(null_host_error.str(), ContainsSubstring("host is NULL"));

	REQUIRE(fers_enable_vita49_udp_output(context.get(), "", 4991) == 1);
	api_test::ApiString const empty_host_error = api_test::lastError();
	REQUIRE(empty_host_error.get() != nullptr);
	REQUIRE_THAT(empty_host_error.str(), ContainsSubstring("host must be non-empty"));

	REQUIRE(fers_enable_vita49_udp_output(context.get(), "127.0.0.1", 0) == 1);
	api_test::ApiString const port_error = api_test::lastError();
	REQUIRE(port_error.get() != nullptr);
	REQUIRE_THAT(port_error.str(), ContainsSubstring("port must be in the range 1..65535"));

	REQUIRE(fers_set_vita49_fullscale(context.get(), 0.0) == 1);
	api_test::ApiString const fullscale_error = api_test::lastError();
	REQUIRE(fullscale_error.get() != nullptr);
	REQUIRE_THAT(fullscale_error.str(), ContainsSubstring("fullscale"));

	REQUIRE(fers_set_vita49_fullscale(context.get(), std::numeric_limits<double>::infinity()) == 1);
	api_test::ApiString const infinite_fullscale_error = api_test::lastError();
	REQUIRE(infinite_fullscale_error.get() != nullptr);
	REQUIRE_THAT(infinite_fullscale_error.str(), ContainsSubstring("positive and finite"));

	REQUIRE(fers_set_vita49_epoch_unix_nanoseconds(context.get(), 4294967296000000000ULL) == 1);
	api_test::ApiString const epoch_error = api_test::lastError();
	REQUIRE(epoch_error.get() != nullptr);
	REQUIRE_THAT(epoch_error.str(), ContainsSubstring("32-bit UTC seconds"));

	REQUIRE(fers_set_vita49_max_udp_payload(context.get(), 0) == 1);
	api_test::ApiString const payload_error = api_test::lastError();
	REQUIRE(payload_error.get() != nullptr);
	REQUIRE_THAT(payload_error.str(), ContainsSubstring("max UDP payload"));

	REQUIRE(fers_set_vita49_queue_depth(context.get(), 0) == 1);
	api_test::ApiString const queue_error = api_test::lastError();
	REQUIRE(queue_error.get() != nullptr);
	REQUIRE_THAT(queue_error.str(), ContainsSubstring("queue depth"));

	REQUIRE(fers_set_vita49_packet_trace_enabled(nullptr, 0) == -1);
	api_test::ApiString const trace_error = api_test::lastError();
	REQUIRE(trace_error.get() != nullptr);
	REQUIRE_THAT(trace_error.str(), ContainsSubstring("context is NULL"));

	REQUIRE(fers_set_vita49_packet_trace_enabled(context.get(), 0) == 0);
	REQUIRE(fers_set_vita49_packet_trace_enabled(context.get(), 1) == 0);
}

TEST_CASE("API run simulation rejects VITA49 mode without fullscale", "[api][runtime][vita49]")
{
	api_test::ParamGuard const guard;
	api_test::clearLastError();
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::minimalScenarioXml("API VITA Missing Fullscale");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);
	REQUIRE(fers_enable_vita49_udp_output(context.get(), "127.0.0.1", 4991) == 0);

	REQUIRE(fers_run_simulation(context.get(), nullptr, nullptr) == 1);
	api_test::ApiString const error = api_test::lastError();
	REQUIRE(error.get() != nullptr);
	REQUIRE_THAT(error.str(), ContainsSubstring("VITA49 fullscale must be a positive finite value"));
}

TEST_CASE("API HDF5 reset clears stale VITA49 output mode", "[api][runtime][vita49]")
{
	api_test::ParamGuard const guard;
	api_test::clearLastError();
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::minimalScenarioXml("API HDF5 Reset");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);
	REQUIRE(fers_enable_vita49_udp_output(context.get(), "127.0.0.1", 4991) == 0);
	REQUIRE(fers_use_hdf5_output(context.get()) == 0);

	REQUIRE(fers_run_simulation(context.get(), nullptr, nullptr) == 0);

	api_test::ApiString const metadata_json(fers_get_last_output_metadata_json(context.get()));
	REQUIRE(metadata_json.get() != nullptr);
	const auto metadata = api_test::json::parse(metadata_json.str());
	CHECK_FALSE(metadata.contains("vita49"));
}

TEST_CASE("API extended simulation reports cooperative cancellation with metadata", "[api][runtime]")
{
	api_test::ParamGuard const guard;
	api_test::clearLastError();
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = api_test::minimalScenarioXml("API Cancel Runtime");
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);

	std::atomic_bool cancel{true};
	REQUIRE(fers_run_simulation_ex(context.get(), nullptr, nullptr, requestCancel, &cancel, nullptr, nullptr) == 2);

	api_test::ApiString const metadata_json(fers_get_last_output_metadata_json(context.get()));
	REQUIRE(metadata_json.get() != nullptr);
	const auto metadata = api_test::json::parse(metadata_json.str());
	CHECK(metadata.at("simulation_name").get<std::string>() == "API Cancel Runtime");
}

TEST_CASE("API run simulation accepts a minimal valid scenario", "[api][runtime]")
{
	api_test::ParamGuard const guard;
	api_test::clearLastError();
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	const auto out_dir = api_test::uniqueTempPath("api_out_dir");
	std::filesystem::create_directories(out_dir);
	api_test::ScopedPath const dir_guard(out_dir);

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

	api_test::ApiString const error = api_test::lastError();
	REQUIRE(error.get() == nullptr);
	REQUIRE(std::filesystem::exists(output_path));

	api_test::ApiString const metadata_json(fers_get_last_output_metadata_json(context.get()));
	REQUIRE(metadata_json.get() != nullptr);
	const auto metadata = api_test::json::parse(metadata_json.str());
	CHECK_FALSE(metadata.contains("vita49"));
}

TEST_CASE("VITA49 metadata section records runtime output config", "[api][runtime][vita49]")
{
	const core::Vita49OutputConfig config{.host = "127.0.0.1",
										  .port = 4991,
										  .adc_fullscale = 2.5,
										  .queue_depth = 17,
										  .epoch_unix_nanoseconds = 1700000000123456789ULL,
										  .max_udp_payload = 1200};
	core::OutputMetadata output_metadata;
	output_metadata.vita49 = core::vita49MetadataFromConfig(config);
	core::Vita49StreamMetadata stream;
	stream.receiver_id = 10;
	stream.receiver_name = "Rx";
	stream.stream_id = 0x12345678u;
	stream.mode = "pulsed";
	stream.sample_rate = 100000.0;
	stream.reference_frequency = 10000000.0;
	stream.packets_emitted = 2930;
	stream.samples_emitted = 1000000;
	stream.context_packet_count = 12;
	stream.first_sample_time = 0.0;
	stream.end_sample_time = 10.0;
	stream.first_timestamp =
		core::Vita49Timestamp{.integer_seconds = 1700000000, .fractional_picoseconds = 123456789000};
	stream.end_timestamp = core::Vita49Timestamp{.integer_seconds = 1700000010, .fractional_picoseconds = 123456789000};
	output_metadata.vita49->streams.push_back(stream);

	const auto metadata = api_test::json::parse(core::outputMetadataToJsonString(output_metadata));
	REQUIRE(metadata.contains("vita49"));
	const auto& vita49 = metadata.at("vita49");
	CHECK(vita49.at("endpoint").get<std::string>() == "127.0.0.1:4991");
	CHECK(vita49.at("endpoint_host").get<std::string>() == "127.0.0.1");
	CHECK(vita49.at("endpoint_port").get<std::uint16_t>() == 4991);
	CHECK(vita49.at("epoch_unix_nanoseconds").get<std::string>() == "1700000000123456789");
	CHECK(vita49.at("class_id").get<std::string>() == "0xFA52530001000101");
	CHECK(vita49.at("adc_fullscale").get<double>() == 2.5);
	CHECK(vita49.at("max_udp_payload").get<std::uint16_t>() == 1200);
	CHECK(vita49.at("queue_depth").get<std::uint32_t>() == 17);
	REQUIRE(vita49.at("streams").is_array());
	REQUIRE(vita49.at("streams").size() == 1u);
	const auto& json_stream = vita49.at("streams").front();
	CHECK(json_stream.at("mode").get<std::string>() == "pulsed");
	CHECK(json_stream.at("samples_emitted").get<std::uint64_t>() == 1000000ULL);
	CHECK(json_stream.at("first_sample_time").get<double>() == 0.0);
	CHECK(json_stream.at("end_sample_time").get<double>() == 10.0);
	CHECK(json_stream.at("first_timestamp").at("integer_seconds").get<std::uint32_t>() == 1700000000u);
	CHECK(json_stream.at("end_timestamp").at("integer_seconds").get<std::uint32_t>() == 1700000010u);
	CHECK_FALSE(json_stream.contains("first_timestamp_unix_ps"));
	CHECK_FALSE(json_stream.contains("last_timestamp_unix_ps"));
}

TEST_CASE("API VITA49 completion waits for wall-clock stream drain", "[api][runtime][vita49]")
{
	api_test::ParamGuard const guard;
	api_test::clearLastError();
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	const std::string xml = vita49DrainTimingScenarioXml();
	REQUIRE(fers_load_scenario_from_xml_string(context.get(), xml.c_str(), 0) == 0);
	REQUIRE(fers_enable_vita49_udp_output(context.get(), "127.0.0.1", 4991) == 0);
	REQUIRE(fers_set_vita49_fullscale(context.get(), 1.0) == 0);

	CallbackState state;
	const auto start = std::chrono::steady_clock::now();
	REQUIRE(fers_run_simulation(context.get(), recordProgress, &state) == 0);

	requireVita49DrainCompletionTiming(state, start);
}

TEST_CASE("API run simulation invokes progress callbacks with caller user data", "[api][runtime]")
{
	api_test::ParamGuard const guard;
	api_test::clearLastError();
	api_test::Context const context;
	REQUIRE(context.get() != nullptr);

	const auto out_dir = api_test::uniqueTempPath("api_out_dir_cb");
	std::filesystem::create_directories(out_dir);
	api_test::ScopedPath const dir_guard(out_dir);

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
