#pragma once

#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <libfers/api.h>
#include <nlohmann/json.hpp>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>

#include "core/parameters.h"

namespace api_test
{
	using json = nlohmann::json;

	struct ParamGuard
	{
		params::Parameters saved;

		ParamGuard() : saved(params::params) {}

		~ParamGuard() { params::params = saved; }
	};

	struct Context
	{
		fers_context_t* ptr;

		Context() : ptr(fers_context_create()) {}

		explicit Context(fers_context_t* raw) : ptr(raw) {}

		Context(const Context&) = delete;
		Context& operator=(const Context&) = delete;

		Context(Context&& other) noexcept : ptr(std::exchange(other.ptr, nullptr)) {}

		Context& operator=(Context&& other) noexcept
		{
			if (this != &other)
			{
				fers_context_destroy(ptr);
				ptr = std::exchange(other.ptr, nullptr);
			}
			return *this;
		}

		~Context() { fers_context_destroy(ptr); }

		[[nodiscard]] fers_context_t* get() const noexcept { return ptr; }

		explicit operator bool() const noexcept { return ptr != nullptr; }
	};

	struct ApiString
	{
		char* ptr;

		ApiString() : ptr(nullptr) {}

		explicit ApiString(char* raw) : ptr(raw) {}

		ApiString(const ApiString&) = delete;
		ApiString& operator=(const ApiString&) = delete;

		ApiString(ApiString&& other) noexcept : ptr(std::exchange(other.ptr, nullptr)) {}

		ApiString& operator=(ApiString&& other) noexcept
		{
			if (this != &other)
			{
				fers_free_string(ptr);
				ptr = std::exchange(other.ptr, nullptr);
			}
			return *this;
		}

		~ApiString() { fers_free_string(ptr); }

		[[nodiscard]] char* get() const noexcept { return ptr; }

		[[nodiscard]] std::string_view view() const noexcept
		{
			return ptr ? std::string_view(ptr) : std::string_view{};
		}

		[[nodiscard]] std::string str() const { return ptr ? std::string(ptr) : std::string{}; }
	};

	struct MotionPath
	{
		fers_interpolated_path_t* ptr;

		explicit MotionPath(fers_interpolated_path_t* raw) : ptr(raw) {}

		MotionPath(const MotionPath&) = delete;
		MotionPath& operator=(const MotionPath&) = delete;

		MotionPath(MotionPath&& other) noexcept : ptr(std::exchange(other.ptr, nullptr)) {}

		MotionPath& operator=(MotionPath&& other) noexcept
		{
			if (this != &other)
			{
				fers_free_interpolated_motion_path(ptr);
				ptr = std::exchange(other.ptr, nullptr);
			}
			return *this;
		}

		~MotionPath() { fers_free_interpolated_motion_path(ptr); }

		[[nodiscard]] fers_interpolated_path_t* get() const noexcept { return ptr; }
	};

	struct RotationPath
	{
		fers_interpolated_rotation_path_t* ptr;

		explicit RotationPath(fers_interpolated_rotation_path_t* raw) : ptr(raw) {}

		RotationPath(const RotationPath&) = delete;
		RotationPath& operator=(const RotationPath&) = delete;

		RotationPath(RotationPath&& other) noexcept : ptr(std::exchange(other.ptr, nullptr)) {}

		RotationPath& operator=(RotationPath&& other) noexcept
		{
			if (this != &other)
			{
				fers_free_interpolated_rotation_path(ptr);
				ptr = std::exchange(other.ptr, nullptr);
			}
			return *this;
		}

		~RotationPath() { fers_free_interpolated_rotation_path(ptr); }

		[[nodiscard]] fers_interpolated_rotation_path_t* get() const noexcept { return ptr; }
	};

	struct AntennaPattern
	{
		fers_antenna_pattern_data_t* ptr;

		explicit AntennaPattern(fers_antenna_pattern_data_t* raw) : ptr(raw) {}

		AntennaPattern(const AntennaPattern&) = delete;
		AntennaPattern& operator=(const AntennaPattern&) = delete;

		AntennaPattern(AntennaPattern&& other) noexcept : ptr(std::exchange(other.ptr, nullptr)) {}

		AntennaPattern& operator=(AntennaPattern&& other) noexcept
		{
			if (this != &other)
			{
				fers_free_antenna_pattern_data(ptr);
				ptr = std::exchange(other.ptr, nullptr);
			}
			return *this;
		}

		~AntennaPattern() { fers_free_antenna_pattern_data(ptr); }

		[[nodiscard]] fers_antenna_pattern_data_t* get() const noexcept { return ptr; }
	};

	struct PreviewLinks
	{
		fers_visual_link_list_t* ptr;

		explicit PreviewLinks(fers_visual_link_list_t* raw) : ptr(raw) {}

		PreviewLinks(const PreviewLinks&) = delete;
		PreviewLinks& operator=(const PreviewLinks&) = delete;

		PreviewLinks(PreviewLinks&& other) noexcept : ptr(std::exchange(other.ptr, nullptr)) {}

		PreviewLinks& operator=(PreviewLinks&& other) noexcept
		{
			if (this != &other)
			{
				fers_free_preview_links(ptr);
				ptr = std::exchange(other.ptr, nullptr);
			}
			return *this;
		}

		~PreviewLinks() { fers_free_preview_links(ptr); }

		[[nodiscard]] fers_visual_link_list_t* get() const noexcept { return ptr; }
	};

	struct ScopedPath
	{
		std::filesystem::path path;

		explicit ScopedPath(std::filesystem::path p) : path(std::move(p)) {}

		ScopedPath(const ScopedPath&) = delete;
		ScopedPath& operator=(const ScopedPath&) = delete;

		ScopedPath(ScopedPath&& other) noexcept : path(std::move(other.path)) { other.path.clear(); }

		ScopedPath& operator=(ScopedPath&& other) noexcept
		{
			if (this != &other)
			{
				cleanup();
				path = std::move(other.path);
				other.path.clear();
			}
			return *this;
		}

		~ScopedPath() { cleanup(); }

		void cleanup() const noexcept
		{
			if (path.empty())
			{
				return;
			}
			std::error_code ec;
			std::filesystem::remove_all(path, ec);
		}
	};

	inline ApiString lastError() { return ApiString(fers_get_last_error_message()); }

	inline void clearLastError()
	{
		Context context;
		REQUIRE(context.get() != nullptr);
	}

	inline std::filesystem::path uniqueTempPath(const std::string& prefix, const std::string& extension = "")
	{
		return std::filesystem::temp_directory_path() /
			(prefix + "_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + extension);
	}

	inline std::string pathString(const std::filesystem::path& path) { return path.string(); }

	inline void writeTextFile(const std::filesystem::path& path, const std::string& content)
	{
		std::ofstream out(path, std::ios::binary);
		REQUIRE(out.is_open());
		out << content;
	}

	inline std::string readTextFile(const std::filesystem::path& path)
	{
		std::ifstream in(path, std::ios::binary);
		REQUIRE(in.is_open());
		return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
	}

	inline std::string malformedXml() { return R"(<simulation><parameters><starttime>0</starttime></parameters>)"; }

	inline std::string malformedJson() { return "{\"simulation\":"; }

	inline std::string minimalScenarioXml(const std::string& simulationName = "API Minimal Scenario")
	{
		return "<simulation name=\"" + simulationName +
			"\">\n"
			"  <parameters>\n"
			"    <starttime>0</starttime>\n"
			"    <endtime>1</endtime>\n"
			"    <rate>1000</rate>\n"
			"    <origin latitude=\"-33.957652\" longitude=\"18.4611991\" altitude=\"111.01\"/>\n"
			"    <coordinatesystem frame=\"ENU\"/>\n"
			"  </parameters>\n"
			"  <waveform name=\"api_min_wave\">\n"
			"    <power>1000</power>\n"
			"    <carrier_frequency>1000000000</carrier_frequency>\n"
			"    <cw/>\n"
			"  </waveform>\n"
			"  <timing name=\"api_clock\">\n"
			"    <frequency>1000000</frequency>\n"
			"  </timing>\n"
			"  <antenna name=\"api_iso\" pattern=\"isotropic\"/>\n"
			"  <platform name=\"api_sensor\">\n"
			"    <motionpath interpolation=\"static\">\n"
			"      <positionwaypoint>\n"
			"        <x>0</x>\n"
			"        <y>0</y>\n"
			"        <altitude>100</altitude>\n"
			"        <time>0</time>\n"
			"      </positionwaypoint>\n"
			"    </motionpath>\n"
			"    <rotationpath interpolation=\"static\">\n"
			"      <rotationwaypoint>\n"
			"        <azimuth>0</azimuth>\n"
			"        <elevation>0</elevation>\n"
			"        <time>0</time>\n"
			"      </rotationwaypoint>\n"
			"    </rotationpath>\n"
			"    <transmitter name=\"api_tx\" waveform=\"api_min_wave\" antenna=\"api_iso\" timing=\"api_clock\">\n"
			"      <cw_mode/>\n"
			"    </transmitter>\n"
			"  </platform>\n"
			"</simulation>\n";
	}

	inline std::string previewScenarioXml(const std::string& simulationName = "API Preview Scenario")
	{
		return "<simulation name=\"" + simulationName +
			"\">\n"
			"  <parameters>\n"
			"    <starttime>0</starttime>\n"
			"    <endtime>0.1</endtime>\n"
			"    <rate>1000000</rate>\n"
			"    <origin latitude=\"-33.957652\" longitude=\"18.4611991\" altitude=\"111.01\"/>\n"
			"    <coordinatesystem frame=\"ENU\"/>\n"
			"  </parameters>\n"
			"  <waveform name=\"api_preview_wave\">\n"
			"    <power>1000000</power>\n"
			"    <carrier_frequency>1000000000</carrier_frequency>\n"
			"    <cw/>\n"
			"  </waveform>\n"
			"  <timing name=\"api_preview_clock\">\n"
			"    <frequency>1000000</frequency>\n"
			"  </timing>\n"
			"  <antenna name=\"api_preview_iso\" pattern=\"isotropic\"/>\n"
			"  <platform name=\"api_tx_platform\">\n"
			"    <motionpath interpolation=\"static\">\n"
			"      <positionwaypoint>\n"
			"        <x>0</x>\n"
			"        <y>0</y>\n"
			"        <altitude>0</altitude>\n"
			"        <time>0</time>\n"
			"      </positionwaypoint>\n"
			"    </motionpath>\n"
			"    <rotationpath interpolation=\"static\">\n"
			"      <rotationwaypoint>\n"
			"        <azimuth>0</azimuth>\n"
			"        <elevation>0</elevation>\n"
			"        <time>0</time>\n"
			"      </rotationwaypoint>\n"
			"    </rotationpath>\n"
			"    <transmitter name=\"api_preview_tx\" waveform=\"api_preview_wave\" antenna=\"api_preview_iso\" "
			"timing=\"api_preview_clock\">\n"
			"      <cw_mode/>\n"
			"    </transmitter>\n"
			"  </platform>\n"
			"  <platform name=\"api_rx_platform\">\n"
			"    <motionpath interpolation=\"static\">\n"
			"      <positionwaypoint>\n"
			"        <x>200</x>\n"
			"        <y>0</y>\n"
			"        <altitude>0</altitude>\n"
			"        <time>0</time>\n"
			"      </positionwaypoint>\n"
			"    </motionpath>\n"
			"    <rotationpath interpolation=\"static\">\n"
			"      <rotationwaypoint>\n"
			"        <azimuth>180</azimuth>\n"
			"        <elevation>0</elevation>\n"
			"        <time>0</time>\n"
			"      </rotationwaypoint>\n"
			"    </rotationpath>\n"
			"    <receiver name=\"api_preview_rx\" antenna=\"api_preview_iso\" timing=\"api_preview_clock\">\n"
			"      <cw_mode/>\n"
			"      <noise_temp>290</noise_temp>\n"
			"    </receiver>\n"
			"  </platform>\n"
			"  <platform name=\"api_target_platform\">\n"
			"    <motionpath interpolation=\"static\">\n"
			"      <positionwaypoint>\n"
			"        <x>100</x>\n"
			"        <y>50</y>\n"
			"        <altitude>0</altitude>\n"
			"        <time>0</time>\n"
			"      </positionwaypoint>\n"
			"    </motionpath>\n"
			"    <rotationpath interpolation=\"static\">\n"
			"      <rotationwaypoint>\n"
			"        <azimuth>0</azimuth>\n"
			"        <elevation>0</elevation>\n"
			"        <time>0</time>\n"
			"      </rotationwaypoint>\n"
			"    </rotationpath>\n"
			"    <target name=\"api_preview_target\">\n"
			"      <rcs type=\"isotropic\">\n"
			"        <value>100</value>\n"
			"      </rcs>\n"
			"      <model type=\"constant\"/>\n"
			"    </target>\n"
			"  </platform>\n"
			"</simulation>\n";
	}

	inline std::string
	monostaticPreviewScenarioXml(const std::string& simulationName = "API Monostatic Preview Scenario")
	{
		return "<simulation name=\"" + simulationName +
			"\">\n"
			"  <parameters>\n"
			"    <starttime>0</starttime>\n"
			"    <endtime>0.1</endtime>\n"
			"    <rate>1000000</rate>\n"
			"    <origin latitude=\"-33.957652\" longitude=\"18.4611991\" altitude=\"111.01\"/>\n"
			"    <coordinatesystem frame=\"ENU\"/>\n"
			"  </parameters>\n"
			"  <waveform name=\"api_preview_wave\">\n"
			"    <power>1000000</power>\n"
			"    <carrier_frequency>1000000000</carrier_frequency>\n"
			"    <cw/>\n"
			"  </waveform>\n"
			"  <timing name=\"api_preview_clock\">\n"
			"    <frequency>1000000</frequency>\n"
			"  </timing>\n"
			"  <antenna name=\"api_preview_iso\" pattern=\"isotropic\"/>\n"
			"  <platform name=\"api_mono_platform\">\n"
			"    <motionpath interpolation=\"static\">\n"
			"      <positionwaypoint>\n"
			"        <x>0</x>\n"
			"        <y>0</y>\n"
			"        <altitude>0</altitude>\n"
			"        <time>0</time>\n"
			"      </positionwaypoint>\n"
			"    </motionpath>\n"
			"    <rotationpath interpolation=\"static\">\n"
			"      <rotationwaypoint>\n"
			"        <azimuth>0</azimuth>\n"
			"        <elevation>0</elevation>\n"
			"        <time>0</time>\n"
			"      </rotationwaypoint>\n"
			"    </rotationpath>\n"
			"    <monostatic name=\"api_preview_mono\" waveform=\"api_preview_wave\" antenna=\"api_preview_iso\" "
			"timing=\"api_preview_clock\">\n"
			"      <cw_mode/>\n"
			"      <noise_temp>290</noise_temp>\n"
			"    </monostatic>\n"
			"  </platform>\n"
			"  <platform name=\"api_target_platform\">\n"
			"    <motionpath interpolation=\"static\">\n"
			"      <positionwaypoint>\n"
			"        <x>100</x>\n"
			"        <y>0</y>\n"
			"        <altitude>0</altitude>\n"
			"        <time>0</time>\n"
			"      </positionwaypoint>\n"
			"    </motionpath>\n"
			"    <rotationpath interpolation=\"static\">\n"
			"      <rotationwaypoint>\n"
			"        <azimuth>0</azimuth>\n"
			"        <elevation>0</elevation>\n"
			"        <time>0</time>\n"
			"      </rotationwaypoint>\n"
			"    </rotationpath>\n"
			"    <target name=\"api_preview_target\">\n"
			"      <rcs type=\"isotropic\">\n"
			"        <value>100</value>\n"
			"      </rcs>\n"
			"      <model type=\"constant\"/>\n"
			"    </target>\n"
			"  </platform>\n"
			"</simulation>\n";
	}

	inline ApiString scenarioAsJson(fers_context_t* context) { return ApiString(fers_get_scenario_as_json(context)); }

	inline ApiString scenarioAsXml(fers_context_t* context) { return ApiString(fers_get_scenario_as_xml(context)); }

	inline json parseScenarioJson(fers_context_t* context)
	{
		ApiString json_text = scenarioAsJson(context);
		REQUIRE(json_text.get() != nullptr);
		return json::parse(json_text.str());
	}

	inline std::uint64_t parseId(const json& value)
	{
		if (value.is_string())
		{
			return std::stoull(value.get<std::string>());
		}
		if (value.is_number_unsigned())
		{
			return value.get<std::uint64_t>();
		}
		if (value.is_number_integer())
		{
			return static_cast<std::uint64_t>(value.get<std::int64_t>());
		}
		FAIL("ID value is not string or integer");
		return 0;
	}
}
