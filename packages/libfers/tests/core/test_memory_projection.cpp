#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>

#include "antenna/antenna_factory.h"
#include "core/logging.h"
#include "core/memory_projection.h"
#include "core/parameters.h"
#include "core/world.h"
#include "math/coord.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/transmitter.h"
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

	struct LogLevelGuard
	{
		explicit LogLevelGuard(const logging::Level level) { logging::logger.setLevel(level); }
		~LogLevelGuard() { logging::logger.setLevel(logging::Level::DEBUG); }
	};

	struct CerrCapture
	{
		std::ostringstream buffer;
		std::streambuf* old{nullptr};
		CerrCapture() { old = std::cerr.rdbuf(buffer.rdbuf()); }
		~CerrCapture() { std::cerr.rdbuf(old); }
		[[nodiscard]] std::string str() const { return buffer.str(); }
	};

	void setupPlatform(radar::Platform& platform)
	{
		platform.getMotionPath()->addCoord(math::Coord{math::Vec3{0.0, 0.0, 0.0}, 0.0});
		platform.getMotionPath()->finalize();
		platform.getRotationPath()->addCoord(math::RotationCoord{0.0, 0.0, 0.0});
		platform.getRotationPath()->finalize();
	}

	std::shared_ptr<timing::Timing> makeTiming(const timing::PrototypeTiming& prototype, const SimId id)
	{
		auto timing = std::make_shared<timing::Timing>("clock_instance", 42, id);
		timing->initializeModel(&prototype);
		return timing;
	}

	std::unique_ptr<core::World> makeProjectionWorld()
	{
		auto world = std::make_unique<core::World>();

		auto platform = std::make_unique<radar::Platform>("Platform", 100);
		setupPlatform(*platform);
		auto* platform_ptr = platform.get();

		auto timing_proto = std::make_unique<timing::PrototypeTiming>("Clock", 200);
		timing_proto->setFrequency(1.0e6);
		timing_proto->setPhaseOffset(0.25);
		const auto timing = makeTiming(*timing_proto, 201);

		auto antenna = std::make_unique<antenna::Isotropic>("Iso", 300);
		auto* antenna_ptr = antenna.get();

		auto wave = std::make_unique<fers_signal::RadarSignal>("Wave", 1.0, 1.0e9, 1.0,
															   std::make_unique<fers_signal::CwSignal>(), 400);
		auto* wave_ptr = wave.get();

		auto tx = std::make_unique<radar::Transmitter>(platform_ptr, "Tx", radar::OperationMode::CW_MODE, 500);
		tx->setTiming(timing);
		tx->setAntenna(antenna_ptr);
		tx->setSignal(wave_ptr);

		auto streaming_rx =
			std::make_unique<radar::Receiver>(platform_ptr, "StreamingRx", 11, radar::OperationMode::CW_MODE, 600);
		streaming_rx->setTiming(timing);
		streaming_rx->setAntenna(antenna_ptr);
		streaming_rx->prepareStreamingData(20);

		auto pulsed_rx =
			std::make_unique<radar::Receiver>(platform_ptr, "PulsedRx", 12, radar::OperationMode::PULSED_MODE, 700);
		pulsed_rx->setTiming(timing);
		pulsed_rx->setAntenna(antenna_ptr);
		pulsed_rx->setWindowProperties(0.2, 2.0, 0.0);

		world->add(std::move(platform));
		world->add(std::move(timing_proto));
		world->add(std::move(antenna));
		world->add(std::move(wave));
		world->add(std::move(tx));
		world->add(std::move(streaming_rx));
		world->add(std::move(pulsed_rx));
		return world;
	}
}

TEST_CASE("Simulation memory projection totals core categories", "[core][memory_projection]")
{
	ParamGuard guard;
	params::setTime(0.0, 1.0);
	params::setRate(10.0);
	params::setOversampleRatio(2);

	const auto world = makeProjectionWorld();
	const auto projection = core::projectSimulationMemory(*world);

	REQUIRE(projection.streaming_sample_count == 20);
	REQUIRE(projection.phase_noise_sample_count == 21);
	REQUIRE(projection.phase_noise_timing_count == 1);
	REQUIRE(projection.enabled_phase_noise_timing_count == 1);
	REQUIRE(projection.streaming_receiver_count == 1);
	REQUIRE(projection.pulsed_receiver_count == 1);
	REQUIRE(projection.pulsed_window_count == 3);
	REQUIRE(projection.rendered_hdf5_sample_count == 16);

	REQUIRE(projection.phase_noise_lookup.bytes == 21 * sizeof(RealType));
	REQUIRE(projection.streaming_iq_buffers.bytes == 20 * sizeof(ComplexType));
	REQUIRE(projection.allocated_streaming_iq_buffers.bytes == 20 * sizeof(ComplexType));
	REQUIRE(projection.rendered_hdf5_payload.bytes == 16 * 2 * sizeof(RealType));

	const auto json = nlohmann::json::parse(core::memoryProjectionToJsonString(projection));
	REQUIRE(json["phase_noise_lookups"]["bytes"] == projection.phase_noise_lookup.bytes);
	REQUIRE(json["streaming_iq_buffers"]["bytes"] == projection.streaming_iq_buffers.bytes);
	REQUIRE(json["rendered_hdf5_dataset_payload"]["bytes"] == projection.rendered_hdf5_payload.bytes);
	REQUIRE(json["resident_baseline"].contains("bytes"));
}

TEST_CASE("Simulation memory projection log names required categories", "[core][memory_projection][logging]")
{
	ParamGuard guard;
	params::setTime(0.0, 1.0);
	params::setRate(10.0);
	params::setOversampleRatio(2);

	const auto world = makeProjectionWorld();
	LogLevelGuard level_guard(logging::Level::DEBUG);
	CerrCapture capture;

	core::logSimulationMemoryProjection(*world);

	const std::string output = capture.str();
	REQUIRE_THAT(output, ContainsSubstring("Projected simulation footprint"));
	REQUIRE_THAT(output, ContainsSubstring("phase_noise_lookup_memory"));
	REQUIRE_THAT(output, ContainsSubstring("streaming_iq_buffer_memory"));
	REQUIRE_THAT(output, ContainsSubstring("rendered_hdf5_dataset_payload"));
	REQUIRE_THAT(output, ContainsSubstring("resident_baseline"));
}
