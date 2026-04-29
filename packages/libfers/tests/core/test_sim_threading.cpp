#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/logging.h"
#include "core/parameters.h"
#include "core/sim_events.h"
#include "core/sim_threading.h"
#include "core/thread_pool.h"
#include "core/world.h"
#include "math/coord.h"
#include "math/geometry_ops.h"
#include "math/path.h"
#include "math/rotation_path.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "signal/radar_signal.h"
#include "simulation/channel_model.h"
#include "timing/prototype_timing.h"
#include "timing/timing.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace
{
	/**
	 * @struct ParamGuard
	 * @brief RAII guard to ensure global parameters are restored after each test.
	 */
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

	/// Initializes a platform with constant position and zero rotation.
	void setupStaticPlatform(radar::Platform& platform, const math::Vec3& position)
	{
		platform.getMotionPath()->addCoord(math::Coord{position, 0.0});
		platform.getMotionPath()->finalize();
		platform.getRotationPath()->addCoord(math::RotationCoord{0.0, 0.0, 0.0});
		platform.getRotationPath()->finalize();
	}

	/// Wraps a phase delta into the [-pi, pi] interval.
	RealType unwrapDelta(RealType delta)
	{
		while (delta > PI)
		{
			delta -= 2.0 * PI;
		}
		while (delta < -PI)
		{
			delta += 2.0 * PI;
		}
		return delta;
	}

	/**
	 * @brief Helper to create a fully configured World with basic physics objects.
	 * @return A unique pointer to a populated World.
	 */
	std::unique_ptr<core::World> createPhysicsWorld()
	{
		auto world = std::make_unique<core::World>();

		// Setup Platforms: Tx at origin, Rx at (100, 0, 0), Target at (50, 0, 0)
		auto tx_plat = std::make_unique<radar::Platform>("TxPlat", 10);
		tx_plat->getMotionPath()->addCoord(math::Coord{math::Vec3(0.0, 0.0, 0.0), 0.0});
		tx_plat->getMotionPath()->addCoord(math::Coord{math::Vec3(0.0, 0.0, 0.0), 100.0});
		tx_plat->getMotionPath()->finalize();
		tx_plat->getRotationPath()->addCoord(math::RotationCoord(0.0, 0.0, 0.0));
		tx_plat->getRotationPath()->addCoord(math::RotationCoord(0.0, 0.0, 100.0));
		tx_plat->getRotationPath()->finalize();

		auto rx_plat = std::make_unique<radar::Platform>("RxPlat", 11);
		rx_plat->getMotionPath()->addCoord(math::Coord{math::Vec3(100.0, 0.0, 0.0), 0.0});
		rx_plat->getMotionPath()->addCoord(math::Coord{math::Vec3(100.0, 0.0, 0.0), 100.0});
		rx_plat->getMotionPath()->finalize();
		rx_plat->getRotationPath()->addCoord(math::RotationCoord(0.0, 0.0, 0.0));
		rx_plat->getRotationPath()->addCoord(math::RotationCoord(0.0, 0.0, 100.0));
		rx_plat->getRotationPath()->finalize();

		auto tgt_plat = std::make_unique<radar::Platform>("TgtPlat", 12);
		tgt_plat->getMotionPath()->addCoord(math::Coord{math::Vec3(50.0, 0.0, 0.0), 0.0});
		tgt_plat->getMotionPath()->addCoord(math::Coord{math::Vec3(50.0, 0.0, 0.0), 100.0});
		tgt_plat->getMotionPath()->finalize();
		tgt_plat->getRotationPath()->addCoord(math::RotationCoord(0.0, 0.0, 0.0));
		tgt_plat->getRotationPath()->addCoord(math::RotationCoord(0.0, 0.0, 100.0));
		tgt_plat->getRotationPath()->finalize();

		// Setup Timing (0 offset for coherent math)
		auto proto_timing = std::make_unique<timing::PrototypeTiming>("Clock", 1);
		proto_timing->setFrequency(1e9); // 1 GHz carrier

		auto timing = std::make_shared<timing::Timing>("ClockInstance", 42);
		timing->initializeModel(proto_timing.get());

		// Setup Antenna (Isotropic, Gain = 1.0)
		auto antenna = std::make_unique<antenna::Isotropic>("IsoAnt", 2);

		// Setup Signal (1 Watt CW)
		auto signal = std::make_unique<fers_signal::RadarSignal>("CWWave", 1.0, 1e9, 1.0,
																 std::make_unique<fers_signal::CwSignal>(), 3);

		// Setup Transmitter
		auto tx = std::make_unique<radar::Transmitter>(tx_plat.get(), "Tx", radar::OperationMode::CW_MODE, 4);
		tx->setTiming(timing);
		tx->setAntenna(antenna.get());
		tx->setSignal(signal.get());
		tx->setPrf(1000.0);

		// Setup Receiver
		auto rx = std::make_unique<radar::Receiver>(rx_plat.get(), "Rx", 42, radar::OperationMode::CW_MODE, 5);
		rx->setTiming(timing);
		rx->setAntenna(antenna.get());
		rx->setWindowProperties(0.001, 1000.0, 0.0);

		// Setup Target (RCS = 1.0 m^2)
		auto target = std::make_unique<radar::IsoTarget>(tgt_plat.get(), "Tgt", 1.0, 42, 6);

		world->add(std::move(tx_plat));
		world->add(std::move(rx_plat));
		world->add(std::move(tgt_plat));
		world->add(std::move(proto_timing));
		world->add(std::move(antenna));
		world->add(std::move(signal));
		world->add(std::move(tx));
		world->add(std::move(rx));
		world->add(std::move(target));

		return world;
	}

	std::unique_ptr<core::World> createFmcwLoggingWorld(const std::optional<std::size_t> chirp_count,
														std::vector<radar::SchedulePeriod> schedule = {})
	{
		auto world = std::make_unique<core::World>();

		auto platform = std::make_unique<radar::Platform>("FmcwPlatform", 100);
		auto* platform_ptr = platform.get();
		auto timing_proto = std::make_unique<timing::PrototypeTiming>("FmcwClock", 101);
		timing_proto->setFrequency(1.0e6);
		auto timing = std::make_shared<timing::Timing>("FmcwClockInstance", 42);
		timing->initializeModel(timing_proto.get());
		auto antenna = std::make_unique<antenna::Isotropic>("FmcwAntenna", 102);
		auto* antenna_ptr = antenna.get();
		auto fmcw = std::make_unique<fers_signal::FmcwChirpSignal>(1.0e6, 0.04, 0.1, 0.0, chirp_count);
		auto wave = std::make_unique<fers_signal::RadarSignal>("FmcwWave", 10.0, 1.0e9, 0.04, std::move(fmcw), 103);
		auto* wave_ptr = wave.get();
		auto transmitter =
			std::make_unique<radar::Transmitter>(platform_ptr, "FmcwTx", radar::OperationMode::FMCW_MODE, 104);
		transmitter->setTiming(timing);
		transmitter->setAntenna(antenna_ptr);
		transmitter->setSignal(wave_ptr);
		if (!schedule.empty())
		{
			transmitter->setSchedule(std::move(schedule));
		}

		world->add(std::move(platform));
		world->add(std::move(timing_proto));
		world->add(std::move(antenna));
		world->add(std::move(wave));
		world->add(std::move(transmitter));
		return world;
	}

	simulation::CwPhaseNoiseLookup makeLookup(const std::shared_ptr<timing::Timing>& timing)
	{
		const std::vector<std::shared_ptr<timing::Timing>> timings = {timing};
		return simulation::CwPhaseNoiseLookup::build(timings, params::startTime(), params::endTime());
	}
}

TEST_CASE("ProgressReporter safely wraps and calls callback", "[core][threading]")
{
	int call_count = 0;
	std::string last_msg;
	int last_current = 0;
	int last_total = 0;

	auto cb = [&](const std::string& msg, int current, int total)
	{
		call_count++;
		last_msg = msg;
		last_current = current;
		last_total = total;
	};

	core::ProgressReporter reporter(cb);

	SECTION("Valid callback is executed")
	{
		reporter.report("Test", 50, 100);
		REQUIRE(call_count == 1);
		REQUIRE(last_msg == "Test");
		REQUIRE(last_current == 50);
		REQUIRE(last_total == 100);
	}

	SECTION("Null callback does not crash")
	{
		core::ProgressReporter null_reporter(nullptr);
		REQUIRE_NOTHROW(null_reporter.report("Test", 0, 0));
	}
}

TEST_CASE("makeActiveSource caches streaming scalars and clips FMCW chirp count", "[core][threading][fmcw]")
{
	ParamGuard guard;

	radar::Platform platform("TxPlatform");
	antenna::Isotropic antenna("Iso");
	auto timing = std::make_shared<timing::Timing>("Clock", 42);

	radar::Transmitter tx(&platform, "FmcwTx", radar::OperationMode::FMCW_MODE, 101);
	tx.setAntenna(&antenna);
	tx.setTiming(timing);
	auto fmcw_signal = std::make_unique<fers_signal::FmcwChirpSignal>(20.0e6, 100.0e-6, 250.0e-6, 1.0e6, 1000);
	fers_signal::RadarSignal wave("FmcwWave", 16.0, 10.0e9, 100.0e-6, std::move(fmcw_signal), 301);
	tx.setSignal(&wave);

	const core::ActiveStreamingSource source = core::makeActiveSource(&tx, 5.0, 65.0);

	REQUIRE(source.transmitter == &tx);
	REQUIRE(source.is_fmcw);
	REQUIRE(source.fmcw == wave.getFmcwChirpSignal());
	REQUIRE_THAT(source.segment_end, WithinAbs(5.25, 1.0e-12));
	REQUIRE_THAT(source.carrier_freq, WithinAbs(10.0e9, 1.0e-6));
	REQUIRE_THAT(source.amplitude, WithinAbs(4.0, 1.0e-12));
	REQUIRE_THAT(source.chirp_duration, WithinAbs(100.0e-6, 1.0e-15));
	REQUIRE_THAT(source.chirp_period, WithinAbs(250.0e-6, 1.0e-15));
	REQUIRE_THAT(source.chirp_rate, WithinAbs(200.0e9, 1.0e-3));
	REQUIRE_THAT(source.signed_chirp_rate, WithinAbs(200.0e9, 1.0e-3));
	REQUIRE_THAT(source.start_freq_off, WithinAbs(1.0e6, 1.0e-9));
	REQUIRE_THAT(source.two_pi_f0, WithinAbs(2.0 * PI * 1.0e6, 1.0e-9));
	REQUIRE_THAT(source.s_pi_alpha, WithinAbs(PI * 200.0e9, 1.0e-3));
	REQUIRE(source.chirp_count.has_value());
	REQUIRE(*source.chirp_count == std::size_t{1000});
}

TEST_CASE("makeActiveSource caches signed FMCW down-chirp coefficient", "[core][threading][fmcw]")
{
	ParamGuard guard;

	radar::Platform platform("TxPlatform");
	antenna::Isotropic antenna("Iso");
	auto timing = std::make_shared<timing::Timing>("Clock", 42);

	radar::Transmitter tx(&platform, "FmcwTx", radar::OperationMode::FMCW_MODE, 101);
	tx.setAntenna(&antenna);
	tx.setTiming(timing);
	auto fmcw_signal = std::make_unique<fers_signal::FmcwChirpSignal>(20.0e6, 100.0e-6, 250.0e-6, 1.0e6, std::nullopt,
																	  fers_signal::FmcwChirpDirection::Down);
	fers_signal::RadarSignal wave("FmcwWave", 16.0, 10.0e9, 100.0e-6, std::move(fmcw_signal), 301);
	tx.setSignal(&wave);

	const core::ActiveStreamingSource source = core::makeActiveSource(&tx, 5.0, 65.0);

	REQUIRE(source.is_fmcw);
	REQUIRE_THAT(source.chirp_rate, WithinAbs(200.0e9, 1.0e-3));
	REQUIRE_THAT(source.signed_chirp_rate, WithinAbs(-200.0e9, 1.0e-3));
	REQUIRE_THAT(source.s_pi_alpha, WithinAbs(-PI * 200.0e9, 1.0e-3));
}

TEST_CASE("SimulationEngine handles streaming state events", "[core][threading]")
{
	ParamGuard guard;
	auto world = createPhysicsWorld();
	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr, ".");

	auto* tx = world->getTransmitters().front().get();
	auto* rx = world->getReceivers().front().get();

	SECTION("Tx streaming start adds to active list")
	{
		engine.handleTxStreamingStart(core::makeActiveSource(tx, params::startTime(), params::endTime()));
		REQUIRE(world->getSimulationState().active_streaming_transmitters.size() == 1);
		REQUIRE(world->getSimulationState().active_streaming_transmitters[0].transmitter == tx);
	}

	SECTION("Tx streaming end keeps source as in-flight candidate")
	{
		engine.handleTxStreamingStart(core::makeActiveSource(tx, params::startTime(), params::endTime()));
		engine.handleTxStreamingEnd(tx);
		REQUIRE(world->getSimulationState().active_streaming_transmitters.size() == 1);
	}

	SECTION("Rx streaming start sets active")
	{
		REQUIRE_FALSE(rx->isActive());
		engine.handleRxStreamingStart(rx);
		REQUIRE(rx->isActive());
	}

	SECTION("Rx streaming end clears active")
	{
		engine.handleRxStreamingStart(rx);
		engine.handleRxStreamingEnd(rx);
		REQUIRE_FALSE(rx->isActive());
	}
}

TEST_CASE("SimulationEngine handles Pulsed Window events", "[core][threading]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	auto world = createPhysicsWorld();
	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr, ".");

	// Reconfigure as pulsed
	auto pulsed_rx_plat = std::make_unique<radar::Platform>("PulsedRxPlat", 20);
	pulsed_rx_plat->getMotionPath()->addCoord(math::Coord{math::Vec3(0.0, 0.0, 0.0), 0.0});
	pulsed_rx_plat->getMotionPath()->addCoord(math::Coord{math::Vec3(0.0, 0.0, 0.0), 100.0});
	pulsed_rx_plat->getMotionPath()->finalize();
	pulsed_rx_plat->getRotationPath()->addCoord(math::RotationCoord(0.0, 0.0, 0.0));
	pulsed_rx_plat->getRotationPath()->addCoord(math::RotationCoord(0.0, 0.0, 100.0));
	pulsed_rx_plat->getRotationPath()->finalize();

	auto pulsed_rx =
		std::make_unique<radar::Receiver>(pulsed_rx_plat.get(), "PulsedRx", 42, radar::OperationMode::PULSED_MODE, 21);
	pulsed_rx->setWindowProperties(0.001, 100.0, 0.0);

	auto proto_timing = std::make_unique<timing::PrototypeTiming>("RxClock", 22);
	auto timing = std::make_shared<timing::Timing>("RxClockInst", 99);
	timing->initializeModel(proto_timing.get());
	pulsed_rx->setTiming(timing);

	auto* rx_ptr = pulsed_rx.get();
	world->add(std::move(pulsed_rx_plat));
	world->add(std::move(proto_timing));
	world->add(std::move(pulsed_rx));

	SECTION("Window Start sets active and schedules End")
	{
		engine.handleRxPulsedWindowStart(rx_ptr, 1.0);
		REQUIRE(rx_ptr->isActive());

		auto& queue = world->getEventQueue();
		REQUIRE(queue.size() == 1);
		REQUIRE(queue.top().type == core::EventType::RX_PULSED_WINDOW_END);
		REQUIRE_THAT(queue.top().timestamp, WithinAbs(1.0 + 0.001, 1e-9));
	}

	SECTION("Window End clears active, enqueues job, and schedules next Start")
	{
		rx_ptr->setActive(true);
		engine.handleRxPulsedWindowEnd(rx_ptr, 1.001);
		REQUIRE_FALSE(rx_ptr->isActive());

		// Check event queue for next start
		auto& queue = world->getEventQueue();
		REQUIRE(queue.size() == 1);
		REQUIRE(queue.top().type == core::EventType::RX_PULSED_WINDOW_START);
		// Next theoretical = 1.001 - 0.001 + (1/100) = 1.01
		REQUIRE_THAT(queue.top().timestamp, WithinAbs(1.01, 1e-9));

		// Check finalizer job queue
		core::RenderingJob job;
		REQUIRE(rx_ptr->waitAndDequeueFinalizerJob(job));
		REQUIRE_THAT(job.duration, WithinAbs(0.001, 1e-9));
	}
}

TEST_CASE("SimulationEngine handles Tx Pulsed Start and routes responses", "[core][threading]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	auto world = createPhysicsWorld();
	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr, ".");

	auto* tx = world->getTransmitters().front().get();
	auto* cw_rx = world->getReceivers().front().get();

	// Add a pulsed receiver to test routing to inbox
	auto pulsed_rx_plat = std::make_unique<radar::Platform>("PulsedRxPlat", 30);
	pulsed_rx_plat->getMotionPath()->addCoord(math::Coord{math::Vec3(200.0, 0.0, 0.0), 0.0});
	pulsed_rx_plat->getMotionPath()->addCoord(math::Coord{math::Vec3(200.0, 0.0, 0.0), 100.0});
	pulsed_rx_plat->getMotionPath()->finalize();
	pulsed_rx_plat->getRotationPath()->addCoord(math::RotationCoord(0.0, 0.0, 0.0));
	pulsed_rx_plat->getRotationPath()->addCoord(math::RotationCoord(0.0, 0.0, 100.0));
	pulsed_rx_plat->getRotationPath()->finalize();

	auto pulsed_rx =
		std::make_unique<radar::Receiver>(pulsed_rx_plat.get(), "PulsedRx", 42, radar::OperationMode::PULSED_MODE, 31);
	pulsed_rx->setTiming(tx->getTiming());
	pulsed_rx->setAntenna(tx->getAntenna());
	auto* pulsed_rx_ptr = pulsed_rx.get();

	world->add(std::move(pulsed_rx_plat));
	world->add(std::move(pulsed_rx));

	SECTION("Calculates responses and routes to correct receiver bin")
	{
		engine.handleTxPulsedStart(tx, 1.0);

		// CW Receiver should get interference logs (1 direct + 1 reflected)
		REQUIRE(cw_rx->getPulsedInterferenceLog().size() == 2);

		// Pulsed Receiver should get inbox responses (1 direct + 1 reflected)
		auto inbox = pulsed_rx_ptr->drainInbox();
		REQUIRE(inbox.size() == 2);

		// Check that next pulse was scheduled
		auto& queue = world->getEventQueue();
		REQUIRE(queue.size() == 1);
		REQUIRE(queue.top().type == core::EventType::TX_PULSED_START);
		REQUIRE_THAT(queue.top().timestamp, WithinAbs(1.0 + 1.0 / 1000.0, 1e-9));
	}
}

TEST_CASE("SimulationEngine calculates mathematically correct CW physics", "[core][threading]")
{
	ParamGuard guard;
	params::setC(3e8); // Override c to 3e8 for clean math

	auto world = createPhysicsWorld();
	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr, ".");

	auto* tx = world->getTransmitters().front().get();
	auto* rx = world->getReceivers().front().get();
	SECTION("Calculates exact Direct and Reflected Path complex envelope")
	{
		// Math Verification:
		// c = 3e8, f = 1e9 -> lambda = 0.3m
		// Tx=(0,0,0), Rx=(100,0,0), Tgt=(50,0,0)
		//
		// DIRECT PATH:
		// d = 100m.
		// Pr/Pt = (Gt * Gr * lambda^2) / (16 * pi^2 * d^2)
		// Pr/Pt = (1 * 1 * 0.09) / (16 * pi^2 * 10000) = 0.09 / (160000 * pi^2)
		// Amplitude = sqrt(Pr/Pt) = 0.3 / (400 * pi) = 0.0002387324146
		// Delay = 100 / 3e8 = 3.33333333e-7 s
		// Phase = -2 * pi * f * delay = -2 * pi * 1e9 * (100/3e8) = -2000*pi/3 rad
		//
		// REFLECTED PATH:
		// Rtx = 50m, Rrx = 50m. RCS = 1.0
		// Pr/Pt = (Gt * Gr * RCS * lambda^2) / (64 * pi^3 * Rtx^2 * Rrx^2)
		// Pr/Pt = (1 * 1 * 1 * 0.09) / (64 * pi^3 * 2500 * 2500) = 0.09 / (400,000,000 * pi^3)
		// Amplitude = sqrt(Pr/Pt) = 0.3 / (20000 * pi * sqrt(pi)) = 0.00000268689
		// Delay = (50+50) / 3e8 = 100 / 3e8 = 3.33333333e-7 s
		// Phase = -2000*pi/3 rad (Same as direct path since total distance is the same)

		const RealType expected_direct_amp = 0.3 / (400.0 * PI);
		const RealType expected_refl_amp = 0.3 / (20000.0 * PI * std::sqrt(PI));
		const RealType expected_phase = -2000.0 * PI / 3.0;

		ComplexType expected_direct = std::polar(expected_direct_amp, expected_phase);
		ComplexType expected_refl = std::polar(expected_refl_amp, expected_phase);
		ComplexType expected_total = expected_direct + expected_refl;

		ComplexType actual_total = simulation::calculateDirectPathContribution(tx, rx, 0.0) +
			simulation::calculateReflectedPathContribution(tx, rx, world->getTargets().front().get(), 0.0);

		REQUIRE_THAT(actual_total.real(), WithinAbs(expected_total.real(), 1e-12));
		REQUIRE_THAT(actual_total.imag(), WithinAbs(expected_total.imag(), 1e-12));
	}
}

TEST_CASE("SimulationEngine processStreamingPhysics steps through time and updates buffers", "[core][threading]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	auto world = createPhysicsWorld();
	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr, ".");

	auto* tx = world->getTransmitters().front().get();
	auto* rx = world->getReceivers().front().get();

	// Prepare receiver buffer
	rx->prepareStreamingData(10);
	rx->setActive(true);

	// Set state
	world->getSimulationState().t_current = 0.0;
	engine.handleTxStreamingStart(core::makeActiveSource(tx, params::startTime(), params::endTime()));

	SECTION("Processes correct number of samples based on dt")
	{
		// dt = 1/1000 = 0.001s. Processing up to t=0.0025 should hit indices 0, 1, 2.
		engine.processStreamingPhysics(0.0025);

		const auto& buffer = rx->getStreamingData();
		// Sample 0 arrives before the propagation delay; samples 1 and 2 contain retarded-time energy.
		REQUIRE(std::abs(buffer[0]) == 0.0);
		REQUIRE(std::abs(buffer[1]) > 0.0);
		REQUIRE(std::abs(buffer[2]) > 0.0);

		// Sample 3 should be untouched (0.0)
		REQUIRE(std::abs(buffer[3]) == 0.0);
	}
}

TEST_CASE("SimulationEngine keeps streaming source through propagation tail after transmit end", "[core][threading]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 0.4);
	params::setC(1000.0);

	auto world = createPhysicsWorld();
	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr, ".");

	auto* tx = world->getTransmitters().front().get();
	auto* rx = world->getReceivers().front().get();

	rx->prepareStreamingData(400);
	rx->setActive(true);

	const core::ActiveStreamingSource source = core::makeActiveSource(tx, 0.0, 0.2);
	engine.handleTxStreamingStart(source);
	engine.handleTxStreamingEnd(tx);
	REQUIRE(world->getSimulationState().active_streaming_transmitters.size() == 1);

	world->getSimulationState().t_current = 0.2;
	engine.processStreamingPhysics(0.205);

	bool saw_tail_energy = false;
	for (std::size_t index = 200; index < 205; ++index)
	{
		saw_tail_energy = saw_tail_energy || std::abs(rx->getStreamingData()[index]) > 0.0;
	}
	REQUIRE(saw_tail_energy);

	world->getSimulationState().t_current = 0.31;
	engine.processStreamingPhysics(0.315);

	for (std::size_t index = 310; index < 315; ++index)
	{
		REQUIRE(std::abs(rx->getStreamingData()[index]) == 0.0);
	}
}

TEST_CASE("SimulationEngine processStreamingPhysics handles active streaming receiver without streaming transmitters",
		  "[core][threading]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	auto world = createPhysicsWorld();
	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr, ".");

	auto* rx = world->getReceivers().front().get();
	rx->prepareStreamingData(10);
	rx->setActive(true);
	world->getSimulationState().t_current = 0.0;

	REQUIRE(world->getSimulationState().active_streaming_transmitters.empty());
	REQUIRE_NOTHROW(engine.processStreamingPhysics(0.0025));

	const auto& buffer = rx->getStreamingData();
	for (const auto& sample : buffer)
	{
		REQUIRE(sample.real() == 0.0);
		REQUIRE(sample.imag() == 0.0);
	}
}

TEST_CASE("SimulationEngine processStreamingPhysics uses buffered shared timing for streaming samples",
		  "[core][threading]")
{
	ParamGuard guard;
	params::setRate(10.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 1.0);
	params::setC(10.0);

	auto world = std::make_unique<core::World>();

	auto tx_plat = std::make_unique<radar::Platform>("TxPlat", 10);
	tx_plat->getMotionPath()->addCoord(math::Coord{math::Vec3{0.0, 0.0, 0.0}, 0.0});
	tx_plat->getMotionPath()->finalize();
	tx_plat->getRotationPath()->addCoord(math::RotationCoord{0.0, 0.0, 0.0});
	tx_plat->getRotationPath()->finalize();

	auto rx_plat = std::make_unique<radar::Platform>("RxPlat", 11);
	rx_plat->getMotionPath()->addCoord(math::Coord{math::Vec3{1.5, 0.0, 0.0}, 0.0});
	rx_plat->getMotionPath()->finalize();
	rx_plat->getRotationPath()->addCoord(math::RotationCoord{0.0, 0.0, 0.0});
	rx_plat->getRotationPath()->finalize();

	auto proto_timing = std::make_unique<timing::PrototypeTiming>("Clock", 1);
	proto_timing->setFrequency(1.0);
	proto_timing->setFreqOffset(1.0);
	auto timing = std::make_shared<timing::Timing>("ClockInstance", 42, 1);
	timing->initializeModel(proto_timing.get());
	const auto lookup = makeLookup(timing);

	auto antenna = std::make_unique<antenna::Isotropic>("IsoAnt", 2);
	auto signal = std::make_unique<fers_signal::RadarSignal>("CWWave", 1.0, 1.0, 1.0,
															 std::make_unique<fers_signal::CwSignal>(), 3);

	auto tx = std::make_unique<radar::Transmitter>(tx_plat.get(), "Tx", radar::OperationMode::CW_MODE, 4);
	tx->setTiming(timing);
	tx->setAntenna(antenna.get());
	tx->setSignal(signal.get());

	auto rx = std::make_unique<radar::Receiver>(rx_plat.get(), "Rx", 42, radar::OperationMode::CW_MODE, 5);
	rx->setTiming(timing);
	rx->setAntenna(antenna.get());
	rx->prepareStreamingData(10);
	rx->setActive(true);

	auto* tx_ptr = tx.get();
	auto* rx_ptr = rx.get();

	world->add(std::move(tx_plat));
	world->add(std::move(rx_plat));
	world->add(std::move(proto_timing));
	world->add(std::move(antenna));
	world->add(std::move(signal));
	world->add(std::move(tx));
	world->add(std::move(rx));

	world->getSimulationState().t_current = 0.0;

	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr, ".");
	engine.handleTxStreamingStart(core::makeActiveSource(tx_ptr, params::startTime(), params::endTime()));
	engine.processStreamingPhysics(0.6);

	const ComplexType expected = simulation::calculateDirectPathContribution(tx_ptr, rx_ptr, 0.5, &lookup);
	const ComplexType actual = rx_ptr->getStreamingData()[5];

	REQUIRE_THAT(actual.real(), WithinAbs(expected.real(), 1e-6));
	REQUIRE_THAT(actual.imag(), WithinAbs(expected.imag(), 1e-6));
}

TEST_CASE("SimulationEngine phase-noise lookup covers pre-start retarded streaming emissions", "[core][threading]")
{
	ParamGuard guard;
	params::setRate(10.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 1.0);
	params::setC(10.0);

	auto world = std::make_unique<core::World>();
	auto tx_plat = std::make_unique<radar::Platform>("TxPlat", 10);
	setupStaticPlatform(*tx_plat, math::Vec3{0.0, 0.0, 0.0});
	auto rx_plat = std::make_unique<radar::Platform>("RxPlat", 11);
	setupStaticPlatform(*rx_plat, math::Vec3{1.5, 0.0, 0.0});

	auto proto_timing = std::make_unique<timing::PrototypeTiming>("Clock", 1);
	proto_timing->setFrequency(1.0);
	proto_timing->setFreqOffset(1.0);
	auto timing = std::make_shared<timing::Timing>("ClockInstance", 42, 1);
	timing->initializeModel(proto_timing.get());
	const std::vector<std::shared_ptr<timing::Timing>> timings = {timing};
	const auto expected_lookup = simulation::CwPhaseNoiseLookup::build(timings, -0.2, params::endTime());

	auto antenna = std::make_unique<antenna::Isotropic>("IsoAnt", 2);
	auto signal = std::make_unique<fers_signal::RadarSignal>("CWWave", 1.0, 1.0, 1.0,
															 std::make_unique<fers_signal::CwSignal>(), 3);
	auto tx = std::make_unique<radar::Transmitter>(tx_plat.get(), "Tx", radar::OperationMode::CW_MODE, 4);
	tx->setTiming(timing);
	tx->setAntenna(antenna.get());
	tx->setSignal(signal.get());
	auto rx = std::make_unique<radar::Receiver>(rx_plat.get(), "Rx", 42, radar::OperationMode::CW_MODE, 5);
	rx->setTiming(timing);
	rx->setAntenna(antenna.get());
	rx->prepareStreamingData(10);
	rx->setActive(true);

	auto* tx_ptr = tx.get();
	auto* rx_ptr = rx.get();
	world->add(std::move(tx_plat));
	world->add(std::move(rx_plat));
	world->add(std::move(proto_timing));
	world->add(std::move(antenna));
	world->add(std::move(signal));
	world->add(std::move(tx));
	world->add(std::move(rx));
	world->getSimulationState().t_current = 0.0;

	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr, ".");
	engine.handleTxStreamingStart(core::makeActiveSource(tx_ptr, -0.2, params::endTime()));
	engine.processStreamingPhysics(0.2);

	const ComplexType expected = simulation::calculateDirectPathContribution(tx_ptr, rx_ptr, 0.1, &expected_lookup);
	const ComplexType actual = rx_ptr->getStreamingData()[1];
	REQUIRE_THAT(actual.real(), WithinAbs(expected.real(), 1.0e-6));
	REQUIRE_THAT(actual.imag(), WithinAbs(expected.imag(), 1.0e-6));
}

TEST_CASE("SimulationEngine native FMCW dechirp produces positive stationary-target beat",
		  "[core][threading][fmcw][dechirp]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(1.0e6);
	params::setOversampleRatio(1);
	params::setSimSamplingRate(1.0e6);
	params::setTime(0.0, 1.0e-3);

	const RealType target_range = 150.0;
	const RealType chirp_bandwidth = 1.0e6;
	const RealType chirp_duration = 1.0e-3;
	const RealType chirp_rate = chirp_bandwidth / chirp_duration;
	const RealType tau = (2.0 * target_range) / params::c();

	auto world = std::make_unique<core::World>();
	auto radar_platform = std::make_unique<radar::Platform>("RadarPlatform", 10);
	setupStaticPlatform(*radar_platform, math::Vec3{0.0, 0.0, 0.0});
	auto target_platform = std::make_unique<radar::Platform>("TargetPlatform", 11);
	setupStaticPlatform(*target_platform, math::Vec3{target_range, 0.0, 0.0});
	auto timing_proto = std::make_unique<timing::PrototypeTiming>("Clock", 1);
	timing_proto->setFrequency(10.0e6);
	auto timing = std::make_shared<timing::Timing>("ClockInstance", 42, 1);
	timing->initializeModel(timing_proto.get());
	auto antenna = std::make_unique<antenna::Isotropic>("IsoAnt", 2);
	auto fmcw = std::make_unique<fers_signal::FmcwChirpSignal>(chirp_bandwidth, chirp_duration, chirp_duration);
	auto signal =
		std::make_unique<fers_signal::RadarSignal>("FmcwWave", 1.0, 10.0e6, chirp_duration, std::move(fmcw), 3);

	auto tx = std::make_unique<radar::Transmitter>(radar_platform.get(), "Tx", radar::OperationMode::FMCW_MODE, 4);
	tx->setTiming(timing);
	tx->setAntenna(antenna.get());
	tx->setSignal(signal.get());

	auto rx = std::make_unique<radar::Receiver>(radar_platform.get(), "Rx", 42, radar::OperationMode::FMCW_MODE, 5);
	rx->setTiming(timing);
	rx->setAntenna(antenna.get());
	rx->setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
	rx->setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
	rx->setDechirpMode(radar::Receiver::DechirpMode::Ideal);
	radar::Receiver::DechirpReference reference;
	reference.source = radar::Receiver::DechirpReferenceSource::Attached;
	rx->setDechirpReference(reference);
	rx->prepareStreamingData(1000);
	rx->setActive(true);

	auto target = radar::createIsoTarget(target_platform.get(), "Target", 1.0, 7, 6);
	auto* tx_ptr = tx.get();
	auto* rx_ptr = rx.get();
	tx->setAttached(rx_ptr);
	rx->setAttached(tx_ptr);

	world->add(std::move(radar_platform));
	world->add(std::move(target_platform));
	world->add(std::move(timing_proto));
	world->add(std::move(antenna));
	world->add(std::move(signal));
	world->add(std::move(tx));
	world->add(std::move(rx));
	world->add(std::move(target));
	world->resolveReceiverDechirpReferences();
	world->getSimulationState().t_current = 0.0;

	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr, ".");
	engine.handleTxStreamingStart(core::makeActiveSource(tx_ptr, params::startTime(), params::endTime()));
	engine.processStreamingPhysics(0.0007);

	const RealType dt = 1.0 / params::simSamplingRate();
	const std::size_t first_index = 20;
	const std::size_t sample_count = 500;
	RealType unwrapped_span = 0.0;
	RealType previous_phase = std::arg(rx_ptr->getStreamingData()[first_index]);
	for (std::size_t i = 1; i < sample_count; ++i)
	{
		const RealType phase = std::arg(rx_ptr->getStreamingData()[first_index + i]);
		unwrapped_span += unwrapDelta(phase - previous_phase);
		previous_phase = phase;
	}

	const RealType measured_beat_hz = unwrapped_span / (2.0 * PI * dt * static_cast<RealType>(sample_count - 1));
	REQUIRE_THAT(measured_beat_hz, WithinRel(chirp_rate * tau, 1.0e-3));
}

TEST_CASE("SimulationEngine physical FMCW dechirp keeps timing decorrelation absent from ideal mode",
		  "[core][threading][fmcw][dechirp]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(2.0e6);
	params::setOversampleRatio(1);
	params::setSimSamplingRate(2.0e6);
	params::setTime(0.0, 1.0e-3);

	const RealType target_range = 150.0;
	const RealType chirp_duration = 1.0e-3;
	const RealType freq_offset = 250.0e3;
	const RealType tau = (2.0 * target_range) / params::c();
	const std::size_t sample_index = 80;
	const RealType sample_time = static_cast<RealType>(sample_index) / params::simSamplingRate();

	const auto run_sample = [&](const radar::Receiver::DechirpMode mode)
	{
		auto world = std::make_unique<core::World>();
		auto radar_platform = std::make_unique<radar::Platform>("RadarPlatform", 10);
		setupStaticPlatform(*radar_platform, math::Vec3{0.0, 0.0, 0.0});
		auto target_platform = std::make_unique<radar::Platform>("TargetPlatform", 11);
		setupStaticPlatform(*target_platform, math::Vec3{target_range, 0.0, 0.0});
		auto timing_proto = std::make_unique<timing::PrototypeTiming>("Clock", 1);
		timing_proto->setFrequency(10.0e6);
		timing_proto->setFreqOffset(freq_offset);
		auto timing = std::make_shared<timing::Timing>("ClockInstance", 42, 1);
		timing->initializeModel(timing_proto.get());
		auto antenna = std::make_unique<antenna::Isotropic>("IsoAnt", 2);
		auto fmcw = std::make_unique<fers_signal::FmcwChirpSignal>(1.0e6, chirp_duration, chirp_duration);
		auto signal =
			std::make_unique<fers_signal::RadarSignal>("FmcwWave", 1.0, 10.0e6, chirp_duration, std::move(fmcw), 3);

		auto tx = std::make_unique<radar::Transmitter>(radar_platform.get(), "Tx", radar::OperationMode::FMCW_MODE, 4);
		tx->setTiming(timing);
		tx->setAntenna(antenna.get());
		tx->setSignal(signal.get());
		auto rx = std::make_unique<radar::Receiver>(radar_platform.get(), "Rx", 42, radar::OperationMode::FMCW_MODE, 5);
		rx->setTiming(timing);
		rx->setAntenna(antenna.get());
		rx->setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
		rx->setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
		rx->setDechirpMode(mode);
		radar::Receiver::DechirpReference reference;
		reference.source = radar::Receiver::DechirpReferenceSource::Attached;
		rx->setDechirpReference(reference);
		rx->prepareStreamingData(sample_index + 2);
		rx->setActive(true);

		auto target = radar::createIsoTarget(target_platform.get(), "Target", 1.0, 7, 6);
		auto* tx_ptr = tx.get();
		auto* rx_ptr = rx.get();
		tx->setAttached(rx_ptr);
		rx->setAttached(tx_ptr);

		world->add(std::move(radar_platform));
		world->add(std::move(target_platform));
		world->add(std::move(timing_proto));
		world->add(std::move(antenna));
		world->add(std::move(signal));
		world->add(std::move(tx));
		world->add(std::move(rx));
		world->add(std::move(target));
		world->resolveReceiverDechirpReferences();
		world->getSimulationState().t_current = 0.0;

		pool::ThreadPool pool(1);
		core::SimulationEngine engine(world.get(), pool, nullptr, ".");
		engine.handleTxStreamingStart(core::makeActiveSource(tx_ptr, params::startTime(), params::endTime()));
		engine.processStreamingPhysics(sample_time + 1.0 / params::simSamplingRate());
		return rx_ptr->getStreamingData()[sample_index];
	};

	const ComplexType physical = run_sample(radar::Receiver::DechirpMode::Physical);
	const ComplexType ideal = run_sample(radar::Receiver::DechirpMode::Ideal);
	const RealType measured_delta = std::arg(physical * std::conj(ideal));
	const RealType expected_delta = 2.0 * PI * freq_offset * tau;

	REQUIRE_THAT(measured_delta, WithinAbs(expected_delta, 1.0e-6));
}

TEST_CASE("SimulationEngine runEventDrivenSim executes full loop", "[core][threading]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 0.005); // Very short simulation

	auto world = createPhysicsWorld();
	pool::ThreadPool pool(1);

	// Add an event to ensure the loop runs
	auto* tx = world->getTransmitters().front().get();
	world->getEventQueue().push({0.001, core::EventType::TX_STREAMING_START, tx});
	world->getEventQueue().push({0.004, core::EventType::TX_STREAMING_END, tx});

	int progress_updates = 0;
	auto progress_cb = [&](const std::string&, int, int) { progress_updates++; };

	SECTION("Simulation runs to completion and shuts down cleanly")
	{
		REQUIRE_NOTHROW(core::runEventDrivenSim(world.get(), pool, progress_cb, "."));

		// Should have received at least initialization and completion updates
		REQUIRE(progress_updates >= 2);

		// Simulation time should have advanced to the last event
		REQUIRE_THAT(world->getSimulationState().t_current, WithinAbs(0.004, 1e-9));
	}
}

TEST_CASE("SimulationEngine reports progress while processing long streaming spans", "[core][threading]")
{
	ParamGuard guard;
	params::setRate(100.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 1.0);

	auto world = createPhysicsWorld();
	auto* rx = world->getReceivers().front().get();
	rx->prepareStreamingData(100);
	rx->setActive(true);

	pool::ThreadPool pool(1);
	std::vector<int> progress_updates;
	std::vector<std::string> messages;
	auto reporter = std::make_shared<core::ProgressReporter>(
		[&](const std::string& message, const int current, int)
		{
			messages.push_back(message);
			progress_updates.push_back(current);
		});

	core::SimulationEngine engine(world.get(), pool, reporter, ".");
	engine.processStreamingPhysics(1.0);

	REQUIRE_FALSE(progress_updates.empty());
	bool saw_intermediate_progress = false;
	for (const int progress : progress_updates)
	{
		if (progress > 0 && progress < 100)
		{
			saw_intermediate_progress = true;
			break;
		}
	}
	REQUIRE(saw_intermediate_progress);
	REQUIRE_THAT(messages.front(), ContainsSubstring("Simulating..."));
}

TEST_CASE("SimulationEngine logs FMCW derived chirp counts at startup", "[core][threading][fmcw]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);

	SECTION("unbounded always-on transmitter logs total started chirps over simulation span")
	{
		params::setTime(0.0, 0.26);
		auto world = createFmcwLoggingWorld(std::nullopt);
		pool::ThreadPool pool(1);
		LogLevelGuard log_level(logging::Level::INFO);
		CerrCapture capture;

		core::SimulationEngine engine(world.get(), pool, nullptr, ".");
		engine.run();

		const std::string output = capture.str();
		REQUIRE_THAT(output, ContainsSubstring("shape=linear up"));
		REQUIRE_THAT(output, ContainsSubstring("chirp_count=unbounded"));
		REQUIRE_THAT(output, ContainsSubstring("total_chirp_count=3"));
	}

	SECTION("scheduled capped transmitter logs per-segment and cumulative started chirps")
	{
		params::setTime(0.05, 0.45);
		auto world = createFmcwLoggingWorld(std::size_t{3}, {{-0.05, 0.25}, {0.3, 0.36}});
		pool::ThreadPool pool(1);
		LogLevelGuard log_level(logging::Level::INFO);
		CerrCapture capture;

		core::SimulationEngine engine(world.get(), pool, nullptr, ".");
		engine.run();

		const std::string output = capture.str();
		REQUIRE_THAT(output, ContainsSubstring("shape=linear up"));
		REQUIRE_THAT(output, ContainsSubstring("chirp_count=3"));
		REQUIRE_THAT(output, ContainsSubstring("segment_chirp_count=2"));
		REQUIRE_THAT(output, ContainsSubstring("total_chirp_count=2"));
		REQUIRE_THAT(output, ContainsSubstring("segment_chirp_count=1"));
		REQUIRE_THAT(output, ContainsSubstring("total_chirp_count=3"));
	}
}

TEST_CASE("SimulationEngine handles Pulsed receiver finalizer thread lifecycle", "[core][threading]")
{
	// This test covers:
	// 1. initializeFinalizers() -> _finalizer_threads.emplace_back(...)
	// 6. shutdown() -> else if (PULSED_MODE) { enqueue shutdown_job }
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 0.005);

	auto world = createPhysicsWorld();

	// Add a pulsed receiver to trigger the thread creation and shutdown logic
	auto rx_plat = std::make_unique<radar::Platform>("PulsedRxPlat", 100);
	rx_plat->getMotionPath()->addCoord(math::Coord{math::Vec3(0, 0, 0), 0});
	rx_plat->getMotionPath()->finalize();
	rx_plat->getRotationPath()->addCoord(math::RotationCoord(0, 0, 0));
	rx_plat->getRotationPath()->finalize();

	auto rx = std::make_unique<radar::Receiver>(rx_plat.get(), "PulsedRx", 42, radar::OperationMode::PULSED_MODE, 101);
	auto proto_timing = std::make_unique<timing::PrototypeTiming>("RxClock", 102);
	auto timing = std::make_shared<timing::Timing>("RxClockInst", 99);
	timing->initializeModel(proto_timing.get());
	rx->setTiming(timing);
	rx->setWindowProperties(0.001, 1000.0, 0.0);

	world->add(std::move(rx_plat));
	world->add(std::move(proto_timing));
	world->add(std::move(rx));

	pool::ThreadPool pool(1);

	// Running the engine will initialize the finalizer thread and then shut it down.
	// If the shutdown logic fails, the std::jthread will hang indefinitely waiting for a job.
	REQUIRE_NOTHROW(core::runEventDrivenSim(world.get(), pool, nullptr, "."));
}

TEST_CASE("SimulationEngine processStreamingPhysics exits early if t_event <= t_current", "[core][threading]")
{
	// This test covers:
	// 2. processStreamingPhysics(t_event) -> if (t_event <= t_current) { return; }
	ParamGuard guard;
	auto world = createPhysicsWorld();
	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr, ".");

	auto* rx = world->getReceivers().front().get();
	rx->prepareStreamingData(10);
	rx->setActive(true);

	// Advance the simulation state time
	world->getSimulationState().t_current = 1.0;

	// Call with t_event < t_current
	REQUIRE_NOTHROW(engine.processStreamingPhysics(0.5));

	// Call with t_event == t_current
	REQUIRE_NOTHROW(engine.processStreamingPhysics(1.0));

	// The buffer should be completely untouched (all zeros) because the method returned early
	for (const auto& sample : rx->getStreamingData())
	{
		REQUIRE(sample.real() == 0.0);
		REQUIRE(sample.imag() == 0.0);
	}
}

TEST_CASE("SimulationEngine processEvent dispatches all event types correctly", "[core][threading]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	auto world = createPhysicsWorld();
	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr, ".");

	auto* tx = world->getTransmitters().front().get();
	auto* rx = world->getReceivers().front().get();

	// Test TX_STREAMING_START
	engine.processEvent({1.0, core::EventType::TX_STREAMING_START, tx});
	REQUIRE(world->getSimulationState().active_streaming_transmitters.size() == 1);

	// Test TX_STREAMING_END
	engine.processEvent({2.0, core::EventType::TX_STREAMING_END, tx});
	REQUIRE(world->getSimulationState().active_streaming_transmitters.size() == 1);

	// Test RX_STREAMING_START
	engine.processEvent({3.0, core::EventType::RX_STREAMING_START, rx});
	REQUIRE(rx->isActive());

	// Test RX_STREAMING_END
	engine.processEvent({4.0, core::EventType::RX_STREAMING_END, rx});
	REQUIRE_FALSE(rx->isActive());

	// Test TX_PULSED_START
	auto queue_size_before = world->getEventQueue().size();
	engine.processEvent({5.0, core::EventType::TX_PULSED_START, tx});
	// Should schedule next pulse, increasing the queue size
	REQUIRE(world->getEventQueue().size() > queue_size_before);

	// Test RX_PULSED_WINDOW_START
	queue_size_before = world->getEventQueue().size();
	engine.processEvent({6.0, core::EventType::RX_PULSED_WINDOW_START, rx});
	REQUIRE(rx->isActive());
	// Should schedule the end of the window
	REQUIRE(world->getEventQueue().size() > queue_size_before);

	// Test RX_PULSED_WINDOW_END
	queue_size_before = world->getEventQueue().size();
	engine.processEvent({7.0, core::EventType::RX_PULSED_WINDOW_END, rx});
	REQUIRE_FALSE(rx->isActive());
	// Should schedule the next window start
	REQUIRE(world->getEventQueue().size() > queue_size_before);
}

TEST_CASE("SimulationEngine routeResponse handles null responses safely", "[core][threading]")
{
	// This test covers:
	// 4. routeResponse(...) -> if (!response) { return; }
	ParamGuard guard;
	auto world = std::make_unique<core::World>();

	// Put Tx and Rx on the EXACT SAME platform to force calculateResponse to return nullptr
	// (Simulation skips direct path calculations for co-located far-field antennas).
	auto plat = std::make_unique<radar::Platform>("SharedPlat", 200);
	plat->getMotionPath()->addCoord(math::Coord{math::Vec3(0, 0, 0), 0});
	plat->getMotionPath()->finalize();
	plat->getRotationPath()->addCoord(math::RotationCoord(0, 0, 0));
	plat->getRotationPath()->finalize();

	auto proto_timing = std::make_unique<timing::PrototypeTiming>("Clock", 201);
	auto timing = std::make_shared<timing::Timing>("ClockInst", 99);
	timing->initializeModel(proto_timing.get());

	auto antenna = std::make_unique<antenna::Isotropic>("IsoAnt", 202);
	auto signal = std::make_unique<fers_signal::RadarSignal>("CWWave", 1.0, 1e9, 1.0,
															 std::make_unique<fers_signal::CwSignal>(), 203);

	auto tx = std::make_unique<radar::Transmitter>(plat.get(), "Tx", radar::OperationMode::PULSED_MODE, 204);
	tx->setTiming(timing);
	tx->setAntenna(antenna.get());
	tx->setSignal(signal.get());
	tx->setPrf(1000.0);

	auto rx = std::make_unique<radar::Receiver>(plat.get(), "Rx", 42, radar::OperationMode::PULSED_MODE, 205);
	rx->setTiming(timing);
	rx->setAntenna(antenna.get());

	auto* tx_ptr = tx.get();
	auto* rx_ptr = rx.get();

	world->add(std::move(plat));
	world->add(std::move(proto_timing));
	world->add(std::move(antenna));
	world->add(std::move(signal));
	world->add(std::move(tx));
	world->add(std::move(rx));

	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr, ".");

	// This will call calculateResponse, which returns nullptr because they share a platform,
	// and then pass it to routeResponse. It should return early without crashing.
	REQUIRE_NOTHROW(engine.handleTxPulsedStart(tx_ptr, 0.0));

	// Ensure nothing was added to the inbox
	REQUIRE(rx_ptr->drainInbox().empty());
}

TEST_CASE("SimulationEngine updateProgress safely handles null reporter", "[core][threading]")
{
	// This test covers:
	// 5. updateProgress() -> if (!_reporter) { return; }
	ParamGuard guard;
	params::setTime(0.0, 1.0);
	auto world = createPhysicsWorld();
	pool::ThreadPool pool(1);

	// Add an event to ensure the main loop runs at least once, triggering updateProgress()
	auto* tx = world->getTransmitters().front().get();
	world->getEventQueue().push({0.5, core::EventType::TX_STREAMING_START, tx});

	// Passing nullptr for the progress callback ensures the internal _reporter is null
	REQUIRE_NOTHROW(core::runEventDrivenSim(world.get(), pool, nullptr, "."));
}

// TODO: A deeply synchronized test for the `shutdown()` thread joining logic is currently omitted.
// The `shutdown()` method enqueues jobs to `std::jthread`s which join automatically on destruction.
// Testing the precise thread joining deterministic behavior requires complex synchronization hooks
// inside the finalizer loops which breaks the "no-mocking" rule for this suite. The side-effects
// (enqueuing the shutdown job) are implicitly tested in the `runEventDrivenSim` integration test above.
