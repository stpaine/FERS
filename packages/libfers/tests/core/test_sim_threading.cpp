#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <complex>
#include <memory>
#include <string>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/config.h"
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
#include "timing/prototype_timing.h"
#include "timing/timing.h"

using Catch::Matchers::WithinAbs;

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

TEST_CASE("SimulationEngine handles CW state events", "[core][threading]")
{
	ParamGuard guard;
	auto world = createPhysicsWorld();
	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr);

	auto* tx = world->getTransmitters().front().get();
	auto* rx = world->getReceivers().front().get();

	SECTION("Tx CW Start adds to active list")
	{
		engine.handleTxCwStart(tx);
		REQUIRE(world->getSimulationState().active_cw_transmitters.size() == 1);
		REQUIRE(world->getSimulationState().active_cw_transmitters[0] == tx);
	}

	SECTION("Tx CW End removes from active list")
	{
		engine.handleTxCwStart(tx);
		engine.handleTxCwEnd(tx);
		REQUIRE(world->getSimulationState().active_cw_transmitters.empty());
	}

	SECTION("Rx CW Start sets active")
	{
		REQUIRE_FALSE(rx->isActive());
		engine.handleRxCwStart(rx);
		REQUIRE(rx->isActive());
	}

	SECTION("Rx CW End clears active")
	{
		engine.handleRxCwStart(rx);
		engine.handleRxCwEnd(rx);
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
	core::SimulationEngine engine(world.get(), pool, nullptr);

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
	core::SimulationEngine engine(world.get(), pool, nullptr);

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
	core::SimulationEngine engine(world.get(), pool, nullptr);

	auto* tx = world->getTransmitters().front().get();
	auto* rx = world->getReceivers().front().get();
	std::vector<radar::Transmitter*> active_tx = {tx};

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

		ComplexType actual_total = engine.calculateCwSample(rx, 0.0, active_tx);

		REQUIRE_THAT(actual_total.real(), WithinAbs(expected_total.real(), 1e-12));
		REQUIRE_THAT(actual_total.imag(), WithinAbs(expected_total.imag(), 1e-12));
	}
}

TEST_CASE("SimulationEngine processCwPhysics steps through time and updates buffers", "[core][threading]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);
	params::setTime(0.0, 10.0);

	auto world = createPhysicsWorld();
	pool::ThreadPool pool(1);
	core::SimulationEngine engine(world.get(), pool, nullptr);

	auto* tx = world->getTransmitters().front().get();
	auto* rx = world->getReceivers().front().get();

	// Prepare receiver buffer
	rx->prepareCwData(10);
	rx->setActive(true);

	// Set state
	world->getSimulationState().t_current = 0.0;
	world->getSimulationState().active_cw_transmitters.push_back(tx);

	SECTION("Processes correct number of samples based on dt")
	{
		// dt = 1/1000 = 0.001s. Processing up to t=0.0025 should hit indices 0, 1, 2.
		engine.processCwPhysics(0.0025);

		const auto& buffer = rx->getCwData();
		// Samples 0, 1, 2 should be non-zero (populated with physics)
		REQUIRE(std::abs(buffer[0]) > 0.0);
		REQUIRE(std::abs(buffer[1]) > 0.0);
		REQUIRE(std::abs(buffer[2]) > 0.0);

		// Sample 3 should be untouched (0.0)
		REQUIRE(std::abs(buffer[3]) == 0.0);
	}
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
	world->getEventQueue().push({0.001, core::EventType::TX_CW_START, tx});
	world->getEventQueue().push({0.004, core::EventType::TX_CW_END, tx});

	int progress_updates = 0;
	auto progress_cb = [&](const std::string&, int, int) { progress_updates++; };

	SECTION("Simulation runs to completion and shuts down cleanly")
	{
		REQUIRE_NOTHROW(core::runEventDrivenSim(world.get(), pool, progress_cb));

		// Should have received at least initialization and completion updates
		REQUIRE(progress_updates >= 2);

		// Simulation time should have advanced to the last event
		REQUIRE_THAT(world->getSimulationState().t_current, WithinAbs(0.004, 1e-9));
	}
}

// TODO: A deeply synchronized test for the `shutdown()` thread joining logic is currently omitted.
// The `shutdown()` method enqueues jobs to `std::jthread`s which join automatically on destruction.
// Testing the precise thread joining deterministic behavior requires complex synchronization hooks
// inside the finalizer loops which breaks the "no-mocking" rule for this suite. The side-effects
// (enqueuing the shutdown job) are implicitly tested in the `runEventDrivenSim` integration test above.
