#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <memory>
#include <vector>

#include "core/parameters.h"
#include "radar/platform.h"
#include "radar/transmitter.h"
#include "signal/radar_signal.h"

using Catch::Matchers::WithinAbs;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
	};
}

TEST_CASE("Transmitter basic accessors and signal setters", "[radar][transmitter]")
{
	radar::Platform platform("TxPlatform");
	radar::Transmitter tx(&platform, "TxA", radar::OperationMode::PULSED_MODE, 5555);

	REQUIRE(tx.getMode() == radar::OperationMode::PULSED_MODE);
	REQUIRE(tx.getId() == 5555);
	REQUIRE(tx.getSignal() == nullptr);

	tx.setMode(radar::OperationMode::CW_MODE);
	REQUIRE(tx.getMode() == radar::OperationMode::CW_MODE);

	auto signal = std::make_unique<fers_signal::CwSignal>();
	fers_signal::RadarSignal radar_signal("Sig", 10.0, 1.0e9, 1.0, std::move(signal));

	tx.setSignal(&radar_signal);
	REQUIRE(tx.getSignal() == &radar_signal);

	tx.setWave(nullptr);
	REQUIRE(tx.getSignal() == nullptr);
}

TEST_CASE("Transmitter setPrf quantizes to sample rate", "[radar][transmitter]")
{
	ParamGuard guard;
	params::setRate(1000.0);
	params::setOversampleRatio(1);

	radar::Platform platform("TxPlatform");
	radar::Transmitter tx(&platform, "TxA", radar::OperationMode::CW_MODE);

	tx.setPrf(333.0);

	const RealType expected = 1.0 / (std::floor(1000.0 / 333.0) / 1000.0);
	REQUIRE_THAT(tx.getPrf(), WithinAbs(expected, 1e-12));
}

TEST_CASE("Transmitter schedule resolves next pulse time", "[radar][transmitter]")
{
	radar::Platform platform("TxPlatform");
	radar::Transmitter tx(&platform, "TxA", radar::OperationMode::PULSED_MODE);

	SECTION("No schedule means always on")
	{
		const auto next = tx.getNextPulseTime(2.5);
		REQUIRE(next.has_value());
		REQUIRE_THAT(*next, WithinAbs(2.5, 1e-12));
	}

	SECTION("Schedule enforces active windows")
	{
		std::vector<radar::SchedulePeriod> schedule = {{1.0, 2.0}, {4.0, 5.0}};
		tx.setSchedule(schedule);

		REQUIRE(tx.getSchedule().size() == 2);

		const auto inside = tx.getNextPulseTime(1.5);
		REQUIRE(inside.has_value());
		REQUIRE_THAT(*inside, WithinAbs(1.5, 1e-12));

		const auto before = tx.getNextPulseTime(3.0);
		REQUIRE(before.has_value());
		REQUIRE_THAT(*before, WithinAbs(4.0, 1e-12));

		const auto after = tx.getNextPulseTime(6.0);
		REQUIRE_FALSE(after.has_value());
	}
}
