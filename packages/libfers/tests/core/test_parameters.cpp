#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <string>

#include "core/logging.h"
#include "core/parameters.h"

using Catch::Matchers::ContainsSubstring;
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

TEST_CASE("Parameters default values are consistent", "[core][parameters]")
{
	ParamGuard guard;
	params::params.reset();

	REQUIRE_THAT(params::c(), WithinAbs(params::Parameters::DEFAULT_C, 0.0));
	REQUIRE_THAT(params::boltzmannK(), WithinAbs(params::Parameters::DEFAULT_BOLTZMANN_K, 0.0));
	REQUIRE_THAT(params::startTime(), WithinAbs(0.0, 0.0));
	REQUIRE_THAT(params::endTime(), WithinAbs(0.0, 0.0));
	REQUIRE_THAT(params::simSamplingRate(), WithinAbs(1000.0, 0.0));
	REQUIRE_THAT(params::rate(), WithinAbs(0.0, 0.0));
	REQUIRE(params::randomSeed() == 0u);
	REQUIRE(params::adcBits() == 0u);
	REQUIRE(params::renderFilterLength() == 33u);
	REQUIRE(params::renderThreads() == 1u);
	REQUIRE(params::oversampleRatio() == 1u);
	REQUIRE(params::coordinateFrame() == params::CoordinateFrame::ENU);
	REQUIRE(params::rotationAngleUnit() == params::RotationAngleUnit::Degrees);
	REQUIRE(params::utmZone() == 0);
	REQUIRE(params::utmNorthHemisphere());
}

TEST_CASE("Parameters setters update getters", "[core][parameters]")
{
	ParamGuard guard;
	params::params.reset();

	params::setC(300000000.0);
	REQUIRE_THAT(params::c(), WithinAbs(300000000.0, 1e-9));

	params::setTime(1.25, 9.5);
	REQUIRE_THAT(params::startTime(), WithinAbs(1.25, 1e-12));
	REQUIRE_THAT(params::endTime(), WithinAbs(9.5, 1e-12));

	params::setSimSamplingRate(2048.0);
	REQUIRE_THAT(params::simSamplingRate(), WithinAbs(2048.0, 1e-12));

	params::setRate(512.0);
	REQUIRE_THAT(params::rate(), WithinAbs(512.0, 1e-12));

	params::setRandomSeed(42);
	REQUIRE(params::randomSeed() == 42u);

	params::setAdcBits(12);
	REQUIRE(params::adcBits() == 12u);

	params::setOversampleRatio(4);
	REQUIRE(params::oversampleRatio() == 4u);

	params::setRotationAngleUnit(params::RotationAngleUnit::Radians);
	REQUIRE(params::rotationAngleUnit() == params::RotationAngleUnit::Radians);
}

TEST_CASE("Parameters setters validate inputs", "[core][parameters]")
{
	ParamGuard guard;
	params::params.reset();

	REQUIRE_THROWS_AS(params::setRate(0.0), std::runtime_error);
	REQUIRE_THROWS_AS(params::setRate(-1.0), std::runtime_error);

	REQUIRE_THROWS_AS(params::setOversampleRatio(0), std::runtime_error);
	REQUIRE_THROWS_WITH(params::setOversampleRatio(9), ContainsSubstring("Oversampling ratios > 8 are not supported"));
}

TEST_CASE("Parameters origin and coordinate settings", "[core][parameters]")
{
	ParamGuard guard;
	params::params.reset();

	params::setOrigin(1.5, 2.5, 3.5);
	REQUIRE_THAT(params::originLatitude(), WithinAbs(1.5, 1e-12));
	REQUIRE_THAT(params::originLongitude(), WithinAbs(2.5, 1e-12));
	REQUIRE_THAT(params::originAltitude(), WithinAbs(3.5, 1e-12));

	params::setCoordinateSystem(params::CoordinateFrame::UTM, 33, false);
	REQUIRE(params::coordinateFrame() == params::CoordinateFrame::UTM);
	REQUIRE(params::utmZone() == 33);
	REQUIRE_FALSE(params::utmNorthHemisphere());
}

TEST_CASE("Parameters setThreads returns expected", "[core][parameters]")
{
	ParamGuard guard;
	params::params.reset();

	const auto error = params::setThreads(0);
	REQUIRE_FALSE(error.has_value());

	const auto ok = params::setThreads(4);
	REQUIRE(ok.has_value());
	REQUIRE(params::renderThreads() == 4u);
}

TEST_CASE("Parameters reset restores defaults", "[core][parameters]")
{
	ParamGuard guard;
	params::setC(1.0);
	params::setTime(5.0, 10.0);
	params::setSimSamplingRate(500.0);
	params::setRate(100.0);
	params::setRandomSeed(7);
	params::setAdcBits(8);
	params::setOversampleRatio(2);
	params::setOrigin(9.0, 8.0, 7.0);
	params::setCoordinateSystem(params::CoordinateFrame::ECEF, 12, false);

	params::params.reset();

	REQUIRE_THAT(params::c(), WithinAbs(params::Parameters::DEFAULT_C, 0.0));
	REQUIRE_THAT(params::startTime(), WithinAbs(0.0, 0.0));
	REQUIRE_THAT(params::endTime(), WithinAbs(0.0, 0.0));
	REQUIRE_THAT(params::simSamplingRate(), WithinAbs(1000.0, 0.0));
	REQUIRE_THAT(params::rate(), WithinAbs(0.0, 0.0));
	REQUIRE(params::randomSeed() == 0u);
	REQUIRE(params::adcBits() == 0u);
	REQUIRE(params::oversampleRatio() == 1u);
	REQUIRE(params::coordinateFrame() == params::CoordinateFrame::ENU);
	REQUIRE(params::rotationAngleUnit() == params::RotationAngleUnit::Degrees);
	REQUIRE(params::utmZone() == 0);
	REQUIRE(params::utmNorthHemisphere());
}
