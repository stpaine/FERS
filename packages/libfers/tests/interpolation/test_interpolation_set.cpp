#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/config.h"
#include "interpolation/interpolation_point.h"
#include "interpolation/interpolation_set.h"

using Catch::Matchers::WithinAbs;

TEST_CASE("InterpPoint defaults to zeros", "[interpolation][point]")
{
	interp::InterpPoint point;
	REQUIRE_THAT(point.power, WithinAbs(0.0, 0.0));
	REQUIRE_THAT(point.time, WithinAbs(0.0, 0.0));
	REQUIRE_THAT(point.delay, WithinAbs(0.0, 0.0));
	REQUIRE_THAT(point.phase, WithinAbs(0.0, 0.0));
}

TEST_CASE("InterpSetData returns nullopt when empty", "[interpolation][set]")
{
	interp::InterpSetData data;
	REQUIRE_FALSE(data.value(0.0).has_value());
	REQUIRE_THAT(data.max(), WithinAbs(0.0, 0.0));
}

TEST_CASE("InterpSetData clamps before first sample", "[interpolation][set]")
{
	interp::InterpSetData data;
	data.insertSample(5.0, 10.0);
	data.insertSample(6.0, 20.0);

	auto value = data.value(1.0);
	REQUIRE(value.has_value());
	REQUIRE_THAT(*value, WithinAbs(10.0, 1e-12));
}

TEST_CASE("InterpSetData clamps after last sample", "[interpolation][set]")
{
	interp::InterpSetData data;
	data.insertSample(1.0, -5.0);
	data.insertSample(3.0, -15.0);

	auto value = data.value(10.0);
	REQUIRE(value.has_value());
	REQUIRE_THAT(*value, WithinAbs(-15.0, 1e-12));
}

TEST_CASE("InterpSetData returns exact sample", "[interpolation][set]")
{
	interp::InterpSetData data;
	data.insertSample(2.0, 4.0);
	data.insertSample(4.0, 8.0);

	auto value = data.value(2.0);
	REQUIRE(value.has_value());
	REQUIRE_THAT(*value, WithinAbs(4.0, 1e-12));
}

TEST_CASE("InterpSetData linearly interpolates", "[interpolation][set]")
{
	interp::InterpSetData data;
	data.insertSample(0.0, 0.0);
	data.insertSample(10.0, 20.0);

	auto mid = data.value(2.5);
	REQUIRE(mid.has_value());
	REQUIRE_THAT(*mid, WithinAbs(5.0, 1e-12));
}

TEST_CASE("InterpSetData max uses absolute magnitude", "[interpolation][set]")
{
	interp::InterpSetData data;
	data.insertSample(0.0, -7.0);
	data.insertSample(1.0, 5.0);
	REQUIRE_THAT(data.max(), WithinAbs(7.0, 1e-12));
}

TEST_CASE("InterpSetData divide scales values", "[interpolation][set]")
{
	interp::InterpSetData data;
	data.insertSample(0.0, -4.0);
	data.insertSample(2.0, 6.0);

	data.divide(2.0);
	REQUIRE_THAT(data.max(), WithinAbs(3.0, 1e-12));

	auto value = data.value(2.0);
	REQUIRE(value.has_value());
	REQUIRE_THAT(*value, WithinAbs(3.0, 1e-12));
}

TEST_CASE("InterpSetData divide rejects zero", "[interpolation][set]")
{
	interp::InterpSetData data;
	data.insertSample(0.0, 1.0);
	REQUIRE_THROWS_AS(data.divide(0.0), std::invalid_argument);
}

TEST_CASE("InterpSet wrapper supports float interpolation", "[interpolation][set]")
{
	interp::InterpSet set;
	set.insertSample(0.0f, 0.0f);
	set.insertSample(1.0f, 1.0f);

	auto value = set.getValueAt(0.5f);
	REQUIRE(value.has_value());
	REQUIRE_THAT(*value, WithinAbs(0.5f, 1e-6f));

	set.divide(2.0f);
	auto scaled = set.getValueAt(0.5f);
	REQUIRE(scaled.has_value());
	REQUIRE_THAT(*scaled, WithinAbs(0.25f, 1e-6f));
}

TEST_CASE("InterpSet wrapper exposes max", "[interpolation][set]")
{
	interp::InterpSet set;
	set.insertSample(0.0, -2.0);
	set.insertSample(1.0, 1.0);
	REQUIRE_THAT(set.getMax(), WithinAbs(2.0, 1e-12));
}
