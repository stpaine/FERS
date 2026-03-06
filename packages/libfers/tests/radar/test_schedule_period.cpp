#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <string>
#include <vector>

#include "core/parameters.h"
#include "radar/schedule_period.h"

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

TEST_CASE("SchedulePeriod filters invalid and out-of-bounds periods", "[radar][schedule_period]")
{
	ParamGuard guard;
	params::setTime(0.0, 10.0);

	std::vector<radar::SchedulePeriod> raw = {{-1.0, -0.5}, {12.0, 15.0}, {5.0, 5.0},
											  {7.0, 3.0},	{1.0, 2.0},	  {2.0, 3.0}};

	const auto processed = radar::processRawSchedule(std::move(raw), "TestRadar", false, 1.0);

	REQUIRE(processed.size() == 1);
	REQUIRE_THAT(processed.front().start, WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(processed.front().end, WithinAbs(3.0, 1e-12));
}

TEST_CASE("SchedulePeriod sorts and merges overlaps and adjacency", "[radar][schedule_period]")
{
	ParamGuard guard;
	params::setTime(0.0, 100.0);

	std::vector<radar::SchedulePeriod> raw = {{5.0, 6.0},	{1.0, 4.0},	  {4.0, 5.0},
											  {10.0, 12.0}, {11.0, 20.0}, {30.0, 31.0}};

	const auto processed = radar::processRawSchedule(std::move(raw), "MergeRadar", false, 1.0);

	REQUIRE(processed.size() == 3);
	REQUIRE_THAT(processed[0].start, WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(processed[0].end, WithinAbs(6.0, 1e-12));
	REQUIRE_THAT(processed[1].start, WithinAbs(10.0, 1e-12));
	REQUIRE_THAT(processed[1].end, WithinAbs(20.0, 1e-12));
	REQUIRE_THAT(processed[2].start, WithinAbs(30.0, 1e-12));
	REQUIRE_THAT(processed[2].end, WithinAbs(31.0, 1e-12));
}

TEST_CASE("SchedulePeriod keeps periods that intersect simulation bounds", "[radar][schedule_period]")
{
	ParamGuard guard;
	params::setTime(5.0, 15.0);

	std::vector<radar::SchedulePeriod> raw = {{0.0, 6.0}, {5.0, 6.0}, {14.0, 16.0}, {20.0, 25.0}};

	const auto processed = radar::processRawSchedule(std::move(raw), "BoundsRadar", false, 1.0);

	REQUIRE(processed.size() == 2);
	REQUIRE_THAT(processed[0].start, WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(processed[0].end, WithinAbs(6.0, 1e-12));
	REQUIRE_THAT(processed[1].start, WithinAbs(14.0, 1e-12));
	REQUIRE_THAT(processed[1].end, WithinAbs(16.0, 1e-12));
}

TEST_CASE("SchedulePeriod honors PRI checks without altering periods", "[radar][schedule_period]")
{
	ParamGuard guard;
	params::setTime(0.0, 10.0);

	std::vector<radar::SchedulePeriod> raw = {{1.0, 1.2}, {2.0, 3.0}};

	const auto processed = radar::processRawSchedule(raw, "PulsedRadar", true, 0.5);

	REQUIRE(processed.size() == 2);
	REQUIRE_THAT(processed[0].start, WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(processed[0].end, WithinAbs(1.2, 1e-12));
	REQUIRE_THAT(processed[1].start, WithinAbs(2.0, 1e-12));
	REQUIRE_THAT(processed[1].end, WithinAbs(3.0, 1e-12));
}

TEST_CASE("SchedulePeriod handles empty input", "[radar][schedule_period]")
{
	ParamGuard guard;
	params::setTime(0.0, 10.0);

	const auto processed = radar::processRawSchedule({}, "EmptyRadar", false, 1.0);
	REQUIRE(processed.empty());
}

TEST_CASE("SchedulePeriod returns empty after filtering all invalid periods", "[radar][schedule_period]")
{
	ParamGuard guard;
	params::setTime(0.0, 10.0);

	std::vector<radar::SchedulePeriod> raw = {{5.0, 5.0}, {7.0, 3.0}, {-5.0, -1.0}, {12.0, 20.0}};

	const auto processed = radar::processRawSchedule(std::move(raw), "FilteredRadar", false, 1.0);
	REQUIRE(processed.empty());
}
