#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <numbers>

#include "core/parameters.h"
#include "interpolation/interpolation_filter.h"

using Catch::Matchers::WithinAbs;

TEST_CASE("InterpFilter sinc definition", "[interpolation][filter]")
{
	REQUIRE_THAT(interp::InterpFilter::sinc(0.0), WithinAbs(1.0, 0.0));
	const RealType pi_value = std::numbers::pi_v<RealType>;
	const RealType expected = std::sin(0.5 * pi_value) / (0.5 * pi_value);
	REQUIRE_THAT(interp::InterpFilter::sinc(0.5), WithinAbs(expected, 1e-12));
}

TEST_CASE("InterpFilter Kaiser window symmetry", "[interpolation][filter]")
{
	const auto length = params::renderFilterLength();
	const RealType alpha = std::floor(length / 2.0);
	interp::InterpFilter& filter = interp::InterpFilter::getInstance();

	auto left = filter.kaiserWinCompute(alpha - 0.25);
	auto right = filter.kaiserWinCompute(alpha + 0.25);
	REQUIRE(left.has_value());
	REQUIRE(right.has_value());
	REQUIRE_THAT(*left, WithinAbs(*right, 1e-12));
}

TEST_CASE("InterpFilter Kaiser window normalization", "[interpolation][filter]")
{
	const auto length = params::renderFilterLength();
	const RealType alpha = std::floor(length / 2.0);
	interp::InterpFilter& filter = interp::InterpFilter::getInstance();

	auto center = filter.kaiserWinCompute(alpha);
	REQUIRE(center.has_value());
	REQUIRE_THAT(*center, WithinAbs(1.0, 1e-9));
}

TEST_CASE("InterpFilter Kaiser window clamps outside support", "[interpolation][filter]")
{
	const auto length = params::renderFilterLength();
	const RealType alpha = std::floor(length / 2.0);
	interp::InterpFilter& filter = interp::InterpFilter::getInstance();

	auto below = filter.kaiserWinCompute(-0.1);
	REQUIRE(below.has_value());
	REQUIRE_THAT(*below, WithinAbs(0.0, 0.0));

	auto above = filter.kaiserWinCompute(2.0 * alpha + 0.1);
	REQUIRE(above.has_value());
	REQUIRE_THAT(*above, WithinAbs(0.0, 0.0));
}

TEST_CASE("InterpFilter windowed sinc at zero", "[interpolation][filter]")
{
	interp::InterpFilter& filter = interp::InterpFilter::getInstance();
	auto value = filter.interpFilter(0.0);
	REQUIRE(value.has_value());
	REQUIRE_THAT(*value, WithinAbs(1.0, 1e-9));
}

TEST_CASE("InterpFilter windowed sinc outside support", "[interpolation][filter]")
{
	const auto length = params::renderFilterLength();
	const RealType alpha = std::floor(length / 2.0);
	interp::InterpFilter& filter = interp::InterpFilter::getInstance();

	auto value = filter.interpFilter(-alpha - 0.5);
	REQUIRE(value.has_value());
	REQUIRE_THAT(*value, WithinAbs(0.0, 0.0));
}

TEST_CASE("InterpFilter getFilter returns length and center", "[interpolation][filter]")
{
	const auto length = params::renderFilterLength();
	const RealType alpha = std::floor(length / 2.0);
	interp::InterpFilter& filter = interp::InterpFilter::getInstance();

	auto span = filter.getFilter(0.0);
	REQUIRE(span.size() == length);
	REQUIRE_THAT(span[static_cast<size_t>(alpha)], WithinAbs(1.0, 1e-6));
}

TEST_CASE("InterpFilter getFilter rejects out of range", "[interpolation][filter]")
{
	interp::InterpFilter& filter = interp::InterpFilter::getInstance();
	REQUIRE_THROWS_AS(filter.getFilter(1.1), std::runtime_error);
	REQUIRE_THROWS_AS(filter.getFilter(-1.1), std::runtime_error);
}

TEST_CASE("InterpFilter getFilter keeps the upper boundary in range", "[interpolation][filter]")
{
	interp::InterpFilter& filter = interp::InterpFilter::getInstance();

	const auto last_filter = filter.getFilter(1.0);
	const auto near_last_filter = filter.getFilter(0.999);

	REQUIRE(last_filter.size() == params::renderFilterLength());
	REQUIRE(near_last_filter.size() == last_filter.size());
	REQUIRE(last_filter.data() >= near_last_filter.data());
}
