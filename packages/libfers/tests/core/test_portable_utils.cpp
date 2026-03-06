#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <limits>

#include "core/logging.h"
#include "core/portable_utils.h"

using Catch::Matchers::WithinRel;

TEST_CASE("besselJ1 matches standard library j1", "[core][portable]")
{
	const RealType x1 = 0.0;
	const RealType x2 = 1.0;
	const RealType x3 = -2.5;

	REQUIRE_THAT(core::besselJ1(x1), WithinRel(j1(x1), 1e-12));
	REQUIRE_THAT(core::besselJ1(x2), WithinRel(j1(x2), 1e-12));
	REQUIRE_THAT(core::besselJ1(x3), WithinRel(j1(x3), 1e-12));
}

TEST_CASE("countProcessors returns a valid count", "[core][portable]")
{
	const unsigned count = core::countProcessors();
	REQUIRE(count >= 1u);
}

TEST_CASE("countProcessors fallback logging not directly testable", "[core][portable]")
{
	// TODO: Cannot reliably force std::thread::hardware_concurrency() to return 0 to hit the error logging branch.
	REQUIRE(true);
}
