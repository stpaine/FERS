#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <limits>
#include <type_traits>

#include "core/config.h"

using Catch::Matchers::WithinAbs;

TEST_CASE("Config defines RealType and ComplexType", "[core][config]")
{
	REQUIRE(std::is_same_v<RealType, double>);
	REQUIRE(std::is_same_v<ComplexType, std::complex<RealType>>);
}

TEST_CASE("Config defines PI and EPSILON constants", "[core][config]")
{
	REQUIRE_THAT(PI, WithinAbs(std::numbers::pi_v<RealType>, 1e-15));
	REQUIRE_THAT(EPSILON, WithinAbs(std::numeric_limits<RealType>::epsilon(), 0.0));
	REQUIRE(EPSILON > 0.0);
}
