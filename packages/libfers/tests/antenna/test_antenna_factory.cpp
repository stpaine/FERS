#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/portable_utils.h"
#include "core/sim_id.h"
#include "math/geometry_ops.h"

using Catch::Matchers::WithinAbs;

namespace
{
	math::SVec3 unitDirection(const RealType azimuth, const RealType elevation) { return {1.0, azimuth, elevation}; }
}

TEST_CASE("Isotropic antenna gain equals efficiency factor", "[antenna]")
{
	antenna::Isotropic antenna("iso", 42);
	antenna.setEfficiencyFactor(0.7);

	const auto gain = antenna.getGain(unitDirection(0.2, -0.1), unitDirection(0.0, 0.0), 1.0);
	REQUIRE_THAT(gain, WithinAbs(0.7, 1e-12));
	REQUIRE(antenna.getId() == 42);
	REQUIRE(antenna.getName() == "iso");
}

TEST_CASE("Antenna ids default to antenna type", "[antenna]")
{
	antenna::Isotropic antenna("iso-default");
	REQUIRE(SimIdGenerator::getType(antenna.getId()) == ObjectType::Antenna);
}

TEST_CASE("Antenna base noise temperature defaults to zero", "[antenna]")
{
	antenna::Isotropic antenna("iso");
	REQUIRE_THAT(antenna.getNoiseTemperature(unitDirection(0.3, -0.2)), WithinAbs(0.0, 0.0));
}

TEST_CASE("Sinc antenna peak and null match physics", "[antenna]")
{
	const RealType alpha = 2.0;
	const RealType beta = 2.0;
	const RealType gamma = 2.0;
	antenna::Sinc antenna("sinc", alpha, beta, gamma);
	antenna.setEfficiencyFactor(0.5);
	REQUIRE(antenna.getAlpha() == alpha);
	REQUIRE(antenna.getBeta() == beta);
	REQUIRE(antenna.getGamma() == gamma);

	const auto ref = unitDirection(0.0, 0.0);
	const auto on_axis = unitDirection(0.0, 0.0);
	const auto gain_on_axis = antenna.getGain(on_axis, ref, 1.0);
	REQUIRE_THAT(gain_on_axis, WithinAbs(alpha * 0.5, 1e-12));

	const auto null_angle = unitDirection(PI / 2.0, 0.0);
	const auto gain_null = antenna.getGain(null_angle, ref, 1.0);
	REQUIRE_THAT(gain_null, WithinAbs(0.0, 1e-12));
}

TEST_CASE("Gaussian antenna gain follows expected Gaussian surface", "[antenna]")
{
	const RealType azscale = 0.5;
	const RealType elscale = 1.5;
	antenna::Gaussian antenna("gauss", azscale, elscale);

	const auto ref = unitDirection(0.0, 0.0);
	const auto angle = unitDirection(0.2, 0.1);
	const RealType expected = std::exp(-0.2 * 0.2 * azscale) * std::exp(-0.1 * 0.1 * elscale);

	const auto gain = antenna.getGain(angle, ref, 1.0);
	REQUIRE_THAT(gain, WithinAbs(expected, 1e-12));
	REQUIRE(antenna.getAzimuthScale() == azscale);
	REQUIRE(antenna.getElevationScale() == elscale);
}

TEST_CASE("Sinc antenna gain is symmetric about boresight", "[antenna]")
{
	antenna::Sinc antenna("sinc", 3.0, 1.7, 2.0);
	antenna.setEfficiencyFactor(0.6);

	const auto ref = unitDirection(0.0, 0.0);
	const auto left = unitDirection(-0.2, 0.0);
	const auto right = unitDirection(0.2, 0.0);

	REQUIRE_THAT(antenna.getGain(left, ref, 1.0), WithinAbs(antenna.getGain(right, ref, 1.0), 1e-12));
}

TEST_CASE("Gaussian antenna gain is symmetric about boresight", "[antenna]")
{
	antenna::Gaussian antenna("gauss", 0.75, 1.25);

	const auto ref = unitDirection(0.0, 0.0);
	const auto positive = unitDirection(0.2, -0.15);
	const auto negative = unitDirection(-0.2, 0.15);

	REQUIRE_THAT(antenna.getGain(positive, ref, 1.0), WithinAbs(antenna.getGain(negative, ref, 1.0), 1e-12));
}

TEST_CASE("Gaussian antenna currently ignores efficiency factor", "[antenna]")
{
	const RealType azscale = 0.5;
	const RealType elscale = 1.5;
	antenna::Gaussian antenna("gauss", azscale, elscale);
	antenna.setEfficiencyFactor(0.2);

	const auto ref = unitDirection(0.0, 0.0);
	const auto angle = unitDirection(0.2, 0.1);
	const RealType expected = std::exp(-0.2 * 0.2 * azscale) * std::exp(-0.1 * 0.1 * elscale);

	REQUIRE_THAT(antenna.getGain(angle, ref, 1.0), WithinAbs(expected, 1e-12));
}

TEST_CASE("Square horn gain matches aperture physics", "[antenna]")
{
	const RealType dimension = 0.4;
	const RealType wavelength = 0.1;
	antenna::SquareHorn antenna("horn", dimension);
	antenna.setEfficiencyFactor(0.9);
	REQUIRE(antenna.getDimension() == dimension);

	const auto ref = unitDirection(0.0, 0.0);
	const auto on_axis = unitDirection(0.0, 0.0);
	const RealType ge = 4.0 * PI * dimension * dimension / (wavelength * wavelength);
	const RealType expected_on_axis = ge * 0.9;

	const auto gain_on_axis = antenna.getGain(on_axis, ref, wavelength);
	REQUIRE_THAT(gain_on_axis, WithinAbs(expected_on_axis, 1e-10));

	const auto off_axis = unitDirection(0.15, 0.0);
	const RealType x = PI * dimension * std::sin(0.15) / wavelength;
	const RealType expected_off_axis = ge * std::pow(std::sin(x) / x, 2) * 0.9;
	const auto gain_off_axis = antenna.getGain(off_axis, ref, wavelength);
	REQUIRE_THAT(gain_off_axis, WithinAbs(expected_off_axis, 1e-10));
	REQUIRE(gain_off_axis < gain_on_axis);
}

TEST_CASE("Parabolic gain matches Bessel-based model", "[antenna]")
{
	const RealType diameter = 1.2;
	const RealType wavelength = 0.3;
	antenna::Parabolic antenna("parabola", diameter);
	antenna.setEfficiencyFactor(0.8);

	const auto ref = unitDirection(0.0, 0.0);
	const auto on_axis = unitDirection(0.0, 0.0);
	const RealType ge = std::pow(PI * diameter / wavelength, 2);
	const RealType expected_on_axis = ge * 4.0 * 0.8;
	const auto gain_on_axis = antenna.getGain(on_axis, ref, wavelength);

	REQUIRE_THAT(gain_on_axis, WithinAbs(expected_on_axis, 1e-8));
	REQUIRE(antenna.getDiameter() == diameter);

	const auto off_axis = unitDirection(0.1, 0.0);
	const RealType x = PI * diameter * std::sin(0.1) / wavelength;
	const RealType expected_off_axis = ge * std::pow(2.0 * (core::besselJ1(x) / x), 2) * 0.8;
	const auto gain_off_axis = antenna.getGain(off_axis, ref, wavelength);
	REQUIRE_THAT(gain_off_axis, WithinAbs(expected_off_axis, 1e-8));
	REQUIRE(gain_off_axis < gain_on_axis);
}

TEST_CASE("Square horn efficiency scales gain", "[antenna]")
{
	const auto ref = unitDirection(0.0, 0.0);
	const auto angle = unitDirection(0.1, -0.05);
	constexpr RealType wavelength = 0.2;

	antenna::SquareHorn low_efficiency("horn-low", 0.4);
	low_efficiency.setEfficiencyFactor(0.25);

	antenna::SquareHorn high_efficiency("horn-high", 0.4);
	high_efficiency.setEfficiencyFactor(0.75);

	const RealType low_gain = low_efficiency.getGain(angle, ref, wavelength);
	const RealType high_gain = high_efficiency.getGain(angle, ref, wavelength);

	REQUIRE_THAT(high_gain, WithinAbs(low_gain * 3.0, 1e-10));
}

TEST_CASE("Parabolic efficiency scales gain", "[antenna]")
{
	const auto ref = unitDirection(0.0, 0.0);
	const auto angle = unitDirection(0.08, 0.0);
	constexpr RealType wavelength = 0.25;

	antenna::Parabolic low_efficiency("dish-low", 1.2);
	low_efficiency.setEfficiencyFactor(0.4);

	antenna::Parabolic high_efficiency("dish-high", 1.2);
	high_efficiency.setEfficiencyFactor(0.9);

	const RealType low_gain = low_efficiency.getGain(angle, ref, wavelength);
	const RealType high_gain = high_efficiency.getGain(angle, ref, wavelength);

	REQUIRE_THAT(high_gain, WithinAbs(low_gain * (0.9 / 0.4), 1e-8));
}

TEST_CASE("Efficiency factor accepts values greater than unity", "[antenna]")
{
	antenna::Isotropic antenna("iso");
	antenna.setEfficiencyFactor(1.25);
	REQUIRE_THAT(antenna.getEfficiencyFactor(), WithinAbs(1.25, 0.0));
}
