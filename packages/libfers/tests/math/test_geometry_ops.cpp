#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "math/geometry_ops.h"

using namespace math;
using Catch::Matchers::WithinAbs;

constexpr RealType PI_VAL = std::numbers::pi_v<RealType>;

TEST_CASE("Vec3 Construction and Basic Access", "[math][geometry][vec3]")
{
	Vec3 v_default;
	REQUIRE(v_default.x == 0.0);
	REQUIRE(v_default.y == 0.0);
	REQUIRE(v_default.z == 0.0);

	Vec3 v_param(1.0, -2.5, 3.14);
	REQUIRE(v_param.x == 1.0);
	REQUIRE(v_param.y == -2.5);
	REQUIRE(v_param.z == 3.14);
}

TEST_CASE("Vec3 Arithmetic Operators", "[math][geometry][vec3]")
{
	Vec3 a(1.0, 2.0, 3.0);
	Vec3 b(4.0, 5.0, 6.0);

	SECTION("Addition")
	{
		Vec3 c = a + b;
		REQUIRE_THAT(c.x, WithinAbs(5.0, 1e-9));
		a += b;
		REQUIRE_THAT(a.x, WithinAbs(5.0, 1e-9));
	}

	SECTION("Subtraction")
	{
		Vec3 c = b - a;
		REQUIRE_THAT(c.x, WithinAbs(3.0, 1e-9));
		b -= a;
		REQUIRE_THAT(b.x, WithinAbs(3.0, 1e-9));
	}

	SECTION("Unary Negation")
	{
		Vec3 neg = -a;
		REQUIRE_THAT(neg.x, WithinAbs(-1.0, 1e-9));
	}

	SECTION("Multiplication & Division")
	{
		Vec3 c = a * b;
		REQUIRE_THAT(c.x, WithinAbs(4.0, 1e-9));
		a *= b;
		REQUIRE_THAT(a.x, WithinAbs(4.0, 1e-9));

		Vec3 d = b / a;
		REQUIRE_THAT(d.x, WithinAbs(1.0, 1e-9));
	}
}

TEST_CASE("Vec3 Scalar Operations", "[math][geometry][vec3]")
{
	Vec3 v(2.0, 4.0, 8.0);
	RealType s = 2.0;

	SECTION("Multiplication")
	{
		Vec3 res1 = v * s;
		REQUIRE_THAT(res1.x, WithinAbs(4.0, 1e-9));
		v *= s;
		REQUIRE_THAT(v.y, WithinAbs(8.0, 1e-9));
	}

	SECTION("Division")
	{
		Vec3 res1 = v / s;
		REQUIRE_THAT(res1.x, WithinAbs(1.0, 1e-9));
		Vec3 res2 = s / v;
		REQUIRE_THAT(res2.x, WithinAbs(1.0, 1e-9));
		v /= s;
		REQUIRE_THAT(v.z, WithinAbs(4.0, 1e-9));
	}
}

TEST_CASE("Vec3 Advanced Math", "[math][geometry][vec3]")
{
	Vec3 v(3.0, 4.0, 0.0);

	SECTION("Length") { REQUIRE_THAT(v.length(), WithinAbs(5.0, 1e-9)); }

	SECTION("Dot Product")
	{
		Vec3 a(1, 0, 0);
		Vec3 b(0, 1, 0);
		Vec3 c(2, 2, 0);
		REQUIRE_THAT(dotProduct(a, b), WithinAbs(0.0, 1e-9)); // Orthogonal
		REQUIRE_THAT(dotProduct(a, c), WithinAbs(2.0, 1e-9));
	}
}

TEST_CASE("Matrix3 Operations", "[math][geometry][matrix]")
{
	Matrix3 m;
	m.elements = {2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
	Vec3 v(1.0, 2.0, 3.0);
	v *= m;
	REQUIRE_THAT(v.x, WithinAbs(2.0, 1e-9));
}

TEST_CASE("SVec3 (Spherical) Operations", "[math][geometry][svec3]")
{
	SECTION("Conversion from Vec3")
	{
		Vec3 vx(10.0, 0.0, 0.0);
		SVec3 svx(vx);
		REQUIRE_THAT(svx.length, WithinAbs(10.0, 1e-9));
		REQUIRE_THAT(svx.azimuth, WithinAbs(0.0, 1e-9));

		Vec3 vz(0.0, 0.0, 5.0);
		SVec3 svz(vz);
		REQUIRE_THAT(svz.elevation, WithinAbs(PI_VAL / 2.0, 1e-9));

		Vec3 vz_back(svz);
		REQUIRE_THAT(vz_back.z, WithinAbs(5.0, 1e-9));
	}

	SECTION("Scalar Ops")
	{
		SVec3 s(10.0, 1.0, 1.0);
		s *= 2.0;
		REQUIRE_THAT(s.length, WithinAbs(20.0, 1e-9));
		s /= 4.0;
		REQUIRE_THAT(s.length, WithinAbs(5.0, 1e-9));
	}

	SECTION("Addition with Wrapping")
	{
		// Normal addition
		SVec3 s1(1.0, 0.5, 0.0);
		SVec3 s2(1.0, 0.5, 0.0);
		SVec3 sum = s1 + s2;
		REQUIRE_THAT(sum.azimuth, WithinAbs(1.0, 1e-9));

		// fmod(-1.5 * PI, 2 * PI) = -1.5 * PI in C++. We must wrap this to +0.5 * PI.
		SVec3 s3(1.0, -PI_VAL, 0.0);
		SVec3 s4(1.0, -PI_VAL / 2.0, 0.0);
		SVec3 sum_neg = s3 + s4;
		REQUIRE_THAT(sum_neg.azimuth, WithinAbs(PI_VAL / 2.0, 1e-9)); // Physically correct wrapping
	}

	SECTION("Subtraction with Shortest Path Wrapping")
	{
		// Normal subtraction
		SVec3 s1(1.0, 1.0, 0.0);
		SVec3 s2(1.0, 0.5, 0.0);
		SVec3 diff = s1 - s2;
		REQUIRE_THAT(diff.azimuth, WithinAbs(0.5, 1e-9));

		// 10 degrees - 350 degrees = -340 degrees -> should wrap to +20 degrees
		SVec3 s3(1.0, 0.1, 0.0);
		SVec3 s4(1.0, 2 * PI_VAL - 0.1, 0.0);
		SVec3 diff_neg_wrap = s3 - s4;
		REQUIRE_THAT(diff_neg_wrap.azimuth, WithinAbs(0.2, 1e-9));

		// 350 degrees - 10 degrees = 340 degrees -> should wrap to -20 degrees
		SVec3 s5(1.0, 2 * PI_VAL - 0.1, 0.0);
		SVec3 s6(1.0, 0.1, 0.0);
		SVec3 diff_pos_wrap = s5 - s6;
		REQUIRE_THAT(diff_pos_wrap.azimuth, WithinAbs(-0.2, 1e-9)); // Physically shortest path
	}
}
