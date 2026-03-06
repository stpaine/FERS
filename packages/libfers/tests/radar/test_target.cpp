#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <memory>

#include "radar/platform.h"
#include "radar/target.h"

using Catch::Matchers::WithinAbs;

namespace
{
	class FixedRcsModel final : public radar::RcsModel
	{
	public:
		explicit FixedRcsModel(RealType value) : _value(value) {}

		RealType sampleModel() override { return _value; }

	private:
		RealType _value;
	};
}

TEST_CASE("RcsConst returns unity", "[radar][target]")
{
	radar::RcsConst model;
	REQUIRE_THAT(model.sampleModel(), WithinAbs(1.0, 1e-12));
}

TEST_CASE("RcsChiSquare reports configured degrees of freedom", "[radar][target]")
{
	std::mt19937 rng(12345);
	radar::RcsChiSquare model(rng, 2.0);
	REQUIRE_THAT(model.getK(), WithinAbs(2.0, 1e-12));

	const RealType sample = model.sampleModel();
	REQUIRE(sample > 0.0);
}

TEST_CASE("IsoTarget returns constant RCS", "[radar][target]")
{
	radar::Platform platform("TargetPlatform");
	radar::IsoTarget target(&platform, "Iso", 12.5, 42, 9002);

	math::SVec3 in_angle(1.0, 0.0, 0.0);
	math::SVec3 out_angle(1.0, 0.0, 0.0);

	REQUIRE_THAT(target.getConstRcs(), WithinAbs(12.5, 1e-12));
	REQUIRE_THAT(target.getRcs(in_angle, out_angle, 0.0), WithinAbs(12.5, 1e-12));
	REQUIRE(target.getId() == 9002);
}

TEST_CASE("createIsoTarget constructs IsoTarget", "[radar][target]")
{
	radar::Platform platform("TargetPlatform");
	const auto target = radar::createIsoTarget(&platform, "IsoFactory", 9.0, 77, 8080);

	const auto* iso_target = dynamic_cast<const radar::IsoTarget*>(target.get());
	REQUIRE(iso_target != nullptr);
	REQUIRE_THAT(iso_target->getConstRcs(), WithinAbs(9.0, 1e-12));
	REQUIRE(iso_target->getId() == 8080);
}

TEST_CASE("IsoTarget applies fluctuation model", "[radar][target]")
{
	radar::Platform platform("TargetPlatform");
	radar::IsoTarget target(&platform, "Iso", 3.0, 7);

	auto model = std::make_unique<FixedRcsModel>(2.0);
	const auto* model_ptr = model.get();
	target.setFluctuationModel(std::move(model));

	math::SVec3 in_angle(1.0, 0.0, 0.0);
	math::SVec3 out_angle(1.0, 0.0, 0.0);

	REQUIRE(target.getFluctuationModel() == model_ptr);
	REQUIRE_THAT(target.getRcs(in_angle, out_angle, 0.0), WithinAbs(6.0, 1e-12));
}

TEST_CASE("Target RNG is deterministic per seed", "[radar][target]")
{
	radar::Platform platform("TargetPlatform");
	radar::IsoTarget target_a(&platform, "IsoA", 1.0, 1337);
	radar::IsoTarget target_b(&platform, "IsoB", 1.0, 1337);

	const auto first_a = target_a.getRngEngine()();
	const auto first_b = target_b.getRngEngine()();

	REQUIRE(first_a == first_b);
}
