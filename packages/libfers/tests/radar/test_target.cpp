#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <system_error>

#include "core/sim_id.h"
#include "math/coord.h"
#include "radar/platform.h"
#include "radar/target.h"

using Catch::Matchers::WithinAbs;

namespace
{
	math::SVec3 unitDirection(const RealType azimuth, const RealType elevation) { return {1.0, azimuth, elevation}; }

	class FixedRcsModel final : public radar::RcsModel
	{
	public:
		explicit FixedRcsModel(RealType value) : _value(value) {}

		RealType sampleModel() override { return _value; }

	private:
		RealType _value;
	};

	std::string uniqueFileName(const std::string& prefix)
	{
		return prefix + "_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".xml";
	}

	std::filesystem::path tempFilePath(const std::string& filename)
	{
		return std::filesystem::temp_directory_path() / filename;
	}

	void removeIfExists(const std::filesystem::path& path)
	{
		std::error_code ec;
		std::filesystem::remove(path, ec);
	}

	void writeXmlTargetFile(const std::filesystem::path& path, const std::string& body)
	{
		std::ofstream out(path, std::ios::binary);
		REQUIRE(out.is_open());
		out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << body;
	}

	std::string xmlTargetFixture()
	{
		return R"(<target>
		<azimuth>
			<rcssample><angle>0.0</angle><rcs>2.0</rcs></rcssample>
			<rcssample><angle>1.0</angle><rcs>6.0</rcs></rcssample>
			<rcssample><angle>2.0</angle><rcs>10.0</rcs></rcssample>
		</azimuth>
		<elevation>
			<rcssample><angle>0.0</angle><rcs>3.0</rcs></rcssample>
			<rcssample><angle>1.0</angle><rcs>5.0</rcs></rcssample>
			<rcssample><angle>2.0</angle><rcs>7.0</rcs></rcssample>
		</elevation>
	</target>)";
	}

	void setStaticRotation(radar::Platform& platform, const RealType azimuth = 0.0, const RealType elevation = 0.0)
	{
		platform.getRotationPath()->addCoord(math::RotationCoord(azimuth, elevation, 0.0));
		platform.getRotationPath()->finalize();
	}

	void setConstantRotation(radar::Platform& platform, const RealType start_azimuth, const RealType start_elevation,
							 const RealType rate_azimuth, const RealType rate_elevation)
	{
		platform.getRotationPath()->setConstantRate(math::RotationCoord(start_azimuth, start_elevation, 0.0),
													math::RotationCoord(rate_azimuth, rate_elevation, 0.0));
	}
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

TEST_CASE("RcsChiSquare sample mean matches configured mean power", "[radar][target]")
{
	std::mt19937 rng(24680);
	radar::RcsChiSquare model(rng, 4.0);

	constexpr int sample_count = 50000;
	RealType sum = 0.0;
	RealType minimum = model.sampleModel();
	sum += minimum;

	for (int i = 1; i < sample_count; ++i)
	{
		const RealType sample = model.sampleModel();
		minimum = std::min(minimum, sample);
		sum += sample;
	}

	REQUIRE(minimum >= 0.0);
	REQUIRE_THAT(sum / static_cast<RealType>(sample_count), WithinAbs(4.0, 0.15));
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

TEST_CASE("IsoTarget defaults to no fluctuation model and target-typed id", "[radar][target]")
{
	radar::Platform platform("TargetPlatform");
	radar::IsoTarget target(&platform, "Iso", 1.0, 99);

	REQUIRE(target.getFluctuationModel() == nullptr);
	REQUIRE(SimIdGenerator::getType(target.getId()) == ObjectType::Target);
}

TEST_CASE("FileTarget multiplies interpolated azimuth and elevation RCS in the target frame", "[radar][target]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("target_rcs_interp"));
	removeIfExists(path);
	writeXmlTargetFile(path, xmlTargetFixture());

	radar::Platform platform("TargetPlatform");
	setStaticRotation(platform);

	constexpr SimId explicit_id = 4001;
	radar::FileTarget target(&platform, "File", path.string(), 12, explicit_id);

	math::SVec3 in_angle = unitDirection(0.6, 0.2);
	math::SVec3 out_angle = unitDirection(0.4, 0.6);

	const RealType expected_azimuth_rcs = 2.0 + 4.0 * 0.5;
	const RealType expected_elevation_rcs = 3.0 + 2.0 * 0.4;
	const RealType expected_rcs = expected_azimuth_rcs * expected_elevation_rcs;

	REQUIRE(target.getId() == explicit_id);
	REQUIRE(target.getFilename() == path.string());
	REQUIRE_THAT(target.getRcs(in_angle, out_angle, 0.0), WithinAbs(expected_rcs, 1e-12));

	removeIfExists(path);
}

TEST_CASE("FileTarget uses time-dependent platform rotation to look up body-frame RCS", "[radar][target]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("target_rcs_rotation"));
	removeIfExists(path);
	writeXmlTargetFile(path, xmlTargetFixture());

	radar::Platform platform("TargetPlatform");
	setConstantRotation(platform, 0.1, 0.2, 0.2, 0.1);

	radar::FileTarget target(&platform, "File", path.string(), 34, 5002);

	math::SVec3 in_angle = unitDirection(0.4, 0.2);
	math::SVec3 out_angle = unitDirection(0.5, 0.4);

	constexpr RealType time = 2.0;
	const RealType expected_azimuth_rcs = 2.0 + 4.0 * 0.2;
	const RealType expected_elevation_rcs = 3.0 + 2.0 * 0.1;
	const RealType expected_rcs = expected_azimuth_rcs * expected_elevation_rcs;

	REQUIRE_THAT(target.getRcs(in_angle, out_angle, time), WithinAbs(expected_rcs, 1e-12));

	removeIfExists(path);
}

TEST_CASE("FileTarget applies fluctuation model as a multiplicative RCS term", "[radar][target]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("target_rcs_fluctuation"));
	removeIfExists(path);
	writeXmlTargetFile(path, xmlTargetFixture());

	radar::Platform platform("TargetPlatform");
	setStaticRotation(platform);

	radar::FileTarget target(&platform, "File", path.string(), 56, 6003);
	auto model = std::make_unique<FixedRcsModel>(1.5);
	const auto* model_ptr = model.get();
	target.setFluctuationModel(std::move(model));

	math::SVec3 in_angle = unitDirection(1.0, 1.0);
	math::SVec3 out_angle = unitDirection(1.0, 1.0);

	REQUIRE(target.getFluctuationModel() == model_ptr);
	REQUIRE_THAT(target.getRcs(in_angle, out_angle, 0.0), WithinAbs(30.0 * 1.5, 1e-12));

	removeIfExists(path);
}

TEST_CASE("createFileTarget constructs FileTarget with explicit id and filename", "[radar][target]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("target_factory"));
	removeIfExists(path);
	writeXmlTargetFile(path, xmlTargetFixture());

	radar::Platform platform("TargetPlatform");
	setStaticRotation(platform);

	const auto target = radar::createFileTarget(&platform, "FileFactory", path.string(), 78, 7004);
	const auto* file_target = dynamic_cast<const radar::FileTarget*>(target.get());

	REQUIRE(file_target != nullptr);
	REQUIRE(file_target->getId() == 7004);
	REQUIRE(file_target->getFilename() == path.string());

	removeIfExists(path);
}

TEST_CASE("FileTarget skips incomplete sample nodes and uses the valid samples that remain", "[radar][target]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("target_partial_samples"));
	removeIfExists(path);
	writeXmlTargetFile(path,
					   R"(<target>
		<azimuth>
			<rcssample><angle>0.0</angle></rcssample>
			<rcssample><angle>1.0</angle><rcs>4.0</rcs></rcssample>
		</azimuth>
		<elevation>
			<rcssample><rcs>9.0</rcs></rcssample>
			<rcssample><angle>1.0</angle><rcs>3.0</rcs></rcssample>
		</elevation>
	</target>)");

	radar::Platform platform("TargetPlatform");
	setStaticRotation(platform);

	radar::FileTarget target(&platform, "File", path.string(), 90, 8005);
	math::SVec3 in_angle = unitDirection(1.0, 1.0);
	math::SVec3 out_angle = unitDirection(1.0, 1.0);

	REQUIRE_THAT(target.getRcs(in_angle, out_angle, 0.0), WithinAbs(12.0, 1e-12));

	removeIfExists(path);
}

TEST_CASE("FileTarget throws for missing target description file", "[radar][target]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("target_missing"));
	removeIfExists(path);

	radar::Platform platform("TargetPlatform");
	setStaticRotation(platform);

	REQUIRE_THROWS_AS(radar::FileTarget(&platform, "File", path.string(), 111, 9006), std::runtime_error);
}

TEST_CASE("FileTarget throws for malformed target description file", "[radar][target]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("target_malformed"));
	removeIfExists(path);

	{
		std::ofstream out(path, std::ios::binary);
		REQUIRE(out.is_open());
		out << "<target><azimuth><rcssample></target>";
	}

	radar::Platform platform("TargetPlatform");
	setStaticRotation(platform);

	REQUIRE_THROWS_AS(radar::FileTarget(&platform, "File", path.string(), 222, 10007), std::runtime_error);

	removeIfExists(path);
}

TEST_CASE("FileTarget currently throws at lookup time when an RCS axis has no samples", "[radar][target]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("target_missing_axis"));
	removeIfExists(path);
	writeXmlTargetFile(path,
					   R"(<target>
		<azimuth>
			<rcssample><angle>0.0</angle><rcs>2.0</rcs></rcssample>
		</azimuth>
		<elevation></elevation>
	</target>)");

	radar::Platform platform("TargetPlatform");
	setStaticRotation(platform);

	radar::FileTarget target(&platform, "File", path.string(), 333, 11008);
	math::SVec3 in_angle = unitDirection(0.0, 0.0);
	math::SVec3 out_angle = unitDirection(0.0, 0.0);

	REQUIRE_THROWS_AS(target.getRcs(in_angle, out_angle, 0.0), std::runtime_error);

	removeIfExists(path);
}

TEST_CASE("Target exposes initial seed", "[radar][target]")
{
	radar::Platform platform("TargetPlatform");
	radar::IsoTarget target(&platform, "Iso", 1.0, 12345);
	REQUIRE(target.getSeed() == 12345);
}

// TODO: FileTarget accepts XML files with missing/empty azimuth or elevation axes
// during construction and only fails later in getRcs(). Once the source validates
// required sample sets in the constructor, replace the lookup-time failure test
// with a constructor-level validation test.
