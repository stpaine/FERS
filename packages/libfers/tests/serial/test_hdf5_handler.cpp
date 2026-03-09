#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <filesystem>
#include <highfive/highfive.hpp>
#include <string>
#include <vector>

#include "core/config.h"
#include "core/parameters.h"
#include "serial/hdf5_handler.h"

using Catch::Matchers::WithinAbs;

namespace
{
	struct RateGuard
	{
		RealType previous_rate = params::rate();

		explicit RateGuard(RealType rate) { params::setRate(rate); }

		~RateGuard() { params::params.rate = previous_rate; }
	};

	std::string tempFilePath(const std::string& name)
	{
		return (std::filesystem::temp_directory_path() / name).string();
	}

	std::string uniqueFileName(const std::string& prefix)
	{
		return prefix + "_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".h5";
	}

	void removeIfExists(const std::string& path)
	{
		std::error_code ec;
		std::filesystem::remove(path, ec);
	}

	HighFive::DataSet createVectorDataSet(HighFive::Group& group, const std::string& name,
										  const std::vector<double>& data)
	{
		return group.createDataSet<double>(name, HighFive::DataSpace::From(data));
	}
}

TEST_CASE("readPulseData throws for missing file", "[serial][hdf5]")
{
	std::vector<ComplexType> data;
	const std::string missing = tempFilePath("missing_pulse_data.h5");
	removeIfExists(missing);

	REQUIRE_THROWS_AS(serial::readPulseData(missing, data), std::runtime_error);
}

TEST_CASE("readPulseData loads I/Q into complex vector", "[serial][hdf5]")
{
	const std::string path = tempFilePath(uniqueFileName("pulse_data"));
	removeIfExists(path);

	{
		HighFive::File file(path, HighFive::File::Overwrite);
		const std::vector<double> i_values{1.0, -2.0, 3.5};
		const std::vector<double> q_values{0.5, -1.5, 2.25};

		auto i_group = file.createGroup("/I");
		auto q_group = file.createGroup("/Q");
		createVectorDataSet(i_group, "value", i_values).write(i_values);
		createVectorDataSet(q_group, "value", q_values).write(q_values);
	}

	std::vector<ComplexType> data;
	serial::readPulseData(path, data);

	REQUIRE(data.size() == 3u);
	REQUIRE_THAT(data[0].real(), WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(data[0].imag(), WithinAbs(0.5, 1e-12));
	REQUIRE_THAT(data[1].real(), WithinAbs(-2.0, 1e-12));
	REQUIRE_THAT(data[1].imag(), WithinAbs(-1.5, 1e-12));
	REQUIRE_THAT(data[2].real(), WithinAbs(3.5, 1e-12));
	REQUIRE_THAT(data[2].imag(), WithinAbs(2.25, 1e-12));

	removeIfExists(path);
}

TEST_CASE("readPulseData throws when Q size mismatches I", "[serial][hdf5]")
{
	const std::string path = tempFilePath(uniqueFileName("pulse_mismatch"));
	removeIfExists(path);

	{
		HighFive::File file(path, HighFive::File::Overwrite);
		const std::vector<double> i_values{1.0, 2.0, 3.0};
		const std::vector<double> q_values{1.0, 2.0};

		auto i_group = file.createGroup("/I");
		auto q_group = file.createGroup("/Q");
		createVectorDataSet(i_group, "value", i_values).write(i_values);
		createVectorDataSet(q_group, "value", q_values).write(q_values);
	}

	std::vector<ComplexType> data;
	REQUIRE_THROWS_AS(serial::readPulseData(path, data), std::runtime_error);

	removeIfExists(path);
}

TEST_CASE("addChunkToFile writes I/Q datasets with attributes", "[serial][hdf5]")
{
	const std::string path = tempFilePath(uniqueFileName("chunks"));
	removeIfExists(path);

	RateGuard rate_guard(2'000.0);
	const std::vector<ComplexType> samples = {ComplexType(1.25, -0.75), ComplexType(-2.5, 3.0)};
	const RealType time = 1.5;
	const RealType fullscale = 4096.0;

	{
		HighFive::File file(path, HighFive::File::Overwrite);
		serial::addChunkToFile(file, samples, time, fullscale, 7);
	}

	{
		HighFive::File file(path, HighFive::File::ReadOnly);
		const auto i_dataset = file.getDataSet("chunk_000007_I");
		const auto q_dataset = file.getDataSet("chunk_000007_Q");

		std::vector<RealType> i_read;
		std::vector<RealType> q_read;
		i_dataset.read(i_read);
		q_dataset.read(q_read);

		REQUIRE(i_read.size() == 2u);
		REQUIRE(q_read.size() == 2u);
		REQUIRE_THAT(i_read[0], WithinAbs(1.25, 1e-12));
		REQUIRE_THAT(i_read[1], WithinAbs(-2.5, 1e-12));
		REQUIRE_THAT(q_read[0], WithinAbs(-0.75, 1e-12));
		REQUIRE_THAT(q_read[1], WithinAbs(3.0, 1e-12));

		RealType time_attr = 0.0;
		RealType rate_attr = 0.0;
		RealType fullscale_attr = 0.0;
		i_dataset.getAttribute("time").read(time_attr);
		i_dataset.getAttribute("rate").read(rate_attr);
		i_dataset.getAttribute("fullscale").read(fullscale_attr);

		REQUIRE_THAT(time_attr, WithinAbs(time, 1e-12));
		REQUIRE_THAT(rate_attr, WithinAbs(params::rate(), 1e-12));
		REQUIRE_THAT(fullscale_attr, WithinAbs(fullscale, 1e-12));

		q_dataset.getAttribute("time").read(time_attr);
		q_dataset.getAttribute("rate").read(rate_attr);
		q_dataset.getAttribute("fullscale").read(fullscale_attr);

		REQUIRE_THAT(time_attr, WithinAbs(time, 1e-12));
		REQUIRE_THAT(rate_attr, WithinAbs(params::rate(), 1e-12));
		REQUIRE_THAT(fullscale_attr, WithinAbs(fullscale, 1e-12));
	}

	removeIfExists(path);
}

TEST_CASE("addChunkToFile throws when dataset already exists", "[serial][hdf5]")
{
	const std::string path = tempFilePath(uniqueFileName("chunk_exists"));
	removeIfExists(path);

	RateGuard rate_guard(1'000.0);
	const std::vector<ComplexType> samples = {ComplexType(1.0, 1.0)};

	{
		HighFive::File file(path, HighFive::File::Overwrite);
		serial::addChunkToFile(file, samples, 0.0, 1.0, 1);
		REQUIRE_THROWS_AS(serial::addChunkToFile(file, samples, 0.0, 1.0, 1), std::runtime_error);
	}

	removeIfExists(path);
}

TEST_CASE("readPattern returns 2D dataset", "[serial][hdf5]")
{
	const std::string path = tempFilePath(uniqueFileName("pattern"));
	removeIfExists(path);

	{
		HighFive::File file(path, HighFive::File::Overwrite);
		std::vector<std::vector<RealType>> pattern = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
		auto dataset = file.createDataSet<RealType>("pattern", HighFive::DataSpace::From(pattern));
		dataset.write(pattern);
	}

	const auto data = serial::readPattern(path, "pattern");
	REQUIRE(data.size() == 2u);
	REQUIRE(data[0].size() == 3u);
	REQUIRE_THAT(data[0][0], WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(data[0][2], WithinAbs(3.0, 1e-12));
	REQUIRE_THAT(data[1][0], WithinAbs(4.0, 1e-12));
	REQUIRE_THAT(data[1][2], WithinAbs(6.0, 1e-12));

	removeIfExists(path);
}

TEST_CASE("readPattern throws for non-2D dataset", "[serial][hdf5]")
{
	const std::string path = tempFilePath(uniqueFileName("pattern_bad"));
	removeIfExists(path);

	{
		HighFive::File file(path, HighFive::File::Overwrite);
		std::vector<RealType> values = {1.0, 2.0, 3.0};
		auto dataset = file.createDataSet<RealType>("pattern", HighFive::DataSpace::From(values));
		dataset.write(values);
	}

	REQUIRE_THROWS_AS(serial::readPattern(path, "pattern"), std::runtime_error);

	removeIfExists(path);
}

TEST_CASE("readPattern throws when dataset missing", "[serial][hdf5]")
{
	const std::string path = tempFilePath(uniqueFileName("pattern_missing"));
	removeIfExists(path);

	{
		HighFive::File file(path, HighFive::File::Overwrite);
		std::vector<std::vector<RealType>> pattern = {{0.0}};
		auto dataset = file.createDataSet<RealType>("other", HighFive::DataSpace::From(pattern));
		dataset.write(pattern);
	}

	REQUIRE_THROWS_AS(serial::readPattern(path, "pattern"), std::runtime_error);

	removeIfExists(path);
}

// TODO: Add tests that force HighFive to throw during dataset creation or attribute
// creation in addChunkToFile; this would require fault injection or an artificial
// filesystem failure that is not reliably reproducible in unit tests.
