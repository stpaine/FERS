#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <complex>
#include <filesystem>
#include <fstream>
#include <highfive/highfive.hpp>
#include <string>
#include <system_error>
#include <vector>

#include "core/config.h"
#include "core/parameters.h"
#include "core/sim_id.h"
#include "serial/waveform_factory.h"
#include "signal/radar_signal.h"

using Catch::Matchers::WithinAbs;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;

		ParamGuard() : saved(params::params) {}

		~ParamGuard() { params::params = saved; }
	};

	std::string uniqueFileName(const std::string& prefix, const std::string& extension)
	{
		return prefix + "_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + extension;
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

	void writeCsvWaveform(const std::filesystem::path& path, const RealType sampleRate,
						  const std::vector<ComplexType>& samples)
	{
		std::ofstream out(path);
		REQUIRE(out.is_open());

		out << samples.size() << "\n";
		out << sampleRate << "\n";
		for (const auto& sample : samples)
		{
			out << sample << "\n";
		}
	}

	void writeHdf5Waveform(const std::filesystem::path& path, const std::vector<ComplexType>& samples)
	{
		HighFive::File file(path.string(), HighFive::File::Overwrite);

		std::vector<RealType> iValues;
		std::vector<RealType> qValues;
		iValues.reserve(samples.size());
		qValues.reserve(samples.size());

		for (const auto& sample : samples)
		{
			iValues.push_back(sample.real());
			qValues.push_back(sample.imag());
		}

		auto iGroup = file.createGroup("/I");
		auto qGroup = file.createGroup("/Q");
		iGroup.createDataSet<RealType>("value", HighFive::DataSpace::From(iValues)).write(iValues);
		qGroup.createDataSet<RealType>("value", HighFive::DataSpace::From(qValues)).write(qValues);
	}
}

TEST_CASE("Waveform factory loads CSV waveform metadata", "[serial][waveform_factory]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(9999.0);
	params::setOversampleRatio(2);

	const std::filesystem::path path = tempFilePath(uniqueFileName("waveform_factory_csv", ".csv"));
	removeIfExists(path);

	const std::vector<ComplexType> samples = {
		ComplexType(1.0, 0.5),
		ComplexType(-0.25, 1.5),
		ComplexType(0.75, -0.125),
	};
	constexpr RealType csvRate = 20.0;
	constexpr RealType power = 12.5;
	constexpr RealType carrier = 915.0e6;
	constexpr SimId explicitId = 42;

	writeCsvWaveform(path, csvRate, samples);

	const auto waveform = serial::loadWaveformFromFile("wave", path.string(), power, carrier, explicitId);

	REQUIRE(waveform->getName() == "wave");
	REQUIRE_THAT(waveform->getPower(), WithinAbs(power, 1e-12));
	REQUIRE_THAT(waveform->getCarrier(), WithinAbs(carrier, 1e-3));
	REQUIRE(waveform->getId() == explicitId);
	REQUIRE(waveform->getFilename().has_value());
	REQUIRE(*waveform->getFilename() == path.string());
	REQUIRE_THAT(waveform->getRate(), WithinAbs(csvRate * params::oversampleRatio(), 1e-12));
	REQUIRE_THAT(waveform->getLength(), WithinAbs(static_cast<RealType>(samples.size()) / csvRate, 1e-12));

	removeIfExists(path);
}

TEST_CASE("Waveform factory loads HDF5 waveform metadata", "[serial][waveform_factory]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(64.0);

	const std::filesystem::path path = tempFilePath(uniqueFileName("waveform_factory_hdf5", ".h5"));
	removeIfExists(path);

	const std::vector<ComplexType> samples = {
		ComplexType(1.0, -1.0),
		ComplexType(0.5, 0.25),
		ComplexType(-0.75, 2.0),
		ComplexType(0.0, -0.5),
	};
	constexpr RealType power = 3.0;
	constexpr RealType carrier = 2.45e9;
	constexpr SimId explicitId = 77;

	writeHdf5Waveform(path, samples);

	const auto waveform = serial::loadWaveformFromFile("wave", path.string(), power, carrier, explicitId);

	REQUIRE(waveform->getName() == "wave");
	REQUIRE_THAT(waveform->getPower(), WithinAbs(power, 1e-12));
	REQUIRE_THAT(waveform->getCarrier(), WithinAbs(carrier, 1e-3));
	REQUIRE(waveform->getId() == explicitId);
	REQUIRE(waveform->getFilename().has_value());
	REQUIRE(*waveform->getFilename() == path.string());
	REQUIRE_THAT(waveform->getRate(), WithinAbs(params::rate(), 1e-12));
	REQUIRE_THAT(waveform->getLength(), WithinAbs(static_cast<RealType>(samples.size()) / params::rate(), 1e-12));

	removeIfExists(path);
}

TEST_CASE("Waveform factory preserves explicit waveform id", "[serial][waveform_factory]")
{
	ParamGuard guard;
	params::params.reset();

	const std::filesystem::path path = tempFilePath(uniqueFileName("waveform_factory_id", ".csv"));
	removeIfExists(path);

	writeCsvWaveform(path, 8.0, {ComplexType(0.0, 1.0)});

	constexpr SimId explicitId = (static_cast<SimId>(0x1234) << 32) | 0x55AA;
	const auto waveform = serial::loadWaveformFromFile("wave", path.string(), 1.0, 10.0, explicitId);

	REQUIRE(waveform->getId() == explicitId);

	removeIfExists(path);
}

TEST_CASE("Waveform factory throws for unsupported extension", "[serial][waveform_factory]")
{
	REQUIRE_THROWS_AS(serial::loadWaveformFromFile("wave", "unsupported_waveform.txt", 1.0, 2.0), std::runtime_error);
}

TEST_CASE("Waveform factory throws for missing CSV file", "[serial][waveform_factory]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("waveform_factory_missing", ".csv"));
	removeIfExists(path);

	REQUIRE_THROWS_AS(serial::loadWaveformFromFile("wave", path.string(), 1.0, 2.0), std::runtime_error);
}

TEST_CASE("Waveform factory throws for missing HDF5 file", "[serial][waveform_factory]")
{
	ParamGuard guard;
	params::params.reset();
	params::setRate(100.0);

	const std::filesystem::path path = tempFilePath(uniqueFileName("waveform_factory_missing", ".h5"));
	removeIfExists(path);

	REQUIRE_THROWS_AS(serial::loadWaveformFromFile("wave", path.string(), 1.0, 2.0), std::runtime_error);
}

TEST_CASE("Waveform factory throws for incomplete CSV waveform data", "[serial][waveform_factory]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("waveform_factory_incomplete", ".csv"));
	removeIfExists(path);

	{
		std::ofstream out(path);
		REQUIRE(out.is_open());
		out << 4 << "\n";
		out << 16.0 << "\n";
		out << ComplexType(1.0, 0.0) << "\n";
		out << ComplexType(0.0, 1.0) << "\n";
		out << ComplexType(-1.0, 0.5) << "\n";
	}

	REQUIRE_THROWS_AS(serial::loadWaveformFromFile("wave", path.string(), 1.0, 2.0), std::runtime_error);

	removeIfExists(path);
}

TEST_CASE("Waveform factory throws for malformed CSV header", "[serial][waveform_factory]")
{
	const std::filesystem::path path = tempFilePath(uniqueFileName("waveform_factory_bad_header", ".csv"));
	removeIfExists(path);

	{
		std::ofstream out(path);
		REQUIRE(out.is_open());
		out << "not-a-count\n";
		out << "still-not-a-rate\n";
	}

	REQUIRE_THROWS_AS(serial::loadWaveformFromFile("wave", path.string(), 1.0, 2.0), std::runtime_error);

	removeIfExists(path);
}

TEST_CASE("Waveform factory rejects invalid CSV sample counts", "[serial][waveform_factory]")
{
	SECTION("fractional sample count")
	{
		const std::filesystem::path path = tempFilePath(uniqueFileName("waveform_factory_fractional_count", ".csv"));
		removeIfExists(path);

		{
			std::ofstream out(path);
			REQUIRE(out.is_open());
			out << 2.5 << "\n";
			out << 16.0 << "\n";
			out << ComplexType(1.0, 0.0) << "\n";
			out << ComplexType(0.0, 1.0) << "\n";
		}

		REQUIRE_THROWS_AS(serial::loadWaveformFromFile("wave", path.string(), 1.0, 2.0), std::runtime_error);
		removeIfExists(path);
	}

	SECTION("negative sample count")
	{
		const std::filesystem::path path = tempFilePath(uniqueFileName("waveform_factory_negative_count", ".csv"));
		removeIfExists(path);

		{
			std::ofstream out(path);
			REQUIRE(out.is_open());
			out << -1 << "\n";
			out << 16.0 << "\n";
		}

		REQUIRE_THROWS_AS(serial::loadWaveformFromFile("wave", path.string(), 1.0, 2.0), std::runtime_error);
		removeIfExists(path);
	}
}
