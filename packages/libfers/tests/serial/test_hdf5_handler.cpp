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

TEST_CASE("HDF5 writer adds backward-compatible output metadata attributes", "[serial][hdf5]")
{
	const std::string path = tempFilePath(uniqueFileName("metadata"));
	removeIfExists(path);

	RateGuard rate_guard(4'000.0);
	core::OutputFileMetadata file_metadata{.receiver_id = 42,
										   .receiver_name = "RxMetadata",
										   .mode = "pulsed",
										   .path = path,
										   .total_samples = 2,
										   .sample_start = 0,
										   .sample_end_exclusive = 2,
										   .pulse_count = 1,
										   .min_pulse_length_samples = 2,
										   .max_pulse_length_samples = 2,
										   .uniform_pulse_length = true};
	core::PulseChunkMetadata chunk_metadata{.chunk_index = 0,
											.i_dataset = "chunk_000000_I",
											.q_dataset = "chunk_000000_Q",
											.start_time = 1.25,
											.sample_count = 2,
											.sample_start = 0,
											.sample_end_exclusive = 2};
	file_metadata.chunks.push_back(chunk_metadata);

	{
		HighFive::File file(path, HighFive::File::Overwrite);
		serial::addChunkToFile(file, {ComplexType(1.0, 0.0), ComplexType(2.0, 0.5)}, 1.25, 8.0, 0, &chunk_metadata);
		std::scoped_lock lock(serial::hdf5_global_mutex);
		serial::writeOutputFileMetadataAttributes(file, file_metadata);
	}

	{
		HighFive::File file(path, HighFive::File::ReadOnly);
		unsigned schema_version = 0;
		unsigned long long total_samples = 0;
		std::string receiver_name;
		std::string metadata_json;
		file.getAttribute("fers_metadata_schema_version").read(schema_version);
		file.getAttribute("receiver_name").read(receiver_name);
		file.getAttribute("total_samples").read(total_samples);
		file.getAttribute("fers_metadata_json").read(metadata_json);

		REQUIRE(schema_version == 1U);
		REQUIRE(receiver_name == "RxMetadata");
		REQUIRE(total_samples == 2ULL);
		REQUIRE(metadata_json.find("\"receiver_name\": \"RxMetadata\"") != std::string::npos);

		const auto i_dataset = file.getDataSet("chunk_000000_I");
		unsigned chunk_index = 1;
		unsigned long long sample_count = 0;
		i_dataset.getAttribute("chunk_index").read(chunk_index);
		i_dataset.getAttribute("sample_count").read(sample_count);

		REQUIRE(chunk_index == 0U);
		REQUIRE(sample_count == 2ULL);
	}

	removeIfExists(path);
}

TEST_CASE("HDF5 writer exposes FMCW segment metadata without cw_segments JSON alias", "[serial][hdf5]")
{
	const std::string path = tempFilePath(uniqueFileName("fmcw_metadata"));
	removeIfExists(path);

	core::OutputFileMetadata file_metadata{.receiver_id = 7,
										   .receiver_name = "RxFmcw",
										   .mode = "fmcw",
										   .path = path,
										   .total_samples = 20,
										   .sample_start = 0,
										   .sample_end_exclusive = 20,
										   .fmcw = core::FmcwMetadata{.chirp_bandwidth = 100.0,
																	  .chirp_duration = 0.01,
																	  .chirp_period = 0.02,
																	  .chirp_rate = 10'000.0,
																	  .chirp_rate_signed = -10'000.0,
																	  .chirp_direction = "down",
																	  .start_frequency_offset = -50.0,
																	  .chirp_count = 3}};
	file_metadata.streaming_segments.push_back({.start_time = 1.0,
												.end_time = 1.2,
												.sample_count = 10,
												.sample_start = 0,
												.sample_end_exclusive = 10,
												.first_chirp_start_time = 1.0,
												.emitted_chirp_count = 3});
	file_metadata.streaming_segments.push_back({.start_time = 2.0,
												.end_time = 2.1,
												.sample_count = 10,
												.sample_start = 10,
												.sample_end_exclusive = 20,
												.first_chirp_start_time = 2.0,
												.emitted_chirp_count = 3});

	{
		HighFive::File file(path, HighFive::File::Overwrite);
		std::scoped_lock lock(serial::hdf5_global_mutex);
		serial::writeOutputFileMetadataAttributes(file, file_metadata);
	}

	{
		HighFive::File file(path, HighFive::File::ReadOnly);
		unsigned long long segment_count = 0;
		RealType signed_rate = 0.0;
		std::string direction;
		std::vector<RealType> first_chirp_starts;
		std::vector<unsigned long long> emitted_chirp_counts;
		std::string metadata_json;

		file.getAttribute("streaming_segment_count").read(segment_count);
		file.getAttribute("fmcw_chirp_rate_signed").read(signed_rate);
		file.getAttribute("fmcw_chirp_direction").read(direction);
		file.getAttribute("streaming_first_chirp_start_time").read(first_chirp_starts);
		file.getAttribute("streaming_emitted_chirp_count").read(emitted_chirp_counts);
		file.getAttribute("fers_metadata_json").read(metadata_json);

		REQUIRE(segment_count == 2ULL);
		REQUIRE_THAT(signed_rate, WithinAbs(-10'000.0, 1e-12));
		REQUIRE(direction == "down");
		REQUIRE(first_chirp_starts.size() == 2u);
		REQUIRE_THAT(first_chirp_starts[0], WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(first_chirp_starts[1], WithinAbs(2.0, 1e-12));
		REQUIRE(emitted_chirp_counts == std::vector<unsigned long long>{3ULL, 3ULL});
		REQUIRE(metadata_json.find("\"streaming_segments\"") != std::string::npos);
		REQUIRE(metadata_json.find("\"chirp_direction\": \"down\"") != std::string::npos);
		REQUIRE(metadata_json.find("\"chirp_rate_signed\": -10000.0") != std::string::npos);
		REQUIRE(metadata_json.find("\"cw_segments\"") == std::string::npos);
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
