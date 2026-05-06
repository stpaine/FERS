#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <chrono>
#include <filesystem>
#include <highfive/highfive.hpp>
#include <string>
#include <vector>

#include "core/config.h"
#include "core/parameters.h"
#include "processing/finalizer_pipeline.h"
#include "signal/dsp_filters.h"

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
	};

	std::filesystem::path uniqueTempPath(const std::string& prefix)
	{
		return std::filesystem::temp_directory_path() /
			(prefix + "_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".h5");
	}

	void removeIfExists(const std::filesystem::path& path)
	{
		std::error_code ec;
		std::filesystem::remove(path, ec);
	}
}

TEST_CASE("applyDownsamplingAndQuantization reduces to direct normalization when no oversampling is configured",
		  "[processing][finalizer][io]")
{
	ParamGuard guard;
	params::setOversampleRatio(1);
	params::setAdcBits(0);

	std::vector<ComplexType> buffer = {ComplexType{2.0, 0.0}, ComplexType{-1.0, -3.0}};

	const RealType fullscale = processing::pipeline::applyDownsamplingAndQuantization(buffer);

	REQUIRE_THAT(fullscale, WithinAbs(3.0, 1e-12));
	REQUIRE_THAT(buffer[0].real(), WithinAbs(2.0 / 3.0, 1e-12));
	REQUIRE_THAT(buffer[0].imag(), WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(buffer[1].real(), WithinAbs(-1.0 / 3.0, 1e-12));
	REQUIRE_THAT(buffer[1].imag(), WithinAbs(-1.0, 1e-12));
}

TEST_CASE("applyDownsamplingAndQuantization preserves low-frequency oversampled content after decimation",
		  "[processing][finalizer][io]")
{
	ParamGuard guard;
	params::params.filter_length = 8;
	params::setOversampleRatio(2);
	params::setAdcBits(0);

	const size_t sample_count = 64;
	const RealType frequency = 0.05;
	std::vector<ComplexType> baseband(sample_count);
	for (size_t i = 0; i < baseband.size(); ++i)
	{
		const RealType phase = 2.0 * PI * frequency * static_cast<RealType>(i);
		baseband[i] = ComplexType{std::cos(phase), std::sin(phase)};
	}

	std::vector<ComplexType> oversampled(sample_count * params::oversampleRatio());
	fers_signal::upsample(baseband, static_cast<unsigned>(baseband.size()), oversampled);

	const RealType fullscale = processing::pipeline::applyDownsamplingAndQuantization(oversampled);

	REQUIRE(oversampled.size() == baseband.size());
	REQUIRE(fullscale > 0.9);
	REQUIRE(fullscale < 1.2);
	ComplexType correlation{0.0, 0.0};
	RealType output_energy = 0.0;
	RealType reference_energy = 0.0;
	for (size_t i = 10; i < oversampled.size() - 10; ++i)
	{
		correlation += oversampled[i] * std::conj(baseband[i]);
		output_energy += std::norm(oversampled[i]);
		reference_energy += std::norm(baseband[i]);
	}
	const RealType normalized_correlation = std::abs(correlation) / std::sqrt(output_energy * reference_energy);

	REQUIRE(normalized_correlation > 0.98);
}

TEST_CASE("applyDownsamplingAndQuantization fails fast when oversample ratio exceeds fixed-filter limit",
		  "[processing][finalizer][io]")
{
	ParamGuard guard;
	params::params.oversample_ratio = 16;
	params::setAdcBits(0);

	std::vector<ComplexType> buffer = {ComplexType{1.0, 0.0}, ComplexType{0.0, 1.0}, ComplexType{-1.0, 0.5}};

	REQUIRE_THROWS_WITH(processing::pipeline::applyDownsamplingAndQuantization(buffer),
						ContainsSubstring("Oversampling ratios > 8 are not supported"));
}

TEST_CASE("exportStreamingToHdf5 writes physically meaningful datasets and metadata", "[processing][finalizer][io]")
{
	ParamGuard guard;
	params::setRate(2000.0);
	params::setTime(1.25, 3.0);

	const std::filesystem::path path = uniqueTempPath("cw_export");
	removeIfExists(path);

	const std::vector<ComplexType> iq_buffer = {ComplexType{1.5, -0.5}, ComplexType{-2.0, 3.0}};
	processing::pipeline::exportStreamingToHdf5(path.string(), iq_buffer, 4096.0, 9.0e8);

	HighFive::File file(path.string(), HighFive::File::ReadOnly);
	std::vector<RealType> i_data;
	std::vector<RealType> q_data;
	file.getDataSet("I_data").read(i_data);
	file.getDataSet("Q_data").read(q_data);

	RealType sampling_rate = 0.0;
	RealType start_time = 0.0;
	RealType fullscale = 0.0;
	RealType reference_frequency = 0.0;
	file.getAttribute("sampling_rate").read(sampling_rate);
	file.getAttribute("start_time").read(start_time);
	file.getAttribute("fullscale").read(fullscale);
	file.getAttribute("reference_carrier_frequency").read(reference_frequency);

	REQUIRE(i_data.size() == 2u);
	REQUIRE(q_data.size() == 2u);
	REQUIRE_THAT(i_data[0], WithinAbs(1.5, 1e-12));
	REQUIRE_THAT(i_data[1], WithinAbs(-2.0, 1e-12));
	REQUIRE_THAT(q_data[0], WithinAbs(-0.5, 1e-12));
	REQUIRE_THAT(q_data[1], WithinAbs(3.0, 1e-12));
	REQUIRE_THAT(sampling_rate, WithinAbs(params::rate(), 1e-12));
	REQUIRE_THAT(start_time, WithinAbs(params::startTime(), 1e-12));
	REQUIRE_THAT(fullscale, WithinAbs(4096.0, 1e-12));
	REQUIRE_THAT(reference_frequency, WithinAbs(9.0e8, 1e-12));

	removeIfExists(path);
}

TEST_CASE("exportStreamingToHdf5 handles file creation failures without throwing", "[processing][finalizer][io]")
{
	const auto parent = std::filesystem::temp_directory_path() /
		("missing_finalizer_dir_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
	const auto path = parent / "results.h5";

	std::error_code ec;
	std::filesystem::remove_all(parent, ec);

	REQUIRE_NOTHROW(processing::pipeline::exportStreamingToHdf5(path.string(), {}, 1.0, 2.0));
	REQUIRE_FALSE(std::filesystem::exists(path));
}

// TODO: The HDF5 error branch is observable only indirectly because exportStreamingToHdf5
// logs failures instead of surfacing the exception details to the caller.
