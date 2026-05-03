#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <filesystem>
#include <format>
#include <highfive/highfive.hpp>
#include <memory>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/parameters.h"
#include "core/rendering_job.h"
#include "core/sim_threading.h"
#include "math/coord.h"
#include "processing/finalizer.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/response.h"
#include "signal/if_resampler.h"
#include "signal/radar_signal.h"
#include "timing/prototype_timing.h"
#include "timing/timing.h"

using Catch::Matchers::WithinAbs;

namespace
{
	struct ParamGuard
	{
		params::Parameters saved;
		ParamGuard() : saved(params::params) {}
		~ParamGuard() { params::params = saved; }
	};

	struct FixedSignal final : public fers_signal::Signal
	{
		std::vector<ComplexType> data;

		std::vector<ComplexType> render(const std::vector<interp::InterpPoint>&, unsigned& size,
										RealType) const override
		{
			size = static_cast<unsigned>(data.size());
			return data;
		}
	};

	std::filesystem::path resultPath(const std::filesystem::path& dir, const std::string& receiver_name)
	{
		return dir / (receiver_name + "_results.h5");
	}

	void setupPlatform(radar::Platform& platform, const math::Vec3& position)
	{
		platform.getMotionPath()->addCoord(math::Coord{position, 0.0});
		platform.getMotionPath()->finalize();
		platform.getRotationPath()->addCoord(math::RotationCoord{0.0, 0.0, 0.0});
		platform.getRotationPath()->finalize();
	}

	std::string uniqueName(const std::string& prefix)
	{
		return prefix + "_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
	}

	void removeIfExists(const std::filesystem::path& path)
	{
		std::error_code ec;
		std::filesystem::remove(path, ec);
	}

	struct OwnedTiming
	{
		std::unique_ptr<timing::PrototypeTiming> prototype;
		std::shared_ptr<timing::Timing> timing;
	};

	OwnedTiming makeQuietTiming(const std::string& name, unsigned seed, RealType frequency)
	{
		OwnedTiming owned{};
		owned.prototype = std::make_unique<timing::PrototypeTiming>(name);
		owned.prototype->setFrequency(frequency);
		owned.timing = std::make_shared<timing::Timing>(name, seed);
		owned.timing->initializeModel(owned.prototype.get());
		return owned;
	}

	OwnedTiming makePhaseTiming(const std::string& name, unsigned seed, RealType frequency, RealType phase_offset,
								bool sync_on_pulse)
	{
		OwnedTiming owned{};
		owned.prototype = std::make_unique<timing::PrototypeTiming>(name);
		owned.prototype->setFrequency(frequency);
		owned.prototype->setPhaseOffset(phase_offset);
		if (sync_on_pulse)
		{
			owned.prototype->setSyncOnPulse();
		}
		owned.timing = std::make_shared<timing::Timing>(name, seed);
		owned.timing->initializeModel(owned.prototype.get());
		return owned;
	}

	struct FmcwTxFixture
	{
		radar::Platform platform;
		std::unique_ptr<fers_signal::RadarSignal> wave;
		radar::Transmitter transmitter;

		FmcwTxFixture(const std::string& name, SimId tx_id, SimId waveform_id, RealType chirp_bandwidth,
					  RealType chirp_duration, RealType chirp_period, RealType start_frequency_offset,
					  std::optional<std::size_t> chirp_count) :
			platform(name + "Platform"),
			wave(std::make_unique<fers_signal::RadarSignal>(
				name + "Wave", 1.0, 10.0e9, chirp_duration,
				std::make_unique<fers_signal::FmcwChirpSignal>(chirp_bandwidth, chirp_duration, chirp_period,
															   start_frequency_offset, chirp_count),
				waveform_id)),
			transmitter(&platform, name, radar::OperationMode::FMCW_MODE, tx_id)
		{
			transmitter.setSignal(wave.get());
		}
	};

	std::unique_ptr<serial::Response>
	makeFixedResponse(const radar::Transmitter* transmitter,
					  std::vector<std::unique_ptr<fers_signal::RadarSignal>>& wave_store,
					  const std::vector<ComplexType>& samples, RealType sample_rate, RealType start_time)
	{
		auto signal = std::make_unique<FixedSignal>();
		signal->data = samples;
		signal->load(samples, static_cast<unsigned>(samples.size()), sample_rate);

		auto wave = std::make_unique<fers_signal::RadarSignal>(
			"wave", 1.0, 1.0e9, static_cast<RealType>(samples.size()) / sample_rate, std::move(signal));
		const auto* wave_ptr = wave.get();
		wave_store.push_back(std::move(wave));

		auto response = std::make_unique<serial::Response>(wave_ptr, transmitter);
		response->addInterpPoint({1.0, start_time, 0.0, 0.0});
		response->addInterpPoint({1.0, start_time + static_cast<RealType>(samples.size() - 1) / sample_rate, 0.0, 0.0});
		return response;
	}

	std::vector<RealType> readDataset(const HighFive::File& file, const std::string& name)
	{
		std::vector<RealType> values;
		file.getDataSet(name).read(values);
		return values;
	}

	nlohmann::json readMetadataJson(const std::filesystem::path& output_path)
	{
		HighFive::File file(output_path.string(), HighFive::File::ReadOnly);
		std::string metadata_json;
		file.getAttribute("fers_metadata_json").read(metadata_json);
		return nlohmann::json::parse(metadata_json);
	}

	struct ProgressCall
	{
		std::string message;
		int current;
		int total;
	};
}

TEST_CASE("finalizeStreamingReceiver exits cleanly when no streaming samples were collected", "[processing][finalizer]")
{
	ParamGuard guard;
	params::setRate(4.0);
	params::setOversampleRatio(1);

	const std::string receiver_name = uniqueName("cw_empty");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("cw_empty_dir");
	std::filesystem::create_directories(out_dir);
	const auto output_path = resultPath(out_dir, receiver_name);
	removeIfExists(output_path);

	radar::Platform platform("RxPlatform");
	radar::Receiver receiver(&platform, receiver_name, 55, radar::OperationMode::CW_MODE);
	auto timing_owner = makeQuietTiming("quiet", 11, 77.0);
	receiver.setTiming(timing_owner.timing);

	std::vector<ProgressCall> progress_calls;
	auto reporter =
		std::make_shared<core::ProgressReporter>([&progress_calls](const std::string& msg, int current, int total)
												 { progress_calls.push_back({msg, current, total}); });

	processing::finalizeStreamingReceiver(&receiver, nullptr, reporter, out_dir.string());

	REQUIRE_FALSE(std::filesystem::exists(output_path));
	REQUIRE(progress_calls.empty());

	std::filesystem::remove_all(out_dir);
}

TEST_CASE("finalizeStreamingReceiver adds logged pulsed interference without rotating the summed streaming buffer and "
		  "exports HDF5",
		  "[processing][finalizer]")
{
	ParamGuard guard;
	params::setTime(0.0, 1.0);
	params::setRate(4.0);
	params::setOversampleRatio(1);
	params::setAdcBits(0);

	const std::string receiver_name = uniqueName("cw_finalize");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("cw_finalize_dir");
	std::filesystem::create_directories(out_dir);
	const auto output_path = resultPath(out_dir, receiver_name);
	removeIfExists(output_path);

	radar::Platform rx_platform("RxPlatform");
	radar::Receiver receiver(&rx_platform, receiver_name, 56, radar::OperationMode::CW_MODE);
	auto timing_owner = makePhaseTiming("phase_clk", 22, 77.0, PI / 2.0, false);
	receiver.setTiming(timing_owner.timing);
	receiver.setNoiseTemperature(0.0);
	receiver.prepareStreamingData(4);

	radar::Platform tx_platform("TxPlatform");
	radar::Transmitter transmitter(&tx_platform, "TxA", radar::OperationMode::PULSED_MODE, 601);

	std::vector<std::unique_ptr<fers_signal::RadarSignal>> wave_store;
	receiver.addInterferenceToLog(
		makeFixedResponse(&transmitter, wave_store, {ComplexType{1.0, 0.0}, ComplexType{1.0, 0.0}}, 4.0, 0.25));

	std::vector<ProgressCall> progress_calls;
	auto reporter =
		std::make_shared<core::ProgressReporter>([&progress_calls](const std::string& msg, int current, int total)
												 { progress_calls.push_back({msg, current, total}); });

	processing::finalizeStreamingReceiver(&receiver, nullptr, reporter, out_dir.string());

	{
		HighFive::File file(output_path.string(), HighFive::File::ReadOnly);
		const auto i_data = readDataset(file, "I_data");
		const auto q_data = readDataset(file, "Q_data");

		RealType fullscale = 0.0;
		RealType reference_frequency = 0.0;
		file.getAttribute("fullscale").read(fullscale);
		file.getAttribute("reference_carrier_frequency").read(reference_frequency);

		REQUIRE(i_data.size() == 4u);
		REQUIRE(q_data.size() == 4u);
		REQUIRE_THAT(i_data[0], WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(i_data[1], WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(i_data[2], WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(i_data[3], WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(q_data[0], WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(q_data[1], WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(q_data[2], WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(q_data[3], WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(fullscale, WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(reference_frequency, WithinAbs(77.0, 1e-12));
	}

	const std::vector<int> expected_progress = {0, 25, 50, 75, 100};
	REQUIRE(progress_calls.size() == expected_progress.size());
	REQUIRE(progress_calls.front().message == std::format("Finalizing CW Receiver {}", receiver_name));
	for (size_t i = 0; i < expected_progress.size(); ++i)
	{
		REQUIRE(progress_calls[i].current == expected_progress[i]);
	}

	std::filesystem::remove_all(out_dir);
}

TEST_CASE("finalizeStreamingReceiver accepts oversampled pulsed interference for non-dechirped CW receivers",
		  "[processing][finalizer][interference]")
{
	ParamGuard guard;
	params::setTime(0.0, 1.0);
	params::setRate(4.0);
	params::setOversampleRatio(2);
	params::setAdcBits(0);
	params::params.filter_length = 8;

	const std::string receiver_name = uniqueName("cw_finalize_oversampled");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("cw_finalize_oversampled_dir");
	std::filesystem::create_directories(out_dir);
	const auto output_path = resultPath(out_dir, receiver_name);
	removeIfExists(output_path);

	radar::Platform rx_platform("RxPlatform");
	radar::Receiver receiver(&rx_platform, receiver_name, 59, radar::OperationMode::CW_MODE);
	auto timing_owner = makeQuietTiming("oversampled_clk", 25, 77.0);
	receiver.setTiming(timing_owner.timing);
	receiver.setNoiseTemperature(0.0);
	receiver.prepareStreamingData(8);

	radar::Platform tx_platform("TxPlatform");
	radar::Transmitter transmitter(&tx_platform, "TxOversampled", radar::OperationMode::PULSED_MODE, 602);

	std::vector<std::unique_ptr<fers_signal::RadarSignal>> wave_store;
	receiver.addInterferenceToLog(
		makeFixedResponse(&transmitter, wave_store,
						  {ComplexType{1.0, 0.0}, ComplexType{1.0, 0.0}, ComplexType{1.0, 0.0}, ComplexType{1.0, 0.0}},
						  params::rate(), 0.25));

	processing::finalizeStreamingReceiver(&receiver, nullptr, nullptr, out_dir.string());

	{
		HighFive::File file(output_path.string(), HighFive::File::ReadOnly);
		const auto i_data = readDataset(file, "I_data");
		const auto q_data = readDataset(file, "Q_data");

		RealType sampling_rate = 0.0;
		RealType fullscale = 0.0;
		file.getAttribute("sampling_rate").read(sampling_rate);
		file.getAttribute("fullscale").read(fullscale);

		REQUIRE(i_data.size() == 4u);
		REQUIRE(q_data.size() == 4u);
		REQUIRE_THAT(sampling_rate, WithinAbs(params::rate(), 1e-12));
		REQUIRE(fullscale > 0.0);
		REQUIRE(
			std::any_of(i_data.begin(), i_data.end(), [](const RealType value) { return std::abs(value) > 1e-12; }));
	}

	std::filesystem::remove_all(out_dir);
}

TEST_CASE("finalizeStreamingReceiver labels FMCW receiver progress distinctly", "[processing][finalizer][fmcw]")
{
	ParamGuard guard;
	params::setTime(0.0, 1.0);
	params::setRate(4.0);
	params::setOversampleRatio(1);
	params::setAdcBits(0);

	const std::string receiver_name = uniqueName("fmcw_finalize");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("fmcw_finalize_dir");
	std::filesystem::create_directories(out_dir);
	const auto output_path = resultPath(out_dir, receiver_name);
	removeIfExists(output_path);

	radar::Platform rx_platform("RxPlatform");
	radar::Receiver receiver(&rx_platform, receiver_name, 57, radar::OperationMode::FMCW_MODE);
	auto timing_owner = makeQuietTiming("fmcw_clk", 23, 77.0);
	receiver.setTiming(timing_owner.timing);
	receiver.setNoiseTemperature(0.0);
	receiver.prepareStreamingData(4);

	std::vector<ProgressCall> progress_calls;
	auto reporter =
		std::make_shared<core::ProgressReporter>([&progress_calls](const std::string& msg, int current, int total)
												 { progress_calls.push_back({msg, current, total}); });

	processing::finalizeStreamingReceiver(&receiver, nullptr, reporter, out_dir.string());

	REQUIRE(std::filesystem::exists(output_path));
	REQUIRE(progress_calls.size() == 5u);
	REQUIRE(progress_calls.front().message == std::format("Finalizing FMCW Receiver {}", receiver_name));
	REQUIRE(progress_calls.back().message == std::format("Finalized {}", receiver_name));

	const auto metadata = readMetadataJson(output_path);
	REQUIRE(metadata.at("mode") == "fmcw");
	REQUIRE(metadata.at("fmcw_sources").empty());
	REQUIRE_FALSE(metadata.contains("fmcw"));

	std::filesystem::remove_all(out_dir);
}

TEST_CASE("finalizeStreamingReceiver records FMCW source metadata for detached receivers",
		  "[processing][finalizer][fmcw][metadata]")
{
	ParamGuard guard;
	params::setTime(0.0, 0.01);
	params::setRate(1'000.0);
	params::setOversampleRatio(1);
	params::setAdcBits(0);

	const std::string receiver_name = uniqueName("fmcw_detached");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("fmcw_detached_dir");
	std::filesystem::create_directories(out_dir);
	const auto output_path = resultPath(out_dir, receiver_name);
	removeIfExists(output_path);

	radar::Platform rx_platform("RxPlatform");
	radar::Receiver receiver(&rx_platform, receiver_name, 58, radar::OperationMode::FMCW_MODE);
	auto timing_owner = makeQuietTiming("detached_clk", 24, 77.0);
	receiver.setTiming(timing_owner.timing);
	receiver.setNoiseTemperature(0.0);
	receiver.prepareStreamingData(10);

	FmcwTxFixture source_fixture("DetachedTx", 901, 902, 200.0, 0.001, 0.002, 5.0, std::size_t{4});
	const auto source = core::makeActiveSource(&source_fixture.transmitter, 0.0, params::endTime());

	processing::finalizeStreamingReceiver(&receiver, nullptr, nullptr, out_dir.string(), nullptr, {source});

	{
		HighFive::File file(output_path.string(), HighFive::File::ReadOnly);
		unsigned long long source_count = 0;
		RealType scalar_bandwidth = 0.0;
		file.getAttribute("fmcw_source_count").read(source_count);
		file.getAttribute("fmcw_chirp_bandwidth").read(scalar_bandwidth);
		REQUIRE(source_count == 1ULL);
		REQUIRE_THAT(scalar_bandwidth, WithinAbs(200.0, 1e-12));
	}

	const auto metadata = readMetadataJson(output_path);
	REQUIRE(metadata.at("mode") == "fmcw");
	REQUIRE(metadata.contains("fmcw"));
	REQUIRE(metadata.at("fmcw_sources").size() == 1u);

	const auto& source_json = metadata.at("fmcw_sources").front();
	REQUIRE(source_json.at("transmitter_id") == 901);
	REQUIRE(source_json.at("transmitter_name") == "DetachedTx");
	REQUIRE(source_json.at("waveform_id") == 902);
	REQUIRE(source_json.at("waveform_name") == "DetachedTxWave");
	REQUIRE_THAT(source_json.at("chirp_bandwidth").get<RealType>(), WithinAbs(200.0, 1e-12));
	REQUIRE_THAT(source_json.at("chirp_duration").get<RealType>(), WithinAbs(0.001, 1e-12));
	REQUIRE_THAT(source_json.at("chirp_period").get<RealType>(), WithinAbs(0.002, 1e-12));
	REQUIRE(source_json.at("chirp_direction") == "up");
	REQUIRE_THAT(source_json.at("chirp_rate_signed").get<RealType>(), WithinAbs(200'000.0, 1e-12));
	REQUIRE(source_json.at("chirp_count") == 4);
	REQUIRE(source_json.at("segments").size() == 1u);
	REQUIRE_THAT(source_json.at("segments").front().at("first_chirp_start_time").get<RealType>(),
				 WithinAbs(0.0, 1e-12));
	REQUIRE(source_json.at("segments").front().at("emitted_chirp_count") == 4);

	const auto& streaming_segment = metadata.at("streaming_segments").front();
	REQUIRE_THAT(streaming_segment.at("first_chirp_start_time").get<RealType>(), WithinAbs(0.0, 1e-12));
	REQUIRE(streaming_segment.at("emitted_chirp_count") == 4);

	std::filesystem::remove_all(out_dir);
}

TEST_CASE("finalizeStreamingReceiver writes IF-rate FMCW metadata", "[processing][finalizer][fmcw][if]")
{
	ParamGuard guard;
	params::setTime(0.0, 1.0);
	params::setRate(256.0);
	params::setOversampleRatio(1);
	params::setAdcBits(0);

	const std::string receiver_name = uniqueName("fmcw_if_finalize");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("fmcw_if_finalize_dir");
	std::filesystem::create_directories(out_dir);
	const auto output_path = resultPath(out_dir, receiver_name);
	removeIfExists(output_path);

	radar::Platform rx_platform("IfRxPlatform");
	radar::Receiver receiver(&rx_platform, receiver_name, 60, radar::OperationMode::FMCW_MODE);
	auto timing_owner = makeQuietTiming("if_clk", 26, 77.0);
	receiver.setTiming(timing_owner.timing);
	receiver.setNoiseTemperature(0.0);
	receiver.setDechirpMode(radar::Receiver::DechirpMode::Physical);
	receiver.setDechirpReference({.source = radar::Receiver::DechirpReferenceSource::Attached,
								  .name = "",
								  .transmitter_name = "",
								  .waveform_name = ""});
	receiver.setFmcwIfChainRequest(
		{.sample_rate_hz = 64.0, .filter_bandwidth_hz = 16.0, .filter_transition_width_hz = 8.0});
	FmcwTxFixture source_fixture("IfTx", 1001, 1002, 1.0, 1.0, 1.0, 0.0, std::size_t{1});
	const auto source = core::makeActiveSource(&source_fixture.transmitter, params::startTime(), params::endTime());
	receiver.setResolvedDechirpSources({source});

	const fers_signal::FmcwIfResamplerRequest request{.input_sample_rate_hz = 256.0,
													  .output_sample_rate_hz = 64.0,
													  .filter_bandwidth_hz = 16.0,
													  .filter_transition_width_hz = 8.0};
	const auto plan = fers_signal::planFmcwIfResampler(request);
	receiver.initializeFmcwIfResampling(plan);
	const std::vector<ComplexType> high_rate_iq(4096, ComplexType{1.0, 0.0});
	receiver.consumeFmcwIfBlock(high_rate_iq);

	processing::finalizeStreamingReceiver(&receiver, nullptr, nullptr, out_dir.string(), nullptr, {source});

	{
		HighFive::File file(output_path.string(), HighFive::File::ReadOnly);
		RealType sampling_rate = 0.0;
		unsigned long long total_samples = 0;
		bool decimation_enabled = false;
		bool legacy_full_rate = true;
		unsigned long long numerator = 0;
		unsigned long long denominator = 0;
		RealType if_sample_rate = 0.0;
		RealType group_delay = 0.0;
		std::string metadata_json;

		file.getAttribute("sampling_rate").read(sampling_rate);
		file.getAttribute("total_samples").read(total_samples);
		file.getAttribute("fmcw_if_decimation_enabled").read(decimation_enabled);
		file.getAttribute("fmcw_if_legacy_full_rate").read(legacy_full_rate);
		file.getAttribute("fmcw_if_resample_numerator").read(numerator);
		file.getAttribute("fmcw_if_resample_denominator").read(denominator);
		file.getAttribute("fmcw_if_sample_rate").read(if_sample_rate);
		file.getAttribute("fmcw_if_filter_group_delay_seconds").read(group_delay);
		file.getAttribute("fers_metadata_json").read(metadata_json);

		REQUIRE_THAT(sampling_rate, WithinAbs(64.0, 1e-12));
		REQUIRE(total_samples == 64ULL);
		REQUIRE(decimation_enabled);
		REQUIRE_FALSE(legacy_full_rate);
		REQUIRE(numerator == 1ULL);
		REQUIRE(denominator == 4ULL);
		REQUIRE_THAT(if_sample_rate, WithinAbs(64.0, 1e-12));
		REQUIRE_THAT(group_delay, WithinAbs(plan.group_delay_seconds, 1e-12));

		const auto metadata = nlohmann::json::parse(metadata_json);
		REQUIRE(metadata.at("fmcw_if_decimation_enabled").get<bool>());
		REQUIRE_FALSE(metadata.at("fmcw_if_legacy_full_rate").get<bool>());
		REQUIRE_THAT(metadata.at("fmcw_if_sample_rate").get<RealType>(), WithinAbs(64.0, 1e-12));
		REQUIRE(metadata.at("fmcw_if_resample_numerator") == 1);
		REQUIRE(metadata.at("fmcw_if_resample_denominator") == 4);
		REQUIRE(metadata.at("fmcw_if_group_delay_compensated").get<bool>());
		REQUIRE(metadata.at("streaming_segments").front().at("sample_count") == 64);
	}

	std::filesystem::remove_all(out_dir);
}

TEST_CASE("finalizeStreamingReceiver keeps multiple FMCW sources unambiguous",
		  "[processing][finalizer][fmcw][metadata]")
{
	ParamGuard guard;
	params::setTime(0.0, 0.01);
	params::setRate(1'000.0);
	params::setOversampleRatio(1);
	params::setAdcBits(0);

	const std::string receiver_name = uniqueName("fmcw_multi");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("fmcw_multi_dir");
	std::filesystem::create_directories(out_dir);
	const auto output_path = resultPath(out_dir, receiver_name);
	removeIfExists(output_path);

	radar::Platform rx_platform("RxPlatform");
	radar::Receiver receiver(&rx_platform, receiver_name, 59, radar::OperationMode::FMCW_MODE);
	auto timing_owner = makeQuietTiming("multi_clk", 25, 77.0);
	receiver.setTiming(timing_owner.timing);
	receiver.setNoiseTemperature(0.0);
	receiver.prepareStreamingData(10);

	FmcwTxFixture first_source("FirstTx", 911, 912, 200.0, 0.001, 0.002, 0.0, std::size_t{4});
	FmcwTxFixture second_source("SecondTx", 921, 922, 300.0, 0.0015, 0.003, 10.0, std::size_t{2});

	processing::finalizeStreamingReceiver(&receiver, nullptr, nullptr, out_dir.string(), nullptr,
										  {core::makeActiveSource(&first_source.transmitter, 0.0, params::endTime()),
										   core::makeActiveSource(&second_source.transmitter, 0.0, params::endTime())});

	{
		HighFive::File file(output_path.string(), HighFive::File::ReadOnly);
		unsigned long long source_count = 0;
		file.getAttribute("fmcw_source_count").read(source_count);
		REQUIRE(source_count == 2ULL);
	}

	const auto metadata = readMetadataJson(output_path);
	REQUIRE(metadata.at("mode") == "fmcw");
	REQUIRE_FALSE(metadata.contains("fmcw"));
	REQUIRE(metadata.at("fmcw_sources").size() == 2u);
	REQUIRE(metadata.at("fmcw_sources")[0].at("transmitter_id") == 911);
	REQUIRE(metadata.at("fmcw_sources")[1].at("transmitter_id") == 921);
	REQUIRE(metadata.at("fmcw_sources")[0].at("segments").front().at("emitted_chirp_count") == 4);
	REQUIRE(metadata.at("fmcw_sources")[1].at("segments").front().at("emitted_chirp_count") == 2);

	std::filesystem::remove_all(out_dir);
}

TEST_CASE("runPulsedFinalizer writes jittered chunks and emits completion progress", "[processing][finalizer]")
{
	ParamGuard guard;
	params::setRate(8.0);
	params::setOversampleRatio(1);
	params::setAdcBits(0);

	const std::string receiver_name = uniqueName("pulsed_finalize");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("pulsed_finalize_dir");
	std::filesystem::create_directories(out_dir);
	const auto output_path = resultPath(out_dir, receiver_name);
	removeIfExists(output_path);

	radar::Platform rx_platform("RxPlatform");
	setupPlatform(rx_platform, math::Vec3{0.0, 0.0, 0.0});
	radar::Receiver receiver(&rx_platform, receiver_name, 57, radar::OperationMode::PULSED_MODE);
	antenna::Isotropic antenna("iso");
	receiver.setAntenna(&antenna);
	auto timing_owner = makePhaseTiming("pulse_clk", 33, 2.0, PI / 2.0, true);
	receiver.setTiming(timing_owner.timing);
	receiver.setNoiseTemperature(0.0);
	receiver.setWindowProperties(0.5, 1.0, 0.125);

	radar::Platform tx_platform("TxPlatform");
	radar::Transmitter transmitter(&tx_platform, "TxA", radar::OperationMode::PULSED_MODE, 701);
	std::vector<std::unique_ptr<fers_signal::RadarSignal>> wave_store;

	core::RenderingJob first_job{};
	first_job.ideal_start_time = 0.0;
	first_job.duration = 0.5;
	first_job.responses.push_back(
		makeFixedResponse(&transmitter, wave_store, {ComplexType{1.0, 0.0}, ComplexType{1.0, 0.0}}, 8.0, 0.125));

	core::RenderingJob second_job{};
	second_job.ideal_start_time = 1.0;
	second_job.duration = 0.5;
	second_job.responses.push_back(
		makeFixedResponse(&transmitter, wave_store, {ComplexType{1.0, 0.0}, ComplexType{1.0, 0.0}}, 8.0, 1.125));

	std::vector<std::unique_ptr<radar::Target>> targets;
	std::vector<ProgressCall> progress_calls;
	auto reporter =
		std::make_shared<core::ProgressReporter>([&progress_calls](const std::string& msg, int current, int total)
												 { progress_calls.push_back({msg, current, total}); });

	std::jthread worker(processing::runPulsedFinalizer, &receiver, &targets, reporter, out_dir.string(), nullptr);
	receiver.enqueueFinalizerJob(std::move(first_job));
	receiver.enqueueFinalizerJob(std::move(second_job));
	core::RenderingJob shutdown_job{};
	shutdown_job.duration = -1.0;
	receiver.enqueueFinalizerJob(std::move(shutdown_job));
	worker.join();

	{
		HighFive::File file(output_path.string(), HighFive::File::ReadOnly);
		const auto i_chunk_0 = readDataset(file, "chunk_000000_I");
		const auto q_chunk_0 = readDataset(file, "chunk_000000_Q");
		const auto i_chunk_1 = readDataset(file, "chunk_000001_I");
		const auto q_chunk_1 = readDataset(file, "chunk_000001_Q");

		RealType time_attr_0 = 0.0;
		RealType time_attr_1 = 0.0;
		file.getDataSet("chunk_000000_I").getAttribute("time").read(time_attr_0);
		file.getDataSet("chunk_000001_I").getAttribute("time").read(time_attr_1);

		REQUIRE(i_chunk_0.size() == 4u);
		REQUIRE(q_chunk_0.size() == 4u);
		REQUIRE(i_chunk_1.size() == 4u);
		REQUIRE(q_chunk_1.size() == 4u);

		for (const auto& sample : i_chunk_0)
		{
			REQUIRE_THAT(sample, WithinAbs(0.0, 1e-12));
		}
		for (const auto& sample : i_chunk_1)
		{
			REQUIRE_THAT(sample, WithinAbs(0.0, 1e-12));
		}
		REQUIRE_THAT(q_chunk_0[0], WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(q_chunk_0[1], WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(q_chunk_0[2], WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(q_chunk_0[3], WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(q_chunk_1[0], WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(q_chunk_1[1], WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(q_chunk_1[2], WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(q_chunk_1[3], WithinAbs(0.0, 1e-12));
		REQUIRE_THAT(time_attr_0, WithinAbs(0.125, 1e-12));
		REQUIRE_THAT(time_attr_1, WithinAbs(1.125, 1e-12));
	}

	REQUIRE(std::ranges::any_of(
		progress_calls, [&receiver_name](const ProgressCall& call)
		{ return call.message == "Finished Exporting " + receiver_name && call.current == 100 && call.total == 100; }));

	std::filesystem::remove_all(out_dir);
}

// TODO: The null-check branches after Timing::clone in finalizer.cpp are not reachable
// through the public Timing API because clone either returns a valid object or throws.
