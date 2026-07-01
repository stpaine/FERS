#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <highfive/highfive.hpp>
#include <memory>
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
#include "serial/hdf5_output_sink.h"
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
		ParamGuard(const ParamGuard&) = delete;
		ParamGuard& operator=(const ParamGuard&) = delete;
		ParamGuard(ParamGuard&&) = delete;
		ParamGuard& operator=(ParamGuard&&) = delete;
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

	struct SfcwTxFixture
	{
		radar::Platform platform;
		std::unique_ptr<fers_signal::RadarSignal> wave;
		radar::Transmitter transmitter;

		SfcwTxFixture(const std::string& name, SimId tx_id, SimId waveform_id, RealType start_frequency_offset,
					  RealType step_size, std::size_t step_count, RealType dwell_time, RealType step_period,
					  std::optional<std::size_t> sweep_count) :
			platform(name + "Platform"),
			wave(std::make_unique<fers_signal::RadarSignal>(
				name + "Wave", 1.0, 10.0e9, dwell_time,
				std::make_unique<fers_signal::SteppedFrequencySignal>(start_frequency_offset, step_size, step_count,
																	  dwell_time, step_period, sweep_count),
				waveform_id)),
			transmitter(&platform, name, radar::OperationMode::SFCW_MODE, tx_id)
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

	std::unique_ptr<core::ReceiverOutputSink> makeTestHdf5Sink(const std::filesystem::path& out_dir)
	{
		auto sink = serial::makeHdf5OutputSink(out_dir.string());
		sink->initializeRun(core::OutputConfig{}, "sim");
		return sink;
	}

	struct ProgressCall
	{
		std::string message;
		int current;
		int total;
	};

	void requireJitteredPulsedChunks(const std::filesystem::path& output_path)
	{
		HighFive::File const file(output_path.string(), HighFive::File::ReadOnly);
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

	[[nodiscard]] bool hasCompletionProgress(const std::vector<ProgressCall>& progress_calls,
											 const std::string& receiver_name)
	{
		return std::ranges::any_of(progress_calls,
								   [&receiver_name](const ProgressCall& call)
								   {
									   return call.message == "Finished Exporting " + receiver_name &&
										   call.current == 100 && call.total == 100;
								   });
	}

	struct CapturedSinkBlock
	{
		core::ReceiverStreamDescriptor stream;
		RealType first_sample_time = 0.0;
		RealType sample_rate = 0.0;
		std::uint64_t sample_start = 0;
		std::vector<ComplexType> samples;
	};

	class CapturingOutputSink final : public core::ReceiverOutputSink
	{
	public:
		void initializeRun(const core::OutputConfig&, std::string) override {}

		std::uint32_t registerStream(const core::ReceiverStreamDescriptor& stream) override
		{
			registered_streams.push_back(stream);
			return next_stream_id++;
		}

		void openStream(const std::uint32_t stream_id, const RealType first_sample_time) override
		{
			opened_streams.push_back(stream_id);
			open_times.push_back(first_sample_time);
		}

		void submitBlock(const core::ReceiverSampleBlock& block) override
		{
			blocks.push_back(
				CapturedSinkBlock{.stream = block.stream,
								  .first_sample_time = block.first_sample_time,
								  .sample_rate = block.sample_rate,
								  .sample_start = block.sample_start,
								  .samples = std::vector<ComplexType>(block.samples.begin(), block.samples.end())});
		}

		void emitContextHeartbeat(RealType) override {}

		void closeStream(const std::uint32_t stream_id) override { closed_streams.push_back(stream_id); }

		core::OutputStats finalize() override { return {}; }

		std::uint32_t next_stream_id = 100;
		std::vector<core::ReceiverStreamDescriptor> registered_streams;
		std::vector<std::uint32_t> opened_streams;
		std::vector<RealType> open_times;
		std::vector<CapturedSinkBlock> blocks;
		std::vector<std::uint32_t> closed_streams;
	};
}

TEST_CASE("Hdf5OutputSink writes streaming blocks through the receiver output contract",
		  "[processing][finalizer][hdf5]")
{
	ParamGuard const guard;
	params::setAdcBits(0);

	const std::string receiver_name = uniqueName("hdf5_sink");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("hdf5_sink_dir");
	std::filesystem::create_directories(out_dir);
	const auto output_path = resultPath(out_dir, receiver_name);
	removeIfExists(output_path);

	serial::Hdf5OutputSink sink(out_dir.string());
	sink.initializeRun(core::OutputConfig{}, "sim");

	const core::ReceiverStreamDescriptor stream{.receiver_id = 42,
												.receiver_name = receiver_name,
												.mode = "cw",
												.sample_rate = 2.0,
												.reference_frequency = 1.0e9,
												.coordinate = {},
												.initial_platform_state = {},
												.fmcw = {}};
	const std::vector<ComplexType> samples{ComplexType{1.0, 0.0}, ComplexType{0.5, -0.5}};
	const auto stream_id = sink.registerStream(stream);
	sink.openStream(stream_id, 0.0);
	sink.submitBlock(core::ReceiverSampleBlock{
		.stream = stream, .first_sample_time = 0.0, .sample_rate = 2.0, .samples = samples, .sample_start = 0});
	sink.closeStream(stream_id);
	sink.finalize();

	{
		HighFive::File const file(output_path.string(), HighFive::File::ReadOnly);
		const auto i_data = readDataset(file, "I_data");
		const auto q_data = readDataset(file, "Q_data");

		REQUIRE(i_data.size() == 2u);
		REQUIRE(q_data.size() == 2u);
		REQUIRE_THAT(i_data[0], WithinAbs(1.0, 1e-12));
		REQUIRE_THAT(i_data[1], WithinAbs(0.5, 1e-12));
		REQUIRE_THAT(q_data[1], WithinAbs(-0.5, 1e-12));
	}

	std::filesystem::remove_all(out_dir);
}

TEST_CASE("buildReceiverSampleBlock captures receiver identity, timing, and sample basis",
		  "[processing][finalizer][vita49]")
{
	ParamGuard const guard;
	params::setAdcBits(12);
	params::setOrigin(-33.5, 18.25, 123.0);
	params::setCoordinateSystem(params::CoordinateFrame::UTM, 34, false);

	radar::Platform platform("RxPlatform");
	setupPlatform(platform, math::Vec3{100.0, 200.0, 300.0});
	radar::Receiver receiver(&platform, "RxDescriptor", 55, radar::OperationMode::FMCW_MODE, 4242);
	auto timing_owner = makeQuietTiming("descriptor_clk", 11, 10.5e9);
	receiver.setTiming(timing_owner.timing);
	receiver.setDechirpMode(radar::Receiver::DechirpMode::Physical);
	FmcwTxFixture tx("DescriptorTx", 6001, 7001, 40.0e6, 2.0e-3, 3.0e-3, -1.0e6, 7);
	setupPlatform(tx.platform, math::Vec3{0.0, 0.0, 0.0});
	std::vector<core::ActiveStreamingSource> streaming_sources{core::makeActiveSource(&tx.transmitter, 0.0, 1.0)};

	const std::vector<ComplexType> samples = {ComplexType{1.0, 0.0}, ComplexType{0.0, 1.0}};
	const auto block =
		processing::buildReceiverSampleBlock(&receiver, 1.25, 2.0e6, samples, 17, streaming_sources, nullptr);

	REQUIRE(block.stream.receiver_id == 4242);
	REQUIRE(block.stream.receiver_name == "RxDescriptor");
	REQUIRE(block.stream.mode == "fmcw");
	REQUIRE(block.stream.dechirped);
	REQUIRE(block.stream.adc_bits == 12u);
	REQUIRE(block.stream.coordinate.frame == "UTM");
	REQUIRE_THAT(block.stream.coordinate.origin_latitude, WithinAbs(-33.5, 1e-12));
	REQUIRE(block.stream.coordinate.utm_zone == 34);
	REQUIRE_FALSE(block.stream.coordinate.utm_north_hemisphere);
	REQUIRE(block.stream.initial_platform_state.platform_id == platform.getId());
	REQUIRE(block.stream.initial_platform_state.platform_name == "RxPlatform");
	REQUIRE_THAT(block.stream.initial_platform_state.position_x, WithinAbs(100.0, 1e-12));
	REQUIRE_THAT(block.stream.initial_platform_state.position_y, WithinAbs(200.0, 1e-12));
	REQUIRE_THAT(block.stream.initial_platform_state.position_z, WithinAbs(300.0, 1e-12));
	REQUIRE(block.stream.fmcw.present);
	REQUIRE(block.stream.fmcw.waveform_shape == "linear");
	REQUIRE(block.stream.fmcw.sweep_direction == "up");
	REQUIRE_THAT(block.stream.fmcw.chirp_bandwidth, WithinAbs(40.0e6, 1e-6));
	REQUIRE_THAT(block.stream.fmcw.chirp_duration, WithinAbs(2.0e-3, 1e-12));
	REQUIRE_THAT(block.stream.fmcw.chirp_period, WithinAbs(3.0e-3, 1e-12));
	REQUIRE_THAT(block.stream.fmcw.chirp_rate, WithinAbs(20.0e9, 1e-3));
	REQUIRE(block.stream.fmcw.chirp_count.has_value());
	REQUIRE(block.stream.fmcw.chirp_count.value_or(0u) == 7u);
	REQUIRE_THAT(block.stream.sample_rate, WithinAbs(2.0e6, 1e-12));
	REQUIRE_THAT(block.stream.reference_frequency, WithinAbs(10.5e9, 1e-6));
	REQUIRE_THAT(block.first_sample_time, WithinAbs(1.25, 1e-12));
	REQUIRE(block.sample_start == 17u);
	REQUIRE(block.samples.size() == samples.size());
}

TEST_CASE("buildReceiverStreamDescriptor keeps CW metadata isolated from active FMCW sources",
		  "[processing][finalizer][vita49]")
{
	ParamGuard const guard;
	params::setTime(0.0, 1.0);
	params::setRate(1'000.0);
	params::setOversampleRatio(1);

	radar::Platform rx_platform("CwRxPlatform");
	radar::Receiver receiver(&rx_platform, "CwRx", 61, radar::OperationMode::CW_MODE);
	auto timing_owner = makeQuietTiming("cw_metadata_clk", 27, 5.0e9);
	receiver.setTiming(timing_owner.timing);

	FmcwTxFixture const source_fixture("MixedFmcwTx", 1101, 1102, 200.0, 0.001, 0.002, 0.0, std::size_t{4});
	const std::vector sources = {core::makeActiveSource(&source_fixture.transmitter, 0.0, params::endTime())};

	const auto descriptor = processing::buildReceiverStreamDescriptor(&receiver, params::rate(), sources);

	REQUIRE(descriptor.mode == "cw");
	REQUIRE(descriptor.cw.present);
	REQUIRE_FALSE(descriptor.fmcw.present);
	REQUIRE_THAT(descriptor.cw.carrier_frequency, WithinAbs(5.0e9, 1e-6));
}

TEST_CASE("buildReceiverStreamDescriptor records SFCW waveform metadata", "[processing][finalizer][sfcw][vita49]")
{
	ParamGuard const guard;
	params::setTime(0.0, 0.01);
	params::setRate(10'000.0);
	params::setOversampleRatio(1);

	radar::Platform rx_platform("SfcwRxPlatform");
	radar::Receiver receiver(&rx_platform, "SfcwRx", 62, radar::OperationMode::SFCW_MODE);
	auto timing_owner = makeQuietTiming("sfcw_metadata_clk", 28, 9.6e9);
	receiver.setTiming(timing_owner.timing);

	SfcwTxFixture const source_fixture("DescriptorSfcwTx", 1201, 1202, -1.0e6, 2.0e5, 8, 1.0e-4, 2.0e-4,
									   std::size_t{2});
	const std::vector sources = {core::makeActiveSource(&source_fixture.transmitter, 0.0, params::endTime())};

	const auto descriptor = processing::buildReceiverStreamDescriptor(&receiver, params::rate(), sources);

	REQUIRE(descriptor.mode == "sfcw");
	REQUIRE(descriptor.sfcw.present);
	REQUIRE_FALSE(descriptor.fmcw.present);
	REQUIRE_FALSE(descriptor.cw.present);
	REQUIRE(descriptor.sfcw.waveform_id == 1202);
	REQUIRE(descriptor.sfcw.step_count == 8u);
	REQUIRE_THAT(descriptor.sfcw.start_frequency_offset, WithinAbs(-1.0e6, 1e-9));
	REQUIRE_THAT(descriptor.sfcw.step_size, WithinAbs(2.0e5, 1e-9));
	REQUIRE_THAT(descriptor.sfcw.dwell_time, WithinAbs(1.0e-4, 1e-12));
	REQUIRE_THAT(descriptor.sfcw.step_period, WithinAbs(2.0e-4, 1e-12));
	REQUIRE_THAT(descriptor.sfcw.sweep_period, WithinAbs(1.6e-3, 1e-12));
	REQUIRE(descriptor.sfcw.sweep_count.has_value());
	REQUIRE(descriptor.sfcw.sweep_count.value_or(0u) == 2u);
	REQUIRE_THAT(descriptor.sfcw.first_frequency, WithinAbs(9.999e9, 1e-3));
	REQUIRE_THAT(descriptor.sfcw.last_frequency, WithinAbs(10.0004e9, 1e-3));
	REQUIRE_THAT(descriptor.sfcw.effective_bandwidth, WithinAbs(1.6e6, 1e-6));
	REQUIRE_THAT(descriptor.sfcw.range_resolution, WithinAbs(params::c() / (2.0 * 1.6e6), 1e-9));
	REQUIRE_THAT(descriptor.sfcw.unambiguous_range, WithinAbs(params::c() / (2.0 * 2.0e5), 1e-9));
}

TEST_CASE("buildStreamingOutputMetadata records FMCW source metadata for detached receivers",
		  "[processing][finalizer][fmcw][metadata]")
{
	ParamGuard const guard;
	params::setTime(0.0, 0.01);
	params::setRate(1'000.0);
	params::setOversampleRatio(1);
	params::setAdcBits(0);

	const std::string receiver_name = uniqueName("fmcw_detached");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("fmcw_detached_dir");
	const auto output_path = resultPath(out_dir, receiver_name);

	radar::Platform rx_platform("RxPlatform");
	radar::Receiver receiver(&rx_platform, receiver_name, 58, radar::OperationMode::FMCW_MODE);
	auto timing_owner = makeQuietTiming("detached_clk", 24, 77.0);
	receiver.setTiming(timing_owner.timing);
	receiver.setNoiseTemperature(0.0);

	FmcwTxFixture const source_fixture("DetachedTx", 901, 902, 200.0, 0.001, 0.002, 5.0, std::size_t{4});
	const auto source = core::makeActiveSource(&source_fixture.transmitter, 0.0, params::endTime());

	const auto metadata =
		processing::buildStreamingOutputMetadata(&receiver, output_path.string(), 10, {source}, params::rate());

	REQUIRE(metadata.mode == "fmcw");
	REQUIRE(metadata.fmcw.has_value());
	REQUIRE(metadata.fmcw_sources.size() == 1u);

	const auto& source_metadata = metadata.fmcw_sources.front();
	REQUIRE(source_metadata.transmitter_id == 901);
	REQUIRE(source_metadata.transmitter_name == "DetachedTx");
	REQUIRE(source_metadata.waveform_id == 902);
	REQUIRE(source_metadata.waveform_name == "DetachedTxWave");
	REQUIRE_THAT(source_metadata.waveform.chirp_bandwidth, WithinAbs(200.0, 1e-12));
	REQUIRE_THAT(source_metadata.waveform.chirp_duration, WithinAbs(0.001, 1e-12));
	REQUIRE_THAT(source_metadata.waveform.chirp_period, WithinAbs(0.002, 1e-12));
	REQUIRE(source_metadata.waveform.chirp_direction == "up");
	REQUIRE_THAT(source_metadata.waveform.chirp_rate_signed, WithinAbs(200'000.0, 1e-12));
	REQUIRE(source_metadata.waveform.chirp_count == std::optional<std::uint64_t>{4});
	REQUIRE(source_metadata.segments.size() == 1u);
	REQUIRE(source_metadata.segments.front().first_chirp_start_time.has_value());
	REQUIRE_THAT(source_metadata.segments.front().first_chirp_start_time.value_or(1.0), WithinAbs(0.0, 1e-12));
	REQUIRE(source_metadata.segments.front().emitted_chirp_count == std::optional<std::uint64_t>{4});

	const auto& streaming_segment = metadata.streaming_segments.front();
	REQUIRE(streaming_segment.first_chirp_start_time.has_value());
	REQUIRE_THAT(streaming_segment.first_chirp_start_time.value_or(1.0), WithinAbs(0.0, 1e-12));
	REQUIRE(streaming_segment.emitted_chirp_count == std::optional<std::uint64_t>{4});
}

TEST_CASE("buildStreamingOutputMetadata records SFCW source metadata for detached receivers",
		  "[processing][finalizer][sfcw][metadata]")
{
	ParamGuard const guard;
	params::setTime(0.0, 0.01);
	params::setRate(10'000.0);
	params::setOversampleRatio(1);
	params::setAdcBits(0);

	const std::string receiver_name = uniqueName("sfcw_detached");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("sfcw_detached_dir");
	const auto output_path = resultPath(out_dir, receiver_name);

	radar::Platform rx_platform("SfcwRxPlatform");
	radar::Receiver receiver(&rx_platform, receiver_name, 58, radar::OperationMode::SFCW_MODE);
	auto timing_owner = makeQuietTiming("sfcw_detached_clk", 24, 77.0);
	receiver.setTiming(timing_owner.timing);
	receiver.setNoiseTemperature(0.0);

	SfcwTxFixture const source_fixture("DetachedSfcwTx", 1301, 1302, 0.0, 1.0e5, 4, 1.0e-4, 2.0e-4, std::size_t{2});
	const auto source = core::makeActiveSource(&source_fixture.transmitter, 0.0, params::endTime());

	const auto metadata =
		processing::buildStreamingOutputMetadata(&receiver, output_path.string(), 20, {source}, params::rate());

	REQUIRE(metadata.mode == "sfcw");
	REQUIRE(metadata.sfcw.has_value());
	REQUIRE(metadata.sfcw_sources.size() == 1u);
	REQUIRE_FALSE(metadata.fmcw.has_value());

	const auto& source_metadata = metadata.sfcw_sources.front();
	REQUIRE(source_metadata.transmitter_id == 1301);
	REQUIRE(source_metadata.transmitter_name == "DetachedSfcwTx");
	REQUIRE(source_metadata.waveform_id == 1302);
	REQUIRE(source_metadata.waveform_name == "DetachedSfcwTxWave");
	REQUIRE_THAT(source_metadata.waveform.step_size, WithinAbs(1.0e5, 1e-9));
	REQUIRE(source_metadata.waveform.step_count == 4u);
	REQUIRE_THAT(source_metadata.waveform.dwell_time, WithinAbs(1.0e-4, 1e-12));
	REQUIRE_THAT(source_metadata.waveform.step_period, WithinAbs(2.0e-4, 1e-12));
	REQUIRE(source_metadata.waveform.sweep_count == std::optional<std::uint64_t>{2});
	REQUIRE_THAT(source_metadata.waveform.effective_bandwidth, WithinAbs(4.0e5, 1e-9));
	REQUIRE_THAT(source_metadata.waveform.range_resolution, WithinAbs(params::c() / (2.0 * 4.0e5), 1e-9));
	REQUIRE_THAT(source_metadata.waveform.unambiguous_range, WithinAbs(params::c() / (2.0 * 1.0e5), 1e-9));
	REQUIRE(source_metadata.segments.size() == 1u);
	REQUIRE(source_metadata.segments.front().first_step_start_time.has_value());
	REQUIRE_THAT(source_metadata.segments.front().first_step_start_time.value_or(1.0), WithinAbs(0.0, 1e-12));
	REQUIRE(source_metadata.segments.front().emitted_step_count == std::optional<std::uint64_t>{8});

	const auto& streaming_segment = metadata.streaming_segments.front();
	REQUIRE(streaming_segment.first_sfcw_step_start_time.has_value());
	REQUIRE_THAT(streaming_segment.first_sfcw_step_start_time.value_or(1.0), WithinAbs(0.0, 1e-12));
	REQUIRE(streaming_segment.emitted_sfcw_step_count == std::optional<std::uint64_t>{8});
}

TEST_CASE("buildStreamingOutputMetadata writes IF-rate FMCW metadata", "[processing][finalizer][fmcw][if]")
{
	ParamGuard const guard;
	params::setTime(0.0, 1.0);
	params::setRate(256.0);
	params::setOversampleRatio(1);
	params::setAdcBits(0);

	const std::string receiver_name = uniqueName("fmcw_if_finalize");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("fmcw_if_finalize_dir");
	const auto output_path = resultPath(out_dir, receiver_name);

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
	FmcwTxFixture const source_fixture("IfTx", 1001, 1002, 1.0, 1.0, 1.0, 0.0, std::size_t{1});
	const auto source = core::makeActiveSource(&source_fixture.transmitter, params::startTime(), params::endTime());
	receiver.setResolvedDechirpSources({source});

	const fers_signal::FmcwIfResamplerRequest request{.input_sample_rate_hz = 256.0,
													  .output_sample_rate_hz = 64.0,
													  .filter_bandwidth_hz = 16.0,
													  .filter_transition_width_hz = 8.0};
	const auto plan = fers_signal::planFmcwIfResampler(request);
	receiver.initializeFmcwIfResampling(plan);
	const auto metadata = processing::buildStreamingOutputMetadata(&receiver, output_path.string(), 64, {source}, 64.0);

	REQUIRE_THAT(metadata.sampling_rate, WithinAbs(64.0, 1e-12));
	REQUIRE(metadata.total_samples == 64u);
	REQUIRE(metadata.fmcw_if_decimation_enabled);
	REQUIRE_FALSE(metadata.fmcw_if_legacy_full_rate);
	REQUIRE(metadata.fmcw_if_resample_numerator == 1u);
	REQUIRE(metadata.fmcw_if_resample_denominator == 4u);
	REQUIRE(metadata.fmcw_if_sample_rate == std::optional<RealType>{64.0});
	REQUIRE(metadata.fmcw_if_filter_group_delay_seconds.has_value());
	REQUIRE_THAT(metadata.fmcw_if_filter_group_delay_seconds.value_or(0.0), WithinAbs(plan.group_delay_seconds, 1e-12));
	REQUIRE(metadata.fmcw_if_group_delay_compensated);
	REQUIRE(metadata.streaming_segments.front().sample_count == 64u);
}

TEST_CASE("buildStreamingOutputMetadata keeps multiple FMCW sources unambiguous",
		  "[processing][finalizer][fmcw][metadata]")
{
	ParamGuard const guard;
	params::setTime(0.0, 0.01);
	params::setRate(1'000.0);
	params::setOversampleRatio(1);
	params::setAdcBits(0);

	const std::string receiver_name = uniqueName("fmcw_multi");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("fmcw_multi_dir");
	const auto output_path = resultPath(out_dir, receiver_name);

	radar::Platform rx_platform("RxPlatform");
	radar::Receiver receiver(&rx_platform, receiver_name, 59, radar::OperationMode::FMCW_MODE);
	auto timing_owner = makeQuietTiming("multi_clk", 25, 77.0);
	receiver.setTiming(timing_owner.timing);
	receiver.setNoiseTemperature(0.0);

	FmcwTxFixture const first_source("FirstTx", 911, 912, 200.0, 0.001, 0.002, 0.0, std::size_t{4});
	FmcwTxFixture const second_source("SecondTx", 921, 922, 300.0, 0.0015, 0.003, 10.0, std::size_t{2});

	const std::vector sources = {core::makeActiveSource(&first_source.transmitter, 0.0, params::endTime()),
								 core::makeActiveSource(&second_source.transmitter, 0.0, params::endTime())};
	const auto metadata =
		processing::buildStreamingOutputMetadata(&receiver, output_path.string(), 10, sources, params::rate());

	REQUIRE(metadata.mode == "fmcw");
	REQUIRE_FALSE(metadata.fmcw.has_value());
	REQUIRE(metadata.fmcw_sources.size() == 2u);
	REQUIRE(metadata.fmcw_sources[0].transmitter_id == 911);
	REQUIRE(metadata.fmcw_sources[1].transmitter_id == 921);
	REQUIRE(metadata.fmcw_sources[0].segments.front().emitted_chirp_count == std::optional<std::uint64_t>{4});
	REQUIRE(metadata.fmcw_sources[1].segments.front().emitted_chirp_count == std::optional<std::uint64_t>{2});
}

TEST_CASE("runPulsedFinalizer writes jittered chunks and emits completion progress", "[processing][finalizer]")
{
	ParamGuard const guard;
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
	radar::Transmitter const transmitter(&tx_platform, "TxA", radar::OperationMode::PULSED_MODE, 701);
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

	auto hdf5_sink = makeTestHdf5Sink(out_dir);
	std::jthread worker(processing::runPulsedFinalizer, &receiver, &targets, reporter, out_dir.string(), nullptr,
						hdf5_sink.get());
	receiver.enqueueFinalizerJob(std::move(first_job));
	receiver.enqueueFinalizerJob(std::move(second_job));
	core::RenderingJob shutdown_job{};
	shutdown_job.duration = -1.0;
	receiver.enqueueFinalizerJob(std::move(shutdown_job));
	worker.join();
	hdf5_sink->finalize();

	requireJitteredPulsedChunks(output_path);
	REQUIRE(hasCompletionProgress(progress_calls, receiver_name));

	std::filesystem::remove_all(out_dir);
}

TEST_CASE("runPulsedFinalizer routes completed acquisition windows to an output sink without HDF5",
		  "[processing][finalizer][vita49]")
{
	ParamGuard const guard;
	params::setTime(0.0, 1.0);
	params::setRate(4.0);
	params::setOversampleRatio(1);
	params::setAdcBits(0);

	const std::string receiver_name = uniqueName("pulsed_sink");
	const auto out_dir = std::filesystem::temp_directory_path() / uniqueName("pulsed_sink_dir");
	std::filesystem::create_directories(out_dir);
	const auto output_path = resultPath(out_dir, receiver_name);
	removeIfExists(output_path);

	radar::Platform rx_platform("RxPlatform");
	setupPlatform(rx_platform, math::Vec3{0.0, 0.0, 0.0});
	radar::Receiver receiver(&rx_platform, receiver_name, 57, radar::OperationMode::PULSED_MODE);
	antenna::Isotropic antenna("iso");
	receiver.setAntenna(&antenna);
	auto timing_owner = makeQuietTiming("pulsed_sink_clk", 33, 2.0);
	receiver.setTiming(timing_owner.timing);
	receiver.setNoiseTemperature(0.0);
	receiver.setWindowProperties(0.5, 1.0, 0.0);

	radar::Platform tx_platform("TxPlatform");
	radar::Transmitter const transmitter(&tx_platform, "TxA", radar::OperationMode::PULSED_MODE, 701);
	std::vector<std::unique_ptr<fers_signal::RadarSignal>> wave_store;

	core::RenderingJob job{};
	job.ideal_start_time = 0.0;
	job.duration = 0.5;
	job.responses.push_back(
		makeFixedResponse(&transmitter, wave_store, {ComplexType{1.0, 0.0}, ComplexType{2.0, 0.0}}, 4.0, 0.0));

	std::vector<std::unique_ptr<radar::Target>> targets;
	CapturingOutputSink sink;

	std::jthread worker(processing::runPulsedFinalizer, &receiver, &targets, nullptr, out_dir.string(), nullptr, &sink);
	receiver.enqueueFinalizerJob(std::move(job));
	core::RenderingJob shutdown_job{};
	shutdown_job.duration = -1.0;
	receiver.enqueueFinalizerJob(std::move(shutdown_job));
	worker.join();

	REQUIRE_FALSE(std::filesystem::exists(output_path));
	REQUIRE(sink.registered_streams.size() == 1u);
	REQUIRE(sink.opened_streams.size() == 1u);
	REQUIRE(sink.closed_streams.size() == 1u);
	REQUIRE(sink.blocks.size() == 1u);

	const auto& block = sink.blocks.front();
	REQUIRE(block.stream.mode == "pulsed");
	REQUIRE_THAT(block.first_sample_time, WithinAbs(0.0, 1e-12));
	REQUIRE_THAT(block.sample_rate, WithinAbs(params::rate(), 1e-12));
	REQUIRE(block.sample_start == 0u);
	REQUIRE(block.samples.size() == 2u);
	REQUIRE_THAT(block.samples[0].real(), WithinAbs(1.0, 1e-12));
	REQUIRE_THAT(block.samples[1].real(), WithinAbs(2.0, 1e-12));

	std::filesystem::remove_all(out_dir);
}

// TODO: The null-check branches after Timing::clone in finalizer.cpp are not reachable
// through the public Timing API because clone either returns a valid object or throws.
