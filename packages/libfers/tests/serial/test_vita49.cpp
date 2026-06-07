// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).

#include <algorithm>
#include <atomic>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <future>
#include <iterator>
#include <memory>
#include <mutex>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "core/parameters.h"
#include "serial/vita49/paced_sender.h"
#include "serial/vita49/stream_registry.h"
#include "serial/vita49/udp_sender.h"
#include "serial/vita49/vita49_context_builder.h"
#include "serial/vita49/vita49_output_sink.h"
#include "serial/vita49/vita49_packetizer.h"
#include "serial/vita49/vita49_serializer.h"

#ifndef _WIN32
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>
#endif

using Catch::Matchers::ContainsSubstring;

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

	[[nodiscard]] std::uint16_t readU16(const std::vector<std::uint8_t>& bytes, const std::size_t offset)
	{
		return static_cast<std::uint16_t>((static_cast<std::uint16_t>(bytes.at(offset)) << 8u) |
										  static_cast<std::uint16_t>(bytes.at(offset + 1u)));
	}

	[[nodiscard]] std::uint32_t readU32(const std::vector<std::uint8_t>& bytes, const std::size_t offset)
	{
		return (static_cast<std::uint32_t>(bytes.at(offset)) << 24u) |
			(static_cast<std::uint32_t>(bytes.at(offset + 1u)) << 16u) |
			(static_cast<std::uint32_t>(bytes.at(offset + 2u)) << 8u) |
			static_cast<std::uint32_t>(bytes.at(offset + 3u));
	}

	[[nodiscard]] std::uint64_t readU64(const std::vector<std::uint8_t>& bytes, const std::size_t offset)
	{
		return (static_cast<std::uint64_t>(readU32(bytes, offset)) << 32u) | readU32(bytes, offset + 4u);
	}

	[[nodiscard]] double readF64(const std::vector<std::uint8_t>& bytes, const std::size_t offset)
	{
		return std::bit_cast<double>(readU64(bytes, offset));
	}

#ifndef _WIN32
	[[nodiscard]] sockaddr* asSockaddr(sockaddr_in& address) noexcept
	{
		return static_cast<sockaddr*>(static_cast<void*>(&address));
	}
#endif

	[[nodiscard]] std::string readAsciiMetadata(const std::vector<std::uint8_t>& bytes, std::size_t offset)
	{
		std::string value;
		while (offset < bytes.size() && bytes.at(offset) != 0)
		{
			++offset;
			value.push_back(static_cast<char>(bytes.at(offset - 1u)));
		}
		return value;
	}

	[[nodiscard]] core::ReceiverStreamDescriptor fmcwContextStreamDescriptor()
	{
		return core::ReceiverStreamDescriptor{.receiver_id = 0x1122334455667788ull,
											  .receiver_name = "rx1",
											  .mode = "fmcw",
											  .sample_rate = 2.5e6,
											  .reference_frequency = 9.6e9,
											  .if_offset = -1.25e6,
											  .bandwidth = 5.0e6,
											  .dechirped = true,
											  .if_resampled = false,
											  .adc_bits = 14,
											  .coordinate = {.frame = "UTM",
															 .origin_latitude = -33.9,
															 .origin_longitude = 18.4,
															 .origin_altitude = 111.0,
															 .utm_zone = 34,
															 .utm_north_hemisphere = false},
											  .initial_platform_state = {.platform_id = 0x99,
																		 .platform_name = "rx-platform",
																		 .position_x = 1.0,
																		 .position_y = 2.0,
																		 .position_z = 3.0,
																		 .velocity_x = 4.0,
																		 .velocity_y = 5.0,
																		 .velocity_z = 6.0,
																		 .azimuth = 0.25,
																		 .elevation = -0.5},
											  .fmcw = {.present = true,
													   .waveform_shape = "linear",
													   .chirp_bandwidth = 20.0e6,
													   .chirp_duration = 1.0e-3,
													   .chirp_period = 2.0e-3,
													   .chirp_rate = 20.0e9,
													   .chirp_rate_signed = -20.0e9,
													   .sweep_direction = "down",
													   .start_frequency_offset = 1.5e6,
													   .chirp_count = 64,
													   .dechirp_mode = "physical",
													   .dechirp_reference_source = "transmitter",
													   .dechirp_reference_transmitter_id = 0x1234,
													   .dechirp_reference_transmitter_name = "tx-lo",
													   .dechirp_reference_waveform_id = 0,
													   .dechirp_reference_waveform_name = ""}};
	}

	[[nodiscard]] std::vector<std::uint8_t> serializedFmcwContextPacket()
	{
		using namespace serial::vita49;
		const auto context = Vita49ContextBuilder::build(ContextBuildRequest{.stream = fmcwContextStreamDescriptor(),
																			 .stream_id = 0x55667788u,
																			 .simulation_name = "sim",
																			 .adc_fullscale = 2.0,
																			 .timestamp = Timestamp{100u, 250u},
																			 .packet_count = 3,
																			 .sample_loss = true,
																			 .stream_open = true});
		return Vita49Serializer::serializeContext(context);
	}

	void requireContextPacketHeader(const std::vector<std::uint8_t>& bytes)
	{
		using namespace serial::vita49;
		REQUIRE(bytes.size() % 4u == 0u);
		REQUIRE((readU32(bytes, 0) & 0xFFFF0000u) == 0x4A630000u);
		REQUIRE((readU32(bytes, 0) & 0x0000FFFFu) == bytes.size() / 4u);
		REQUIRE(readU32(bytes, 4) == 0x55667788u);
		REQUIRE(readU64(bytes, 8) == kFersVrtIqClassId);
		REQUIRE(readU32(bytes, 16) == 100u);
		REQUIRE(readU64(bytes, 20) == 250u);
		REQUIRE(readU32(bytes, 28) == kFersContextCif0);
	}

	void requireContextPacketFields(const std::vector<std::uint8_t>& bytes)
	{
		using namespace serial::vita49;
		REQUIRE((readU32(bytes, 32) & TrailerSampleLoss) == TrailerSampleLoss);
		REQUIRE(readU64(bytes, 36) == makeComplexInt16PayloadFormat());
		REQUIRE(readF64(bytes, 44) == 2.5e6);
		REQUIRE(readF64(bytes, 52) == 9.6e9);
		REQUIRE(readF64(bytes, 60) == -1.25e6);
		REQUIRE(readF64(bytes, 68) == 5.0e6);
		REQUIRE(readF64(bytes, 76) == 2.0);
		REQUIRE(readU64(bytes, 84) == 0x1122334455667788ull);
	}

	void requireAsciiMetadataPadding(const std::vector<std::uint8_t>& bytes, const std::size_t ascii_offset)
	{
		const auto null_pos = std::find(bytes.begin() + static_cast<std::ptrdiff_t>(ascii_offset), bytes.end(), 0);
		REQUIRE(null_pos != bytes.end());
		for (auto it = bytes.begin() + static_cast<std::ptrdiff_t>(ascii_offset); it != null_pos; ++it)
		{
			REQUIRE(*it <= 0x7Fu);
		}
		for (auto it = std::next(null_pos); it != bytes.end(); ++it)
		{
			REQUIRE(*it == 0u);
		}
	}

	void requireContextReceiverMetadata(const nlohmann::json& metadata)
	{
		using namespace serial::vita49;
		REQUIRE(metadata.at("schema") == "fers-vita49-context-v1");
		REQUIRE(metadata.at("simulation_name") == "sim");
		REQUIRE(metadata.at("receiver").at("id").get<std::uint64_t>() == 0x1122334455667788ull);
		REQUIRE(metadata.at("receiver").at("name") == "rx1");
		REQUIRE(metadata.at("receiver").at("mode") == "fmcw");
		REQUIRE(metadata.at("receiver").at("adc_bits") == 14u);
		const auto context_flags = metadata.at("receiver").at("context_flags").get<std::uint32_t>();
		REQUIRE(
			(context_flags &
			 (ContextFlagDechirped | ContextFlagSampleLoss | ContextFlagStreamOpen | ContextFlagFmcwMetadataPresent)) ==
			(ContextFlagDechirped | ContextFlagSampleLoss | ContextFlagStreamOpen | ContextFlagFmcwMetadataPresent));
	}

	void requireContextCoordinateMetadata(const nlohmann::json& metadata)
	{
		REQUIRE(metadata.at("coordinate_frame").at("frame") == "UTM");
		REQUIRE(metadata.at("coordinate_frame").at("origin").at("latitude") == -33.9);
		REQUIRE(metadata.at("coordinate_frame").at("origin").at("longitude") == 18.4);
		REQUIRE(metadata.at("coordinate_frame").at("origin").at("altitude") == 111.0);
		REQUIRE(metadata.at("coordinate_frame").at("utm_zone") == 34);
		REQUIRE_FALSE(metadata.at("coordinate_frame").at("utm_north_hemisphere").get<bool>());
	}

	void requireContextPlatformMetadata(const nlohmann::json& metadata)
	{
		REQUIRE(metadata.at("initial_platform_state").at("platform_id") == 0x99u);
		REQUIRE(metadata.at("initial_platform_state").at("platform_name") == "rx-platform");
		REQUIRE(metadata.at("initial_platform_state").at("position_m").at("x") == 1.0);
		REQUIRE(metadata.at("initial_platform_state").at("position_m").at("y") == 2.0);
		REQUIRE(metadata.at("initial_platform_state").at("position_m").at("z") == 3.0);
		REQUIRE(metadata.at("initial_platform_state").at("velocity_mps").at("x") == 4.0);
		REQUIRE(metadata.at("initial_platform_state").at("velocity_mps").at("y") == 5.0);
		REQUIRE(metadata.at("initial_platform_state").at("velocity_mps").at("z") == 6.0);
		REQUIRE(metadata.at("initial_platform_state").at("rotation_rad").at("azimuth") == 0.25);
		REQUIRE(metadata.at("initial_platform_state").at("rotation_rad").at("elevation") == -0.5);
	}

	void requireContextFmcwMetadata(const nlohmann::json& metadata)
	{
		REQUIRE(metadata.at("fmcw").at("present").get<bool>());
		REQUIRE(metadata.at("fmcw").at("waveform_shape") == "linear");
		REQUIRE(metadata.at("fmcw").at("sweep_direction") == "down");
		REQUIRE(metadata.at("fmcw").at("chirp_bandwidth_hz") == 20.0e6);
		REQUIRE(metadata.at("fmcw").at("chirp_duration_s") == 1.0e-3);
		REQUIRE(metadata.at("fmcw").at("chirp_period_s") == 2.0e-3);
		REQUIRE(metadata.at("fmcw").at("chirp_rate_hz_per_s") == 20.0e9);
		REQUIRE(metadata.at("fmcw").at("chirp_rate_signed_hz_per_s") == -20.0e9);
		REQUIRE(metadata.at("fmcw").at("start_frequency_offset_hz") == 1.5e6);
		REQUIRE(metadata.at("fmcw").at("chirp_count") == 64u);
		REQUIRE(metadata.at("fmcw").at("dechirp_mode") == "physical");
		REQUIRE(metadata.at("fmcw").at("dechirp_reference_source") == "transmitter");
		REQUIRE(metadata.at("fmcw").at("dechirp_reference_transmitter_id") == 0x1234u);
		REQUIRE(metadata.at("fmcw").at("dechirp_reference_transmitter_name") == "tx-lo");
		REQUIRE(metadata.at("waveform").at("kind") == "fmcw");
		REQUIRE(metadata.at("waveform").at("metadata_ref") == "fmcw");
	}

	[[nodiscard]] core::ReceiverStreamDescriptor basicCwStreamDescriptor()
	{
		return core::ReceiverStreamDescriptor{.receiver_id = 77,
											  .receiver_name = "rx",
											  .mode = "cw",
											  .sample_rate = 1.0,
											  .reference_frequency = 1.0e9,
											  .coordinate = {},
											  .initial_platform_state = {},
											  .fmcw = {}};
	}

	void requirePacketsFitPayload(const std::vector<serial::vita49::SerializedPacket>& packets,
								  const std::size_t max_bytes)
	{
		for (const auto& packet : packets)
		{
			REQUIRE(packet.bytes.size() <= max_bytes);
			REQUIRE((readU32(packet.bytes, 0) & 0x0000FFFFu) == packet.bytes.size() / 4u);
		}
	}

	[[nodiscard]] serial::vita49::SerializedPacket testSerializedPacket(const std::uint8_t payload_byte,
																		const std::uint64_t sample_count)
	{
		return serial::vita49::SerializedPacket{.bytes = {0, 0, 0, payload_byte},
												.stream_id = 9,
												.sample_count = sample_count,
												.first_sample_time = 0.0,
												.data_packet = true,
												.timestamp = serial::vita49::Timestamp{}};
	}

	[[nodiscard]] std::future<serial::vita49::EnqueueResult>
	enqueueAsync(serial::vita49::PacedSender& sender, const serial::vita49::SerializedPacket& packet,
				 std::atomic_bool& enqueue_finished)
	{
		return std::async(std::launch::async,
						  [&sender, packet, &enqueue_finished]
						  {
							  auto result = sender.enqueue(packet);
							  enqueue_finished.store(true);
							  return result;
						  });
	}

	void requireLiveTelemetryStats(const core::OutputStats& final_stats,
								   const std::vector<core::OutputStats>& stats_snapshots)
	{
		REQUIRE_FALSE(stats_snapshots.empty());
		REQUIRE(final_stats.streams.size() == 1u);
		const auto& stream_stats = final_stats.streams.front();
		REQUIRE(stream_stats.first_sample_time.has_value());
		REQUIRE(stream_stats.end_sample_time.has_value());
		REQUIRE(stream_stats.first_timestamp.has_value());
		REQUIRE(stream_stats.end_timestamp.has_value());

		const auto first_timestamp = stream_stats.first_timestamp.value_or(core::Vita49Timestamp{});
		const auto end_timestamp = stream_stats.end_timestamp.value_or(core::Vita49Timestamp{});
		CHECK(stream_stats.first_sample_time.value_or(0.0) == 0.0);
		CHECK(stream_stats.end_sample_time.value_or(0.0) == 0.5);
		CHECK(first_timestamp.integer_seconds == 1700000000u);
		CHECK(first_timestamp.fractional_picoseconds == 0u);
		CHECK(end_timestamp.integer_seconds == 1700000000u);
		CHECK(end_timestamp.fractional_picoseconds == 500000000000u);
		CHECK(stats_snapshots.back().streams.size() == 1u);
		CHECK(stats_snapshots.back().streams.front().receiver_id == 7u);
	}

	void requireLivePacketTraces(const std::vector<core::ReceiverOutputPacketTrace>& packet_traces)
	{
		CHECK(std::ranges::any_of(packet_traces, [](const auto& trace) { return trace.event == "context"; }));
		CHECK(std::ranges::any_of(packet_traces, [](const auto& trace) { return trace.event == "data"; }));
		CHECK(std::ranges::all_of(packet_traces, [](const auto& trace) { return trace.sequence > 0u; }));
		const auto data_trace =
			std::ranges::find_if(packet_traces, [](const auto& trace) { return trace.event == "data"; });
		REQUIRE(data_trace != packet_traces.end());
		REQUIRE(data_trace->timestamp.has_value());
		const auto data_timestamp = data_trace->timestamp.value_or(core::Vita49Timestamp{});
		CHECK(data_timestamp.integer_seconds == 1700000000u);
		CHECK(data_timestamp.fractional_picoseconds == 0u);
	}

	class RecordingSender final : public serial::vita49::DatagramSender
	{
	public:
		void open(const std::string&, std::uint16_t) override { opened = true; }
		void send(const std::span<const std::uint8_t> bytes) override { sent.emplace_back(bytes.begin(), bytes.end()); }
		void close() noexcept override { closed = true; }

		bool opened = false;
		bool closed = false;
		std::vector<std::vector<std::uint8_t>> sent;
	};

	class BlockingFirstSendRecordingSender final : public serial::vita49::DatagramSender
	{
	public:
		void open(const std::string&, std::uint16_t) override { opened = true; }

		void send(const std::span<const std::uint8_t> bytes) override
		{
			std::unique_lock lock(mutex);
			sent.emplace_back(bytes.begin(), bytes.end());
			if (!first_send_entered)
			{
				first_send_entered = true;
				cv.notify_all();
				cv.wait(lock, [this] { return release_first_send; });
			}
		}

		void close() noexcept override { closed = true; }

		bool waitUntilFirstSendEntered()
		{
			std::unique_lock lock(mutex);
			return cv.wait_for(lock, std::chrono::seconds(1), [this] { return first_send_entered; });
		}

		void releaseFirstSend()
		{
			{
				std::scoped_lock const lock(mutex);
				release_first_send = true;
			}
			cv.notify_all();
		}

		bool opened = false;
		bool closed = false;
		std::vector<std::vector<std::uint8_t>> sent;

	private:
		std::mutex mutex;
		std::condition_variable cv;
		bool first_send_entered = false;
		bool release_first_send = false;
	};

	class TimedRecordingSender final : public serial::vita49::DatagramSender
	{
	public:
		void open(const std::string&, std::uint16_t) override {}
		void send(const std::span<const std::uint8_t> bytes) override
		{
			send_times.push_back(std::chrono::steady_clock::now());
			sent.emplace_back(bytes.begin(), bytes.end());
		}
		void close() noexcept override {}

		std::vector<std::chrono::steady_clock::time_point> send_times;
		std::vector<std::vector<std::uint8_t>> sent;
	};

	class ThrowingSender final : public serial::vita49::DatagramSender
	{
	public:
		void open(const std::string&, std::uint16_t) override {}
		void send(std::span<const std::uint8_t>) override { throw std::runtime_error("send failed"); }
		void close() noexcept override { closed = true; }

		bool closed = false;
	};

	class FailFirstSignalSender final : public serial::vita49::DatagramSender
	{
	public:
		void open(const std::string&, std::uint16_t) override { opened = true; }
		void send(const std::span<const std::uint8_t> bytes) override
		{
			if (!failed_signal_send && bytes.size() <= 64u)
			{
				failed_signal_send = true;
				throw std::runtime_error("signal send failed");
			}
			sent.emplace_back(bytes.begin(), bytes.end());
		}
		void close() noexcept override { closed = true; }

		bool opened = false;
		bool closed = false;
		bool failed_signal_send = false;
		std::vector<std::vector<std::uint8_t>> sent;
	};
}

TEST_CASE("VITA signal data packet serializes golden header and IQ bytes", "[serial][vita49]")
{
	using namespace serial::vita49;

	const SignalDataPacket packet{
		.stream_id = 0x10203040u,
		.class_id = kFersVrtIqClassId,
		.timestamp = Timestamp{.integer_seconds = 0x01020304u, .fractional_picoseconds = 0x0102030405060708ull},
		.packet_count = 0x0Au,
		.iq_interleaved = {1, -2, 0x1234, static_cast<std::int16_t>(-0x1234)},
		.trailer = makeTrailer(true, true, true, true, false)};

	const auto bytes = Vita49Serializer::serializeSignalData(packet);

	REQUIRE(bytes.size() == 40u);
	REQUIRE(readU32(bytes, 0) == 0x1E6A000Au);
	REQUIRE(readU32(bytes, 4) == 0x10203040u);
	REQUIRE(readU64(bytes, 8) == kFersVrtIqClassId);
	REQUIRE(readU32(bytes, 16) == 0x01020304u);
	REQUIRE(readU64(bytes, 20) == 0x0102030405060708ull);
	REQUIRE(readU16(bytes, 28) == 0x0001u);
	REQUIRE(readU16(bytes, 30) == 0xFFFEu);
	REQUIRE(readU16(bytes, 32) == 0x1234u);
	REQUIRE(readU16(bytes, 34) == 0xEDCCu);
	REQUIRE(readU32(bytes, 36) ==
			(TrailerEnableValidData | TrailerEnableCalibratedTime | TrailerEnableReferenceLock |
			 TrailerEnableOverRange | TrailerEnableSampleLoss | TrailerValidData | TrailerCalibratedTime |
			 TrailerReferenceLock | TrailerOverRange));
}

TEST_CASE("VITA context packet serializes CIF and fields in profile order", "[serial][vita49]")
{
	using namespace serial::vita49;

	const auto bytes = serializedFmcwContextPacket();

	requireContextPacketHeader(bytes);
	requireContextPacketFields(bytes);
	constexpr std::size_t ascii_offset = 92u;
	requireAsciiMetadataPadding(bytes, ascii_offset);

	const auto metadata = nlohmann::json::parse(readAsciiMetadata(bytes, ascii_offset));
	requireContextReceiverMetadata(metadata);
	requireContextCoordinateMetadata(metadata);
	requireContextPlatformMetadata(metadata);
	requireContextFmcwMetadata(metadata);
}

TEST_CASE("VITA context metadata follows receiver mode when stale mode blocks are present", "[serial][vita49]")
{
	using namespace serial::vita49;

	const core::ReceiverStreamDescriptor stream{.receiver_id = 6,
												.receiver_name = "cw-rx",
												.mode = "cw",
												.sample_rate = 4.0,
												.reference_frequency = 5.0e9,
												.coordinate = {},
												.initial_platform_state = {},
												.cw = {.present = true, .carrier_frequency = 5.0e9},
												.fmcw = {.present = true, .waveform_shape = "linear"}};
	const auto context = Vita49ContextBuilder::build(ContextBuildRequest{.stream = stream,
																		 .stream_id = 0x01020304u,
																		 .simulation_name = "sim",
																		 .adc_fullscale = 1.0,
																		 .timestamp = Timestamp{10u, 0u},
																		 .packet_count = 0});
	const auto bytes = Vita49Serializer::serializeContext(context);
	const auto metadata = nlohmann::json::parse(readAsciiMetadata(bytes, 92u));
	const auto context_flags = metadata.at("receiver").at("context_flags").get<std::uint32_t>();

	REQUIRE((context_flags & ContextFlagCwMetadataPresent) == ContextFlagCwMetadataPresent);
	REQUIRE((context_flags & ContextFlagFmcwMetadataPresent) == 0u);
	REQUIRE(metadata.at("waveform").at("kind") == "cw");
	REQUIRE(metadata.at("waveform").at("metadata_ref") == "cw");
	REQUIRE(metadata.contains("cw"));
	REQUIRE_FALSE(metadata.contains("fmcw"));
}

TEST_CASE("VITA context metadata describes pulsed and CW waveform modes", "[serial][vita49]")
{
	using namespace serial::vita49;

	const core::ReceiverStreamDescriptor pulsed_stream{.receiver_id = 5,
													   .receiver_name = "pulse-rx",
													   .mode = "pulsed",
													   .sample_rate = 1.0e6,
													   .reference_frequency = 10.0e9,
													   .adc_bits = 12,
													   .coordinate = {},
													   .initial_platform_state = {},
													   .pulsed = {.present = true,
																  .window_length = 50.0e-6,
																  .window_prf = 1000.0,
																  .window_skip = 5.0e-6,
																  .window_count = 4,
																  .waveform_id = 99,
																  .waveform_name = "pulse",
																  .carrier_frequency = 10.0e9,
																  .power = 1200.0,
																  .pulse_width = 2.0e-6,
																  .native_sample_rate = 20.0e6,
																  .native_sample_count = 40},
													   .cw = {},
													   .fmcw = {}};
	const auto pulsed_context = Vita49ContextBuilder::build(ContextBuildRequest{.stream = pulsed_stream,
																				.stream_id = 0x11111111u,
																				.simulation_name = "mixed",
																				.adc_fullscale = 1.0,
																				.timestamp = Timestamp{10u, 0u}});
	const auto pulsed_bytes = Vita49Serializer::serializeContext(pulsed_context);
	const auto pulsed_metadata = nlohmann::json::parse(readAsciiMetadata(pulsed_bytes, 92u));
	REQUIRE((pulsed_metadata.at("receiver").at("context_flags").get<std::uint32_t>() &
			 ContextFlagPulsedMetadataPresent) == ContextFlagPulsedMetadataPresent);
	REQUIRE(pulsed_metadata.at("waveform").at("kind") == "pulsed");
	REQUIRE(pulsed_metadata.at("waveform").at("metadata_ref") == "pulsed");
	REQUIRE(pulsed_metadata.at("pulsed").at("window_prf_hz") == 1000.0);
	REQUIRE(pulsed_metadata.at("pulsed").at("pri_s") == 0.001);
	REQUIRE(pulsed_metadata.at("pulsed").at("pulse_width_s") == 2.0e-6);

	const core::ReceiverStreamDescriptor cw_stream{
		.receiver_id = 6,
		.receiver_name = "cw-rx",
		.mode = "cw",
		.sample_rate = 1.0e6,
		.reference_frequency = 5.0e9,
		.adc_bits = 12,
		.coordinate = {},
		.initial_platform_state = {},
		.pulsed = {},
		.cw = {.present = true, .waveform_id = 100, .waveform_name = "cw", .carrier_frequency = 5.0e9, .power = 25.0},
		.fmcw = {}};
	const auto cw_context = Vita49ContextBuilder::build(ContextBuildRequest{.stream = cw_stream,
																			.stream_id = 0x22222222u,
																			.simulation_name = "mixed",
																			.adc_fullscale = 1.0,
																			.timestamp = Timestamp{10u, 0u}});
	const auto cw_bytes = Vita49Serializer::serializeContext(cw_context);
	const auto cw_metadata = nlohmann::json::parse(readAsciiMetadata(cw_bytes, 92u));
	REQUIRE((cw_metadata.at("receiver").at("context_flags").get<std::uint32_t>() & ContextFlagCwMetadataPresent) ==
			ContextFlagCwMetadataPresent);
	REQUIRE(cw_metadata.at("waveform").at("kind") == "cw");
	REQUIRE(cw_metadata.at("waveform").at("metadata_ref") == "cw");
	REQUIRE(cw_metadata.at("cw").at("carrier_hz") == 5.0e9);
	REQUIRE(cw_metadata.at("cw").at("power_w") == 25.0);
}

TEST_CASE("VITA context serializer rejects unsupported CIF0", "[serial][vita49]")
{
	using namespace serial::vita49;

	const core::ReceiverStreamDescriptor stream{.receiver_id = 1,
												.receiver_name = "rx",
												.mode = "cw",
												.sample_rate = 1.0,
												.reference_frequency = 1.0e9,
												.coordinate = {},
												.initial_platform_state = {},
												.fmcw = {}};
	auto context = Vita49ContextBuilder::build(ContextBuildRequest{.stream = stream,
																   .stream_id = 0x01020304u,
																   .simulation_name = "sim",
																   .adc_fullscale = 1.0,
																   .timestamp = Timestamp{1u, 0u}});
	context.cif0 &= ~Cif0AsciiMetadata;

	REQUIRE_THROWS_WITH(Vita49Serializer::serializeContext(context),
						ContainsSubstring("Unsupported VITA 49.2 context CIF0"));
}

TEST_CASE("VITA packetizer uses first sample timestamp, packet cap, and big-endian payload", "[serial][vita49]")
{
	using namespace serial::vita49;

	Vita49Packetizer const packetizer(1'700'000'000'123'456'789ull, 1.0, 1400);
	REQUIRE(packetizer.maxComplexSamplesPerPacket() == 342u);

	std::vector<ComplexType> samples(1000, ComplexType(0.5, -0.5));
	const core::ReceiverSampleBlock block{.stream = basicCwStreamDescriptor(),
										  .first_sample_time = 1.25,
										  .sample_rate = 1.0,
										  .samples = samples,
										  .sample_start = 0};
	PacketCountSequencer counts;
	const auto result = packetizer.packetize(block, 0xABCDEF01u, counts, false);

	REQUIRE(result.packets.size() == 3u);
	REQUIRE(result.samples_emitted == 1000u);
	requirePacketsFitPayload(result.packets, 1400u);
	REQUIRE(readU32(result.packets.front().bytes, 16) == 1'700'000'001u);
	REQUIRE(readU64(result.packets.front().bytes, 20) == 373'456'789'000ull);
	REQUIRE(readU16(result.packets.front().bytes, 28) == 0x4000u);
	REQUIRE(readU16(result.packets.front().bytes, 30) == 0xC000u);
}

TEST_CASE("VITA packet count rolls over modulo 16", "[serial][vita49]")
{
	serial::vita49::PacketCountSequencer counts;
	counts.setNext(14);
	REQUIRE(counts.next() == 14u);
	REQUIRE(counts.next() == 15u);
	REQUIRE(counts.next() == 0u);
	REQUIRE(counts.next() == 1u);
}

TEST_CASE("VITA stream registry allocates stable non-truncated stream IDs", "[serial][vita49]")
{
	serial::vita49::StreamRegistry registry;
	const core::ReceiverStreamDescriptor a{.receiver_id = 0x00030000FFFFFFFFull,
										   .receiver_name = "rx-a",
										   .mode = "cw",
										   .coordinate = {},
										   .initial_platform_state = {},
										   .fmcw = {}};
	const core::ReceiverStreamDescriptor b{.receiver_id = 0x00030001FFFFFFFFull,
										   .receiver_name = "rx-b",
										   .mode = "cw",
										   .coordinate = {},
										   .initial_platform_state = {},
										   .fmcw = {}};
	const core::ReceiverStreamDescriptor a_pulsed{.receiver_id = a.receiver_id,
												  .receiver_name = "rx-a",
												  .mode = "pulsed",
												  .coordinate = {},
												  .initial_platform_state = {},
												  .fmcw = {}};

	const auto stream_a = registry.registerStream(a);
	const auto stream_a_again = registry.registerStream(a);
	const auto stream_b = registry.registerStream(b);
	const auto stream_a_pulsed = registry.registerStream(a_pulsed);

	REQUIRE(stream_a != 0u);
	REQUIRE(stream_a == stream_a_again);
	REQUIRE(stream_a != stream_b);
	REQUIRE(stream_a != stream_a_pulsed);
	REQUIRE(stream_a != static_cast<std::uint32_t>(a.receiver_id));
}

TEST_CASE("VITA paced sender blocks at queue depth instead of dropping", "[serial][vita49]")
{
	using namespace serial::vita49;
	using namespace std::chrono_literals;

	auto recording = std::make_unique<BlockingFirstSendRecordingSender>();
	auto* recording_raw = recording.get();
	PacedSender sender(std::move(recording), 1);
	sender.open("127.0.0.1", 1);
	sender.start(0.0);

	const auto first = testSerializedPacket(1, 10);
	const auto second = testSerializedPacket(2, 12);

	const auto first_result = sender.enqueue(first);
	REQUIRE(first_result.enqueued);
	REQUIRE_FALSE(first_result.dropped);
	REQUIRE(recording_raw->waitUntilFirstSendEntered());
	std::atomic_bool second_enqueue_finished = false;
	auto second_result_future = enqueueAsync(sender, second, second_enqueue_finished);
	std::this_thread::sleep_for(20ms);
	const auto second_finished_before_release = second_enqueue_finished.load();
	recording_raw->releaseFirstSend();

	REQUIRE(recording_raw->opened);
	REQUIRE_FALSE(second_finished_before_release);
	const auto second_result = second_result_future.get();
	REQUIRE(second_result.enqueued);
	REQUIRE_FALSE(second_result.dropped);
	sender.stop();

	REQUIRE(recording_raw->sent.size() == 2u);
	REQUIRE(sender.sentPacketCount(9) == 2u);
	REQUIRE(sender.droppedDataPacketCount(9) == 0u);
	REQUIRE(sender.droppedSampleCount(9) == 0u);
}

TEST_CASE("VITA packetizer encodes explicit sample loss flags", "[serial][vita49]")
{
	using namespace serial::vita49;
	Vita49Packetizer const packetizer(1'700'000'000'000'000'000ull, 1.0, 1400);
	std::vector<ComplexType> samples{ComplexType(0.0, 0.0)};
	const core::ReceiverStreamDescriptor stream{.receiver_id = 9,
												.receiver_name = "rx",
												.mode = "cw",
												.sample_rate = 1.0,
												.coordinate = {},
												.initial_platform_state = {},
												.fmcw = {}};
	const core::ReceiverSampleBlock block{.stream = stream, .sample_rate = 1.0, .samples = samples};
	PacketCountSequencer counts;
	const auto packets = packetizer.packetize(block, 9, counts, true);

	REQUIRE((readU32(packets.packets.front().bytes, packets.packets.front().bytes.size() - 4u) & TrailerSampleLoss) ==
			TrailerSampleLoss);

	const auto context = Vita49ContextBuilder::build(ContextBuildRequest{.stream = stream,
																		 .stream_id = 9,
																		 .simulation_name = "sim",
																		 .adc_fullscale = 1.0,
																		 .timestamp = Timestamp{10u, 0u},
																		 .packet_count = counts.next(),
																		 .sample_loss = true});
	const auto context_bytes = Vita49Serializer::serializeContext(context);
	REQUIRE((readU32(context_bytes, 32) & TrailerSampleLoss) == TrailerSampleLoss);
	const auto context_metadata = nlohmann::json::parse(readAsciiMetadata(context_bytes, 92u));
	REQUIRE((context_metadata.at("receiver").at("context_flags").get<std::uint32_t>() & ContextFlagSampleLoss) ==
			ContextFlagSampleLoss);
}

TEST_CASE("VITA paced sender drains due packets on stop", "[serial][vita49]")
{
	using namespace serial::vita49;

	auto recording = std::make_unique<RecordingSender>();
	auto* recording_raw = recording.get();
	PacedSender sender(std::move(recording), 4);
	sender.open("127.0.0.1", 1);
	sender.start(1000.0);

	const SerializedPacket packet{.bytes = {0, 0, 0, 1},
								  .stream_id = 5,
								  .sample_count = 1,
								  .first_sample_time = 1000.0,
								  .data_packet = true,
								  .timestamp = Timestamp{}};
	REQUIRE(sender.enqueue(packet).enqueued);
	sender.stop();

	REQUIRE(recording_raw->sent.size() == 1u);
	REQUIRE(sender.sentPacketCount(5) == 1u);
}

TEST_CASE("VITA paced sender stop drains future packets at wall-clock pace", "[serial][vita49]")
{
	using namespace serial::vita49;

	auto recording = std::make_unique<RecordingSender>();
	auto* recording_raw = recording.get();
	PacedSender sender(std::move(recording), 4);
	sender.open("127.0.0.1", 1);
	sender.start(0.0);

	const SerializedPacket packet{.bytes = {0, 0, 0, 1},
								  .stream_id = 5,
								  .sample_count = 1,
								  .first_sample_time = 0.08,
								  .data_packet = true,
								  .timestamp = Timestamp{}};
	REQUIRE(sender.enqueue(packet).enqueued);
	const SerializedPacket context_packet{.bytes = {0, 0, 0, 2},
										  .stream_id = 5,
										  .sample_count = 0,
										  .first_sample_time = 0.08,
										  .data_packet = false,
										  .context_packet = true,
										  .timestamp = Timestamp{}};
	REQUIRE(sender.enqueue(context_packet).enqueued);
	const auto start = std::chrono::steady_clock::now();
	sender.stop();
	const auto elapsed = std::chrono::steady_clock::now() - start;

	REQUIRE(elapsed >= std::chrono::milliseconds(50));
	REQUIRE(elapsed < std::chrono::milliseconds(500));
	REQUIRE(recording_raw->sent.size() == 2u);
	REQUIRE(sender.sentPacketCount(5) == 2u);
	REQUIRE(sender.droppedDataPacketCount(5) == 0u);
	REQUIRE(sender.droppedContextPacketCount(5) == 0u);
	REQUIRE(sender.droppedSampleCount(5) == 0u);
}

TEST_CASE("VITA paced sender never sends packets before scheduled wall-clock time", "[serial][vita49]")
{
	using namespace serial::vita49;
	using namespace std::chrono_literals;

	auto recording = std::make_unique<TimedRecordingSender>();
	auto* recording_raw = recording.get();
	PacedSender sender(std::move(recording), 4);
	sender.open("127.0.0.1", 1);
	const auto wall_start = std::chrono::steady_clock::now();
	sender.start(0.0);

	const SerializedPacket first{.bytes = {0, 0, 0, 1},
								 .stream_id = 5,
								 .sample_count = 1,
								 .first_sample_time = 0.04,
								 .data_packet = true,
								 .timestamp = Timestamp{}};
	const SerializedPacket second{.bytes = {0, 0, 0, 2},
								  .stream_id = 5,
								  .sample_count = 1,
								  .first_sample_time = 0.09,
								  .data_packet = true,
								  .timestamp = Timestamp{}};
	REQUIRE(sender.enqueue(first).enqueued);
	REQUIRE(sender.enqueue(second).enqueued);
	sender.stop();

	REQUIRE(recording_raw->send_times.size() == 2u);
	CHECK(recording_raw->send_times.at(0) >= wall_start + 40ms);
	CHECK(recording_raw->send_times.at(1) >= wall_start + 90ms);
	CHECK(recording_raw->sent.size() == 2u);
}

TEST_CASE("VITA paced sender catches datagram send failures", "[serial][vita49]")
{
	using namespace serial::vita49;

	auto throwing = std::make_unique<ThrowingSender>();
	auto* throwing_raw = throwing.get();
	PacedSender sender(std::move(throwing), 4);
	sender.open("127.0.0.1", 1);
	sender.start(0.0);

	const SerializedPacket packet{.bytes = {0, 0, 0, 1},
								  .stream_id = 7,
								  .sample_count = 1,
								  .first_sample_time = 0.0,
								  .data_packet = true,
								  .timestamp = Timestamp{}};
	REQUIRE(sender.enqueue(packet).enqueued);
	REQUIRE_NOTHROW(sender.stop());
	REQUIRE(throwing_raw->closed);
	REQUIRE(sender.sendFailureCount(7) == 1u);
	REQUIRE(sender.droppedDataPacketCount(7) == 1u);
	REQUIRE(sender.droppedSampleCount(7) == 1u);
}

TEST_CASE("VITA output sink emits context, signal data, and stats through injected sender", "[serial][vita49]")
{
	using namespace serial::vita49;

	auto recording = std::make_unique<RecordingSender>();
	auto* recording_raw = recording.get();
	Vita49OutputSink sink(std::move(recording));
	const core::OutputConfig config{.mode = core::OutputMode::Vita49Udp,
									.vita49 = {.host = "127.0.0.1",
											   .port = 1,
											   .adc_fullscale = 1.0,
											   .queue_depth = 8,
											   .epoch_unix_nanoseconds = 1'700'000'000'000'000'000ull,
											   .max_udp_payload = 1400}};
	sink.initializeRun(config, "sim");

	const core::ReceiverStreamDescriptor stream{.receiver_id = 42,
												.receiver_name = "rx",
												.mode = "cw",
												.sample_rate = 2.0,
												.reference_frequency = 9.0e8,
												.coordinate = {},
												.initial_platform_state = {},
												.fmcw = {}};
	std::vector<ComplexType> samples{ComplexType(0.25, -0.25), ComplexType(2.0, 0.0)};
	const core::ReceiverSampleBlock block{
		.stream = stream, .first_sample_time = 0.0, .sample_rate = 2.0, .samples = samples, .sample_start = 0};

	sink.submitBlock(block);
	const auto stats = sink.finalize();

	REQUIRE(recording_raw->opened);
	REQUIRE(recording_raw->closed);
	REQUIRE(recording_raw->sent.size() >= 3u);
	REQUIRE(stats.streams.size() == 1u);
	CHECK(stats.streams.front().mode == "cw");
	CHECK(stats.streams.front().samples_emitted == 2u);
	CHECK(stats.streams.front().packets_emitted == 1u);
	CHECK(stats.streams.front().context_packets >= 2u);
	CHECK(stats.streams.front().over_range_count == 1u);
}

TEST_CASE("VITA output sink keeps mixed receiver modes as distinct streams", "[serial][vita49]")
{
	using namespace serial::vita49;

	auto recording = std::make_unique<RecordingSender>();
	Vita49OutputSink sink(std::move(recording));
	const core::OutputConfig config{.mode = core::OutputMode::Vita49Udp,
									.vita49 = {.host = "127.0.0.1",
											   .port = 1,
											   .adc_fullscale = 1.0,
											   .queue_depth = 16,
											   .epoch_unix_nanoseconds = 1'700'000'000'000'000'000ull,
											   .max_udp_payload = 1400}};
	sink.initializeRun(config, "mixed");

	const core::ReceiverStreamDescriptor pulsed{.receiver_id = 42,
												.receiver_name = "rx",
												.mode = "pulsed",
												.sample_rate = 2.0,
												.reference_frequency = 9.0e8,
												.coordinate = {},
												.initial_platform_state = {},
												.pulsed = {.present = true, .window_length = 1.0, .window_prf = 2.0},
												.cw = {},
												.fmcw = {}};
	const core::ReceiverStreamDescriptor cw{.receiver_id = 42,
											.receiver_name = "rx",
											.mode = "cw",
											.sample_rate = 2.0,
											.reference_frequency = 9.0e8,
											.coordinate = {},
											.initial_platform_state = {},
											.pulsed = {},
											.cw = {.present = true, .carrier_frequency = 9.0e8},
											.fmcw = {}};
	const core::ReceiverStreamDescriptor fmcw{.receiver_id = 43,
											  .receiver_name = "rx-fmcw",
											  .mode = "fmcw",
											  .sample_rate = 2.0,
											  .reference_frequency = 9.0e8,
											  .coordinate = {},
											  .initial_platform_state = {},
											  .pulsed = {},
											  .cw = {},
											  .fmcw = {.present = true,
													   .waveform_shape = "linear",
													   .chirp_bandwidth = 1.0e6,
													   .chirp_duration = 1.0e-3,
													   .chirp_period = 2.0e-3,
													   .chirp_rate = 1.0e9}};
	std::vector<ComplexType> samples{ComplexType(0.25, -0.25)};
	sink.submitBlock(
		core::ReceiverSampleBlock{.stream = pulsed, .first_sample_time = 0.0, .sample_rate = 2.0, .samples = samples});
	sink.submitBlock(
		core::ReceiverSampleBlock{.stream = cw, .first_sample_time = 0.0, .sample_rate = 2.0, .samples = samples});
	sink.submitBlock(
		core::ReceiverSampleBlock{.stream = fmcw, .first_sample_time = 0.0, .sample_rate = 2.0, .samples = samples});

	const auto stats = sink.finalize();

	REQUIRE(stats.streams.size() == 3u);
	std::vector<std::string> modes;
	std::vector<std::uint32_t> stream_ids;
	for (const auto& stream : stats.streams)
	{
		modes.push_back(stream.mode);
		stream_ids.push_back(stream.stream_id);
		CHECK(stream.samples_emitted == 1u);
		CHECK(stream.packets_emitted == 1u);
	}
	std::ranges::sort(modes);
	std::ranges::sort(stream_ids);
	REQUIRE(modes == std::vector<std::string>{"cw", "fmcw", "pulsed"});
	REQUIRE(std::ranges::unique(stream_ids).begin() == stream_ids.end());
}

TEST_CASE("VITA output sink reports sample loss only for send failures", "[serial][vita49]")
{
	using namespace serial::vita49;

	auto failing = std::make_unique<FailFirstSignalSender>();
	auto* failing_raw = failing.get();
	Vita49OutputSink sink(std::move(failing));
	const core::OutputConfig config{.mode = core::OutputMode::Vita49Udp,
									.vita49 = {.host = "127.0.0.1",
											   .port = 1,
											   .adc_fullscale = 1.0,
											   .queue_depth = 8,
											   .epoch_unix_nanoseconds = 1'700'000'000'000'000'000ull,
											   .max_udp_payload = 1400}};
	sink.initializeRun(config, "sim");

	const core::ReceiverStreamDescriptor stream{.receiver_id = 43,
												.receiver_name = "rx-loss",
												.mode = "cw",
												.sample_rate = 1.0,
												.reference_frequency = 9.0e8,
												.coordinate = {},
												.initial_platform_state = {},
												.fmcw = {}};
	std::vector<ComplexType> samples{ComplexType(0.25, -0.25)};
	const core::ReceiverSampleBlock block{
		.stream = stream, .first_sample_time = 0.0, .sample_rate = 1.0, .samples = samples, .sample_start = 0};

	sink.submitBlock(block);
	const auto stats = sink.finalize();

	REQUIRE(failing_raw->opened);
	REQUIRE(failing_raw->closed);
	REQUIRE(failing_raw->failed_signal_send);
	REQUIRE(stats.streams.size() == 1u);
	CHECK(stats.streams.front().packets_emitted == 0u);
	CHECK(stats.streams.front().samples_emitted == 0u);
	CHECK(stats.streams.front().packets_dropped == 1u);
	CHECK(stats.streams.front().samples_dropped == 1u);
	CHECK(stats.streams.front().context_packets >= 2u);
	CHECK(std::ranges::any_of(failing_raw->sent, [](const auto& bytes)
							  { return bytes.size() > 64u && (readU32(bytes, 32) & TrailerSampleLoss) != 0u; }));
}

TEST_CASE("VITA output sink emits live telemetry snapshots and packet traces", "[serial][vita49][telemetry]")
{
	using namespace serial::vita49;

	auto recording = std::make_unique<RecordingSender>();
	std::vector<core::OutputStats> stats_snapshots;
	std::vector<core::ReceiverOutputPacketTrace> packet_traces;
	Vita49OutputSink sink(
		std::move(recording),
		[&](const std::optional<core::OutputStats>& stats, std::span<const core::ReceiverOutputPacketTrace> packets)
		{
			if (stats.has_value())
			{
				stats_snapshots.push_back(*stats);
			}
			packet_traces.insert(packet_traces.end(), packets.begin(), packets.end());
		});
	const core::OutputConfig config{.mode = core::OutputMode::Vita49Udp,
									.vita49 = {.host = "127.0.0.1",
											   .port = 1,
											   .adc_fullscale = 1.0,
											   .queue_depth = 8,
											   .epoch_unix_nanoseconds = 1'700'000'000'000'000'000ull,
											   .max_udp_payload = 1400}};
	sink.initializeRun(config, "sim");

	const core::ReceiverStreamDescriptor stream{.receiver_id = 7,
												.receiver_name = "rx-live",
												.mode = "cw",
												.sample_rate = 4.0,
												.reference_frequency = 1.0e9,
												.coordinate = {},
												.initial_platform_state = {},
												.fmcw = {}};
	std::vector<ComplexType> samples{ComplexType(0.1, 0.2), ComplexType(0.3, 0.4)};
	const core::ReceiverSampleBlock block{
		.stream = stream, .first_sample_time = 0.0, .sample_rate = 4.0, .samples = samples, .sample_start = 0};

	sink.submitBlock(block);
	const auto final_stats = sink.finalize();

	requireLiveTelemetryStats(final_stats, stats_snapshots);
	requireLivePacketTraces(packet_traces);
}

TEST_CASE("VITA output sink batches packet trace telemetry", "[serial][vita49][telemetry]")
{
	using namespace serial::vita49;

	auto recording = std::make_unique<RecordingSender>();
	std::size_t packet_callback_count = 0;
	std::size_t largest_packet_batch = 0;
	std::vector<core::ReceiverOutputPacketTrace> packet_traces;
	Vita49OutputSink sink(
		std::move(recording),
		[&](const std::optional<core::OutputStats>&, std::span<const core::ReceiverOutputPacketTrace> packets)
		{
			if (!packets.empty())
			{
				++packet_callback_count;
				largest_packet_batch = std::max(largest_packet_batch, packets.size());
			}
			packet_traces.insert(packet_traces.end(), packets.begin(), packets.end());
		});
	const core::OutputConfig config{.mode = core::OutputMode::Vita49Udp,
									.vita49 = {.host = "127.0.0.1",
											   .port = 1,
											   .adc_fullscale = 1.0,
											   .queue_depth = 128,
											   .epoch_unix_nanoseconds = 1'700'000'000'000'000'000ull,
											   .max_udp_payload = 1400}};
	sink.initializeRun(config, "sim");

	const core::ReceiverStreamDescriptor stream{.receiver_id = 8,
												.receiver_name = "rx-batched",
												.mode = "cw",
												.sample_rate = 1.0e9,
												.reference_frequency = 1.0e9,
												.coordinate = {},
												.initial_platform_state = {},
												.fmcw = {}};
	std::vector<ComplexType> samples(22'000, ComplexType(0.1, 0.2));
	const core::ReceiverSampleBlock block{
		.stream = stream, .first_sample_time = 0.0, .sample_rate = 1.0e9, .samples = samples, .sample_start = 0};

	sink.submitBlock(block);
	(void)sink.finalize();

	const auto data_trace_count = std::ranges::count_if(packet_traces, [](const auto& trace)
														{ return trace.event == "data" && trace.data_packet; });
	CHECK(data_trace_count > 64);
	CHECK(largest_packet_batch >= 64u);
	CHECK(packet_callback_count < static_cast<std::size_t>(data_trace_count));
	CHECK(std::ranges::all_of(packet_traces, [](const auto& trace) { return trace.sequence > 0u; }));
}

TEST_CASE("VITA output sink can disable packet trace telemetry", "[serial][vita49][telemetry]")
{
	using namespace serial::vita49;

	auto recording = std::make_unique<RecordingSender>();
	std::vector<core::OutputStats> stats_snapshots;
	std::vector<core::ReceiverOutputPacketTrace> packet_traces;
	Vita49OutputSink sink(
		std::move(recording),
		[&](const std::optional<core::OutputStats>& stats, std::span<const core::ReceiverOutputPacketTrace> packets)
		{
			if (stats.has_value())
			{
				stats_snapshots.push_back(*stats);
			}
			packet_traces.insert(packet_traces.end(), packets.begin(), packets.end());
		});
	const core::OutputConfig config{.mode = core::OutputMode::Vita49Udp,
									.vita49 = {.host = "127.0.0.1",
											   .port = 1,
											   .adc_fullscale = 1.0,
											   .queue_depth = 8,
											   .epoch_unix_nanoseconds = 1'700'000'000'000'000'000ull,
											   .max_udp_payload = 1400,
											   .packet_trace_enabled = false}};
	sink.initializeRun(config, "sim");

	const core::ReceiverStreamDescriptor stream{.receiver_id = 9,
												.receiver_name = "rx-no-trace",
												.mode = "cw",
												.sample_rate = 4.0,
												.reference_frequency = 1.0e9,
												.coordinate = {},
												.initial_platform_state = {},
												.fmcw = {}};
	std::vector<ComplexType> samples{ComplexType(0.1, 0.2), ComplexType(0.3, 0.4)};
	const core::ReceiverSampleBlock block{
		.stream = stream, .first_sample_time = 0.0, .sample_rate = 4.0, .samples = samples, .sample_start = 0};

	sink.submitBlock(block);
	const auto final_stats = sink.finalize();

	REQUIRE_FALSE(stats_snapshots.empty());
	REQUIRE(final_stats.streams.size() == 1u);
	CHECK(stats_snapshots.back().streams.size() == 1u);
	CHECK(stats_snapshots.back().streams.front().receiver_id == 9u);
	CHECK(packet_traces.empty());
}

TEST_CASE("VITA output sink starts pacing at simulation start time", "[serial][vita49]")
{
	using namespace serial::vita49;

	ParamGuard const guard;
	params::setTime(2.0, 3.0);

	auto recording = std::make_unique<RecordingSender>();
	auto* recording_raw = recording.get();
	Vita49OutputSink sink(std::move(recording));
	const core::OutputConfig config{.mode = core::OutputMode::Vita49Udp,
									.vita49 = {.host = "127.0.0.1",
											   .port = 1,
											   .adc_fullscale = 1.0,
											   .queue_depth = 8,
											   .epoch_unix_nanoseconds = 1'700'000'000'000'000'000ull,
											   .max_udp_payload = 1400}};
	sink.initializeRun(config, "sim");

	const core::ReceiverStreamDescriptor stream{.receiver_id = 42,
												.receiver_name = "rx",
												.mode = "cw",
												.sample_rate = 1.0,
												.reference_frequency = 9.0e8,
												.coordinate = {},
												.initial_platform_state = {},
												.fmcw = {}};
	std::vector<ComplexType> samples{ComplexType(0.25, -0.25)};
	const core::ReceiverSampleBlock block{
		.stream = stream, .first_sample_time = 2.0, .sample_rate = 1.0, .samples = samples, .sample_start = 0};

	const auto start = std::chrono::steady_clock::now();
	sink.submitBlock(block);
	(void)sink.finalize();
	const auto elapsed = std::chrono::steady_clock::now() - start;

	REQUIRE(recording_raw->sent.size() >= 3u);
	REQUIRE(elapsed < std::chrono::milliseconds(500));
}

TEST_CASE("VITA output sink drains future packets before final stats", "[serial][vita49]")
{
	using namespace serial::vita49;

	ParamGuard const guard;
	params::setTime(0.0, 1.0);

	auto recording = std::make_unique<RecordingSender>();
	auto* recording_raw = recording.get();
	Vita49OutputSink sink(std::move(recording));
	const core::OutputConfig config{.mode = core::OutputMode::Vita49Udp,
									.vita49 = {.host = "127.0.0.1",
											   .port = 1,
											   .adc_fullscale = 1.0,
											   .queue_depth = 8,
											   .epoch_unix_nanoseconds = 1'700'000'000'000'000'000ull,
											   .max_udp_payload = 1400}};
	sink.initializeRun(config, "sim");

	const core::ReceiverStreamDescriptor stream{.receiver_id = 42,
												.receiver_name = "rx",
												.mode = "cw",
												.sample_rate = 1.0,
												.reference_frequency = 9.0e8,
												.coordinate = {},
												.initial_platform_state = {},
												.fmcw = {}};
	std::vector<ComplexType> samples{ComplexType(0.25, -0.25)};
	const core::ReceiverSampleBlock block{
		.stream = stream, .first_sample_time = 0.08, .sample_rate = 1.0, .samples = samples, .sample_start = 0};

	const auto start = std::chrono::steady_clock::now();
	sink.submitBlock(block);
	const auto stats = sink.finalize();
	const auto elapsed = std::chrono::steady_clock::now() - start;

	REQUIRE(elapsed >= std::chrono::milliseconds(50));
	REQUIRE(recording_raw->sent.size() >= 3u);
	REQUIRE(stats.streams.size() == 1u);
	CHECK(stats.streams.front().packets_emitted == 1u);
	CHECK(stats.streams.front().samples_emitted == 1u);
	CHECK(stats.streams.front().context_packets >= 2u);
	CHECK(stats.streams.front().packets_dropped == 0u);
	CHECK(stats.streams.front().samples_dropped == 0u);
}

#ifndef _WIN32
TEST_CASE("VITA UDP sender loopback carries stream ID and packet count", "[serial][vita49][udp]")
{
	using namespace serial::vita49;

	const int receiver = ::socket(AF_INET, SOCK_DGRAM, 0);
	if (receiver < 0)
	{
		SKIP("loopback socket unavailable");
	}

	sockaddr_in bind_addr{};
	bind_addr.sin_family = AF_INET;
	bind_addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
	bind_addr.sin_port = 0;
	if (::bind(receiver, asSockaddr(bind_addr), sizeof(bind_addr)) != 0)
	{
		::close(receiver);
		SKIP("loopback bind unavailable");
	}

	socklen_t addr_len = sizeof(bind_addr);
	REQUIRE(::getsockname(receiver, asSockaddr(bind_addr), &addr_len) == 0);

	timeval timeout{};
	timeout.tv_sec = 1;
	REQUIRE(::setsockopt(receiver, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof(timeout)) == 0);

	const SignalDataPacket packet{.stream_id = 0x01020304u,
								  .timestamp = Timestamp{1u, 2u},
								  .packet_count = 15u,
								  .iq_interleaved = {1, 2},
								  .trailer = makeTrailer(true, true, true, false, false)};
	const auto bytes = Vita49Serializer::serializeSignalData(packet);

	UdpSender sender;
	sender.open("127.0.0.1", ntohs(bind_addr.sin_port));
	sender.send(bytes);

	std::vector<std::uint8_t> received(256);
	const auto received_count = ::recv(receiver, received.data(), received.size(), 0);
	::close(receiver);
	REQUIRE(received_count == static_cast<ssize_t>(bytes.size()));
	received.resize(static_cast<std::size_t>(received_count));

	REQUIRE(readU32(received, 4) == 0x01020304u);
	REQUIRE(((readU32(received, 0) >> 16u) & 0x0Fu) == 15u);
}
#endif
