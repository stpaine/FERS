// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include "core/config.h"
#include "core/receiver_output.h"

namespace serial::vita49
{
	constexpr std::uint16_t kDefaultMaxUdpPayloadBytes = 1400;
	constexpr std::uint32_t kSignalDataFixedBytes = 32;

	// Internal placeholder until FERS obtains/finalizes an assigned VRT profile
	// identifier. The OUI is deliberately private and must be documented in the ICD.
	constexpr std::uint32_t kFersInternalOui = 0xFA5253u;
	constexpr std::uint16_t kFersInformationClassIqStream = 0x0001u;
	constexpr std::uint16_t kFersPacketClassV1 = 0x0001u;
	constexpr std::uint64_t kFersVrtIqClassId = (static_cast<std::uint64_t>(kFersInternalOui) << 40u) |
		(static_cast<std::uint64_t>(kFersInformationClassIqStream) << 24u) |
		(static_cast<std::uint64_t>(kFersPacketClassV1) << 8u) | 0x01u;

	enum class PacketType : std::uint8_t
	{
		SignalDataWithStreamId = 0x1,
		Context = 0x4,
	};

	enum class IntegerTimestampMode : std::uint8_t
	{
		None = 0,
		Utc = 1,
		Gps = 2,
		Other = 3,
	};

	enum class FractionalTimestampMode : std::uint8_t
	{
		None = 0,
		SampleCount = 1,
		RealTimePicoseconds = 2,
		FreeRunning = 3,
	};

	enum TrailerIndicator : std::uint32_t
	{
		TrailerEnableCalibratedTime = 1u << 31u,
		TrailerEnableValidData = 1u << 30u,
		TrailerEnableReferenceLock = 1u << 29u,
		TrailerEnableOverRange = 1u << 28u,
		TrailerEnableSampleLoss = 1u << 27u,
		TrailerCalibratedTime = 1u << 15u,
		TrailerValidData = 1u << 14u,
		TrailerReferenceLock = 1u << 13u,
		TrailerOverRange = 1u << 12u,
		TrailerSampleLoss = 1u << 11u,
	};

	enum ContextIndicator0 : std::uint32_t
	{
		Cif0StateIndicators = 1u << 31u,
		Cif0PayloadFormat = 1u << 30u,
		Cif0SampleRate = 1u << 29u,
		Cif0ReferenceFrequency = 1u << 28u,
		Cif0IfOffset = 1u << 27u,
		Cif0Bandwidth = 1u << 26u,
		Cif0ReferenceLevel = 1u << 25u,
		Cif0DeviceIdentifier = 1u << 24u,
		Cif0AsciiMetadata = 1u << 23u,
	};

	constexpr std::uint32_t kFersContextCif0 = Cif0StateIndicators | Cif0PayloadFormat | Cif0SampleRate |
		Cif0ReferenceFrequency | Cif0IfOffset | Cif0Bandwidth | Cif0ReferenceLevel | Cif0DeviceIdentifier |
		Cif0AsciiMetadata;

	enum ContextFlags : std::uint32_t // NOLINT(performance-enum-size)
	{
		ContextFlagDechirped = 1u << 0u,
		ContextFlagIfResampled = 1u << 1u,
		ContextFlagSampleLoss = 1u << 2u,
		ContextFlagStreamOpen = 1u << 3u,
		ContextFlagStreamClose = 1u << 4u,
		ContextFlagFmcwMetadataPresent = 1u << 5u,
		ContextFlagCwMetadataPresent = 1u << 6u,
		ContextFlagPulsedMetadataPresent = 1u << 7u,
		ContextFlagSfcwMetadataPresent = 1u << 8u,
	};

	struct Timestamp
	{
		std::uint32_t integer_seconds = 0;
		std::uint64_t fractional_picoseconds = 0;
	};

	struct SignalDataPacket
	{
		std::uint32_t stream_id = 0;
		std::uint64_t class_id = kFersVrtIqClassId;
		Timestamp timestamp;
		std::uint8_t packet_count = 0;
		std::vector<std::int16_t> iq_interleaved;
		std::uint32_t trailer = 0;
	};

	struct ContextPacket
	{
		std::uint32_t stream_id = 0;
		std::uint64_t class_id = kFersVrtIqClassId;
		Timestamp timestamp;
		std::uint8_t packet_count = 0;
		std::uint32_t cif0 = kFersContextCif0;
		std::uint32_t state_indicators = 0;
		std::uint64_t payload_format = 0;
		RealType sample_rate = 0.0;
		RealType reference_frequency = 0.0;
		RealType if_offset = 0.0;
		RealType bandwidth = 0.0;
		RealType adc_fullscale = 0.0;
		SimId receiver_id = 0;
		unsigned adc_bits = 0;
		std::uint32_t context_flags = 0;
		std::string receiver_name;
		std::string simulation_name;
		std::string receiver_mode;
		core::ReceiverStreamDescriptor::CoordinateContext coordinate;
		core::ReceiverStreamDescriptor::PlatformState initial_platform_state;
		core::ReceiverStreamDescriptor::PulsedContext pulsed;
		core::ReceiverStreamDescriptor::CwContext cw;
		core::ReceiverStreamDescriptor::FmcwContext fmcw;
		core::ReceiverStreamDescriptor::SfcwContext sfcw;
	};

	struct SerializedPacket
	{
		std::vector<std::uint8_t> bytes;
		std::uint32_t stream_id = 0;
		std::uint64_t sample_count = 0;
		RealType first_sample_time = 0.0;
		bool data_packet = false;
		bool context_packet = false;
		bool over_range = false;
		bool sample_loss = false;
		Timestamp timestamp;
	};

	struct PacketizerResult
	{
		std::vector<SerializedPacket> packets;
		std::uint64_t samples_emitted = 0;
		std::uint64_t over_range_count = 0;
	};

	struct SignalDataSerializationResult
	{
		std::vector<std::uint8_t> bytes;
		std::uint64_t clipped_sample_count = 0;
	};

	class PacketCountSequencer
	{
	public:
		[[nodiscard]] std::uint8_t next() noexcept;
		void setNext(std::uint8_t packet_count) noexcept;
		[[nodiscard]] std::uint8_t peek() const noexcept;

	private:
		std::uint8_t _next = 0;
	};

	[[nodiscard]] std::uint32_t makeHeader(PacketType type, bool class_id_present, bool trailer_present,
										   IntegerTimestampMode tsi, FractionalTimestampMode tsf,
										   std::uint8_t packet_count, std::uint16_t packet_size_words) noexcept;

	[[nodiscard]] Timestamp timestampFromEpoch(std::uint64_t epoch_unix_nanoseconds, RealType sample_time_seconds);
	[[nodiscard]] std::uint32_t makeTrailer(bool valid_data, bool calibrated_time, bool reference_lock, bool over_range,
											bool sample_loss) noexcept;
	[[nodiscard]] std::uint32_t makeContextStateIndicators(bool valid_data, bool calibrated_time, bool reference_lock,
														   bool over_range, bool sample_loss) noexcept;
	[[nodiscard]] std::uint64_t makeComplexInt16PayloadFormat() noexcept;
	[[nodiscard]] std::size_t maxComplexSamplesPerSignalPacket(std::uint16_t max_udp_payload_bytes);
}
