// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

#include "serial/vita49/vita49_types.h"

namespace serial::vita49
{
	class ByteWriter
	{
	public:
		explicit ByteWriter(std::size_t reserve_bytes = 0);

		void writeU16(std::uint16_t value);
		void writeI16(std::int16_t value);
		void writeU32(std::uint32_t value);
		void writeU64(std::uint64_t value);
		void writeF64(RealType value);
		void writeAsciiMetadata(std::string_view value);
		void writeBytes(std::span<const std::uint8_t> bytes);

		[[nodiscard]] const std::vector<std::uint8_t>& bytes() const noexcept;
		[[nodiscard]] std::vector<std::uint8_t> takeBytes() noexcept;

	private:
		std::vector<std::uint8_t> _bytes;
	};

	struct FixedFullscaleSignalDataPacket
	{
		std::uint32_t stream_id;
		std::uint64_t class_id;
		Timestamp timestamp;
		std::uint8_t packet_count;
		bool valid_data;
		bool calibrated_time;
		bool reference_lock;
		bool sample_loss;
		std::span<const ComplexType> samples;
		RealType fullscale;
	};

	class Vita49Serializer
	{
	public:
		[[nodiscard]] static std::vector<std::uint8_t> serializeSignalData(const SignalDataPacket& packet);
		[[nodiscard]] static SignalDataSerializationResult
		serializeSignalDataFixedFullscale(const FixedFullscaleSignalDataPacket& packet);
		[[nodiscard]] static std::vector<std::uint8_t> serializeContext(const ContextPacket& packet);
	};
}
