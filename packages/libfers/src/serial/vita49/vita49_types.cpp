// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "serial/vita49/vita49_types.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace serial::vita49
{
	std::uint8_t PacketCountSequencer::next() noexcept
	{
		const auto value = static_cast<std::uint8_t>(_next & 0x0Fu);
		_next = static_cast<std::uint8_t>((_next + 1u) & 0x0Fu);
		return value;
	}

	void PacketCountSequencer::setNext(const std::uint8_t packet_count) noexcept
	{
		_next = static_cast<std::uint8_t>(packet_count & 0x0Fu);
	}

	std::uint8_t PacketCountSequencer::peek() const noexcept { return static_cast<std::uint8_t>(_next & 0x0Fu); }

	std::uint32_t makeHeader(const PacketType type, const bool class_id_present, const bool trailer_present,
							 const IntegerTimestampMode tsi, const FractionalTimestampMode tsf,
							 const std::uint8_t packet_count, const std::uint16_t packet_size_words) noexcept
	{
		std::uint32_t header = static_cast<std::uint32_t>(type) << 28u;
		if (class_id_present)
		{
			header |= 1u << 27u;
		}
		if (trailer_present)
		{
			header |= 1u << 26u;
		}

		header |= 1u << 25u; // VITA 49.2 format bit.
		header |= static_cast<std::uint32_t>(tsi) << 22u;
		header |= static_cast<std::uint32_t>(tsf) << 20u;
		header |= static_cast<std::uint32_t>(packet_count & 0x0Fu) << 16u;
		header |= packet_size_words;
		return header;
	}

	Timestamp timestampFromEpoch(const std::uint64_t epoch_unix_nanoseconds, const RealType sample_time_seconds)
	{
		if (!std::isfinite(sample_time_seconds))
		{
			throw std::invalid_argument("VITA timestamp sample time must be finite");
		}

		const auto epoch_seconds = static_cast<std::uint64_t>(epoch_unix_nanoseconds / 1'000'000'000ull);
		const auto epoch_nanoseconds = static_cast<std::uint64_t>(epoch_unix_nanoseconds % 1'000'000'000ull);

		const auto offset_seconds = static_cast<long double>(sample_time_seconds);
		const long double offset_floor = std::floor(offset_seconds);
		const auto whole_offset_seconds = static_cast<std::int64_t>(offset_floor);
		long double fractional_seconds =
			(static_cast<long double>(epoch_nanoseconds) / 1.0e9L) + (offset_seconds - offset_floor);

		std::int64_t carry_seconds = 0;
		while (fractional_seconds >= 1.0L)
		{
			fractional_seconds -= 1.0L;
			++carry_seconds;
		}
		while (fractional_seconds < 0.0L)
		{
			fractional_seconds += 1.0L;
			--carry_seconds;
		}

		const auto signed_seconds = static_cast<std::int64_t>(epoch_seconds) + whole_offset_seconds + carry_seconds;
		if (signed_seconds < 0 || signed_seconds > static_cast<std::int64_t>(std::numeric_limits<std::uint32_t>::max()))
		{
			throw std::out_of_range("VITA UTC integer timestamp is outside 32-bit range");
		}

		auto fractional_picoseconds =
			static_cast<std::uint64_t>(std::llround(fractional_seconds * 1'000'000'000'000.0L));
		auto integer_seconds = static_cast<std::uint32_t>(signed_seconds);
		if (fractional_picoseconds >= 1'000'000'000'000ull)
		{
			fractional_picoseconds -= 1'000'000'000'000ull;
			if (integer_seconds == std::numeric_limits<std::uint32_t>::max())
			{
				throw std::out_of_range("VITA UTC integer timestamp overflow");
			}
			++integer_seconds;
		}

		return Timestamp{.integer_seconds = integer_seconds, .fractional_picoseconds = fractional_picoseconds};
	}

	std::uint32_t makeTrailer(const bool valid_data, const bool calibrated_time, const bool reference_lock,
							  const bool over_range, const bool sample_loss) noexcept
	{
		std::uint32_t trailer = TrailerEnableValidData | TrailerEnableCalibratedTime | TrailerEnableReferenceLock |
			TrailerEnableOverRange | TrailerEnableSampleLoss;
		if (valid_data)
		{
			trailer |= TrailerValidData;
		}
		if (calibrated_time)
		{
			trailer |= TrailerCalibratedTime;
		}
		if (reference_lock)
		{
			trailer |= TrailerReferenceLock;
		}
		if (over_range)
		{
			trailer |= TrailerOverRange;
		}
		if (sample_loss)
		{
			trailer |= TrailerSampleLoss;
		}
		return trailer;
	}

	std::uint32_t makeContextStateIndicators(const bool valid_data, const bool calibrated_time,
											 const bool reference_lock, const bool over_range,
											 const bool sample_loss) noexcept
	{
		return makeTrailer(valid_data, calibrated_time, reference_lock, over_range, sample_loss);
	}

	std::uint64_t makeComplexInt16PayloadFormat() noexcept
	{
		// Internal FERS profile payload-format word:
		// [63:60] complex Cartesian, [59:56] two's-complement integer,
		// [55:48] item bits, [47:40] packing bits, [39:32] vector size,
		// [31:16] repeat count, [15:0] reserved/version.
		constexpr std::uint64_t complex_cartesian = 0x2ull;
		constexpr std::uint64_t twos_complement = 0x1ull;
		constexpr std::uint64_t item_bits = 16ull;
		constexpr std::uint64_t packing_bits = 16ull;
		constexpr std::uint64_t vector_size = 2ull;
		constexpr std::uint64_t repeat_count = 1ull;
		constexpr std::uint64_t version = 1ull;
		return (complex_cartesian << 60u) | (twos_complement << 56u) | (item_bits << 48u) | (packing_bits << 40u) |
			(vector_size << 32u) | (repeat_count << 16u) | version;
	}

	std::size_t maxComplexSamplesPerSignalPacket(const std::uint16_t max_udp_payload_bytes)
	{
		if (max_udp_payload_bytes < kSignalDataFixedBytes + sizeof(std::int16_t) * 2u)
		{
			throw std::invalid_argument("VITA max UDP payload cannot fit one complex int16 sample");
		}
		return (max_udp_payload_bytes - kSignalDataFixedBytes) / (sizeof(std::int16_t) * 2u);
	}
}
