// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "serial/vita49/vita49_packetizer.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "serial/vita49/vita49_serializer.h"

namespace serial::vita49
{
	Vita49Packetizer::Vita49Packetizer(const std::uint64_t epoch_unix_nanoseconds, const RealType adc_fullscale,
									   const std::uint16_t max_udp_payload_bytes) :
		_epoch_unix_nanoseconds(epoch_unix_nanoseconds), _adc_fullscale(adc_fullscale),
		_max_udp_payload_bytes(max_udp_payload_bytes),
		_max_complex_samples_per_packet(maxComplexSamplesPerSignalPacket(max_udp_payload_bytes))
	{
		if (!std::isfinite(adc_fullscale) || adc_fullscale <= 0.0)
		{
			throw std::invalid_argument("VITA ADC full-scale must be positive");
		}
	}

	std::uint64_t Vita49Packetizer::epochUnixNanoseconds() const noexcept { return _epoch_unix_nanoseconds; }

	RealType Vita49Packetizer::adcFullscale() const noexcept { return _adc_fullscale; }

	std::uint16_t Vita49Packetizer::maxUdpPayloadBytes() const noexcept { return _max_udp_payload_bytes; }

	std::size_t Vita49Packetizer::maxComplexSamplesPerPacket() const noexcept
	{
		return _max_complex_samples_per_packet;
	}

	PacketizerResult Vita49Packetizer::packetize(const core::ReceiverSampleBlock& block, const std::uint32_t stream_id,
												 PacketCountSequencer& packet_counts,
												 const bool sample_loss_pending) const
	{
		PacketizerResult result;
		if (block.samples.empty())
		{
			return result;
		}

		const RealType sample_rate = block.sample_rate > 0.0 ? block.sample_rate : block.stream.sample_rate;
		if (!std::isfinite(sample_rate) || sample_rate <= 0.0)
		{
			throw std::invalid_argument("VITA sample block must have positive sample rate");
		}

		std::size_t offset = 0;
		bool loss_flag_for_next_packet = sample_loss_pending;
		result.packets.reserve((block.samples.size() + _max_complex_samples_per_packet - 1u) /
							   _max_complex_samples_per_packet);
		while (offset < block.samples.size())
		{
			const auto count = std::min(_max_complex_samples_per_packet, block.samples.size() - offset);
			const RealType packet_first_sample_time =
				block.first_sample_time + static_cast<RealType>(offset) / sample_rate;
			const auto timestamp = timestampFromEpoch(_epoch_unix_nanoseconds, packet_first_sample_time);

			auto serialized = Vita49Serializer::serializeSignalDataFixedFullscale(
				FixedFullscaleSignalDataPacket{.stream_id = stream_id,
											   .class_id = kFersVrtIqClassId,
											   .timestamp = timestamp,
											   .packet_count = packet_counts.next(),
											   .valid_data = block.valid_data,
											   .calibrated_time = block.calibrated_time,
											   .reference_lock = block.reference_lock,
											   .sample_loss = loss_flag_for_next_packet,
											   .samples = block.samples.subspan(offset, count),
											   .fullscale = _adc_fullscale});
			const bool packet_over_range = serialized.clipped_sample_count > 0;
			result.over_range_count += serialized.clipped_sample_count;

			auto& bytes = serialized.bytes;
			if (bytes.size() > _max_udp_payload_bytes)
			{
				throw std::logic_error("VITA signal packet exceeded max UDP payload");
			}

			result.packets.push_back(SerializedPacket{.bytes = std::move(bytes),
													  .stream_id = stream_id,
													  .sample_count = count,
													  .first_sample_time = packet_first_sample_time,
													  .data_packet = true,
													  .context_packet = false,
													  .over_range = packet_over_range,
													  .sample_loss = loss_flag_for_next_packet,
													  .timestamp = timestamp});
			result.samples_emitted += count;
			offset += count;
			loss_flag_for_next_packet = false;
		}

		return result;
	}

	SerializedPacket Vita49Packetizer::makeContextPacket(const ContextPacket& context) const
	{
		auto bytes = Vita49Serializer::serializeContext(context);
		if (bytes.size() > _max_udp_payload_bytes)
		{
			throw std::logic_error("VITA context packet exceeded max UDP payload");
		}
		return SerializedPacket{.bytes = std::move(bytes),
								.stream_id = context.stream_id,
								.sample_count = 0,
								.first_sample_time = 0.0,
								.data_packet = false,
								.context_packet = true,
								.over_range = false,
								.sample_loss = (context.context_flags & ContextFlagSampleLoss) != 0u,
								.timestamp = context.timestamp};
	}

}
