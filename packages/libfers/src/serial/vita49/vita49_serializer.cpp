// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "serial/vita49/vita49_serializer.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <limits>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string_view>

namespace serial::vita49
{
	namespace
	{
		[[nodiscard]] std::int16_t scaleComponentToInt16(const RealType value, const RealType fullscale,
														 bool& clipped) noexcept
		{
			if (value > fullscale)
			{
				clipped = true;
				return std::numeric_limits<std::int16_t>::max();
			}
			if (value < -fullscale)
			{
				clipped = true;
				return std::numeric_limits<std::int16_t>::min();
			}
			if (value == fullscale)
			{
				return std::numeric_limits<std::int16_t>::max();
			}
			if (value == -fullscale)
			{
				return std::numeric_limits<std::int16_t>::min();
			}

			const RealType scaled =
				std::round((value / fullscale) * static_cast<RealType>(std::numeric_limits<std::int16_t>::max()));
			return static_cast<std::int16_t>(
				std::clamp<RealType>(scaled, static_cast<RealType>(std::numeric_limits<std::int16_t>::min()),
									 static_cast<RealType>(std::numeric_limits<std::int16_t>::max())));
		}

		[[nodiscard]] std::uint16_t checkedWordCount(const std::size_t byte_count)
		{
			if (byte_count % sizeof(std::uint32_t) != 0u)
			{
				throw std::logic_error("VITA packet size must be 32-bit aligned");
			}
			const auto words = byte_count / sizeof(std::uint32_t);
			if (words > std::numeric_limits<std::uint16_t>::max())
			{
				throw std::length_error("VITA packet exceeds 16-bit word count");
			}
			return static_cast<std::uint16_t>(words);
		}

		[[nodiscard]] nlohmann::json makeCwMetadataJson(const ContextPacket& packet)
		{
			const RealType carrier_frequency =
				packet.cw.carrier_frequency != 0.0 ? packet.cw.carrier_frequency : packet.reference_frequency;
			return {{"present", packet.cw.present},
					{"waveform_id", packet.cw.waveform_id},
					{"waveform_name", packet.cw.waveform_name},
					{"carrier_hz", carrier_frequency},
					{"power_w", packet.cw.power}};
		}

		[[nodiscard]] nlohmann::json makePulsedMetadataJson(const ContextPacket& packet)
		{
			const RealType carrier_frequency =
				packet.pulsed.carrier_frequency != 0.0 ? packet.pulsed.carrier_frequency : packet.reference_frequency;
			nlohmann::json result = {{"present", packet.pulsed.present},
									 {"waveform_id", packet.pulsed.waveform_id},
									 {"waveform_name", packet.pulsed.waveform_name},
									 {"carrier_hz", carrier_frequency},
									 {"power_w", packet.pulsed.power},
									 {"pulse_width_s", packet.pulsed.pulse_width},
									 {"native_sample_rate_hz", packet.pulsed.native_sample_rate},
									 {"native_sample_count", packet.pulsed.native_sample_count},
									 {"window_length_s", packet.pulsed.window_length},
									 {"window_prf_hz", packet.pulsed.window_prf},
									 {"window_skip_s", packet.pulsed.window_skip},
									 {"window_count", packet.pulsed.window_count}};
			result["pri_s"] = packet.pulsed.window_prf > 0.0 ? nlohmann::json(1.0 / packet.pulsed.window_prf)
															 : nlohmann::json(nullptr);
			return result;
		}

		[[nodiscard]] nlohmann::json makeFmcwMetadataJson(const ContextPacket& packet)
		{
			return {{"present", packet.fmcw.present},
					{"waveform_shape", packet.fmcw.waveform_shape},
					{"chirp_bandwidth_hz", packet.fmcw.chirp_bandwidth},
					{"chirp_duration_s", packet.fmcw.chirp_duration},
					{"chirp_period_s", packet.fmcw.chirp_period},
					{"chirp_rate_hz_per_s", packet.fmcw.chirp_rate},
					{"chirp_rate_signed_hz_per_s", packet.fmcw.chirp_rate_signed},
					{"sweep_direction", packet.fmcw.sweep_direction},
					{"start_frequency_offset_hz", packet.fmcw.start_frequency_offset},
					{"triangle_period_s", packet.fmcw.triangle_period},
					{"chirp_count", packet.fmcw.chirp_count},
					{"triangle_count", packet.fmcw.triangle_count},
					{"dechirp_mode", packet.fmcw.dechirp_mode},
					{"dechirp_reference_source", packet.fmcw.dechirp_reference_source},
					{"dechirp_reference_transmitter_id", packet.fmcw.dechirp_reference_transmitter_id},
					{"dechirp_reference_transmitter_name", packet.fmcw.dechirp_reference_transmitter_name},
					{"dechirp_reference_waveform_id", packet.fmcw.dechirp_reference_waveform_id},
					{"dechirp_reference_waveform_name", packet.fmcw.dechirp_reference_waveform_name}};
		}

		[[nodiscard]] nlohmann::json makeWaveformMetadataJson(const ContextPacket& packet)
		{
			if (packet.receiver_mode == "fmcw")
			{
				return {{"kind", "fmcw"}, {"metadata_ref", "fmcw"}};
			}
			if (packet.receiver_mode == "pulsed")
			{
				return {{"kind", "pulsed"}, {"metadata_ref", "pulsed"}};
			}
			if (packet.receiver_mode == "cw")
			{
				return {{"kind", "cw"}, {"metadata_ref", "cw"}};
			}
			if (packet.fmcw.present)
			{
				return {{"kind", "fmcw"}, {"metadata_ref", "fmcw"}};
			}
			if (packet.pulsed.present)
			{
				return {{"kind", "pulsed"}, {"metadata_ref", "pulsed"}};
			}
			if (packet.cw.present)
			{
				return {{"kind", "cw"}, {"metadata_ref", "cw"}};
			}
			return {{"kind", packet.receiver_mode.empty() ? "unknown" : packet.receiver_mode}};
		}

		[[nodiscard]] bool hasKnownReceiverMode(const ContextPacket& packet) noexcept
		{
			return packet.receiver_mode == "fmcw" || packet.receiver_mode == "pulsed" || packet.receiver_mode == "cw";
		}

		[[nodiscard]] nlohmann::json makeContextMetadataJson(const ContextPacket& packet)
		{
			nlohmann::json metadata{{"schema", "fers-vita49-context-v1"},
									{"simulation_name", packet.simulation_name},
									{"receiver",
									 {{"id", packet.receiver_id},
									  {"name", packet.receiver_name},
									  {"mode", packet.receiver_mode},
									  {"adc_bits", packet.adc_bits},
									  {"context_flags", packet.context_flags}}},
									{"coordinate_frame",
									 {{"frame", packet.coordinate.frame},
									  {"origin",
									   {{"latitude", packet.coordinate.origin_latitude},
										{"longitude", packet.coordinate.origin_longitude},
										{"altitude", packet.coordinate.origin_altitude}}},
									  {"utm_zone", packet.coordinate.utm_zone},
									  {"utm_north_hemisphere", packet.coordinate.utm_north_hemisphere}}},
									{"initial_platform_state",
									 {{"platform_id", packet.initial_platform_state.platform_id},
									  {"platform_name", packet.initial_platform_state.platform_name},
									  {"position_m",
									   {{"x", packet.initial_platform_state.position_x},
										{"y", packet.initial_platform_state.position_y},
										{"z", packet.initial_platform_state.position_z}}},
									  {"velocity_mps",
									   {{"x", packet.initial_platform_state.velocity_x},
										{"y", packet.initial_platform_state.velocity_y},
										{"z", packet.initial_platform_state.velocity_z}}},
									  {"rotation_rad",
									   {{"azimuth", packet.initial_platform_state.azimuth},
										{"elevation", packet.initial_platform_state.elevation}}}}},
									{"waveform", makeWaveformMetadataJson(packet)}};
			const bool known_mode = hasKnownReceiverMode(packet);
			if (packet.receiver_mode == "pulsed" || (!known_mode && packet.pulsed.present))
			{
				metadata["pulsed"] = makePulsedMetadataJson(packet);
			}
			if (packet.receiver_mode == "cw" || (!known_mode && packet.cw.present))
			{
				metadata["cw"] = makeCwMetadataJson(packet);
			}
			if (packet.receiver_mode == "fmcw" || (!known_mode && packet.fmcw.present))
			{
				metadata["fmcw"] = makeFmcwMetadataJson(packet);
			}
			return metadata;
		}
	}

	ByteWriter::ByteWriter(const std::size_t reserve_bytes)
	{
		if (reserve_bytes > 0u)
		{
			_bytes.reserve(reserve_bytes);
		}
	}

	void ByteWriter::writeU16(const std::uint16_t value)
	{
		_bytes.push_back(static_cast<std::uint8_t>((value >> 8u) & 0xFFu));
		_bytes.push_back(static_cast<std::uint8_t>(value & 0xFFu));
	}

	void ByteWriter::writeI16(const std::int16_t value) { writeU16(static_cast<std::uint16_t>(value)); }

	void ByteWriter::writeU32(const std::uint32_t value)
	{
		_bytes.push_back(static_cast<std::uint8_t>((value >> 24u) & 0xFFu));
		_bytes.push_back(static_cast<std::uint8_t>((value >> 16u) & 0xFFu));
		_bytes.push_back(static_cast<std::uint8_t>((value >> 8u) & 0xFFu));
		_bytes.push_back(static_cast<std::uint8_t>(value & 0xFFu));
	}

	void ByteWriter::writeU64(const std::uint64_t value)
	{
		writeU32(static_cast<std::uint32_t>((value >> 32u) & 0xFFFFFFFFull));
		writeU32(static_cast<std::uint32_t>(value & 0xFFFFFFFFull));
	}

	void ByteWriter::writeF64(const RealType value)
	{
		static_assert(sizeof(RealType) == sizeof(std::uint64_t));
		static_assert(std::numeric_limits<RealType>::is_iec559,
					  "VITA F64 context serialization requires IEEE 754 binary64 RealType");
		// VRT context doubles are serialized from their IEEE 754 bit pattern as a
		// big-endian unsigned integer. This assumes the host uses a conventional
		// IEEE representation for RealType.
		const auto bits = std::bit_cast<std::uint64_t>(value);
		writeU64(bits);
	}

	void ByteWriter::writeAsciiMetadata(const std::string_view value)
	{
		if (value.size() >= std::numeric_limits<std::uint32_t>::max())
		{
			throw std::length_error("VITA ASCII metadata field too large");
		}
		for (const auto ch : value)
		{
			if (static_cast<unsigned char>(ch) > 0x7Fu)
			{
				throw std::invalid_argument("VITA ASCII metadata must contain ASCII bytes only");
			}
		}
		_bytes.insert(_bytes.end(), value.begin(), value.end());
		_bytes.push_back(0);
		while (_bytes.size() % sizeof(std::uint32_t) != 0u)
		{
			_bytes.push_back(0);
		}
	}

	void ByteWriter::writeBytes(const std::span<const std::uint8_t> bytes)
	{
		_bytes.insert(_bytes.end(), bytes.begin(), bytes.end());
	}

	const std::vector<std::uint8_t>& ByteWriter::bytes() const noexcept { return _bytes; }

	std::vector<std::uint8_t> ByteWriter::takeBytes() noexcept { return std::move(_bytes); }

	std::vector<std::uint8_t> Vita49Serializer::serializeSignalData(const SignalDataPacket& packet)
	{
		if (packet.iq_interleaved.size() % 2u != 0u)
		{
			throw std::invalid_argument("VITA signal IQ payload must contain I/Q pairs");
		}

		const std::size_t byte_count = kSignalDataFixedBytes + packet.iq_interleaved.size() * sizeof(std::int16_t);
		const auto packet_size_words = checkedWordCount(byte_count);

		ByteWriter writer(byte_count);
		writer.writeU32(makeHeader(PacketType::SignalDataWithStreamId, true, true, IntegerTimestampMode::Utc,
								   FractionalTimestampMode::RealTimePicoseconds, packet.packet_count,
								   packet_size_words));
		writer.writeU32(packet.stream_id);
		writer.writeU64(packet.class_id);
		writer.writeU32(packet.timestamp.integer_seconds);
		writer.writeU64(packet.timestamp.fractional_picoseconds);
		for (const auto item : packet.iq_interleaved)
		{
			writer.writeI16(item);
		}
		writer.writeU32(packet.trailer);

		if (writer.bytes().size() != byte_count)
		{
			throw std::logic_error("VITA signal packet byte count mismatch");
		}
		return writer.takeBytes();
	}

	SignalDataSerializationResult
	Vita49Serializer::serializeSignalDataFixedFullscale(const FixedFullscaleSignalDataPacket& packet)
	{
		if (!std::isfinite(packet.fullscale) || packet.fullscale <= 0.0)
		{
			throw std::invalid_argument("VITA signal full-scale must be positive and finite");
		}

		const std::size_t byte_count = kSignalDataFixedBytes + packet.samples.size() * sizeof(std::int16_t) * 2u;
		const auto packet_size_words = checkedWordCount(byte_count);

		ByteWriter writer(byte_count);
		writer.writeU32(makeHeader(PacketType::SignalDataWithStreamId, true, true, IntegerTimestampMode::Utc,
								   FractionalTimestampMode::RealTimePicoseconds, packet.packet_count,
								   packet_size_words));
		writer.writeU32(packet.stream_id);
		writer.writeU64(packet.class_id);
		writer.writeU32(packet.timestamp.integer_seconds);
		writer.writeU64(packet.timestamp.fractional_picoseconds);

		std::uint64_t clipped_sample_count = 0;
		for (const auto& sample : packet.samples)
		{
			bool clipped = false;
			writer.writeI16(scaleComponentToInt16(sample.real(), packet.fullscale, clipped));
			writer.writeI16(scaleComponentToInt16(sample.imag(), packet.fullscale, clipped));
			if (clipped)
			{
				++clipped_sample_count;
			}
		}
		writer.writeU32(makeTrailer(packet.valid_data, packet.calibrated_time, packet.reference_lock,
									clipped_sample_count > 0, packet.sample_loss));

		if (writer.bytes().size() != byte_count)
		{
			throw std::logic_error("VITA direct signal packet byte count mismatch");
		}
		return SignalDataSerializationResult{.bytes = writer.takeBytes(), .clipped_sample_count = clipped_sample_count};
	}

	std::vector<std::uint8_t> Vita49Serializer::serializeContext(const ContextPacket& packet)
	{
		if (packet.cif0 != kFersContextCif0)
		{
			throw std::invalid_argument(
				"Unsupported VITA 49.2 context CIF0: FERS serializer only supports kFersContextCif0");
		}

		ByteWriter payload;
		payload.writeU32(packet.cif0);
		payload.writeU32(packet.state_indicators);
		payload.writeU64(packet.payload_format);
		payload.writeF64(packet.sample_rate);
		payload.writeF64(packet.reference_frequency);
		payload.writeF64(packet.if_offset);
		payload.writeF64(packet.bandwidth);
		payload.writeF64(packet.adc_fullscale);
		payload.writeU64(packet.receiver_id);
		const auto metadata = makeContextMetadataJson(packet).dump(-1, ' ', true);
		payload.writeAsciiMetadata(metadata);

		const std::size_t byte_count = 4u + 4u + 8u + 4u + 8u + payload.bytes().size();
		const auto packet_size_words = checkedWordCount(byte_count);

		ByteWriter writer(byte_count);
		writer.writeU32(makeHeader(PacketType::Context, true, false, IntegerTimestampMode::Utc,
								   FractionalTimestampMode::RealTimePicoseconds, packet.packet_count,
								   packet_size_words));
		writer.writeU32(packet.stream_id);
		writer.writeU64(packet.class_id);
		writer.writeU32(packet.timestamp.integer_seconds);
		writer.writeU64(packet.timestamp.fractional_picoseconds);
		writer.writeBytes(payload.bytes());
		return writer.takeBytes();
	}
}
