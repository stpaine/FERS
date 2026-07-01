// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "core/config.h"
#include "core/output_config.h"
#include "core/sim_id.h"
#include "core/vita49_timestamp.h"

namespace core
{
	struct OutputFileMetadata;

	struct ReceiverStreamDescriptor
	{
		struct CoordinateContext
		{
			std::string frame = "ENU";
			RealType origin_latitude = 0.0;
			RealType origin_longitude = 0.0;
			RealType origin_altitude = 0.0;
			std::int32_t utm_zone = 0;
			bool utm_north_hemisphere = true;
		};

		struct PlatformState
		{
			SimId platform_id = 0;
			std::string platform_name = {};
			RealType position_x = 0.0;
			RealType position_y = 0.0;
			RealType position_z = 0.0;
			RealType velocity_x = 0.0;
			RealType velocity_y = 0.0;
			RealType velocity_z = 0.0;
			RealType azimuth = 0.0;
			RealType elevation = 0.0;
		};

		struct FmcwContext
		{
			bool present = false;
			std::string waveform_shape = {};
			RealType chirp_bandwidth = 0.0;
			RealType chirp_duration = 0.0;
			RealType chirp_period = 0.0;
			RealType chirp_rate = 0.0;
			RealType chirp_rate_signed = 0.0;
			std::string sweep_direction = {};
			RealType start_frequency_offset = 0.0;
			std::optional<RealType> triangle_period = std::nullopt;
			std::optional<std::uint64_t> chirp_count = std::nullopt;
			std::optional<std::uint64_t> triangle_count = std::nullopt;
			std::string dechirp_mode = "none";
			std::string dechirp_reference_source = "none";
			SimId dechirp_reference_transmitter_id = 0;
			std::string dechirp_reference_transmitter_name = {};
			SimId dechirp_reference_waveform_id = 0;
			std::string dechirp_reference_waveform_name = {};
		};

		struct PulsedContext
		{
			bool present = false;
			RealType window_length = 0.0;
			RealType window_prf = 0.0;
			RealType window_skip = 0.0;
			std::uint64_t window_count = 0;
			SimId waveform_id = 0;
			std::string waveform_name = {};
			RealType carrier_frequency = 0.0;
			RealType power = 0.0;
			RealType pulse_width = 0.0;
			RealType native_sample_rate = 0.0;
			std::uint64_t native_sample_count = 0;
		};

		struct CwContext
		{
			bool present = false;
			SimId waveform_id = 0;
			std::string waveform_name = {};
			RealType carrier_frequency = 0.0;
			RealType power = 0.0;
		};

		struct SfcwContext
		{
			bool present = false;
			SimId waveform_id = 0;
			std::string waveform_name = {};
			RealType carrier_frequency = 0.0;
			RealType power = 0.0;
			RealType start_frequency_offset = 0.0;
			RealType step_size = 0.0;
			std::uint64_t step_count = 0;
			RealType dwell_time = 0.0;
			RealType step_period = 0.0;
			RealType sweep_period = 0.0;
			std::optional<std::uint64_t> sweep_count = std::nullopt;
			RealType first_frequency = 0.0;
			RealType last_frequency = 0.0;
			RealType frequency_span = 0.0;
			RealType effective_bandwidth = 0.0;
			RealType range_resolution = 0.0;
			RealType unambiguous_range = 0.0;
		};

		SimId receiver_id = 0;
		std::string receiver_name = {};
		std::string mode = {};
		RealType sample_rate = 0.0;
		RealType reference_frequency = 0.0;
		RealType if_offset = 0.0;
		RealType bandwidth = 0.0;
		bool dechirped = false;
		bool if_resampled = false;
		unsigned adc_bits = 0;
		CoordinateContext coordinate = {};
		PlatformState initial_platform_state = {};
		PulsedContext pulsed = {};
		CwContext cw = {};
		FmcwContext fmcw = {};
		SfcwContext sfcw = {};
	};

	struct ReceiverSampleBlock
	{
		ReceiverStreamDescriptor stream;
		RealType first_sample_time = 0.0;
		RealType sample_rate = 0.0;
		std::span<const ComplexType> samples;
		std::uint64_t sample_start = 0;
		bool valid_data = true;
		bool calibrated_time = true;
		bool reference_lock = true;
		std::shared_ptr<const OutputFileMetadata> file_metadata = nullptr;
	};

	struct ReceiverStreamStats
	{
		SimId receiver_id = 0;
		std::string receiver_name;
		std::uint32_t stream_id = 0;
		std::string mode = "unknown";
		RealType sample_rate = 0.0;
		RealType reference_frequency = 0.0;
		std::uint64_t packets_emitted = 0;
		std::uint64_t context_packets = 0;
		std::uint64_t samples_emitted = 0;
		std::uint64_t packets_dropped = 0;
		std::uint64_t samples_dropped = 0;
		std::uint64_t over_range_count = 0;
		std::uint64_t late_packet_count = 0;
		std::optional<RealType> first_sample_time = std::nullopt;
		std::optional<RealType> end_sample_time = std::nullopt;
		std::optional<Vita49Timestamp> first_timestamp = std::nullopt;
		std::optional<Vita49Timestamp> end_timestamp = std::nullopt;
	};

	struct OutputStats
	{
		OutputMode mode = OutputMode::Hdf5;
		std::optional<std::uint64_t> epoch_unix_nanoseconds = std::nullopt;
		std::vector<ReceiverStreamStats> streams;
	};

	struct ReceiverOutputPacketTrace
	{
		std::uint64_t sequence = 0;
		std::string event;
		std::uint32_t stream_id = 0;
		std::uint64_t byte_count = 0;
		std::uint64_t sample_count = 0;
		RealType first_sample_time = 0.0;
		std::optional<Vita49Timestamp> timestamp = std::nullopt;
		bool data_packet = false;
		bool context_packet = false;
		bool dropped = false;
		bool over_range = false;
		bool sample_loss = false;
	};

	using ReceiverOutputTelemetryCallback =
		std::function<void(const std::optional<OutputStats>&, std::span<const ReceiverOutputPacketTrace>)>;

	class ReceiverOutputSink
	{
	public:
		virtual ~ReceiverOutputSink() = default;

		virtual void initializeRun(const OutputConfig& config, std::string simulation_name) = 0;
		virtual std::uint32_t registerStream(const ReceiverStreamDescriptor& stream) = 0;
		virtual void openStream(std::uint32_t stream_id, RealType first_sample_time) = 0;
		virtual void submitBlock(const ReceiverSampleBlock& block) = 0;
		virtual void emitContextHeartbeat(RealType simulation_time) = 0;
		virtual void closeStream(std::uint32_t stream_id) = 0;
		virtual OutputStats finalize() = 0;
		[[nodiscard]] virtual OutputStats snapshotStats() const { return {}; }
	};
}
