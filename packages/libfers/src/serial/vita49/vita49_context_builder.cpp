// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "serial/vita49/vita49_context_builder.h"

#include <string_view>

namespace serial::vita49
{
	namespace
	{
		[[nodiscard]] bool modeAllowsMetadata(const std::string& receiver_mode, const std::string_view metadata_mode)
		{
			return receiver_mode.empty() || receiver_mode == metadata_mode;
		}
	}

	ContextPacket Vita49ContextBuilder::build(const ContextBuildRequest& request)
	{
		std::uint32_t flags = 0;
		if (request.stream.dechirped)
		{
			flags |= ContextFlagDechirped;
		}
		if (request.stream.if_resampled)
		{
			flags |= ContextFlagIfResampled;
		}
		if (request.sample_loss)
		{
			flags |= ContextFlagSampleLoss;
		}
		if (request.stream_open)
		{
			flags |= ContextFlagStreamOpen;
		}
		if (request.stream_close)
		{
			flags |= ContextFlagStreamClose;
		}
		if (request.stream.fmcw.present && modeAllowsMetadata(request.stream.mode, "fmcw"))
		{
			flags |= ContextFlagFmcwMetadataPresent;
		}
		if (request.stream.cw.present && modeAllowsMetadata(request.stream.mode, "cw"))
		{
			flags |= ContextFlagCwMetadataPresent;
		}
		if (request.stream.pulsed.present && modeAllowsMetadata(request.stream.mode, "pulsed"))
		{
			flags |= ContextFlagPulsedMetadataPresent;
		}
		if (request.stream.sfcw.present && modeAllowsMetadata(request.stream.mode, "sfcw"))
		{
			flags |= ContextFlagSfcwMetadataPresent;
		}

		return ContextPacket{
			.stream_id = request.stream_id,
			.class_id = kFersVrtIqClassId,
			.timestamp = request.timestamp,
			.packet_count = request.packet_count,
			.cif0 = kFersContextCif0,
			.state_indicators =
				makeContextStateIndicators(request.valid_data, request.calibrated_time, request.reference_lock,
										   request.over_range, request.sample_loss),
			.payload_format = makeComplexInt16PayloadFormat(),
			.sample_rate = request.stream.sample_rate,
			.reference_frequency = request.stream.reference_frequency,
			.if_offset = request.stream.if_offset,
			.bandwidth = request.stream.bandwidth,
			.adc_fullscale = request.adc_fullscale,
			.receiver_id = request.stream.receiver_id,
			.adc_bits = request.stream.adc_bits,
			.context_flags = flags,
			.receiver_name = request.stream.receiver_name,
			.simulation_name = request.simulation_name,
			.receiver_mode = request.stream.mode,
			.coordinate = request.stream.coordinate,
			.initial_platform_state = request.stream.initial_platform_state,
			.pulsed = request.stream.pulsed,
			.cw = request.stream.cw,
			.fmcw = request.stream.fmcw,
			.sfcw = request.stream.sfcw,
		};
	}
}
