// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#pragma once

#include <cstdint>
#include <optional>
#include <string>

#include "core/config.h"

namespace core
{
	enum class OutputMode : std::uint8_t
	{
		Hdf5,
		Vita49Udp
	};

	struct Vita49OutputConfig
	{
		std::string host;
		std::uint16_t port = 0;
		RealType adc_fullscale = 0.0;
		std::uint32_t queue_depth = 1024;
		std::optional<std::uint64_t> epoch_unix_nanoseconds = std::nullopt;
		std::uint16_t max_udp_payload = 1400;
		bool packet_trace_enabled = true;
	};

	struct OutputConfig
	{
		OutputMode mode = OutputMode::Hdf5;
		Vita49OutputConfig vita49;
	};

	[[nodiscard]] inline bool isVita49Enabled(const OutputConfig& config) noexcept
	{
		return config.mode == OutputMode::Vita49Udp;
	}
}
