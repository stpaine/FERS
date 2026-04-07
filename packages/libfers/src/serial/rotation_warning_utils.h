// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "core/config.h"
#include "core/parameters.h"

namespace serial::rotation_warning_utils
{
	enum class ValueKind
	{
		Angle,
		Rate
	};

	enum class Confidence
	{
		None,
		Low,
		Medium,
		High
	};

	enum class WarningSensitivity
	{
		HighConfidence,
		MediumOrHigh,
		Aggressive
	};

	struct InferenceResult
	{
		params::RotationAngleUnit inferred_unit = params::RotationAngleUnit::Degrees;
		Confidence confidence = Confidence::None;
		int degree_score = 0;
		int radian_score = 0;
	};

	inline constexpr WarningSensitivity kWarningSensitivity = WarningSensitivity::MediumOrHigh;

	InferenceResult infer_unit_from_value(RealType value, ValueKind kind) noexcept;

	bool should_warn(Confidence confidence, WarningSensitivity sensitivity) noexcept;

	void clear_captured_warnings() noexcept;

	std::vector<std::string> take_captured_warnings();

	void maybe_warn_about_rotation_value(RealType value, params::RotationAngleUnit declared_unit, ValueKind kind,
										 std::string_view source, std::string_view owner, std::string_view field);
}
