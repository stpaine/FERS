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
	/// Kind of rotation value being inspected for unit mistakes.
	enum class ValueKind
	{
		Angle, ///< Absolute rotation angle.
		Rate ///< Rotation rate.
	};

	/// Confidence level for an inferred rotation unit.
	enum class Confidence
	{
		None, ///< No useful inference could be made.
		Low, ///< Weak evidence for the inferred unit.
		Medium, ///< Moderate evidence for the inferred unit.
		High ///< Strong evidence for the inferred unit.
	};

	/// Minimum confidence threshold for emitting rotation-unit warnings.
	enum class WarningSensitivity
	{
		HighConfidence, ///< Warn only on high-confidence mismatches.
		MediumOrHigh, ///< Warn on medium or high-confidence mismatches.
		Aggressive ///< Warn on any non-low mismatch.
	};

	/// Result of inferring likely units from a rotation value.
	struct InferenceResult
	{
		params::RotationAngleUnit inferred_unit = params::RotationAngleUnit::Degrees; ///< Likely unit for the value.
		Confidence confidence = Confidence::None; ///< Confidence for the inferred unit.
		int degree_score = 0; ///< Heuristic score for degree-like values.
		int radian_score = 0; ///< Heuristic score for radian-like values.
	};

	/// Default warning sensitivity used by parser and serializer warnings.
	inline constexpr WarningSensitivity kWarningSensitivity = WarningSensitivity::MediumOrHigh;

	/// Infers the likely unit of a rotation value.
	InferenceResult infer_unit_from_value(RealType value, ValueKind kind) noexcept;

	/// Returns true when a warning should be emitted for the confidence and sensitivity.
	bool should_warn(Confidence confidence, WarningSensitivity sensitivity) noexcept;

	/// Clears the thread-local captured rotation warnings.
	void clear_captured_warnings() noexcept;

	/// Returns and clears the thread-local captured rotation warnings.
	std::vector<std::string> take_captured_warnings();

	/// Emits or captures a warning when a rotation value likely uses the wrong unit.
	void maybe_warn_about_rotation_value(RealType value, params::RotationAngleUnit declared_unit, ValueKind kind,
										 std::string_view source, std::string_view owner, std::string_view field);
}
