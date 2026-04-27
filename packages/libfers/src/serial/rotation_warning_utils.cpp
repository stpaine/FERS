// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "rotation_warning_utils.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <limits>
#include <numbers>
#include <numeric>
#include <ranges>
#include <string>

#include "core/logging.h"

namespace serial::rotation_warning_utils
{
	namespace
	{
		thread_local std::vector<std::string> captured_warnings; ///< Thread-local warning buffer for API callers.

		constexpr RealType kEpsilon = 1e-9; ///< Near-zero tolerance for rotation heuristics.
		constexpr RealType kMatchTolerance = 1e-3; ///< Match tolerance for common-angle comparisons.
		constexpr RealType kIrrationalLookingTolerance = 5e-3; ///< Tolerance for non-round radian-looking values.

		/// Common degree values used by unit inference.
		constexpr std::array<RealType, 12> common_degree_angles = {0.0,	  30.0,	 45.0,	60.0,  90.0,  120.0,
																   135.0, 150.0, 180.0, 225.0, 270.0, 360.0};

		/// Common radian values used by unit inference.
		constexpr std::array<RealType, 9> common_radian_angles = {0.0,
																  std::numbers::pi_v<RealType> / 6.0,
																  std::numbers::pi_v<RealType> / 4.0,
																  std::numbers::pi_v<RealType> / 3.0,
																  std::numbers::pi_v<RealType> / 2.0,
																  2.0 * std::numbers::pi_v<RealType> / 3.0,
																  3.0 * std::numbers::pi_v<RealType> / 4.0,
																  5.0 * std::numbers::pi_v<RealType> / 6.0,
																  std::numbers::pi_v<RealType>};

		/// Returns true when two values are within the provided tolerance.
		[[nodiscard]] bool is_close(const RealType a, const RealType b,
									const RealType tolerance = kMatchTolerance) noexcept
		{
			return std::abs(a - b) <= tolerance;
		}

		/// Returns true when a value is nearly integral.
		[[nodiscard]] bool is_near_integer(const RealType value) noexcept
		{
			return is_close(value, std::round(value), 1e-6);
		}

		/// Returns true when a value is nearly a tenth increment.
		[[nodiscard]] bool is_near_tenth(const RealType value) noexcept
		{
			return is_close(value * 10.0, std::round(value * 10.0), kIrrationalLookingTolerance * 10.0);
		}

		/// Returns true when an absolute value matches a common-angle candidate.
		[[nodiscard]] bool matches_common_angle(const RealType abs_value,
												const std::ranges::input_range auto& candidates) noexcept
		{
			for (const RealType candidate : candidates)
			{
				if (is_close(abs_value, candidate))
				{
					return true;
				}
			}
			return false;
		}

		/// Returns true when an absolute value matches a simple rational multiple of pi.
		[[nodiscard]] bool matches_simple_pi_ratio(const RealType abs_value) noexcept
		{
			if (abs_value <= kEpsilon)
			{
				return false;
			}

			const RealType pi = std::numbers::pi_v<RealType>;
			for (int numerator = 1; numerator <= 12; ++numerator)
			{
				for (int denominator = 1; denominator <= 12; ++denominator)
				{
					if (std::gcd(numerator, denominator) != 1)
					{
						continue;
					}

					const RealType candidate =
						pi * static_cast<RealType>(numerator) / static_cast<RealType>(denominator);
					if (is_close(abs_value, candidate))
					{
						return true;
					}
				}
			}
			return false;
		}

		/// Finds the closest clean sine/cosine value for a radian input.
		[[nodiscard]] RealType best_clean_trig_distance(const RealType radians) noexcept
		{
			constexpr std::array<RealType, 5> clean_values = {0.0, 0.5, std::numbers::sqrt2_v<RealType> / 2.0,
															  std::numbers::sqrt3_v<RealType> / 2.0, 1.0};

			RealType best = std::numeric_limits<RealType>::max();
			for (const RealType value : {std::abs(std::sin(radians)), std::abs(std::cos(radians))})
			{
				for (const RealType clean : clean_values)
				{
					best = std::min(best, std::abs(value - clean));
				}
			}
			return best;
		}

		/// Returns the external token for a rotation angle unit.
		[[nodiscard]] std::string unit_token(const params::RotationAngleUnit unit)
		{
			return std::string(params::rotationAngleUnitToken(unit));
		}

		/// Returns a display token for an inference confidence.
		[[nodiscard]] std::string confidence_token(const Confidence confidence)
		{
			switch (confidence)
			{
			case Confidence::High:
				return "high";
			case Confidence::Medium:
				return "medium";
			case Confidence::Low:
				return "low";
			case Confidence::None:
			default:
				return "none";
			}
		}
	}

	InferenceResult infer_unit_from_value(const RealType value, const ValueKind kind) noexcept
	{
		const RealType abs_value = std::abs(value);
		const RealType pi = std::numbers::pi_v<RealType>;
		const RealType two_pi = 2.0 * pi;

		InferenceResult result;
		if (abs_value <= kEpsilon)
		{
			return result;
		}

		if (abs_value > 360.0)
		{
			result.degree_score += 8;
		}
		else if (abs_value > 180.0)
		{
			result.degree_score += 7;
		}
		else if (abs_value > two_pi)
		{
			result.degree_score += abs_value < 10.0 ? 3 : 6;
		}
		else if (abs_value > pi)
		{
			result.degree_score += 2;
		}

		const bool degree_famous = matches_common_angle(abs_value, common_degree_angles);
		const bool radian_famous = matches_common_angle(abs_value, common_radian_angles);
		if (degree_famous && !radian_famous)
		{
			result.degree_score += 3;
		}
		else if (radian_famous && !degree_famous)
		{
			result.radian_score += 3;
		}

		if (is_near_integer(abs_value))
		{
			const RealType rounded = std::round(abs_value);
			if ((std::fmod(rounded, 45.0) == 0.0) || (std::fmod(rounded, 30.0) == 0.0) ||
				(std::fmod(rounded, 15.0) == 0.0) || (std::fmod(rounded, 10.0) == 0.0) ||
				(std::fmod(rounded, 5.0) == 0.0))
			{
				result.degree_score += 2;
			}
			else
			{
				result.degree_score += 1;
			}
		}
		else if ((abs_value <= two_pi + 1.0) && !is_near_tenth(abs_value))
		{
			result.radian_score += 1;
		}

		if (matches_simple_pi_ratio(abs_value))
		{
			result.radian_score += 3;
		}

		if ((kind == ValueKind::Angle) && (abs_value <= two_pi + 0.5) &&
			(std::abs(result.degree_score - result.radian_score) <= 1))
		{
			const RealType degree_distance = best_clean_trig_distance(value * pi / 180.0);
			const RealType radian_distance = best_clean_trig_distance(value);
			if (degree_distance + 1e-4 < radian_distance)
			{
				++result.degree_score;
			}
			else if (radian_distance + 1e-4 < degree_distance)
			{
				++result.radian_score;
			}
		}

		if (result.degree_score == result.radian_score)
		{
			result.confidence = Confidence::None;
			return result;
		}

		result.inferred_unit = result.degree_score > result.radian_score ? params::RotationAngleUnit::Degrees
																		 : params::RotationAngleUnit::Radians;
		const int lead = std::abs(result.degree_score - result.radian_score);
		const int max_score = std::max(result.degree_score, result.radian_score);

		if ((lead >= 5) || ((max_score >= 7) && (lead >= 3)))
		{
			result.confidence = Confidence::High;
		}
		else if ((lead >= 3) || ((max_score >= 5) && (lead >= 2)))
		{
			result.confidence = Confidence::Medium;
		}
		else if ((lead >= 2) || (max_score >= 4))
		{
			result.confidence = Confidence::Low;
		}

		return result;
	}

	bool should_warn(const Confidence confidence, const WarningSensitivity sensitivity) noexcept
	{
		switch (sensitivity)
		{
		case WarningSensitivity::HighConfidence:
			return confidence == Confidence::High;
		case WarningSensitivity::MediumOrHigh:
			return (confidence == Confidence::Medium) || (confidence == Confidence::High);
		case WarningSensitivity::Aggressive:
			return confidence != Confidence::None;
		default:
			return false;
		}
	}

	void clear_captured_warnings() noexcept { captured_warnings.clear(); }

	std::vector<std::string> take_captured_warnings()
	{
		std::vector<std::string> warnings = std::move(captured_warnings);
		captured_warnings.clear();
		return warnings;
	}

	void maybe_warn_about_rotation_value(const RealType value, const params::RotationAngleUnit declared_unit,
										 const ValueKind kind, const std::string_view source,
										 const std::string_view owner, const std::string_view field)
	{
		const InferenceResult inference = infer_unit_from_value(value, kind);
		if ((inference.inferred_unit == declared_unit) || !should_warn(inference.confidence, kWarningSensitivity))
		{
			return;
		}

		const std::string message =
			std::format("{} rotation {} '{}' looks like {} but '{}' was declared (confidence: {}, value: {}). "
						"Change rotationangleunit or convert existing values.",
						source, owner, field, unit_token(inference.inferred_unit), unit_token(declared_unit),
						confidence_token(inference.confidence), value);

		if (std::ranges::find(captured_warnings, message) == captured_warnings.end())
		{
			captured_warnings.push_back(message);
		}

		LOG(logging::Level::WARNING, "{}", message);
	}
}
