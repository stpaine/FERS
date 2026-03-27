// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file xml_serializer_utils.h
 * @brief Core utility layer for serializing FERS XML scenario files.
 */

#pragma once

#include <array>
#include <charconv>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "core/parameters.h"
#include "serial/libxml_wrapper.h"

namespace antenna
{
	class Antenna;
}
namespace core
{
	class World;
}
namespace fers_signal
{
	class RadarSignal;
}
namespace math
{
	class Path;
	class RotationPath;
}
namespace radar
{
	class Platform;
	class Receiver;
	class Target;
	class Transmitter;
	struct SchedulePeriod;
}
namespace timing
{
	class PrototypeTiming;
}

namespace serial::xml_serializer_utils
{
	void addChildWithText(const XmlElement& parent, const std::string& name, const std::string& text);

	template <typename T>
	void addChildWithNumber(const XmlElement& parent, const std::string& name, T value)
	{
		if constexpr (std::is_floating_point_v<T>)
		{
			std::array<char, 64> buffer{};
			if (auto [ptr, ec] = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value); ec == std::errc())
			{
				addChildWithText(parent, name, std::string(buffer.data(), ptr - buffer.data()));
			}
			else
			{
				std::stringstream ss;
				ss << std::setprecision(std::numeric_limits<T>::max_digits10) << value;
				addChildWithText(parent, name, ss.str());
			}
		}
		else
		{
			addChildWithText(parent, name, std::to_string(value));
		}
	}

	void setAttributeFromBool(const XmlElement& element, const std::string& name, bool value);

	void serializeSchedule(const std::vector<radar::SchedulePeriod>& schedule, const XmlElement& parent);
	void serializeParameters(const XmlElement& parent, const params::Parameters& p);
	void serializeWaveform(const fers_signal::RadarSignal& waveform, const XmlElement& parent);
	void serializeTiming(const timing::PrototypeTiming& timing, const XmlElement& parent);
	void serializeAntenna(const antenna::Antenna& antenna, const XmlElement& parent);
	void serializeMotionPath(const math::Path& path, const XmlElement& parent);
	void serializeRotation(const math::RotationPath& rotPath, const XmlElement& parent);
	void serializeTransmitter(const radar::Transmitter& tx, const XmlElement& parent);
	void serializeReceiver(const radar::Receiver& rx, const XmlElement& parent);
	void serializeMonostatic(const radar::Transmitter& tx, const radar::Receiver& rx, const XmlElement& parent);
	void serializeTarget(const radar::Target& target, const XmlElement& parent);
	void serializePlatform(const radar::Platform& platform, const core::World& world, const XmlElement& parent);
}
