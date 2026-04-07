// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file xml_serializer_utils.h
 * @brief Core utility layer for serializing FERS XML scenario files.
 *
 * This file provides the internal mechanisms and data structures required to serialize
 * individual simulation objects into their corresponding XML elements. It defines
 * a context-driven serialization approach, converting global simulation state into
 * textual XML data.
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
	/**
	 * @brief Adds a child element with the specified text content.
	 * @param parent The parent XML element.
	 * @param name The name of the child element to create.
	 * @param text The text content to set for the child element.
	 */
	void addChildWithText(const XmlElement& parent, const std::string& name, const std::string& text);

	/**
	 * @brief Adds a child element with the specified numeric content.
	 * @tparam T The numeric type (automatically deduced).
	 * @param parent The parent XML element.
	 * @param name The name of the child element to create.
	 * @param value The numeric value to set for the child element.
	 */
	template <typename T>
	void addChildWithNumber(const XmlElement& parent, const std::string& name, T value)
	{
		if constexpr (std::is_floating_point_v<T>)
		{
			std::array<char, 64> buffer{};
			if (auto [ptr, ec] = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value); ec == std::errc())
			{
				const auto length = static_cast<std::string::size_type>(ptr - buffer.data());
				addChildWithText(parent, name, std::string(buffer.data(), length));
			}
			else
			{
				// This only executes if std::to_chars(...) returns ec != std::errc()
				// which will practically never happen unless the output buffer is too small
				// to hold the formatted value, or the standard library implementation fails
				// to support/format the given floating-point value.
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

	/**
	 * @brief Sets a boolean attribute on an XML element.
	 * @param element The XML element to modify.
	 * @param name The name of the attribute to set.
	 * @param value The boolean value to set.
	 */
	void setAttributeFromBool(const XmlElement& element, const std::string& name, bool value);

	/**
	 * @brief Serializes a schedule (active periods) into a parent XML element.
	 * @param schedule The schedule periods to serialize.
	 * @param parent The parent XML element.
	 */
	void serializeSchedule(const std::vector<radar::SchedulePeriod>& schedule, const XmlElement& parent);

	/**
	 * @brief Serializes a Parameters object into a parent XML element.
	 * @param parent The parent XML element.
	 * @param p The parameter struct to serialize.
	 */
	void serializeParameters(const XmlElement& parent, const params::Parameters& p);

	/**
	 * @brief Serializes a waveform into a parent XML element.
	 * @param waveform The waveform object to serialize.
	 * @param parent The parent XML element.
	 */
	void serializeWaveform(const fers_signal::RadarSignal& waveform, const XmlElement& parent);

	/**
	 * @brief Serializes a timing object into a parent XML element.
	 * @param timing The timing object to serialize.
	 * @param parent The parent XML element.
	 */
	void serializeTiming(const timing::PrototypeTiming& timing, const XmlElement& parent);

	/**
	 * @brief Serializes an antenna into a parent XML element.
	 * @param antenna The antenna object to serialize.
	 * @param parent The parent XML element.
	 */
	void serializeAntenna(const antenna::Antenna& antenna, const XmlElement& parent);

	/**
	 * @brief Serializes a motion path into a parent XML element.
	 * @param path The motion path object to serialize.
	 * @param parent The parent XML element.
	 */
	void serializeMotionPath(const math::Path& path, const XmlElement& parent);

	/**
	 * @brief Serializes a rotation path into a parent XML element.
	 * @param rotPath The rotation path object to serialize.
	 * @param parent The parent XML element.
	 */
	void serializeRotation(const math::RotationPath& rotPath, const XmlElement& parent);

	/**
	 * @brief Serializes a transmitter into a parent XML element.
	 * @param tx The transmitter object to serialize.
	 * @param parent The parent XML element.
	 */
	void serializeTransmitter(const radar::Transmitter& tx, const XmlElement& parent);

	/**
	 * @brief Serializes a receiver into a parent XML element.
	 * @param rx The receiver object to serialize.
	 * @param parent The parent XML element.
	 */
	void serializeReceiver(const radar::Receiver& rx, const XmlElement& parent);

	/**
	 * @brief Serializes a monostatic radar setup containing both a transmitter and receiver.
	 * @param tx The transmitter object.
	 * @param rx The receiver object.
	 * @param parent The parent XML element.
	 */
	void serializeMonostatic(const radar::Transmitter& tx, const radar::Receiver& rx, const XmlElement& parent);

	/**
	 * @brief Serializes a target into a parent XML element.
	 * @param target The target object to serialize.
	 * @param parent The parent XML element.
	 */
	void serializeTarget(const radar::Target& target, const XmlElement& parent);

	/**
	 * @brief Serializes a platform and its attached components into a parent XML element.
	 * @param platform The platform object to serialize.
	 * @param world The simulation world containing global state.
	 * @param parent The parent XML element.
	 */
	void serializePlatform(const radar::Platform& platform, const core::World& world, const XmlElement& parent);
}
