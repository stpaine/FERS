// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2024-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file xml_serializer.cpp
 * @brief Implementation for serializing the simulation world to XML.
 *
 * This file contains the logic to traverse the in-memory C++ object representation
 * of a FERS simulation and build a corresponding XML document. The process involves
 * converting internal data representations (like angles in radians) back to the
 * user-facing format defined by the FERS XML schema (like compass degrees).
 */

#include "xml_serializer.h"

#include <iomanip>
#include <ranges>
#include <sstream>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/parameters.h"
#include "core/sim_id.h"
#include "core/world.h"
#include "libxml_wrapper.h"
#include "math/coord.h"
#include "math/path.h"
#include "math/rotation_path.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "signal/radar_signal.h"
#include "timing/prototype_timing.h"
#include "timing/timing.h"

namespace
{
	// --- Helper Functions ---

	void addChildWithText(const XmlElement& parent, const std::string& name, const std::string& text)
	{
		parent.addChild(name).setText(text);
	}

	template <typename T>
	void addChildWithNumber(const XmlElement& parent, const std::string& name, T value)
	{
		if constexpr (std::is_floating_point_v<T>)
		{
			// `std::to_chars` is used for floating-point types to ensure that the
			// serialization is locale-independent and maintains full precision.
			// This avoids issues where stream-based methods might be affected by
			// the system's locale or might truncate precision.
			std::array<char, 64> buffer{};
			if (auto [ptr, ec] = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value); ec == std::errc())
			{
				addChildWithText(parent, name, std::string(buffer.data(), ptr - buffer.data()));
			}
			else
			{
				// Fallback for the rare case that std::to_chars fails.
				std::stringstream ss;
				ss << std::setprecision(std::numeric_limits<T>::max_digits10) << value;
				addChildWithText(parent, name, ss.str());
			}
		}
		else
		{
			// For integer types, std::to_string is sufficient and clear.
			addChildWithText(parent, name, std::to_string(value));
		}
	}

	void setAttributeFromBool(const XmlElement& element, const std::string& name, const bool value)
	{
		element.setAttribute(name, value ? "true" : "false");
	}

	// --- Component Serialization Functions ---

	void serializeSchedule(const std::vector<radar::SchedulePeriod>& schedule, const XmlElement& parent)
	{
		if (schedule.empty())
		{
			return;
		}
		const XmlElement sched_elem = parent.addChild("schedule");
		for (const auto& period : schedule)
		{
			XmlElement p_elem = sched_elem.addChild("period");
			p_elem.setAttribute("start", std::to_string(period.start));
			p_elem.setAttribute("end", std::to_string(period.end));
		}
	}

	void serializeParameters(const XmlElement& parent)
	{
		addChildWithNumber(parent, "starttime", params::startTime());
		addChildWithNumber(parent, "endtime", params::endTime());
		addChildWithNumber(parent, "rate", params::rate());

		if (params::c() != params::Parameters::DEFAULT_C)
		{
			addChildWithNumber(parent, "c", params::c());
		}
		if (params::simSamplingRate() != 1000.0)
		{
			addChildWithNumber(parent, "simSamplingRate", params::simSamplingRate());
		}
		if (params::params.random_seed)
		{
			addChildWithNumber(parent, "randomseed", *params::params.random_seed);
		}
		if (params::adcBits() != 0)
		{
			addChildWithNumber(parent, "adc_bits", params::adcBits());
		}
		if (params::oversampleRatio() != 1)
		{
			addChildWithNumber(parent, "oversample", params::oversampleRatio());
		}

		const XmlElement origin = parent.addChild("origin");
		origin.setAttribute("latitude", std::to_string(params::originLatitude()));
		origin.setAttribute("longitude", std::to_string(params::originLongitude()));
		origin.setAttribute("altitude", std::to_string(params::originAltitude()));

		const XmlElement cs = parent.addChild("coordinatesystem");
		switch (params::coordinateFrame())
		{
		case params::CoordinateFrame::ENU:
			cs.setAttribute("frame", "ENU");
			break;
		case params::CoordinateFrame::UTM:
			cs.setAttribute("frame", "UTM");
			cs.setAttribute("zone", std::to_string(params::utmZone()));
			cs.setAttribute("hemisphere", params::utmNorthHemisphere() ? "N" : "S");
			break;
		case params::CoordinateFrame::ECEF:
			cs.setAttribute("frame", "ECEF");
			break;
		}
	}

	void serializeWaveform(const fers_signal::RadarSignal& waveform, const XmlElement& parent)
	{
		parent.setAttribute("id", std::to_string(waveform.getId()));
		parent.setAttribute("name", waveform.getName());

		addChildWithNumber(parent, "power", waveform.getPower());
		addChildWithNumber(parent, "carrier_frequency", waveform.getCarrier());

		if (dynamic_cast<const fers_signal::CwSignal*>(waveform.getSignal()))
		{
			(void)parent.addChild("cw"); // Empty element
		}
		else
		{
			const XmlElement pulsed_file = parent.addChild("pulsed_from_file");
			// Gracefully handle missing filename to prevent export crashes.
			// While the schema requires a filename, it is better to export an empty string
			// for work-in-progress scenarios than to crash the application.
			const auto& filename = waveform.getFilename();
			pulsed_file.setAttribute("filename", filename.value_or(""));
		}
	}

	void serializeTiming(const timing::PrototypeTiming& timing, const XmlElement& parent)
	{
		parent.setAttribute("id", std::to_string(timing.getId()));
		parent.setAttribute("name", timing.getName());
		setAttributeFromBool(parent, "synconpulse", timing.getSyncOnPulse());

		addChildWithNumber(parent, "frequency", timing.getFrequency());
		if (const auto val = timing.getFreqOffset())
		{
			addChildWithNumber(parent, "freq_offset", *val);
		}
		if (const auto val = timing.getRandomFreqOffsetStdev())
		{
			addChildWithNumber(parent, "random_freq_offset_stdev", *val);
		}
		if (const auto val = timing.getPhaseOffset())
		{
			addChildWithNumber(parent, "phase_offset", *val);
		}
		if (const auto val = timing.getRandomPhaseOffsetStdev())
		{
			addChildWithNumber(parent, "random_phase_offset_stdev", *val);
		}

		std::vector<RealType> alphas, weights;
		timing.copyAlphas(alphas, weights);
		for (size_t i = 0; i < alphas.size(); ++i)
		{
			XmlElement entry = parent.addChild("noise_entry");
			addChildWithNumber(entry, "alpha", alphas[i]);
			addChildWithNumber(entry, "weight", weights[i]);
		}
	}

	void serializeAntenna(const antenna::Antenna& antenna, const XmlElement& parent)
	{
		parent.setAttribute("id", std::to_string(antenna.getId()));
		parent.setAttribute("name", antenna.getName());

		if (const auto* sinc = dynamic_cast<const antenna::Sinc*>(&antenna))
		{
			parent.setAttribute("pattern", "sinc");
			addChildWithNumber(parent, "alpha", sinc->getAlpha());
			addChildWithNumber(parent, "beta", sinc->getBeta());
			addChildWithNumber(parent, "gamma", sinc->getGamma());
		}
		else if (const auto* gaussian = dynamic_cast<const antenna::Gaussian*>(&antenna))
		{
			parent.setAttribute("pattern", "gaussian");
			addChildWithNumber(parent, "azscale", gaussian->getAzimuthScale());
			addChildWithNumber(parent, "elscale", gaussian->getElevationScale());
		}
		else if (const auto* sh = dynamic_cast<const antenna::SquareHorn*>(&antenna))
		{
			parent.setAttribute("pattern", "squarehorn");
			addChildWithNumber(parent, "diameter", sh->getDimension());
		}
		else if (const auto* parabolic = dynamic_cast<const antenna::Parabolic*>(&antenna))
		{
			parent.setAttribute("pattern", "parabolic");
			addChildWithNumber(parent, "diameter", parabolic->getDiameter());
		}
		else if (const auto* xml_ant = dynamic_cast<const antenna::XmlAntenna*>(&antenna))
		{
			parent.setAttribute("pattern", "xml");
			parent.setAttribute("filename", xml_ant->getFilename());
		}
		else if (const auto* h5_ant = dynamic_cast<const antenna::H5Antenna*>(&antenna))
		{
			parent.setAttribute("pattern", "file");
			parent.setAttribute("filename", h5_ant->getFilename());
		}
		else
		{
			parent.setAttribute("pattern", "isotropic");
		}

		if (antenna.getEfficiencyFactor() != 1.0)
		{
			addChildWithNumber(parent, "efficiency", antenna.getEfficiencyFactor());
		}
	}

	void serializeMotionPath(const math::Path& path, const XmlElement& parent)
	{
		switch (path.getType())
		{
		case math::Path::InterpType::INTERP_STATIC:
			parent.setAttribute("interpolation", "static");
			break;
		case math::Path::InterpType::INTERP_LINEAR:
			parent.setAttribute("interpolation", "linear");
			break;
		case math::Path::InterpType::INTERP_CUBIC:
			parent.setAttribute("interpolation", "cubic");
			break;
		}

		for (const auto& [pos, t] : path.getCoords())
		{
			XmlElement wp_elem = parent.addChild("positionwaypoint");
			addChildWithNumber(wp_elem, "x", pos.x);
			addChildWithNumber(wp_elem, "y", pos.y);
			addChildWithNumber(wp_elem, "altitude", pos.z);
			addChildWithNumber(wp_elem, "time", t);
		}
	}

	void serializeRotation(const math::RotationPath& rotPath, const XmlElement& parent)
	{
		if (rotPath.getType() == math::RotationPath::InterpType::INTERP_CONSTANT)
		{
			const XmlElement fixed_elem = parent.addChild("fixedrotation");
			const auto start = rotPath.getStart();
			const auto rate = rotPath.getRate();

			// Convert internal mathematical angles (radians, CCW from East) back to
			// the XML format's compass degrees (CW from North). This transformation
			// is necessary at the serialization boundary to ensure the output XML
			// conforms to the user-facing FERS schema and is human-readable.
			const RealType start_az_deg = std::fmod(90.0 - start.azimuth * 180.0 / PI + 360.0, 360.0);
			const RealType start_el_deg = start.elevation * 180.0 / PI;
			const RealType rate_az_deg_s = -rate.azimuth * 180.0 / PI; // Invert for CW rotation
			const RealType rate_el_deg_s = rate.elevation * 180.0 / PI;

			addChildWithNumber(fixed_elem, "startazimuth", start_az_deg);
			addChildWithNumber(fixed_elem, "startelevation", start_el_deg);
			addChildWithNumber(fixed_elem, "azimuthrate", rate_az_deg_s);
			addChildWithNumber(fixed_elem, "elevationrate", rate_el_deg_s);
		}
		else
		{
			const XmlElement rot_elem = parent.addChild("rotationpath");
			switch (rotPath.getType())
			{
			case math::RotationPath::InterpType::INTERP_STATIC:
				rot_elem.setAttribute("interpolation", "static");
				break;
			case math::RotationPath::InterpType::INTERP_LINEAR:
				rot_elem.setAttribute("interpolation", "linear");
				break;
			case math::RotationPath::InterpType::INTERP_CUBIC:
				rot_elem.setAttribute("interpolation", "cubic");
				break;
			default:
				break; // Should not happen
			}
			for (const auto& wp : rotPath.getCoords())
			{
				XmlElement wp_elem = rot_elem.addChild("rotationwaypoint");
				// Convert angles back to compass degrees for XML output. This is
				// done to ensure the output conforms to the FERS XML schema.
				const RealType az_deg = std::fmod(90.0 - wp.azimuth * 180.0 / PI + 360.0, 360.0);
				const RealType el_deg = wp.elevation * 180.0 / PI;
				addChildWithNumber(wp_elem, "azimuth", az_deg);
				addChildWithNumber(wp_elem, "elevation", el_deg);
				addChildWithNumber(wp_elem, "time", wp.t);
			}
		}
	}

	void serializeTransmitter(const radar::Transmitter& tx, const XmlElement& parent)
	{
		const XmlElement tx_elem = parent.addChild("transmitter");
		tx_elem.setAttribute("id", std::to_string(tx.getId()));
		tx_elem.setAttribute("name", tx.getName());
		tx_elem.setAttribute("waveform", tx.getSignal() ? std::to_string(tx.getSignal()->getId()) : "");
		tx_elem.setAttribute("antenna", tx.getAntenna() ? std::to_string(tx.getAntenna()->getId()) : "");
		tx_elem.setAttribute("timing", tx.getTiming() ? std::to_string(tx.getTiming()->getId()) : "");

		if (tx.getMode() == radar::OperationMode::PULSED_MODE)
		{
			const XmlElement mode_elem = tx_elem.addChild("pulsed_mode");
			addChildWithNumber(mode_elem, "prf", tx.getPrf());
		}
		else
		{
			(void)tx_elem.addChild("cw_mode");
		}

		serializeSchedule(tx.getSchedule(), tx_elem);
	}

	void serializeReceiver(const radar::Receiver& rx, const XmlElement& parent)
	{
		const XmlElement rx_elem = parent.addChild("receiver");
		rx_elem.setAttribute("id", std::to_string(rx.getId()));
		rx_elem.setAttribute("name", rx.getName());
		rx_elem.setAttribute("antenna", rx.getAntenna() ? std::to_string(rx.getAntenna()->getId()) : "");
		rx_elem.setAttribute("timing", rx.getTiming() ? std::to_string(rx.getTiming()->getId()) : "");
		setAttributeFromBool(rx_elem, "nodirect", rx.checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT));
		setAttributeFromBool(rx_elem, "nopropagationloss", rx.checkFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS));

		if (rx.getMode() == radar::OperationMode::PULSED_MODE)
		{
			const XmlElement mode_elem = rx_elem.addChild("pulsed_mode");
			addChildWithNumber(mode_elem, "prf", rx.getWindowPrf());
			addChildWithNumber(mode_elem, "window_skip", rx.getWindowSkip());
			addChildWithNumber(mode_elem, "window_length", rx.getWindowLength());
		}
		else
		{
			(void)rx_elem.addChild("cw_mode");
		}

		if (rx.getNoiseTemperature() > 0)
		{
			addChildWithNumber(rx_elem, "noise_temp", rx.getNoiseTemperature());
		}

		serializeSchedule(rx.getSchedule(), rx_elem);
	}

	void serializeMonostatic(const radar::Transmitter& tx, const radar::Receiver& rx, const XmlElement& parent)
	{
		const XmlElement mono_elem = parent.addChild("monostatic");
		mono_elem.setAttribute("tx_id", std::to_string(tx.getId()));
		mono_elem.setAttribute("rx_id", std::to_string(rx.getId()));
		mono_elem.setAttribute("name", tx.getName());
		mono_elem.setAttribute("antenna", tx.getAntenna() ? std::to_string(tx.getAntenna()->getId()) : "");
		mono_elem.setAttribute("waveform", tx.getSignal() ? std::to_string(tx.getSignal()->getId()) : "");
		mono_elem.setAttribute("timing", tx.getTiming() ? std::to_string(tx.getTiming()->getId()) : "");
		setAttributeFromBool(mono_elem, "nodirect", rx.checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT));
		setAttributeFromBool(mono_elem, "nopropagationloss", rx.checkFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS));

		if (tx.getMode() == radar::OperationMode::PULSED_MODE)
		{
			const XmlElement mode_elem = mono_elem.addChild("pulsed_mode");
			addChildWithNumber(mode_elem, "prf", tx.getPrf());
			addChildWithNumber(mode_elem, "window_skip", rx.getWindowSkip());
			addChildWithNumber(mode_elem, "window_length", rx.getWindowLength());
		}
		else
		{
			(void)mono_elem.addChild("cw_mode");
		}

		if (rx.getNoiseTemperature() > 0)
		{
			addChildWithNumber(mono_elem, "noise_temp", rx.getNoiseTemperature());
		}

		serializeSchedule(tx.getSchedule(), mono_elem);
	}

	void serializeTarget(const radar::Target& target, const XmlElement& parent)
	{
		const XmlElement target_elem = parent.addChild("target");
		target_elem.setAttribute("id", std::to_string(target.getId()));
		target_elem.setAttribute("name", target.getName());

		const XmlElement rcs_elem = target_elem.addChild("rcs");
		if (const auto* iso = dynamic_cast<const radar::IsoTarget*>(&target))
		{
			rcs_elem.setAttribute("type", "isotropic");
			addChildWithNumber(rcs_elem, "value", iso->getConstRcs());
		}
		else if (const auto* file_target = dynamic_cast<const radar::FileTarget*>(&target))
		{
			rcs_elem.setAttribute("type", "file");
			rcs_elem.setAttribute("filename", file_target->getFilename());
		}

		// Serialize fluctuation model if present (e.g. Swerling/Chi-Square)
		if (const auto* model = target.getFluctuationModel())
		{
			if (const auto* chi = dynamic_cast<const radar::RcsChiSquare*>(model))
			{
				XmlElement model_elem = target_elem.addChild("model");
				model_elem.setAttribute("type", "chisquare");
				addChildWithNumber(model_elem, "k", chi->getK());
			}
		}
	}

	void serializePlatform(const radar::Platform& platform, const core::World& world, const XmlElement& parent)
	{
		parent.setAttribute("id", std::to_string(platform.getId()));
		parent.setAttribute("name", platform.getName());

		const XmlElement motion_elem = parent.addChild("motionpath");
		serializeMotionPath(*platform.getMotionPath(), motion_elem);

		serializeRotation(*platform.getRotationPath(), parent);

		// Transmitters (including Monostatic pairs)
		for (const auto& tx : world.getTransmitters())
		{
			if (tx->getPlatform() == &platform)
			{
				if (tx->getAttached())
				{
					// Serialize as <monostatic> if it is part of a pair
					serializeMonostatic(*tx, *dynamic_cast<const radar::Receiver*>(tx->getAttached()), parent);
				}
				else
				{
					// Serialize as <transmitter>
					serializeTransmitter(*tx, parent);
				}
			}
		}

		// Standalone Receivers
		for (const auto& rx : world.getReceivers())
		{
			// Only serialize receivers that are NOT attached to a transmitter,
			// as those were handled in the loop above.
			if (rx->getPlatform() == &platform && !rx->getAttached())
			{
				serializeReceiver(*rx, parent);
			}
		}

		// Targets
		for (const auto& target : world.getTargets())
		{
			if (target->getPlatform() == &platform)
			{
				serializeTarget(*target, parent);
			}
		}
	}
}

namespace serial
{
	std::string world_to_xml_string(const core::World& world)
	{
		XmlDocument doc;
		xmlNodePtr sim_node = xmlNewNode(nullptr, reinterpret_cast<const xmlChar*>("simulation"));
		XmlElement root(sim_node);
		doc.setRootElement(root);

		if (!params::params.simulation_name.empty())
		{
			root.setAttribute("name", params::params.simulation_name);
		}
		else
		{
			root.setAttribute("name", "FERS Scenario");
		}

		const XmlElement params_elem = root.addChild("parameters");
		serializeParameters(params_elem);

		// Assets (waveforms, timings, antennas) are serialized first. This is
		// necessary because platforms reference these assets by ID. By defining
		// them at the top of the document, we ensure that any XML parser can
		// resolve these references when it later encounters the platform definitions.
		for (const auto& waveform : world.getWaveforms() | std::views::values)
		{
			XmlElement waveform_elem = root.addChild("waveform");
			serializeWaveform(*waveform, waveform_elem);
		}
		for (const auto& timing : world.getTimings() | std::views::values)
		{
			XmlElement timing_elem = root.addChild("timing");
			serializeTiming(*timing, timing_elem);
		}
		for (const auto& antenna : world.getAntennas() | std::views::values)
		{
			XmlElement antenna_elem = root.addChild("antenna");
			serializeAntenna(*antenna, antenna_elem);
		}
		for (const auto& platform : world.getPlatforms())
		{
			XmlElement plat_elem = root.addChild("platform");
			serializePlatform(*platform, world, plat_elem);
		}

		return doc.dumpToString();
	}
}
