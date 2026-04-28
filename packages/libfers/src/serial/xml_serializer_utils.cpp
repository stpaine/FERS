// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "xml_serializer_utils.h"

#include <cmath>

#include "antenna/antenna_factory.h"
#include "core/world.h"
#include "math/coord.h"
#include "math/path.h"
#include "math/rotation_path.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/rotation_angle_utils.h"
#include "signal/radar_signal.h"
#include "timing/prototype_timing.h"
#include "timing/timing.h"

namespace serial::xml_serializer_utils
{
	void addChildWithText(const XmlElement& parent, const std::string& name, const std::string& text)
	{
		parent.addChild(name).setText(text);
	}

	void setAttributeFromBool(const XmlElement& element, const std::string& name, const bool value)
	{
		element.setAttribute(name, value ? "true" : "false");
	}

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

	void serializeParameters(const XmlElement& parent, const params::Parameters& p)
	{
		addChildWithNumber(parent, "starttime", p.start);
		addChildWithNumber(parent, "endtime", p.end);
		addChildWithNumber(parent, "rate", p.rate);

		if (p.c != params::Parameters::DEFAULT_C)
		{
			addChildWithNumber(parent, "c", p.c);
		}
		if (p.sim_sampling_rate != 1000.0)
		{
			addChildWithNumber(parent, "simSamplingRate", p.sim_sampling_rate);
		}
		if (p.random_seed)
		{
			addChildWithNumber(parent, "randomseed", *p.random_seed);
		}
		if (p.adc_bits != 0)
		{
			addChildWithNumber(parent, "adc_bits", p.adc_bits);
		}
		if (p.oversample_ratio != 1)
		{
			addChildWithNumber(parent, "oversample", p.oversample_ratio);
		}
		if (p.rotation_angle_unit != params::RotationAngleUnit::Degrees)
		{
			addChildWithText(parent, "rotationangleunit",
							 std::string(params::rotationAngleUnitToken(p.rotation_angle_unit)));
		}

		const XmlElement origin = parent.addChild("origin");
		origin.setAttribute("latitude", std::to_string(p.origin_latitude));
		origin.setAttribute("longitude", std::to_string(p.origin_longitude));
		origin.setAttribute("altitude", std::to_string(p.origin_altitude));

		const XmlElement cs = parent.addChild("coordinatesystem");
		switch (p.coordinate_frame)
		{
		case params::CoordinateFrame::ENU:
			cs.setAttribute("frame", "ENU");
			break;
		case params::CoordinateFrame::UTM:
			cs.setAttribute("frame", "UTM");
			cs.setAttribute("zone", std::to_string(p.utm_zone));
			cs.setAttribute("hemisphere", p.utm_north_hemisphere ? "N" : "S");
			break;
		case params::CoordinateFrame::ECEF:
			cs.setAttribute("frame", "ECEF");
			break;
		}
	}

	void serializeWaveform(const fers_signal::RadarSignal& waveform, const XmlElement& parent)
	{
		parent.setAttribute("name", waveform.getName());

		addChildWithNumber(parent, "power", waveform.getPower());
		addChildWithNumber(parent, "carrier_frequency", waveform.getCarrier());

		if (dynamic_cast<const fers_signal::CwSignal*>(waveform.getSignal()) != nullptr)
		{
			(void)parent.addChild("cw"); // Empty element
		}
		else if (const auto* fmcw = waveform.getFmcwChirpSignal(); fmcw != nullptr)
		{
			const XmlElement fmcw_elem = parent.addChild("fmcw_linear_chirp");
			fmcw_elem.setAttribute("direction",
								   std::string(fers_signal::fmcwChirpDirectionToken(fmcw->getDirection())));
			addChildWithNumber(fmcw_elem, "chirp_bandwidth", fmcw->getChirpBandwidth());
			addChildWithNumber(fmcw_elem, "chirp_duration", fmcw->getChirpDuration());
			addChildWithNumber(fmcw_elem, "chirp_period", fmcw->getChirpPeriod());
			if (std::abs(fmcw->getStartFrequencyOffset()) > EPSILON)
			{
				addChildWithNumber(fmcw_elem, "start_frequency_offset", fmcw->getStartFrequencyOffset());
			}
			if (fmcw->getChirpCount().has_value())
			{
				addChildWithNumber(fmcw_elem, "chirp_count", static_cast<RealType>(*fmcw->getChirpCount()));
			}
		}
		else
		{
			const XmlElement pulsed_file = parent.addChild("pulsed_from_file");
			const auto& filename = waveform.getFilename();
			pulsed_file.setAttribute("filename", filename.value_or(""));
		}
	}

	void serializeTiming(const timing::PrototypeTiming& timing, const XmlElement& parent)
	{
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
			const auto unit = params::rotationAngleUnit();

			RealType start_az = rotation_angle_utils::internal_azimuth_to_external(start.azimuth, unit);
			if (unit == params::RotationAngleUnit::Degrees)
			{
				start_az = std::fmod(start_az + 360.0, 360.0);
			}
			const RealType start_el = rotation_angle_utils::internal_elevation_to_external(start.elevation, unit);
			const RealType rate_az = rotation_angle_utils::internal_azimuth_rate_to_external(rate.azimuth, unit);
			const RealType rate_el = rotation_angle_utils::internal_elevation_rate_to_external(rate.elevation, unit);

			addChildWithNumber(fixed_elem, "startazimuth", start_az);
			addChildWithNumber(fixed_elem, "startelevation", start_el);
			addChildWithNumber(fixed_elem, "azimuthrate", rate_az);
			addChildWithNumber(fixed_elem, "elevationrate", rate_el);
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
				break;
			}
			const auto unit = params::rotationAngleUnit();
			for (const auto& wp : rotPath.getCoords())
			{
				XmlElement wp_elem = rot_elem.addChild("rotationwaypoint");
				RealType azimuth = rotation_angle_utils::internal_azimuth_to_external(wp.azimuth, unit);
				if (unit == params::RotationAngleUnit::Degrees)
				{
					azimuth = std::fmod(azimuth + 360.0, 360.0);
				}
				const RealType elevation = rotation_angle_utils::internal_elevation_to_external(wp.elevation, unit);
				addChildWithNumber(wp_elem, "azimuth", azimuth);
				addChildWithNumber(wp_elem, "elevation", elevation);
				addChildWithNumber(wp_elem, "time", wp.t);
			}
		}
	}

	void serializeTransmitter(const radar::Transmitter& tx, const XmlElement& parent)
	{
		const XmlElement tx_elem = parent.addChild("transmitter");
		tx_elem.setAttribute("name", tx.getName());
		tx_elem.setAttribute("waveform", (tx.getSignal() != nullptr) ? tx.getSignal()->getName() : "");
		tx_elem.setAttribute("antenna", (tx.getAntenna() != nullptr) ? tx.getAntenna()->getName() : "");
		tx_elem.setAttribute("timing", tx.getTiming() ? tx.getTiming()->getName() : "");

		if (tx.getMode() == radar::OperationMode::PULSED_MODE)
		{
			const XmlElement mode_elem = tx_elem.addChild("pulsed_mode");
			addChildWithNumber(mode_elem, "prf", tx.getPrf());
		}
		else if (tx.getMode() == radar::OperationMode::FMCW_MODE)
		{
			(void)tx_elem.addChild("fmcw_mode");
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
		rx_elem.setAttribute("name", rx.getName());
		rx_elem.setAttribute("antenna", (rx.getAntenna() != nullptr) ? rx.getAntenna()->getName() : "");
		rx_elem.setAttribute("timing", rx.getTiming() ? rx.getTiming()->getName() : "");
		setAttributeFromBool(rx_elem, "nodirect", rx.checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT));
		setAttributeFromBool(rx_elem, "nopropagationloss", rx.checkFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS));

		if (rx.getMode() == radar::OperationMode::PULSED_MODE)
		{
			const XmlElement mode_elem = rx_elem.addChild("pulsed_mode");
			addChildWithNumber(mode_elem, "prf", rx.getWindowPrf());
			addChildWithNumber(mode_elem, "window_skip", rx.getWindowSkip());
			addChildWithNumber(mode_elem, "window_length", rx.getWindowLength());
		}
		else if (rx.getMode() == radar::OperationMode::FMCW_MODE)
		{
			(void)rx_elem.addChild("fmcw_mode");
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
		mono_elem.setAttribute("name", tx.getName());
		mono_elem.setAttribute("antenna", (tx.getAntenna() != nullptr) ? tx.getAntenna()->getName() : "");
		mono_elem.setAttribute("waveform", (tx.getSignal() != nullptr) ? tx.getSignal()->getName() : "");
		mono_elem.setAttribute("timing", tx.getTiming() ? tx.getTiming()->getName() : "");
		setAttributeFromBool(mono_elem, "nodirect", rx.checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT));
		setAttributeFromBool(mono_elem, "nopropagationloss", rx.checkFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS));

		if (tx.getMode() == radar::OperationMode::PULSED_MODE)
		{
			const XmlElement mode_elem = mono_elem.addChild("pulsed_mode");
			addChildWithNumber(mode_elem, "prf", tx.getPrf());
			addChildWithNumber(mode_elem, "window_skip", rx.getWindowSkip());
			addChildWithNumber(mode_elem, "window_length", rx.getWindowLength());
		}
		else if (tx.getMode() == radar::OperationMode::FMCW_MODE)
		{
			(void)mono_elem.addChild("fmcw_mode");
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

		// Serialize fluctuation model if present
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
		parent.setAttribute("name", platform.getName());

		const XmlElement motion_elem = parent.addChild("motionpath");
		serializeMotionPath(*platform.getMotionPath(), motion_elem);

		serializeRotation(*platform.getRotationPath(), parent);

		for (const auto& tx : world.getTransmitters())
		{
			if (tx->getPlatform() == &platform)
			{
				if (tx->getAttached() != nullptr)
				{
					serializeMonostatic(*tx, *dynamic_cast<const radar::Receiver*>(tx->getAttached()), parent);
				}
				else
				{
					serializeTransmitter(*tx, parent);
				}
			}
		}

		for (const auto& rx : world.getReceivers())
		{
			if (rx->getPlatform() == &platform && (rx->getAttached() == nullptr))
			{
				serializeReceiver(*rx, parent);
			}
		}

		for (const auto& target : world.getTargets())
		{
			if (target->getPlatform() == &platform)
			{
				serializeTarget(*target, parent);
			}
		}
	}
}
