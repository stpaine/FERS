// SPDX-License-Identifier: GPL-2.0-only
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).

/**
 * @file json_serializer.cpp
 * @brief Implements JSON serialization and deserialization for FERS objects.
 *
 * This file leverages the `nlohmann/json` library's support for automatic
 * serialization via `to_json` and `from_json` free functions. By placing these
 * functions within the namespaces of the objects they serialize, we enable
 * Argument-Dependent Lookup (ADL). This design choice allows the library to
 * automatically find the correct conversion functions, keeping the serialization
 * logic decoupled from the core object definitions and improving modularity.
 */

#include "serial/json_serializer.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <nlohmann/json.hpp>
#include <random>
#include <unordered_map>

#include "antenna/antenna_factory.h"
#include "core/parameters.h"
#include "core/sim_id.h"
#include "core/world.h"
#include "math/coord.h"
#include "math/path.h"
#include "math/rotation_path.h"
#include "radar/platform.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/fmcw_validation.h"
#include "serial/rotation_angle_utils.h"
#include "serial/rotation_warning_utils.h"
#include "signal/radar_signal.h"
#include "timing/prototype_timing.h"
#include "timing/timing.h"
#include "waveform_factory.h"

// TODO: Add file path validation and error handling as needed.

namespace
{
	using TimingInstanceMap = std::unordered_map<SimId, std::shared_ptr<timing::Timing>>;

	nlohmann::json sim_id_to_json(const SimId id) { return std::to_string(id); }

	SimId parse_json_id(const nlohmann::json& j, const std::string& key, const std::string& owner)
	{
		if (!j.contains(key))
		{
			throw std::runtime_error("Missing required '" + key + "' for " + owner + ".");
		}
		try
		{
			if (j.at(key).is_number_unsigned())
			{
				return j.at(key).get<SimId>();
			}
			if (j.at(key).is_number_integer())
			{
				const auto value = j.at(key).get<long long>();
				if (value < 0)
				{
					throw std::runtime_error("negative id");
				}
				return static_cast<SimId>(value);
			}
			if (j.at(key).is_string())
			{
				const auto str = j.at(key).get<std::string>();
				size_t idx = 0;
				const unsigned long long parsed = std::stoull(str, &idx, 10);
				if (idx != str.size())
				{
					throw std::runtime_error("trailing characters");
				}
				return static_cast<SimId>(parsed);
			}
		}
		catch (const std::exception& e)
		{
			throw std::runtime_error("Invalid '" + key + "' for " + owner + ": " + e.what());
		}
		throw std::runtime_error("Invalid '" + key + "' type for " + owner + ".");
	}

	std::shared_ptr<timing::Timing> resolve_timing_instance(core::World& world, std::mt19937& masterSeeder,
															TimingInstanceMap& timing_instances, const SimId timing_id)
	{
		if (const auto it = timing_instances.find(timing_id); it != timing_instances.end())
		{
			return it->second;
		}

		auto* const timing_proto = world.findTiming(timing_id);
		if (timing_proto == nullptr)
		{
			return nullptr;
		}

		auto timing = std::make_shared<timing::Timing>(timing_proto->getName(), static_cast<unsigned>(masterSeeder()),
													   timing_proto->getId());
		timing->initializeModel(timing_proto);
		timing_instances.emplace(timing_id, timing);
		return timing;
	}

	void throw_json_validation_error(const std::string& message) { throw std::runtime_error(message); }

	void validate_fmcw_waveform(const fers_signal::RadarSignal& wave, const std::string& owner)
	{
		serial::fmcw_validation::validateWaveform(wave, owner, throw_json_validation_error);
	}

	void validate_waveform_mode_match(const fers_signal::RadarSignal& wave, const radar::OperationMode mode,
									  const std::string& owner)
	{
		serial::fmcw_validation::validateWaveformModeMatch(wave, mode, owner, throw_json_validation_error);
	}

	void validate_fmcw_schedule(const std::vector<radar::SchedulePeriod>& schedule,
								const fers_signal::FmcwChirpSignal& fmcw, const std::string& owner)
	{
		serial::fmcw_validation::validateSchedule(schedule, fmcw, owner, throw_json_validation_error);
	}
}

namespace math
{
	void to_json(nlohmann::json& j, const Vec3& v)
	{
		j = {{"x", v.x}, {"y", v.y}, {"z", v.z}};
	} // NOLINT(*-use-internal-linkage)

	void from_json(const nlohmann::json& j, Vec3& v) // NOLINT(*-use-internal-linkage)
	{
		j.at("x").get_to(v.x);
		j.at("y").get_to(v.y);
		j.at("z").get_to(v.z);
	}

	void to_json(nlohmann::json& j, const Coord& c) // NOLINT(*-use-internal-linkage)
	{
		j = {{"time", c.t}, {"x", c.pos.x}, {"y", c.pos.y}, {"altitude", c.pos.z}};
	}

	void from_json(const nlohmann::json& j, Coord& c) // NOLINT(*-use-internal-linkage)
	{
		j.at("time").get_to(c.t);
		j.at("x").get_to(c.pos.x);
		j.at("y").get_to(c.pos.y);
		j.at("altitude").get_to(c.pos.z);
	}

	void to_json(nlohmann::json& j, const RotationCoord& rc) // NOLINT(*-use-internal-linkage)
	{
		const auto unit = params::rotationAngleUnit();
		j = {{"time", rc.t},
			 {"azimuth", serial::rotation_angle_utils::internal_azimuth_to_external(rc.azimuth, unit)},
			 {"elevation", serial::rotation_angle_utils::internal_elevation_to_external(rc.elevation, unit)}};
	}

	void from_json(const nlohmann::json& j, RotationCoord& rc) // NOLINT(*-use-internal-linkage)
	{
		j.at("time").get_to(rc.t);
		const auto external = serial::rotation_angle_utils::external_rotation_to_internal(
			j.at("azimuth").get<RealType>(), j.at("elevation").get<RealType>(), rc.t, params::rotationAngleUnit());
		rc.azimuth = external.azimuth;
		rc.elevation = external.elevation;
	}

	NLOHMANN_JSON_SERIALIZE_ENUM(Path::InterpType,
								 {{Path::InterpType::INTERP_STATIC, "static"},
								  {Path::InterpType::INTERP_LINEAR, "linear"},
								  {Path::InterpType::INTERP_CUBIC, "cubic"}})

	void to_json(nlohmann::json& j, const Path& p) // NOLINT(*-use-internal-linkage)
	{
		j = {{"interpolation", p.getType()}, {"positionwaypoints", p.getCoords()}};
	}

	void from_json(const nlohmann::json& j, Path& p) // NOLINT(*-use-internal-linkage)
	{
		p.setInterp(j.at("interpolation").get<Path::InterpType>());
		for (const auto waypoints = j.at("positionwaypoints").get<std::vector<Coord>>(); const auto& wp : waypoints)
		{
			p.addCoord(wp);
		}
		p.finalize();
	}

	NLOHMANN_JSON_SERIALIZE_ENUM(RotationPath::InterpType,
								 {{RotationPath::InterpType::INTERP_STATIC, "static"},
								  {RotationPath::InterpType::INTERP_CONSTANT,
								   "constant"}, // Not used in xml_parser or UI yet, but for completeness
								  {RotationPath::InterpType::INTERP_LINEAR, "linear"},
								  {RotationPath::InterpType::INTERP_CUBIC, "cubic"}})

	void to_json(nlohmann::json& j, const RotationPath& p) // NOLINT(*-use-internal-linkage)
	{
		j["interpolation"] = p.getType();
		// This logic exists to map the two different rotation definitions from the
		// XML schema (<fixedrotation> and <rotationpath>) into a unified JSON
		// structure that the frontend can more easily handle.
		if (p.getType() == RotationPath::InterpType::INTERP_CONSTANT)
		{
			// A constant-rate rotation path corresponds to the <fixedrotation> XML element.
			// The start and rate values are converted to compass degrees per second.
			// No normalization is applied to preserve negative start angles.
			const auto unit = params::rotationAngleUnit();
			j["startazimuth"] = serial::rotation_angle_utils::internal_azimuth_to_external(p.getStart().azimuth, unit);
			j["startelevation"] =
				serial::rotation_angle_utils::internal_elevation_to_external(p.getStart().elevation, unit);
			j["azimuthrate"] =
				serial::rotation_angle_utils::internal_azimuth_rate_to_external(p.getRate().azimuth, unit);
			j["elevationrate"] =
				serial::rotation_angle_utils::internal_elevation_rate_to_external(p.getRate().elevation, unit);
		}
		else
		{
			j["rotationwaypoints"] = p.getCoords();
		}
	}

	void from_json(const nlohmann::json& j, RotationPath& p) // NOLINT(*-use-internal-linkage)
	{
		p.setInterp(j.at("interpolation").get<RotationPath::InterpType>());
		for (const auto waypoints = j.at("rotationwaypoints").get<std::vector<RotationCoord>>();
			 const auto& wp : waypoints)
		{
			p.addCoord(wp);
		}
		p.finalize();
	}

}

namespace timing
{
	void to_json(nlohmann::json& j, const PrototypeTiming& pt) // NOLINT(*-use-internal-linkage)
	{
		j = nlohmann::json{{"id", sim_id_to_json(pt.getId())},
						   {"name", pt.getName()},
						   {"frequency", pt.getFrequency()},
						   {"synconpulse", pt.getSyncOnPulse()}};

		if (pt.getFreqOffset().has_value())
		{
			j["freq_offset"] = pt.getFreqOffset().value();
		}
		if (pt.getRandomFreqOffsetStdev().has_value())
		{
			j["random_freq_offset_stdev"] = pt.getRandomFreqOffsetStdev().value();
		}
		if (pt.getPhaseOffset().has_value())
		{
			j["phase_offset"] = pt.getPhaseOffset().value();
		}
		if (pt.getRandomPhaseOffsetStdev().has_value())
		{
			j["random_phase_offset_stdev"] = pt.getRandomPhaseOffsetStdev().value();
		}

		std::vector<RealType> alphas;
		std::vector<RealType> weights;
		pt.copyAlphas(alphas, weights);
		if (!alphas.empty())
		{
			nlohmann::json noise_entries = nlohmann::json::array();
			for (size_t i = 0; i < alphas.size(); ++i)
			{
				noise_entries.push_back({{"alpha", alphas[i]}, {"weight", weights[i]}});
			}
			j["noise_entries"] = noise_entries;
		}
	}

	void from_json(const nlohmann::json& j, PrototypeTiming& pt) // NOLINT(*-use-internal-linkage)
	{
		pt.setFrequency(j.at("frequency").get<RealType>());
		if (j.value("synconpulse", false))
		{
			pt.setSyncOnPulse();
		}
		else
		{
			pt.clearSyncOnPulse();
		}

		if (j.contains("freq_offset"))
		{
			pt.setFreqOffset(j.at("freq_offset").get<RealType>());
		}
		else
			pt.clearFreqOffset();
		if (j.contains("random_freq_offset_stdev"))
		{
			pt.setRandomFreqOffsetStdev(j.at("random_freq_offset_stdev").get<RealType>());
		}
		else
			pt.clearRandomFreqOffsetStdev();
		if (j.contains("phase_offset"))
		{
			pt.setPhaseOffset(j.at("phase_offset").get<RealType>());
		}
		else
			pt.clearPhaseOffset();
		if (j.contains("random_phase_offset_stdev"))
		{
			pt.setRandomPhaseOffsetStdev(j.at("random_phase_offset_stdev").get<RealType>());
		}
		else
			pt.clearRandomPhaseOffsetStdev();

		pt.clearNoiseEntries();
		if (j.contains("noise_entries"))
		{
			for (const auto& entry : j.at("noise_entries"))
			{
				pt.setAlpha(entry.at("alpha").get<RealType>(), entry.at("weight").get<RealType>());
			}
		}
	}
}

namespace fers_signal
{
	void to_json(nlohmann::json& j, const RadarSignal& rs) // NOLINT(*-use-internal-linkage)
	{
		j = nlohmann::json{{"id", sim_id_to_json(rs.getId())},
						   {"name", rs.getName()},
						   {"power", rs.getPower()},
						   {"carrier_frequency", rs.getCarrier()}};
		if (dynamic_cast<const CwSignal*>(rs.getSignal()) != nullptr)
		{
			j["cw"] = nlohmann::json::object();
		}
		else if (const auto* fmcw = rs.getFmcwChirpSignal(); fmcw != nullptr)
		{
			j["fmcw_up_chirp"] = {{"chirp_bandwidth", fmcw->getChirpBandwidth()},
								  {"chirp_duration", fmcw->getChirpDuration()},
								  {"chirp_period", fmcw->getChirpPeriod()}};
			if (std::abs(fmcw->getStartFrequencyOffset()) > EPSILON)
			{
				j["fmcw_up_chirp"]["start_frequency_offset"] = fmcw->getStartFrequencyOffset();
			}
			if (fmcw->getChirpCount().has_value())
			{
				j["fmcw_up_chirp"]["chirp_count"] = *fmcw->getChirpCount();
			}
		}
		else
		{
			if (const auto& filename = rs.getFilename(); filename.has_value())
			{
				j["pulsed_from_file"] = {{"filename", *filename}};
			}
			else
			{
				throw std::logic_error("Attempted to serialize a file-based waveform named '" + rs.getName() +
									   "' without a source filename.");
			}
		}
	}

	void from_json(const nlohmann::json& j, std::unique_ptr<RadarSignal>& rs) // NOLINT(*-use-internal-linkage)
	{
		const auto name = j.at("name").get<std::string>();
		const auto id = parse_json_id(j, "id", "waveform");
		const auto power = j.at("power").get<RealType>();
		const auto carrier = j.at("carrier_frequency").get<RealType>();

		if (j.contains("cw"))
		{
			auto cw_signal = std::make_unique<CwSignal>();
			rs = std::make_unique<RadarSignal>(name, power, carrier, params::endTime() - params::startTime(),
											   std::move(cw_signal), id);
		}
		else if (j.contains("fmcw_up_chirp"))
		{
			const auto& fmcw_json = j.at("fmcw_up_chirp");
			std::optional<std::size_t> chirp_count;
			if (fmcw_json.contains("chirp_count"))
			{
				const auto parsed_count = fmcw_json.at("chirp_count").get<long long>();
				if (parsed_count <= 0)
				{
					throw std::runtime_error("Waveform '" + name + "' has an invalid chirp_count.");
				}
				chirp_count = static_cast<std::size_t>(parsed_count);
			}

			auto fmcw_signal = std::make_unique<FmcwChirpSignal>(
				fmcw_json.at("chirp_bandwidth").get<RealType>(), fmcw_json.at("chirp_duration").get<RealType>(),
				fmcw_json.at("chirp_period").get<RealType>(), fmcw_json.value("start_frequency_offset", 0.0),
				chirp_count);
			rs = std::make_unique<RadarSignal>(name, power, carrier, fmcw_signal->getChirpDuration(),
											   std::move(fmcw_signal), id);
			validate_fmcw_waveform(*rs, "Waveform '" + name + "'");
		}
		else if (j.contains("pulsed_from_file"))
		{
			const auto& pulsed_file = j.at("pulsed_from_file");
			const auto filename = pulsed_file.value("filename", "");
			if (filename.empty())
			{
				LOG(logging::Level::WARNING, "Skipping load of file-based waveform '{}': filename is empty.", name);
				return; // rs remains nullptr
			}
			rs = serial::loadWaveformFromFile(name, filename, power, carrier, id);
		}
		else
		{
			throw std::runtime_error("Unsupported waveform type in from_json for '" + name + "'");
		}
	}
}

namespace antenna
{
	void to_json(nlohmann::json& j, const Antenna& a) // NOLINT(*-use-internal-linkage)
	{
		j = {{"id", sim_id_to_json(a.getId())}, {"name", a.getName()}, {"efficiency", a.getEfficiencyFactor()}};

		if (const auto* sinc = dynamic_cast<const Sinc*>(&a))
		{
			j["pattern"] = "sinc";
			j["alpha"] = sinc->getAlpha();
			j["beta"] = sinc->getBeta();
			j["gamma"] = sinc->getGamma();
		}
		else if (const auto* gaussian = dynamic_cast<const Gaussian*>(&a))
		{
			j["pattern"] = "gaussian";
			j["azscale"] = gaussian->getAzimuthScale();
			j["elscale"] = gaussian->getElevationScale();
		}
		else if (const auto* sh = dynamic_cast<const SquareHorn*>(&a))
		{
			j["pattern"] = "squarehorn";
			j["diameter"] = sh->getDimension();
		}
		else if (const auto* parabolic = dynamic_cast<const Parabolic*>(&a))
		{
			j["pattern"] = "parabolic";
			j["diameter"] = parabolic->getDiameter();
		}
		else if (const auto* xml = dynamic_cast<const XmlAntenna*>(&a))
		{
			j["pattern"] = "xml";
			j["filename"] = xml->getFilename();
		}
		else if (const auto* h5 = dynamic_cast<const H5Antenna*>(&a))
		{
			j["pattern"] = "file";
			j["filename"] = h5->getFilename();
		}
		else
		{
			j["pattern"] = "isotropic";
		}
	}

	void from_json(const nlohmann::json& j, std::unique_ptr<Antenna>& ant) // NOLINT(*-use-internal-linkage)
	{
		const auto name = j.at("name").get<std::string>();
		const auto id = parse_json_id(j, "id", "Antenna");
		const auto pattern = j.value("pattern", "isotropic");

		if (pattern == "isotropic")
		{
			ant = std::make_unique<Isotropic>(name, id);
		}
		else if (pattern == "sinc")
		{
			ant = std::make_unique<Sinc>(name, j.at("alpha").get<RealType>(), j.at("beta").get<RealType>(),
										 j.at("gamma").get<RealType>(), id);
		}
		else if (pattern == "gaussian")
		{
			ant =
				std::make_unique<Gaussian>(name, j.at("azscale").get<RealType>(), j.at("elscale").get<RealType>(), id);
		}
		else if (pattern == "squarehorn")
		{
			ant = std::make_unique<SquareHorn>(name, j.at("diameter").get<RealType>(), id);
		}
		else if (pattern == "parabolic")
		{
			ant = std::make_unique<Parabolic>(name, j.at("diameter").get<RealType>(), id);
		}
		else if (pattern == "xml")
		{
			const auto filename = j.value("filename", "");
			if (filename.empty())
			{
				LOG(logging::Level::WARNING, "Skipping load of XML antenna '{}': filename is empty.", name);
				return; // ant remains nullptr
			}
			ant = std::make_unique<XmlAntenna>(name, filename, id);
		}
		else if (pattern == "file")
		{
			const auto filename = j.value("filename", "");
			if (filename.empty())
			{
				LOG(logging::Level::WARNING, "Skipping load of H5 antenna '{}': filename is empty.", name);
				return; // ant remains nullptr
			}
			ant = std::make_unique<H5Antenna>(name, filename, id);
		}
		else
		{
			throw std::runtime_error("Unsupported antenna pattern in from_json: " + pattern);
		}

		ant->setEfficiencyFactor(j.value("efficiency", 1.0));
	}
}

namespace radar
{
	void to_json(nlohmann::json& j, const SchedulePeriod& p)
	{
		j = {{"start", p.start}, {"end", p.end}};
	} // NOLINT(*-use-internal-linkage)

	void from_json(const nlohmann::json& j, SchedulePeriod& p) // NOLINT(*-use-internal-linkage)
	{
		j.at("start").get_to(p.start);
		j.at("end").get_to(p.end);
	}

	void to_json(nlohmann::json& j, const Transmitter& t) // NOLINT(*-use-internal-linkage)
	{
		j = nlohmann::json{{"id", sim_id_to_json(t.getId())},
						   {"name", t.getName()},
						   {"waveform", sim_id_to_json((t.getSignal() != nullptr) ? t.getSignal()->getId() : 0)},
						   {"antenna", sim_id_to_json((t.getAntenna() != nullptr) ? t.getAntenna()->getId() : 0)},
						   {"timing", sim_id_to_json(t.getTiming() ? t.getTiming()->getId() : 0)}};

		if (t.getMode() == OperationMode::PULSED_MODE)
		{
			j["pulsed_mode"] = {{"prf", t.getPrf()}};
		}
		else if (t.getMode() == OperationMode::FMCW_MODE)
		{
			j["fmcw_mode"] = nlohmann::json::object();
		}
		else
		{
			j["cw_mode"] = nlohmann::json::object();
		}
		if (!t.getSchedule().empty())
		{
			j["schedule"] = t.getSchedule();
		}
	}

	void to_json(nlohmann::json& j, const Receiver& r) // NOLINT(*-use-internal-linkage)
	{
		j = nlohmann::json{{"id", sim_id_to_json(r.getId())},
						   {"name", r.getName()},
						   {"noise_temp", r.getNoiseTemperature()},
						   {"antenna", sim_id_to_json((r.getAntenna() != nullptr) ? r.getAntenna()->getId() : 0)},
						   {"timing", sim_id_to_json(r.getTiming() ? r.getTiming()->getId() : 0)},
						   {"nodirect", r.checkFlag(Receiver::RecvFlag::FLAG_NODIRECT)},
						   {"nopropagationloss", r.checkFlag(Receiver::RecvFlag::FLAG_NOPROPLOSS)}};

		if (r.getMode() == OperationMode::PULSED_MODE)
		{
			j["pulsed_mode"] = {
				{"prf", r.getWindowPrf()}, {"window_skip", r.getWindowSkip()}, {"window_length", r.getWindowLength()}};
		}
		else if (r.getMode() == OperationMode::FMCW_MODE)
		{
			j["fmcw_mode"] = nlohmann::json::object();
		}
		else
		{
			j["cw_mode"] = nlohmann::json::object();
		}
		if (!r.getSchedule().empty())
		{
			j["schedule"] = r.getSchedule();
		}
	}

	void to_json(nlohmann::json& j, const Target& t) // NOLINT(*-use-internal-linkage)
	{
		j["id"] = sim_id_to_json(t.getId());
		j["name"] = t.getName();
		nlohmann::json rcs_json;
		if (const auto* iso = dynamic_cast<const IsoTarget*>(&t))
		{
			rcs_json["type"] = "isotropic";
			rcs_json["value"] = iso->getConstRcs();
		}
		else if (const auto* file = dynamic_cast<const FileTarget*>(&t))
		{
			rcs_json["type"] = "file";
			rcs_json["filename"] = file->getFilename();
		}
		j["rcs"] = rcs_json;

		// Serialize the fluctuation model if it exists.
		if (const auto* model_base = t.getFluctuationModel())
		{
			nlohmann::json model_json;
			if (const auto* chi_model = dynamic_cast<const RcsChiSquare*>(model_base))
			{
				model_json["type"] = "chisquare";
				model_json["k"] = chi_model->getK();
			}
			else // Default to constant if it's not a recognized type (e.g., RcsConst)
			{
				model_json["type"] = "constant";
			}
			j["model"] = model_json;
		}
	}

	void to_json(nlohmann::json& j, const Platform& p) // NOLINT(*-use-internal-linkage)
	{
		j = {{"id", sim_id_to_json(p.getId())}, {"name", p.getName()}, {"motionpath", *p.getMotionPath()}};

		if (p.getRotationPath()->getType() == math::RotationPath::InterpType::INTERP_CONSTANT)
		{
			j["fixedrotation"] = *p.getRotationPath();
		}
		else
		{
			j["rotationpath"] = *p.getRotationPath();
		}
	}

}

namespace params
{
	NLOHMANN_JSON_SERIALIZE_ENUM(CoordinateFrame,
								 {{CoordinateFrame::ENU, "ENU"},
								  {CoordinateFrame::UTM, "UTM"},
								  {CoordinateFrame::ECEF, "ECEF"}})
	NLOHMANN_JSON_SERIALIZE_ENUM(RotationAngleUnit,
								 {{RotationAngleUnit::Degrees, "deg"}, {RotationAngleUnit::Radians, "rad"}})

	void to_json(nlohmann::json& j, const Parameters& p) // NOLINT(*-use-internal-linkage)
	{
		j = nlohmann::json{{"starttime", p.start},
						   {"endtime", p.end},
						   {"rate", p.rate},
						   {"c", p.c},
						   {"simSamplingRate", p.sim_sampling_rate},
						   {"adc_bits", p.adc_bits},
						   {"oversample", p.oversample_ratio},
						   {"rotationangleunit", p.rotation_angle_unit}};

		if (p.random_seed.has_value())
		{
			j["randomseed"] = p.random_seed.value();
		}

		j["origin"] = {
			{"latitude", p.origin_latitude}, {"longitude", p.origin_longitude}, {"altitude", p.origin_altitude}};

		j["coordinatesystem"] = {{"frame", p.coordinate_frame}};
		if (p.coordinate_frame == CoordinateFrame::UTM)
		{
			j["coordinatesystem"]["zone"] = p.utm_zone;
			j["coordinatesystem"]["hemisphere"] = p.utm_north_hemisphere ? "N" : "S";
		}
	}

	void from_json(const nlohmann::json& j, Parameters& p) // NOLINT(*-use-internal-linkage)
	{
		p.start = j.at("starttime").get<RealType>();
		p.end = j.at("endtime").get<RealType>();
		p.rate = j.at("rate").get<RealType>();
		p.c = j.value("c", Parameters::DEFAULT_C);
		p.sim_sampling_rate = j.value("simSamplingRate", 1000.0);
		p.adc_bits = j.value("adc_bits", 0u);
		p.oversample_ratio = j.value("oversample", 1u);
		params::validateOversampleRatio(p.oversample_ratio);
		p.rotation_angle_unit = j.value("rotationangleunit", RotationAngleUnit::Degrees);
		p.random_seed = j.value<std::optional<unsigned>>("randomseed", std::nullopt);

		const auto& origin = j.at("origin");
		p.origin_latitude = origin.at("latitude").get<double>();
		p.origin_longitude = origin.at("longitude").get<double>();
		p.origin_altitude = origin.at("altitude").get<double>();

		const auto& cs = j.at("coordinatesystem");
		p.coordinate_frame = cs.at("frame").get<CoordinateFrame>();
		if (p.coordinate_frame == CoordinateFrame::UTM)
		{
			p.utm_zone = cs.at("zone").get<int>();
			p.utm_north_hemisphere = cs.at("hemisphere").get<std::string>() == "N";
		}
	}
}

namespace
{
	nlohmann::json serialize_platform(const radar::Platform* p, const core::World& world)
	{
		nlohmann::json plat_json = *p;

		// Initialize components array to ensure it exists even if empty
		plat_json["components"] = nlohmann::json::array();

		// Add Transmitters and Monostatic Radars
		for (const auto& t : world.getTransmitters())
		{
			if (t->getPlatform() == p)
			{
				if (t->getAttached() != nullptr)
				{
					nlohmann::json monostatic_comp;
					monostatic_comp["name"] = t->getName();
					monostatic_comp["tx_id"] = sim_id_to_json(t->getId());
					monostatic_comp["rx_id"] = sim_id_to_json(t->getAttached()->getId());
					monostatic_comp["waveform"] =
						sim_id_to_json((t->getSignal() != nullptr) ? t->getSignal()->getId() : 0);
					monostatic_comp["antenna"] =
						sim_id_to_json((t->getAntenna() != nullptr) ? t->getAntenna()->getId() : 0);
					monostatic_comp["timing"] = sim_id_to_json(t->getTiming() ? t->getTiming()->getId() : 0);

					if (const auto* recv = dynamic_cast<const radar::Receiver*>(t->getAttached()))
					{
						monostatic_comp["noise_temp"] = recv->getNoiseTemperature();
						monostatic_comp["nodirect"] = recv->checkFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
						monostatic_comp["nopropagationloss"] =
							recv->checkFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);

						if (!t->getSchedule().empty())
						{
							monostatic_comp["schedule"] = t->getSchedule();
						}

						if (t->getMode() == radar::OperationMode::PULSED_MODE)
						{
							monostatic_comp["pulsed_mode"] = {{"prf", t->getPrf()},
															  {"window_skip", recv->getWindowSkip()},
															  {"window_length", recv->getWindowLength()}};
						}
						else
						{
							if (t->getMode() == radar::OperationMode::FMCW_MODE)
							{
								monostatic_comp["fmcw_mode"] = nlohmann::json::object();
							}
							else
							{
								monostatic_comp["cw_mode"] = nlohmann::json::object();
							}
						}
					}
					plat_json["components"].push_back(nlohmann::json{{"monostatic", monostatic_comp}});
				}
				else
				{
					plat_json["components"].push_back(nlohmann::json{{"transmitter", *t}});
				}
			}
		}

		// Add Standalone Receivers
		for (const auto& r : world.getReceivers())
		{
			if (r->getPlatform() == p)
			{
				// This must be a standalone receiver, as monostatic cases were handled above.
				if (r->getAttached() == nullptr)
				{
					plat_json["components"].push_back(nlohmann::json{{"receiver", *r}});
				}
			}
		}

		// Add Targets
		for (const auto& target : world.getTargets())
		{
			if (target->getPlatform() == p)
			{
				plat_json["components"].push_back(nlohmann::json{{"target", *target}});
			}
		}

		return plat_json;
	}

	void parse_parameters(const nlohmann::json& sim, std::mt19937& masterSeeder)
	{
		auto new_params = sim.at("parameters").get<params::Parameters>();

		// If a random seed is present in the incoming JSON, it is used to re-seed
		// the master generator. This is crucial for allowing the UI to control
		// simulation reproducibility.
		if (sim.at("parameters").contains("randomseed"))
		{
			params::params.random_seed = new_params.random_seed;
			if (params::params.random_seed)
			{
				LOG(logging::Level::INFO, "Master seed updated from JSON to: {}", *params::params.random_seed);
				masterSeeder.seed(*params::params.random_seed);
			}
		}

		new_params.random_seed = params::params.random_seed;
		params::params = new_params;
		params::params.simulation_name = sim.value("name", "");
	}

	void parse_assets(const nlohmann::json& sim, core::World& world)
	{
		if (sim.contains("waveforms"))
		{
			for (auto waveforms = sim.at("waveforms").get<std::vector<std::unique_ptr<fers_signal::RadarSignal>>>();
				 auto& waveform : waveforms)
			{
				// Only add valid waveforms. If filename was empty, waveform is nullptr.
				if (waveform)
				{
					world.add(std::move(waveform));
				}
			}
		}

		if (sim.contains("antennas"))
		{
			for (auto antennas = sim.at("antennas").get<std::vector<std::unique_ptr<antenna::Antenna>>>();
				 auto& antenna : antennas)
			{
				// Only add valid antennas.
				if (antenna)
				{
					world.add(std::move(antenna));
				}
			}
		}

		if (sim.contains("timings"))
		{
			for (const auto& timing_json : sim.at("timings"))
			{
				auto name = timing_json.at("name").get<std::string>();
				const auto timing_id = parse_json_id(timing_json, "id", "Timing");
				auto timing_obj = std::make_unique<timing::PrototypeTiming>(name, timing_id);
				timing_json.get_to(*timing_obj);
				world.add(std::move(timing_obj));
			}
		}
	}

	radar::OperationMode parse_mode(const nlohmann::json& comp_json, const std::string& error_context)
	{
		if (comp_json.contains("pulsed_mode"))
		{
			return radar::OperationMode::PULSED_MODE;
		}
		if (comp_json.contains("fmcw_mode"))
		{
			return radar::OperationMode::FMCW_MODE;
		}
		if (comp_json.contains("cw_mode"))
		{
			return radar::OperationMode::CW_MODE;
		}
		throw std::runtime_error(error_context + " must have a 'pulsed_mode', 'cw_mode', or 'fmcw_mode' block.");
	}

	void parse_transmitter(const nlohmann::json& comp_json, radar::Platform* plat, core::World& world,
						   std::mt19937& masterSeeder, TimingInstanceMap& timing_instances)
	{
		// --- Dependency Check ---
		// Validate Waveform and Timing existence before creation to prevent core crashes.
		const auto wave_id = parse_json_id(comp_json, "waveform", "Transmitter");
		const auto timing_id = parse_json_id(comp_json, "timing", "Transmitter");
		const auto antenna_id = parse_json_id(comp_json, "antenna", "Transmitter");

		if (world.findWaveform(wave_id) == nullptr)
		{
			LOG(logging::Level::WARNING, "Skipping Transmitter '{}': Missing or invalid waveform '{}'.",
				comp_json.value("name", "Unnamed"), comp_json.value("waveform", ""));
			return;
		}
		if (world.findTiming(timing_id) == nullptr)
		{
			LOG(logging::Level::WARNING, "Skipping Transmitter '{}': Missing or invalid timing source '{}'.",
				comp_json.value("name", "Unnamed"), comp_json.value("timing", ""));
			return;
		}
		if (world.findAntenna(antenna_id) == nullptr)
		{
			LOG(logging::Level::WARNING, "Skipping Transmitter '{}': Missing or invalid antenna '{}'.",
				comp_json.value("name", "Unnamed"), comp_json.value("antenna", ""));
			return;
		}

		radar::OperationMode mode =
			parse_mode(comp_json, "Transmitter component '" + comp_json.value("name", "Unnamed") + "'");

		const auto trans_id = parse_json_id(comp_json, "id", "Transmitter");
		auto trans = std::make_unique<radar::Transmitter>(plat, comp_json.value("name", "Unnamed"), mode, trans_id);
		if (mode == radar::OperationMode::PULSED_MODE && comp_json.contains("pulsed_mode"))
		{
			trans->setPrf(comp_json.at("pulsed_mode").value("prf", 0.0));
		}

		auto* const waveform = world.findWaveform(wave_id);
		validate_fmcw_waveform(*waveform, "Waveform '" + waveform->getName() + "'");
		validate_waveform_mode_match(*waveform, mode,
									 "Transmitter component '" + comp_json.value("name", "Unnamed") + "'");
		trans->setWave(waveform);
		trans->setAntenna(world.findAntenna(antenna_id));

		if (const auto timing = resolve_timing_instance(world, masterSeeder, timing_instances, timing_id))
		{
			trans->setTiming(timing);
		}

		if (comp_json.contains("schedule"))
		{
			auto raw = comp_json.at("schedule").get<std::vector<radar::SchedulePeriod>>();
			RealType pri = 0.0;
			if (mode == radar::OperationMode::PULSED_MODE)
			{
				pri = 1.0 / trans->getPrf();
			}
			auto schedule = radar::processRawSchedule(std::move(raw), trans->getName(),
													  mode == radar::OperationMode::PULSED_MODE, pri);
			if (const auto* fmcw = waveform->getFmcwChirpSignal(); fmcw != nullptr)
			{
				validate_fmcw_schedule(schedule, *fmcw, "Transmitter component '" + trans->getName() + "'");
			}
			trans->setSchedule(std::move(schedule));
		}

		world.add(std::move(trans));
	}

	void parse_receiver(const nlohmann::json& comp_json, radar::Platform* plat, core::World& world,
						std::mt19937& masterSeeder, TimingInstanceMap& timing_instances)
	{
		// --- Dependency Check ---
		// Receiver strictly requires a Timing source.
		const auto timing_id = parse_json_id(comp_json, "timing", "Receiver");
		const auto antenna_id = parse_json_id(comp_json, "antenna", "Receiver");

		if (world.findTiming(timing_id) == nullptr)
		{
			LOG(logging::Level::WARNING, "Skipping Receiver '{}': Missing or invalid timing source '{}'.",
				comp_json.value("name", "Unnamed"), comp_json.value("timing", ""));
			return;
		}

		if (world.findAntenna(antenna_id) == nullptr)
		{
			LOG(logging::Level::WARNING, "Skipping Receiver '{}': Missing or invalid antenna '{}'.",
				comp_json.value("name", "Unnamed"), comp_json.value("antenna", ""));
			return;
		}

		radar::OperationMode mode =
			parse_mode(comp_json, "Receiver component '" + comp_json.value("name", "Unnamed") + "'");

		const auto recv_id = parse_json_id(comp_json, "id", "Receiver");
		auto recv =
			std::make_unique<radar::Receiver>(plat, comp_json.value("name", "Unnamed"), masterSeeder(), mode, recv_id);
		if (mode == radar::OperationMode::PULSED_MODE && comp_json.contains("pulsed_mode"))
		{
			const auto& mode_json = comp_json.at("pulsed_mode");
			recv->setWindowProperties(mode_json.value("window_length", 0.0), mode_json.value("prf", 0.0),
									  mode_json.value("window_skip", 0.0));
		}

		recv->setNoiseTemperature(comp_json.value("noise_temp", 0.0));

		recv->setAntenna(world.findAntenna(antenna_id));

		if (const auto timing = resolve_timing_instance(world, masterSeeder, timing_instances, timing_id))
		{
			recv->setTiming(timing);
		}

		if (comp_json.value("nodirect", false))
		{
			recv->setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
		}
		if (comp_json.value("nopropagationloss", false))
		{
			recv->setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
		}

		if (comp_json.contains("schedule"))
		{
			auto raw = comp_json.at("schedule").get<std::vector<radar::SchedulePeriod>>();
			RealType pri = 0.0;
			if (mode == radar::OperationMode::PULSED_MODE)
			{
				pri = 1.0 / recv->getWindowPrf();
			}
			recv->setSchedule(radar::processRawSchedule(std::move(raw), recv->getName(),
														mode == radar::OperationMode::PULSED_MODE, pri));
		}

		world.add(std::move(recv));
	}

	void parse_target(const nlohmann::json& comp_json, radar::Platform* plat, core::World& world,
					  std::mt19937& masterSeeder)
	{
		const auto& rcs_json = comp_json.at("rcs");
		const auto rcs_type = rcs_json.at("type").get<std::string>();
		std::unique_ptr<radar::Target> target_obj;

		if (rcs_type == "isotropic")
		{
			const auto target_id = parse_json_id(comp_json, "id", "Target");
			target_obj = radar::createIsoTarget(plat, comp_json.at("name").get<std::string>(),
												rcs_json.at("value").get<RealType>(),
												static_cast<unsigned>(masterSeeder()), target_id);
		}
		else if (rcs_type == "file")
		{
			const auto filename = rcs_json.value("filename", "");
			if (filename.empty())
			{
				LOG(logging::Level::WARNING, "Skipping load of file target '{}': RCS filename is empty.",
					comp_json.value("name", "Unknown"));
				return;
			}
			const auto target_id = parse_json_id(comp_json, "id", "Target");
			target_obj = radar::createFileTarget(plat, comp_json.at("name").get<std::string>(), filename,
												 static_cast<unsigned>(masterSeeder()), target_id);
		}
		else
		{
			throw std::runtime_error("Unsupported target RCS type: " + rcs_type);
		}
		world.add(std::move(target_obj));

		// After creating the target, check for and apply the fluctuation model.
		if (comp_json.contains("model"))
		{
			const auto& model_json = comp_json.at("model");
			const auto model_type = model_json.at("type").get<std::string>();

			if (model_type == "chisquare" || model_type == "gamma")
			{
				auto model = std::make_unique<radar::RcsChiSquare>(world.getTargets().back()->getRngEngine(),
																   model_json.at("k").get<RealType>());
				world.getTargets().back()->setFluctuationModel(std::move(model));
			}
			else if (model_type == "constant")
			{
				world.getTargets().back()->setFluctuationModel(std::make_unique<radar::RcsConst>());
			}
			else
			{
				throw std::runtime_error("Unsupported fluctuation model type: " + model_type);
			}
		}
	}

	void parse_monostatic(const nlohmann::json& comp_json, radar::Platform* plat, core::World& world,
						  std::mt19937& masterSeeder, TimingInstanceMap& timing_instances)
	{
		// This block reconstructs the internal C++ representation of a
		// monostatic radar (a linked Transmitter and Receiver) from the
		// single 'monostatic' component in the JSON.
		// --- Dependency Check ---
		const auto wave_id = parse_json_id(comp_json, "waveform", "Monostatic");
		const auto timing_id = parse_json_id(comp_json, "timing", "Monostatic");
		const auto antenna_id = parse_json_id(comp_json, "antenna", "Monostatic");

		if (world.findWaveform(wave_id) == nullptr)
		{
			LOG(logging::Level::WARNING, "Skipping Monostatic '{}': Missing or invalid waveform '{}'.",
				comp_json.value("name", "Unnamed"), comp_json.value("waveform", ""));
			return;
		}
		if (world.findTiming(timing_id) == nullptr)
		{
			LOG(logging::Level::WARNING, "Skipping Monostatic '{}': Missing or invalid timing source '{}'.",
				comp_json.value("name", "Unnamed"), comp_json.value("timing", ""));
			return;
		}
		if (world.findAntenna(antenna_id) == nullptr)
		{
			LOG(logging::Level::WARNING, "Skipping Monostatic '{}': Missing or invalid antenna '{}'.",
				comp_json.value("name", "Unnamed"), comp_json.value("antenna", ""));
			return;
		}

		radar::OperationMode mode =
			parse_mode(comp_json, "Monostatic component '" + comp_json.value("name", "Unnamed") + "'");

		// Transmitter part
		const auto tx_id = parse_json_id(comp_json, "tx_id", "Monostatic");
		auto trans = std::make_unique<radar::Transmitter>(plat, comp_json.value("name", "Unnamed"), mode, tx_id);
		if (mode == radar::OperationMode::PULSED_MODE && comp_json.contains("pulsed_mode"))
		{
			trans->setPrf(comp_json.at("pulsed_mode").value("prf", 0.0));
		}

		auto* const waveform = world.findWaveform(wave_id);
		validate_fmcw_waveform(*waveform, "Waveform '" + waveform->getName() + "'");
		validate_waveform_mode_match(*waveform, mode,
									 "Monostatic component '" + comp_json.value("name", "Unnamed") + "'");
		trans->setWave(waveform);
		trans->setAntenna(world.findAntenna(antenna_id));
		if (const auto shared_timing = resolve_timing_instance(world, masterSeeder, timing_instances, timing_id))
		{
			trans->setTiming(shared_timing);
		}

		// Receiver part
		const auto rx_id = parse_json_id(comp_json, "rx_id", "Monostatic");
		auto recv =
			std::make_unique<radar::Receiver>(plat, comp_json.value("name", "Unnamed"), masterSeeder(), mode, rx_id);
		if (mode == radar::OperationMode::PULSED_MODE && comp_json.contains("pulsed_mode"))
		{
			const auto& mode_json = comp_json.at("pulsed_mode");
			recv->setWindowProperties(mode_json.value("window_length", 0.0),
									  trans->getPrf(), // Use transmitter's PRF
									  mode_json.value("window_skip", 0.0));
		}
		recv->setNoiseTemperature(comp_json.value("noise_temp", 0.0));

		recv->setAntenna(world.findAntenna(antenna_id));
		if (const auto shared_timing = resolve_timing_instance(world, masterSeeder, timing_instances, timing_id))
		{
			recv->setTiming(shared_timing);
		}

		if (comp_json.value("nodirect", false))
		{
			recv->setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
		}
		if (comp_json.value("nopropagationloss", false))
		{
			recv->setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
		}
		if (comp_json.contains("schedule"))
		{
			auto raw = comp_json.at("schedule").get<std::vector<radar::SchedulePeriod>>();
			RealType pri = 0.0;
			if (mode == radar::OperationMode::PULSED_MODE)
			{
				pri = 1.0 / trans->getPrf();
			}

			// Process once, apply to both
			auto processed_schedule = radar::processRawSchedule(std::move(raw), trans->getName(),
																mode == radar::OperationMode::PULSED_MODE, pri);
			if (const auto* fmcw = waveform->getFmcwChirpSignal(); fmcw != nullptr)
			{
				validate_fmcw_schedule(processed_schedule, *fmcw,
									   "Monostatic component '" + comp_json.value("name", "Unnamed") + "'");
			}

			trans->setSchedule(processed_schedule);
			recv->setSchedule(processed_schedule);
		}

		// Link them and add to world
		trans->setAttached(recv.get());
		recv->setAttached(trans.get());
		world.add(std::move(trans));
		world.add(std::move(recv));
	}

	void parse_platform(const nlohmann::json& plat_json, core::World& world, std::mt19937& masterSeeder,
						TimingInstanceMap& timing_instances)
	{
		auto name = plat_json.at("name").get<std::string>();
		const auto platform_id = parse_json_id(plat_json, "id", "Platform");
		auto plat = std::make_unique<radar::Platform>(name, platform_id);

		serial::update_platform_paths_from_json(plat_json, plat.get());

		// Components - Strict array format
		if (plat_json.contains("components"))
		{
			for (const auto& comp_json_outer : plat_json.at("components"))
			{
				if (comp_json_outer.contains("transmitter"))
				{
					parse_transmitter(comp_json_outer.at("transmitter"), plat.get(), world, masterSeeder,
									  timing_instances);
				}
				else if (comp_json_outer.contains("receiver"))
				{
					parse_receiver(comp_json_outer.at("receiver"), plat.get(), world, masterSeeder, timing_instances);
				}
				else if (comp_json_outer.contains("target"))
				{
					parse_target(comp_json_outer.at("target"), plat.get(), world, masterSeeder);
				}
				else if (comp_json_outer.contains("monostatic"))
				{
					parse_monostatic(comp_json_outer.at("monostatic"), plat.get(), world, masterSeeder,
									 timing_instances);
				}
			}
		}

		world.add(std::move(plat));
	}
}

namespace serial
{
	std::unique_ptr<antenna::Antenna> parse_antenna_from_json(const nlohmann::json& j)
	{
		std::unique_ptr<antenna::Antenna> ant;
		antenna::from_json(j, ant);
		return ant;
	}

	std::unique_ptr<fers_signal::RadarSignal> parse_waveform_from_json(const nlohmann::json& j)
	{
		std::unique_ptr<fers_signal::RadarSignal> wf;
		fers_signal::from_json(j, wf);
		return wf;
	}

	std::unique_ptr<timing::PrototypeTiming> parse_timing_from_json(const nlohmann::json& j, const SimId id)
	{
		auto timing = std::make_unique<timing::PrototypeTiming>(j.at("name").get<std::string>(), id);
		j.get_to(*timing);
		return timing;
	}

	void update_parameters_from_json(const nlohmann::json& j, std::mt19937& masterSeeder)
	{
		nlohmann::json sim;
		sim["parameters"] = j;
		parse_parameters(sim, masterSeeder);
	}

	void update_antenna_from_json(const nlohmann::json& j, antenna::Antenna* ant, core::World& world)
	{
		auto new_pattern = j.value("pattern", "isotropic");
		bool type_changed = false;

		const auto parse_required_antenna = [&j]()
		{
			auto parsed = parse_antenna_from_json(j);
			if (parsed == nullptr)
			{
				const auto name = j.value("name", std::string{});
				const auto pattern = j.value("pattern", "isotropic");
				throw std::runtime_error("Cannot update antenna '" + name + "' to pattern '" + pattern +
										 "' without a filename.");
			}
			return parsed;
		};

		if (new_pattern == "isotropic" && (dynamic_cast<antenna::Isotropic*>(ant) == nullptr))
			type_changed = true;
		else if (new_pattern == "sinc" && (dynamic_cast<antenna::Sinc*>(ant) == nullptr))
			type_changed = true;
		else if (new_pattern == "gaussian" && (dynamic_cast<antenna::Gaussian*>(ant) == nullptr))
			type_changed = true;
		else if (new_pattern == "squarehorn" && (dynamic_cast<antenna::SquareHorn*>(ant) == nullptr))
			type_changed = true;
		else if (new_pattern == "parabolic" && (dynamic_cast<antenna::Parabolic*>(ant) == nullptr))
			type_changed = true;
		else if (new_pattern == "xml" && (dynamic_cast<antenna::XmlAntenna*>(ant) == nullptr))
			type_changed = true;
		else if (new_pattern == "file" && (dynamic_cast<antenna::H5Antenna*>(ant) == nullptr))
			type_changed = true;

		if (type_changed)
		{
			world.replace(parse_required_antenna());
			return;
		}

		ant->setName(j.at("name").get<std::string>());
		ant->setEfficiencyFactor(j.value("efficiency", 1.0));

		if (auto* sinc = dynamic_cast<antenna::Sinc*>(ant))
		{
			sinc->setAlpha(j.value("alpha", 1.0));
			sinc->setBeta(j.value("beta", 1.0));
			sinc->setGamma(j.value("gamma", 2.0));
		}
		else if (auto* gauss = dynamic_cast<antenna::Gaussian*>(ant))
		{
			gauss->setAzimuthScale(j.value("azscale", 1.0));
			gauss->setElevationScale(j.value("elscale", 1.0));
		}
		else if (auto* horn = dynamic_cast<antenna::SquareHorn*>(ant))
		{
			horn->setDimension(j.value("diameter", 0.5));
		}
		else if (auto* para = dynamic_cast<antenna::Parabolic*>(ant))
		{
			para->setDiameter(j.value("diameter", 0.5));
		}
		else if (auto* xml = dynamic_cast<antenna::XmlAntenna*>(ant))
		{
			if (xml->getFilename() != j.value("filename", ""))
			{
				world.replace(parse_required_antenna());
			}
		}
		else if (auto* h5 = dynamic_cast<antenna::H5Antenna*>(ant))
		{
			if (h5->getFilename() != j.value("filename", ""))
			{
				world.replace(parse_required_antenna());
			}
		}
	}

	void update_platform_paths_from_json(const nlohmann::json& j, radar::Platform* plat)
	{
		if (j.contains("motionpath"))
		{
			auto path = std::make_unique<math::Path>();
			j.at("motionpath").get_to(*path);
			plat->setMotionPath(std::move(path));
		}
		if (j.contains("rotationpath"))
		{
			auto rot_path = std::make_unique<math::RotationPath>();
			const auto& rotation_json = j.at("rotationpath");
			rot_path->setInterp(rotation_json.at("interpolation").get<math::RotationPath::InterpType>());
			unsigned waypoint_index = 0;
			for (const auto& waypoint_json : rotation_json.at("rotationwaypoints"))
			{
				const RealType azimuth = waypoint_json.at("azimuth").get<RealType>();
				const RealType elevation = waypoint_json.at("elevation").get<RealType>();
				const RealType time = waypoint_json.at("time").get<RealType>();
				const std::string owner =
					std::format("platform '{}' rotation waypoint {}", plat->getName(), waypoint_index);

				rotation_warning_utils::maybe_warn_about_rotation_value(azimuth, params::rotationAngleUnit(),
																		rotation_warning_utils::ValueKind::Angle,
																		"JSON", owner, "azimuth");
				rotation_warning_utils::maybe_warn_about_rotation_value(elevation, params::rotationAngleUnit(),
																		rotation_warning_utils::ValueKind::Angle,
																		"JSON", owner, "elevation");

				rot_path->addCoord(rotation_angle_utils::external_rotation_to_internal(azimuth, elevation, time,
																					   params::rotationAngleUnit()));
				++waypoint_index;
			}
			rot_path->finalize();
			plat->setRotationPath(std::move(rot_path));
		}
		else if (j.contains("fixedrotation"))
		{
			auto rot_path = std::make_unique<math::RotationPath>();
			const auto& fixed_json = j.at("fixedrotation");
			const RealType start_az_deg = fixed_json.at("startazimuth").get<RealType>();
			const RealType start_el_deg = fixed_json.at("startelevation").get<RealType>();
			const RealType rate_az_deg_s = fixed_json.at("azimuthrate").get<RealType>();
			const RealType rate_el_deg_s = fixed_json.at("elevationrate").get<RealType>();
			const std::string owner = std::format("platform '{}' fixedrotation", plat->getName());

			rotation_warning_utils::maybe_warn_about_rotation_value(start_az_deg, params::rotationAngleUnit(),
																	rotation_warning_utils::ValueKind::Angle, "JSON",
																	owner, "startazimuth");
			rotation_warning_utils::maybe_warn_about_rotation_value(start_el_deg, params::rotationAngleUnit(),
																	rotation_warning_utils::ValueKind::Angle, "JSON",
																	owner, "startelevation");
			rotation_warning_utils::maybe_warn_about_rotation_value(rate_az_deg_s, params::rotationAngleUnit(),
																	rotation_warning_utils::ValueKind::Rate, "JSON",
																	owner, "azimuthrate");
			rotation_warning_utils::maybe_warn_about_rotation_value(rate_el_deg_s, params::rotationAngleUnit(),
																	rotation_warning_utils::ValueKind::Rate, "JSON",
																	owner, "elevationrate");

			const auto start = serial::rotation_angle_utils::external_rotation_to_internal(
				start_az_deg, start_el_deg, 0.0, params::rotationAngleUnit());
			const auto rate = serial::rotation_angle_utils::external_rotation_rate_to_internal(
				rate_az_deg_s, rate_el_deg_s, 0.0, params::rotationAngleUnit());
			rot_path->setConstantRate(start, rate);
			rot_path->finalize();
			plat->setRotationPath(std::move(rot_path));
		}
	}

	void update_transmitter_from_json(const nlohmann::json& j, radar::Transmitter* tx, core::World& world,
									  std::mt19937& /*masterSeeder*/)
	{
		if (j.contains("name"))
			tx->setName(j.at("name").get<std::string>());

		if (j.contains("pulsed_mode"))
		{
			tx->setMode(radar::OperationMode::PULSED_MODE);
			tx->setPrf(j.at("pulsed_mode").value("prf", 0.0));
		}
		else if (j.contains("fmcw_mode"))
		{
			tx->setMode(radar::OperationMode::FMCW_MODE);
		}
		else if (j.contains("cw_mode"))
		{
			tx->setMode(radar::OperationMode::CW_MODE);
		}

		if (j.contains("waveform"))
		{
			auto id = parse_json_id(j, "waveform", "Transmitter");
			auto* wf = world.findWaveform(id);
			if (wf == nullptr)
				throw std::runtime_error("Waveform ID " + std::to_string(id) + " not found.");
			validate_fmcw_waveform(*wf, "Waveform '" + wf->getName() + "'");
			validate_waveform_mode_match(*wf, tx->getMode(), "Transmitter '" + tx->getName() + "'");
			tx->setWave(wf);
		}

		if (j.contains("antenna"))
		{
			auto id = parse_json_id(j, "antenna", "Transmitter");
			auto* ant = world.findAntenna(id);
			if (ant == nullptr)
				throw std::runtime_error("Antenna ID " + std::to_string(id) + " not found.");
			tx->setAntenna(ant);
		}

		if (j.contains("timing"))
		{
			auto timing_id = parse_json_id(j, "timing", "Transmitter");
			if (auto* const timing_proto = world.findTiming(timing_id))
			{
				unsigned seed = tx->getTiming() ? tx->getTiming()->getSeed() : 0;
				auto timing = std::make_shared<timing::Timing>(timing_proto->getName(), seed, timing_proto->getId());
				timing->initializeModel(timing_proto);
				tx->setTiming(timing);
			}
			else
			{
				throw std::runtime_error("Timing ID " + std::to_string(timing_id) + " not found.");
			}
		}
		if (j.contains("schedule"))
		{
			auto raw = j.at("schedule").get<std::vector<radar::SchedulePeriod>>();
			RealType pri = 0.0;
			if (tx->getMode() == radar::OperationMode::PULSED_MODE)
				pri = 1.0 / tx->getPrf();
			auto schedule = radar::processRawSchedule(std::move(raw), tx->getName(),
													  tx->getMode() == radar::OperationMode::PULSED_MODE, pri);
			if (const auto* fmcw = tx->getFmcwSignal(); fmcw != nullptr)
			{
				validate_fmcw_schedule(schedule, *fmcw, "Transmitter '" + tx->getName() + "'");
			}
			tx->setSchedule(std::move(schedule));
		}
		if (tx->getSignal() != nullptr)
		{
			validate_fmcw_waveform(*tx->getSignal(), "Waveform '" + tx->getSignal()->getName() + "'");
			validate_waveform_mode_match(*tx->getSignal(), tx->getMode(), "Transmitter '" + tx->getName() + "'");
		}
	}

	void update_receiver_from_json(const nlohmann::json& j, radar::Receiver* rx, core::World& world,
								   std::mt19937& /*masterSeeder*/)
	{
		if (j.contains("name"))
			rx->setName(j.at("name").get<std::string>());

		if (j.contains("pulsed_mode"))
		{
			rx->setMode(radar::OperationMode::PULSED_MODE);
			const auto& mode_json = j.at("pulsed_mode");
			rx->setWindowProperties(mode_json.value("window_length", 0.0), mode_json.value("prf", 0.0),
									mode_json.value("window_skip", 0.0));
		}
		else if (j.contains("fmcw_mode"))
		{
			rx->setMode(radar::OperationMode::FMCW_MODE);
		}
		else if (j.contains("cw_mode"))
		{
			rx->setMode(radar::OperationMode::CW_MODE);
		}

		if (j.contains("noise_temp"))
			rx->setNoiseTemperature(j.value("noise_temp", 0.0));

		if (j.contains("nodirect"))
		{
			if (j.value("nodirect", false))
				rx->setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
			else
				rx->clearFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
		}
		if (j.contains("nopropagationloss"))
		{
			if (j.value("nopropagationloss", false))
				rx->setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
			else
				rx->clearFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
		}

		if (j.contains("antenna"))
		{
			auto id = parse_json_id(j, "antenna", "Receiver");
			auto* ant = world.findAntenna(id);
			if (ant == nullptr)
				throw std::runtime_error("Antenna ID " + std::to_string(id) + " not found.");
			rx->setAntenna(ant);
		}

		if (j.contains("timing"))
		{
			auto timing_id = parse_json_id(j, "timing", "Receiver");
			if (auto* const timing_proto = world.findTiming(timing_id))
			{
				unsigned seed = rx->getTiming() ? rx->getTiming()->getSeed() : 0;
				auto timing = std::make_shared<timing::Timing>(timing_proto->getName(), seed, timing_proto->getId());
				timing->initializeModel(timing_proto);
				rx->setTiming(timing);
			}
			else
			{
				throw std::runtime_error("Timing ID " + std::to_string(timing_id) + " not found.");
			}
		}
		if (j.contains("schedule"))
		{
			auto raw = j.at("schedule").get<std::vector<radar::SchedulePeriod>>();
			RealType pri = 0.0;
			if (rx->getMode() == radar::OperationMode::PULSED_MODE)
				pri = 1.0 / rx->getWindowPrf();
			rx->setSchedule(radar::processRawSchedule(std::move(raw), rx->getName(),
													  rx->getMode() == radar::OperationMode::PULSED_MODE, pri));
		}
	}

	void update_monostatic_from_json(const nlohmann::json& j, radar::Transmitter* tx, radar::Receiver* rx,
									 core::World& world, std::mt19937& masterSeeder)
	{
		update_transmitter_from_json(j, tx, world, masterSeeder);

		if (j.contains("name"))
			rx->setName(j.at("name").get<std::string>());
		rx->setMode(tx->getMode());
		if (rx->getMode() == radar::OperationMode::PULSED_MODE && j.contains("pulsed_mode"))
		{
			const auto& mode_json = j.at("pulsed_mode");
			rx->setWindowProperties(mode_json.value("window_length", 0.0), tx->getPrf(),
									mode_json.value("window_skip", 0.0));
		}
		if (j.contains("noise_temp"))
			rx->setNoiseTemperature(j.value("noise_temp", 0.0));
		if (j.contains("nodirect"))
		{
			if (j.value("nodirect", false))
				rx->setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
			else
				rx->clearFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
		}
		if (j.contains("nopropagationloss"))
		{
			if (j.value("nopropagationloss", false))
				rx->setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
			else
				rx->clearFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
		}
		if (j.contains("antenna"))
			rx->setAntenna(world.findAntenna(parse_json_id(j, "antenna", "Monostatic")));
		if (j.contains("timing"))
		{
			auto timing_id = parse_json_id(j, "timing", "Monostatic");
			if (auto* const timing_proto = world.findTiming(timing_id))
			{
				unsigned seed = rx->getTiming() ? rx->getTiming()->getSeed() : 0;
				auto shared_timing =
					std::make_shared<timing::Timing>(timing_proto->getName(), seed, timing_proto->getId());
				shared_timing->initializeModel(timing_proto);
				tx->setTiming(shared_timing);
				rx->setTiming(shared_timing);
			}
			else
			{
				throw std::runtime_error("Timing ID " + std::to_string(timing_id) + " not found.");
			}
		}
		if (j.contains("schedule"))
		{
			auto raw = j.at("schedule").get<std::vector<radar::SchedulePeriod>>();
			RealType pri = 0.0;
			if (tx->getMode() == radar::OperationMode::PULSED_MODE)
				pri = 1.0 / tx->getPrf();
			auto processed_schedule = radar::processRawSchedule(
				std::move(raw), tx->getName(), tx->getMode() == radar::OperationMode::PULSED_MODE, pri);
			if (const auto* fmcw = tx->getFmcwSignal(); fmcw != nullptr)
			{
				validate_fmcw_schedule(processed_schedule, *fmcw, "Monostatic '" + tx->getName() + "'");
			}
			tx->setSchedule(processed_schedule);
			rx->setSchedule(processed_schedule);
		}
		if (tx->getSignal() != nullptr)
		{
			validate_fmcw_waveform(*tx->getSignal(), "Waveform '" + tx->getSignal()->getName() + "'");
			validate_waveform_mode_match(*tx->getSignal(), tx->getMode(), "Monostatic '" + tx->getName() + "'");
		}
	}

	void update_target_from_json(const nlohmann::json& j, radar::Target* existing_tgt, core::World& world,
								 std::mt19937& /*masterSeeder*/)
	{
		auto* plat = existing_tgt->getPlatform();
		const auto& rcs_json = j.at("rcs");
		const auto rcs_type = rcs_json.at("type").get<std::string>();
		std::unique_ptr<radar::Target> target_obj;

		const auto target_id = existing_tgt->getId();
		const auto name = j.value("name", existing_tgt->getName());
		unsigned seed = existing_tgt->getSeed();

		if (rcs_type == "isotropic")
		{
			target_obj = radar::createIsoTarget(plat, name, rcs_json.value("value", 1.0), seed, target_id);
		}
		else if (rcs_type == "file")
		{
			const auto filename = rcs_json.value("filename", "");
			target_obj = radar::createFileTarget(plat, name, filename, seed, target_id);
		}
		else
		{
			throw std::runtime_error("Unsupported target RCS type: " + rcs_type);
		}

		if (j.contains("model"))
		{
			const auto& model_json = j.at("model");
			const auto model_type = model_json.at("type").get<std::string>();
			if (model_type == "chisquare" || model_type == "gamma")
			{
				auto model =
					std::make_unique<radar::RcsChiSquare>(target_obj->getRngEngine(), model_json.value("k", 1.0));
				target_obj->setFluctuationModel(std::move(model));
			}
			else if (model_type == "constant")
			{
				target_obj->setFluctuationModel(std::make_unique<radar::RcsConst>());
			}
		}

		world.replace(std::move(target_obj));
	}

	void update_timing_from_json(const nlohmann::json& j, core::World& world, const SimId id)
	{
		auto* existing = world.findTiming(id);
		if (existing == nullptr)
		{
			throw std::runtime_error("Timing ID " + std::to_string(id) + " not found.");
		}

		auto patched = j;
		if (!patched.contains("name"))
		{
			patched["name"] = existing->getName();
		}

		world.replace(parse_timing_from_json(patched, id));
	}

	nlohmann::json world_to_json(const core::World& world)
	{
		nlohmann::json sim_json;

		sim_json["name"] = params::params.simulation_name;
		sim_json["parameters"] = params::params;

		sim_json["waveforms"] = nlohmann::json::array();
		for (const auto& waveform : world.getWaveforms() | std::views::values)
		{
			sim_json["waveforms"].push_back(*waveform);
		}

		sim_json["antennas"] = nlohmann::json::array();
		for (const auto& antenna : world.getAntennas() | std::views::values)
		{
			sim_json["antennas"].push_back(*antenna);
		}

		sim_json["timings"] = nlohmann::json::array();
		for (const auto& timing : world.getTimings() | std::views::values)
		{
			sim_json["timings"].push_back(*timing);
		}

		sim_json["platforms"] = nlohmann::json::array();
		for (const auto& p : world.getPlatforms())
		{
			sim_json["platforms"].push_back(serialize_platform(p.get(), world));
		}

		return {{"simulation", sim_json}};
	}

	void json_to_world(const nlohmann::json& j, core::World& world, std::mt19937& masterSeeder)
	{
		// 1. Clear the existing world state. This function always performs a full
		//    replacement to ensure the C++ state is a perfect mirror of the UI state.
		world.clear();

		const auto& sim = j.at("simulation");

		parse_parameters(sim, masterSeeder);
		parse_assets(sim, world);

		// 3. Restore platforms and their components.
		if (sim.contains("platforms"))
		{
			TimingInstanceMap timing_instances;
			for (const auto& plat_json : sim.at("platforms"))
			{
				parse_platform(plat_json, world, masterSeeder, timing_instances);
			}
		}

		// 4. Finalize world state after all objects are loaded.

		// Prepare CW receiver buffers before starting simulation
		const RealType start_time = params::startTime();
		const RealType end_time = params::endTime();
		const RealType dt_sim = 1.0 / (params::rate() * params::oversampleRatio());
		const auto num_samples = static_cast<size_t>(std::ceil((end_time - start_time) / dt_sim));

		for (const auto& receiver : world.getReceivers())
		{
			if (receiver->getMode() == radar::OperationMode::CW_MODE ||
				receiver->getMode() == radar::OperationMode::FMCW_MODE)
			{
				receiver->prepareStreamingData(num_samples);
			}
		}

		// Schedule initial events after all objects are loaded.
		world.scheduleInitialEvents();
	}
}
