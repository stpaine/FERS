// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "xml_parser_utils.h"

#include <GeographicLib/UTMUPS.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <format>
#include <limits>
#include <optional>
#include <string_view>

#include "antenna/antenna_factory.h"
#include "core/config.h"
#include "core/logging.h"
#include "core/world.h"
#include "fers_xml_dtd.h"
#include "fers_xml_xsd.h"
#include "math/coord.h"
#include "math/geometry_ops.h"
#include "math/path.h"
#include "math/rotation_path.h"
#include "radar/platform.h"
#include "radar/radar_obj.h"
#include "radar/receiver.h"
#include "radar/target.h"
#include "radar/transmitter.h"
#include "serial/fmcw_validation.h"
#include "serial/rotation_angle_utils.h"
#include "serial/rotation_warning_utils.h"
#include "serial/waveform_factory.h"
#include "signal/radar_signal.h"
#include "timing/prototype_timing.h"
#include "timing/timing.h"

namespace fs = std::filesystem;

namespace serial::xml_parser_utils
{
	namespace
	{
		/// Parses the mutually exclusive radar operation mode elements.
		radar::OperationMode parse_mode_elements(const XmlElement& parent, const std::string& owner)
		{
			std::optional<radar::OperationMode> selected_mode;
			const auto select_mode = [&](const char* const element_name, const radar::OperationMode mode)
			{
				if (!parent.childElement(element_name, 0).isValid())
				{
					return;
				}
				if (selected_mode.has_value())
				{
					throw XmlException(owner +
									   " must specify exactly one radar mode (<pulsed_mode>, <cw_mode>, or "
									   "<fmcw_mode>).");
				}
				selected_mode = mode;
			};

			select_mode("pulsed_mode", radar::OperationMode::PULSED_MODE);
			select_mode("cw_mode", radar::OperationMode::CW_MODE);
			select_mode("fmcw_mode", radar::OperationMode::FMCW_MODE);
			if (!selected_mode.has_value())
			{
				throw XmlException(owner +
								   " must specify exactly one radar mode (<pulsed_mode>, <cw_mode>, or <fmcw_mode>).");
			}
			return *selected_mode;
		}

		/// Throws an XML validation exception with the provided message.
		void throw_xml_validation_error(const std::string& message) { throw XmlException(message); }

		/// Validates an FMCW waveform while adapting validation errors to XmlException.
		void validate_fmcw_waveform(const fers_signal::RadarSignal& wave, const std::string& owner)
		{
			serial::fmcw_validation::validateWaveform(wave, owner, throw_xml_validation_error);
		}

		/// Validates waveform/mode compatibility while adapting validation errors to XmlException.
		void validate_waveform_mode_match(const fers_signal::RadarSignal& wave, const radar::OperationMode mode,
										  const std::string& owner)
		{
			serial::fmcw_validation::validateWaveformModeMatch(wave, mode, owner, throw_xml_validation_error);
		}

		/// Validates an FMCW schedule while adapting validation errors to XmlException.
		void validate_fmcw_schedule(const std::vector<radar::SchedulePeriod>& schedule,
									const fers_signal::RadarSignal& wave, const std::string& owner)
		{
			serial::fmcw_validation::validateSchedule(schedule, wave, owner, throw_xml_validation_error);
		}

		/// Draws the next unsigned seed from the master random generator.
		[[nodiscard]] unsigned next_seed(std::mt19937& master_seeder)
		{
			static_assert(std::mt19937::max() <= std::numeric_limits<unsigned>::max(),
						  "std::mt19937 output must fit into unsigned seeds.");
			return static_cast<unsigned>(master_seeder());
		}

		/// Resolves or instantiates a shared timing instance by prototype SimId.
		std::shared_ptr<timing::Timing> resolve_timing_instance(const SimId timing_id, ParserContext& ctx,
																const std::string& owner)
		{
			if (const auto it = ctx.timing_instances.find(timing_id); it != ctx.timing_instances.end())
			{
				return it->second;
			}

			const timing::PrototypeTiming* proto = ctx.world->findTiming(timing_id);
			if (proto == nullptr)
			{
				throw XmlException("Timing ID '" + std::to_string(timing_id) + "' not found for " + owner + "'");
			}

			auto timing_obj =
				std::make_shared<timing::Timing>(proto->getName(), next_seed(*ctx.master_seeder), proto->getId());
			timing_obj->initializeModel(proto);
			ctx.timing_instances.emplace(timing_id, timing_obj);
			return timing_obj;
		}
	}

	RealType get_child_real_type(const XmlElement& element, const std::string& elementName)
	{
		const std::string text = element.childElement(elementName, 0).getText();
		if (text.empty())
		{
			throw XmlException("Element " + elementName + " is empty!");
		}
		return std::stod(text);
	}

	bool get_attribute_bool(const XmlElement& element, const std::string& attributeName, const bool defaultVal)
	{
		const auto attr_value = XmlElement::getOptionalAttribute(element, attributeName);
		if (!attr_value.has_value())
		{
			LOG(logging::Level::DEBUG, "Attribute '{}' not specified. Defaulting to {}.", attributeName, defaultVal);
			return defaultVal;
		}
		if (*attr_value == "true")
		{
			return true;
		}
		if (*attr_value == "false")
		{
			return false;
		}

		LOG(logging::Level::WARNING, "Invalid boolean value '{}' for attribute '{}'. Defaulting to {}.", *attr_value,
			attributeName, defaultVal);
		return defaultVal;
	}

	SimId assign_id_from_attribute(const std::string& owner, ObjectType type)
	{
		const SimId id = SimIdGenerator::instance().generateId(type);
		LOG(logging::Level::TRACE, "Assigned ID {} to {} (generated)", id, owner);
		return id;
	}

	SimId resolve_reference_id(const XmlElement& element, const std::string& attributeName, const std::string& owner,
							   const std::unordered_map<std::string, SimId>& name_map)
	{
		const std::string value = XmlElement::getSafeAttribute(element, attributeName);
		if (value.empty())
		{
			throw XmlException("Missing " + attributeName + " for " + owner + ".");
		}
		const auto it = name_map.find(value);
		if (it != name_map.end())
		{
			return it->second;
		}
		throw XmlException("Unknown " + attributeName + " '" + value + "' for " + owner + ".");
	}

	std::vector<radar::SchedulePeriod> parseSchedule(const XmlElement& parent, const std::string& parentName,
													 const bool isPulsed, const RealType pri)
	{
		std::vector<radar::SchedulePeriod> raw_periods;
		if (const XmlElement schedule_element = parent.childElement("schedule", 0); schedule_element.isValid())
		{
			unsigned p_idx = 0;
			while (true)
			{
				XmlElement period_element = schedule_element.childElement("period", p_idx++);
				if (!period_element.isValid())
				{
					break;
				}
				try
				{
					const RealType start = std::stod(XmlElement::getSafeAttribute(period_element, "start"));
					const RealType end = std::stod(XmlElement::getSafeAttribute(period_element, "end"));
					raw_periods.push_back({start, end});
				}
				catch (const std::exception& e)
				{
					LOG(logging::Level::WARNING, "Failed to parse schedule period for '{}': {}", parentName, e.what());
				}
			}
		}
		return radar::processRawSchedule(std::move(raw_periods), parentName, isPulsed, pri);
	}

	void parseParameters(const XmlElement& parameters, params::Parameters& params_out)
	{
		params_out.start = get_child_real_type(parameters, "starttime");
		params_out.end = get_child_real_type(parameters, "endtime");
		LOG(logging::Level::INFO, "Simulation time set from {:.5f} to {:.5f} seconds", params_out.start,
			params_out.end);

		params_out.rate = get_child_real_type(parameters, "rate");
		if (params_out.rate <= 0)
		{
			throw std::runtime_error("Sampling rate must be > 0");
		}
		LOG(logging::Level::DEBUG, "Sample rate set to: {:.5f}", params_out.rate);

		const auto parse_unsigned_parameter = [&](const std::string_view param_name, const RealType raw_value)
		{
			if (!std::isfinite(raw_value))
			{
				throw XmlException(std::format("Parameter '{}' must be finite.", param_name));
			}
			if (raw_value < 0.0)
			{
				throw XmlException(std::format("Parameter '{}' must be non-negative.", param_name));
			}

			const RealType floored_value = std::floor(raw_value);
			if (floored_value > static_cast<RealType>(std::numeric_limits<unsigned>::max()))
			{
				throw XmlException(std::format("Parameter '{}' exceeds the supported unsigned range.", param_name));
			}

			return static_cast<unsigned>(floored_value);
		};

		auto set_optional_real_parameter = [&](const std::string& param_name, const RealType default_value, auto setter)
		{
			if (!parameters.childElement(param_name, 0).isValid())
			{
				LOG(logging::Level::DEBUG, "Parameter '{}' not specified. Using default value {}.", param_name,
					default_value);
				return;
			}

			setter(get_child_real_type(parameters, param_name));
		};

		auto set_optional_unsigned_parameter =
			[&](const std::string& param_name, const unsigned default_value, auto setter)
		{
			if (!parameters.childElement(param_name, 0).isValid())
			{
				LOG(logging::Level::DEBUG, "Parameter '{}' not specified. Using default value {}.", param_name,
					default_value);
				return;
			}

			setter(parse_unsigned_parameter(param_name, get_child_real_type(parameters, param_name)));
		};

		set_optional_real_parameter("c", params::Parameters::DEFAULT_C,
									[&](const RealType value)
									{
										params_out.c = value;
										LOG(logging::Level::INFO, "Propagation speed (c) set to: {:.5f}", value);
									});

		set_optional_real_parameter("simSamplingRate", 1000.0,
									[&](const RealType value)
									{
										params_out.sim_sampling_rate = value;
										LOG(logging::Level::DEBUG, "Simulation sampling rate set to: {:.5f} Hz", value);
									});

		if (parameters.childElement("randomseed", 0).isValid())
		{
			const auto seed = parse_unsigned_parameter("randomseed", get_child_real_type(parameters, "randomseed"));
			params_out.random_seed = seed;
			LOG(logging::Level::DEBUG, "Random seed set to: {}", seed);
		}

		set_optional_unsigned_parameter("adc_bits", 0,
										[&](const unsigned value)
										{
											params_out.adc_bits = value;
											LOG(logging::Level::DEBUG, "ADC quantization bits set to: {}", value);
										});

		set_optional_unsigned_parameter("oversample", 1,
										[&](const unsigned value)
										{
											params::validateOversampleRatio(value);
											params_out.oversample_ratio = value;
											LOG(logging::Level::DEBUG, "Oversampling enabled with ratio: {}", value);
										});

		try
		{
			const auto unit_token = parameters.childElement("rotationangleunit", 0).getText();
			if (!unit_token.empty())
			{
				if (const auto unit = params::rotationAngleUnitFromToken(unit_token))
				{
					params_out.rotation_angle_unit = *unit;
				}
				else
				{
					throw XmlException("Unsupported rotation angle unit '" + unit_token + "'.");
				}
			}
		}
		catch (const XmlException&)
		{
		}

		bool origin_set = false;
		if (const XmlElement origin_element = parameters.childElement("origin", 0); origin_element.isValid())
		{
			try
			{
				params_out.origin_latitude = std::stod(XmlElement::getSafeAttribute(origin_element, "latitude"));
				params_out.origin_longitude = std::stod(XmlElement::getSafeAttribute(origin_element, "longitude"));
				if (const auto altitude = XmlElement::getOptionalAttribute(origin_element, "altitude"))
				{
					params_out.origin_altitude = std::stod(*altitude);
				}
				else
				{
					params_out.origin_altitude = 0.0;
					LOG(logging::Level::DEBUG, "Origin altitude not specified. Defaulting to 0.");
				}
				origin_set = true;
				LOG(logging::Level::INFO, "Origin set to lat: {}, lon: {}, alt: {}", params_out.origin_latitude,
					params_out.origin_longitude, params_out.origin_altitude);
			}
			catch (const std::exception& e)
			{
				LOG(logging::Level::WARNING, "Could not parse origin from XML, using defaults. Error: {}", e.what());
			}
		}

		if (const XmlElement cs_element = parameters.childElement("coordinatesystem", 0); cs_element.isValid())
		{
			try
			{
				const std::string frame_str = XmlElement::getSafeAttribute(cs_element, "frame");

				if (frame_str == "UTM")
				{
					params_out.coordinate_frame = params::CoordinateFrame::UTM;
					params_out.utm_zone = std::stoi(XmlElement::getSafeAttribute(cs_element, "zone"));
					const std::string hem_str = XmlElement::getSafeAttribute(cs_element, "hemisphere");

					if (params_out.utm_zone < GeographicLib::UTMUPS::MINUTMZONE ||
						params_out.utm_zone > GeographicLib::UTMUPS::MAXUTMZONE)
					{
						throw XmlException("UTM zone " + std::to_string(params_out.utm_zone) +
										   " is invalid; must be in [1, 60].");
					}
					if (hem_str == "N" || hem_str == "n")
					{
						params_out.utm_north_hemisphere = true;
					}
					else if (hem_str == "S" || hem_str == "s")
					{
						params_out.utm_north_hemisphere = false;
					}
					else
					{
						throw XmlException("UTM hemisphere '" + hem_str + "' is invalid; must be 'N' or 'S'.");
					}
					LOG(logging::Level::INFO, "Coordinate system set to UTM, zone {}{}", params_out.utm_zone,
						params_out.utm_north_hemisphere ? 'N' : 'S');
				}
				else if (frame_str == "ECEF")
				{
					params_out.coordinate_frame = params::CoordinateFrame::ECEF;
					LOG(logging::Level::INFO, "Coordinate system set to ECEF.");
				}
				else if (frame_str == "ENU")
				{
					params_out.coordinate_frame = params::CoordinateFrame::ENU;
					if (!origin_set)
					{
						LOG(logging::Level::WARNING,
							"ENU frame specified but no <origin> tag found. Using default origin at UCT.");
					}
					LOG(logging::Level::INFO, "Coordinate system set to ENU local tangent plane.");
				}
				else
				{
					throw XmlException("Unsupported coordinate frame: " + frame_str);
				}
			}
			catch (const std::exception& e)
			{
				LOG(logging::Level::WARNING, "Could not parse <coordinatesystem> from XML: {}. Defaulting to ENU.",
					e.what());
				params_out.coordinate_frame = params::CoordinateFrame::ENU;
				params_out.utm_zone = 0;
				params_out.utm_north_hemisphere = true;
			}
		}
	}

	void parseWaveform(const XmlElement& waveform, ParserContext& ctx)
	{
		const std::string name = XmlElement::getSafeAttribute(waveform, "name");
		const SimId id = assign_id_from_attribute("waveform '" + name + "'", ObjectType::Waveform);

		const auto power = get_child_real_type(waveform, "power");
		const auto carrier = get_child_real_type(waveform, "carrier_frequency");

		if (const XmlElement pulsed_file = waveform.childElement("pulsed_from_file", 0); pulsed_file.isValid())
		{
			const std::string filename_str = XmlElement::getSafeAttribute(pulsed_file, "filename");
			fs::path pulse_path(filename_str);

			if (!fs::exists(pulse_path))
			{
				pulse_path = ctx.base_dir / filename_str;
			}

			// Defer to dependency-injected file loader
			auto wave = ctx.loaders.loadWaveform(name, pulse_path, power, carrier, id);
			ctx.world->add(std::move(wave));
		}
		else if (waveform.childElement("cw", 0).isValid())
		{
			auto cw_signal = std::make_unique<fers_signal::CwSignal>();
			auto wave = std::make_unique<fers_signal::RadarSignal>(
				name, power, carrier, ctx.parameters.end - ctx.parameters.start, std::move(cw_signal), id);
			ctx.world->add(std::move(wave));
		}
		else if (const XmlElement fmcw_element = waveform.childElement("fmcw_linear_chirp", 0); fmcw_element.isValid())
		{
			const auto direction =
				fers_signal::parseFmcwChirpDirection(XmlElement::getSafeAttribute(fmcw_element, "direction"));
			const RealType chirp_bandwidth = get_child_real_type(fmcw_element, "chirp_bandwidth");
			const RealType chirp_duration = get_child_real_type(fmcw_element, "chirp_duration");
			const RealType chirp_period = get_child_real_type(fmcw_element, "chirp_period");

			RealType start_frequency_offset = 0.0;
			if (const auto start_offset = fmcw_element.childElement("start_frequency_offset", 0);
				start_offset.isValid())
			{
				start_frequency_offset = get_child_real_type(fmcw_element, "start_frequency_offset");
			}

			std::optional<std::size_t> chirp_count;
			if (const auto chirp_count_element = fmcw_element.childElement("chirp_count", 0);
				chirp_count_element.isValid())
			{
				const RealType raw_count = get_child_real_type(fmcw_element, "chirp_count");
				if (raw_count <= 0.0 || std::floor(raw_count) != raw_count)
				{
					throw XmlException("Waveform '" + name + "' has an invalid chirp_count.");
				}
				chirp_count = static_cast<std::size_t>(raw_count);
			}

			auto fmcw_signal = std::make_unique<fers_signal::FmcwChirpSignal>(
				chirp_bandwidth, chirp_duration, chirp_period, start_frequency_offset, chirp_count, direction);
			// RadarSignal length is the active chirp duration, not T_rep. The repeat period only spaces chirps.
			auto wave = std::make_unique<fers_signal::RadarSignal>(name, power, carrier, chirp_duration,
																   std::move(fmcw_signal), id);
			validate_fmcw_waveform(*wave, "Waveform '" + name + "'");
			ctx.world->add(std::move(wave));
		}
		else if (const XmlElement fmcw_triangle_element = waveform.childElement("fmcw_triangle", 0);
				 fmcw_triangle_element.isValid())
		{
			const RealType chirp_bandwidth = get_child_real_type(fmcw_triangle_element, "chirp_bandwidth");
			const RealType chirp_duration = get_child_real_type(fmcw_triangle_element, "chirp_duration");

			RealType start_frequency_offset = 0.0;
			if (const auto start_offset = fmcw_triangle_element.childElement("start_frequency_offset", 0);
				start_offset.isValid())
			{
				start_frequency_offset = get_child_real_type(fmcw_triangle_element, "start_frequency_offset");
			}

			std::optional<std::size_t> triangle_count;
			if (const auto triangle_count_element = fmcw_triangle_element.childElement("triangle_count", 0);
				triangle_count_element.isValid())
			{
				const RealType raw_count = get_child_real_type(fmcw_triangle_element, "triangle_count");
				if (raw_count <= 0.0 || std::floor(raw_count) != raw_count)
				{
					throw XmlException("Waveform '" + name + "' has an invalid triangle_count.");
				}
				triangle_count = static_cast<std::size_t>(raw_count);
			}

			auto fmcw_signal = std::make_unique<fers_signal::FmcwTriangleSignal>(
				chirp_bandwidth, chirp_duration, start_frequency_offset, triangle_count);
			auto wave = std::make_unique<fers_signal::RadarSignal>(
				name, power, carrier, fmcw_signal->getTrianglePeriod(), std::move(fmcw_signal), id);
			validate_fmcw_waveform(*wave, "Waveform '" + name + "'");
			ctx.world->add(std::move(wave));
		}
		else
		{
			LOG(logging::Level::FATAL, "Unsupported waveform type for '{}'", name);
			throw XmlException("Unsupported waveform type for '" + name + "'");
		}
	}

	void parseTiming(const XmlElement& timing, ParserContext& ctx)
	{
		const std::string name = XmlElement::getSafeAttribute(timing, "name");
		const SimId id = assign_id_from_attribute("timing '" + name + "'", ObjectType::Timing);
		const RealType freq = get_child_real_type(timing, "frequency");
		auto timing_obj = std::make_unique<timing::PrototypeTiming>(name, id);

		timing_obj->setFrequency(freq);

		unsigned noise_index = 0;
		while (true)
		{
			XmlElement noise_element = timing.childElement("noise_entry", noise_index++);
			if (!noise_element.isValid())
			{
				break;
			}

			timing_obj->setAlpha(get_child_real_type(noise_element, "alpha"),
								 get_child_real_type(noise_element, "weight"));
		}

		const auto set_optional_timing_parameter =
			[&](const std::string& element_name, const std::string& description, auto setter)
		{
			if (!timing.childElement(element_name, 0).isValid())
			{
				LOG(logging::Level::DEBUG, "Clock section '{}' does not specify {}.", name, description);
				return;
			}
			try
			{
				setter(get_child_real_type(timing, element_name));
			}
			catch (const XmlException&)
			{
				LOG(logging::Level::WARNING, "Clock section '{}' has an empty {}. Using default.", name, description);
			}
		};

		set_optional_timing_parameter("freq_offset", "frequency offset",
									  [&](const RealType value) { timing_obj->setFreqOffset(value); });
		set_optional_timing_parameter("random_freq_offset_stdev", "random frequency offset",
									  [&](const RealType value) { timing_obj->setRandomFreqOffsetStdev(value); });
		set_optional_timing_parameter("phase_offset", "phase offset",
									  [&](const RealType value) { timing_obj->setPhaseOffset(value); });
		set_optional_timing_parameter("random_phase_offset_stdev", "random phase offset",
									  [&](const RealType value) { timing_obj->setRandomPhaseOffsetStdev(value); });

		if (get_attribute_bool(timing, "synconpulse", false))
		{
			timing_obj->setSyncOnPulse();
		}

		ctx.world->add(std::move(timing_obj));
	}

	void parseAntenna(const XmlElement& antenna, ParserContext& ctx)
	{
		std::string name = XmlElement::getSafeAttribute(antenna, "name");
		const SimId id = assign_id_from_attribute("antenna '" + name + "'", ObjectType::Antenna);
		const std::string pattern = XmlElement::getSafeAttribute(antenna, "pattern");

		std::unique_ptr<antenna::Antenna> ant;

		LOG(logging::Level::DEBUG, "Adding antenna '{}' with pattern '{}'", name, pattern);
		if (pattern == "isotropic")
		{
			ant = std::make_unique<antenna::Isotropic>(name, id);
		}
		else if (pattern == "sinc")
		{
			ant = std::make_unique<antenna::Sinc>(name, get_child_real_type(antenna, "alpha"),
												  get_child_real_type(antenna, "beta"),
												  get_child_real_type(antenna, "gamma"), id);
		}
		else if (pattern == "gaussian")
		{
			ant = std::make_unique<antenna::Gaussian>(name, get_child_real_type(antenna, "azscale"),
													  get_child_real_type(antenna, "elscale"), id);
		}
		else if (pattern == "squarehorn")
		{
			ant = std::make_unique<antenna::SquareHorn>(name, get_child_real_type(antenna, "diameter"), id);
		}
		else if (pattern == "parabolic")
		{
			ant = std::make_unique<antenna::Parabolic>(name, get_child_real_type(antenna, "diameter"), id);
		}
		else if (pattern == "xml")
		{
			ant = ctx.loaders.loadXmlAntenna(name, XmlElement::getSafeAttribute(antenna, "filename"), id);
		}
		else if (pattern == "file")
		{
			ant = ctx.loaders.loadH5Antenna(name, XmlElement::getSafeAttribute(antenna, "filename"), id);
		}
		else
		{
			LOG(logging::Level::FATAL, "Unsupported antenna pattern: {}", pattern);
			throw XmlException("Unsupported antenna pattern: " + pattern);
		}

		if (!antenna.childElement("efficiency", 0).isValid())
		{
			LOG(logging::Level::DEBUG, "Antenna '{}' does not specify efficiency, assuming unity.", name);
		}
		else
		{
			try
			{
				ant->setEfficiencyFactor(get_child_real_type(antenna, "efficiency"));
			}
			catch (const XmlException&)
			{
				LOG(logging::Level::WARNING, "Antenna '{}' has an empty efficiency, assuming unity.", name);
			}
		}

		ctx.world->add(std::move(ant));
	}

	void parseMotionPath(const XmlElement& motionPath, radar::Platform* platform)
	{
		math::Path* path = platform->getMotionPath();
		if (const auto interp_value = XmlElement::getOptionalAttribute(motionPath, "interpolation"))
		{
			if (*interp_value == "linear")
			{
				path->setInterp(math::Path::InterpType::INTERP_LINEAR);
			}
			else if (*interp_value == "cubic")
			{
				path->setInterp(math::Path::InterpType::INTERP_CUBIC);
			}
			else if (*interp_value == "static")
			{
				path->setInterp(math::Path::InterpType::INTERP_STATIC);
			}
			else
			{
				LOG(logging::Level::ERROR, "Unsupported interpolation type: {} for platform {}. Defaulting to static",
					*interp_value, platform->getName());
				path->setInterp(math::Path::InterpType::INTERP_STATIC);
			}
		}
		else
		{
			LOG(logging::Level::DEBUG,
				"MotionPath interpolation type for platform {} not specified. Defaulting to static.",
				platform->getName());
			path->setInterp(math::Path::InterpType::INTERP_STATIC);
		}

		unsigned waypoint_index = 0;
		while (true)
		{
			XmlElement waypoint = motionPath.childElement("positionwaypoint", waypoint_index);
			if (!waypoint.isValid())
			{
				break;
			}

			try
			{
				math::Coord coord;
				coord.t = get_child_real_type(waypoint, "time");
				coord.pos = math::Vec3(get_child_real_type(waypoint, "x"), get_child_real_type(waypoint, "y"),
									   get_child_real_type(waypoint, "altitude"));
				path->addCoord(coord);
				LOG(logging::Level::TRACE, "Added waypoint {} to motion path for platform {}.", waypoint_index,
					platform->getName());
			}
			catch (const XmlException& e)
			{
				LOG(logging::Level::ERROR, "Failed to add waypoint to motion path. Discarding waypoint. {}", e.what());
			}

			waypoint_index++;
		}
		path->finalize();
	}

	void parseRotationPath(const XmlElement& rotation, radar::Platform* platform, const params::RotationAngleUnit unit)
	{
		math::RotationPath* path = platform->getRotationPath();
		try
		{
			if (const std::string interp = XmlElement::getSafeAttribute(rotation, "interpolation"); interp == "linear")
			{
				path->setInterp(math::RotationPath::InterpType::INTERP_LINEAR);
			}
			else if (interp == "cubic")
			{
				path->setInterp(math::RotationPath::InterpType::INTERP_CUBIC);
			}
			else if (interp == "static")
			{
				path->setInterp(math::RotationPath::InterpType::INTERP_STATIC);
			}
			else
			{
				throw XmlException("Unsupported interpolation type: " + interp);
			}
		}
		catch (XmlException&)
		{
			LOG(logging::Level::ERROR,
				"Failed to set RotationPath interpolation type for platform {}. Defaulting to static",
				platform->getName());
			path->setInterp(math::RotationPath::InterpType::INTERP_STATIC);
		}

		unsigned waypoint_index = 0;
		while (true)
		{
			XmlElement waypoint = rotation.childElement("rotationwaypoint", waypoint_index);
			if (!waypoint.isValid())
			{
				break;
			}

			try
			{
				LOG(logging::Level::TRACE, "Adding waypoint {} to rotation path for platform {}.", waypoint_index,
					platform->getName());

				const RealType az_deg = get_child_real_type(waypoint, "azimuth");
				const RealType el_deg = get_child_real_type(waypoint, "elevation");
				const RealType time = get_child_real_type(waypoint, "time");
				const std::string owner =
					std::format("platform '{}' rotation waypoint {}", platform->getName(), waypoint_index);

				rotation_warning_utils::maybe_warn_about_rotation_value(
					az_deg, unit, rotation_warning_utils::ValueKind::Angle, "XML", owner, "azimuth");
				rotation_warning_utils::maybe_warn_about_rotation_value(
					el_deg, unit, rotation_warning_utils::ValueKind::Angle, "XML", owner, "elevation");

				path->addCoord(rotation_angle_utils::external_rotation_to_internal(az_deg, el_deg, time, unit));
			}
			catch (const XmlException& e)
			{
				LOG(logging::Level::ERROR, "Failed to add waypoint to rotation path. Discarding waypoint. {}",
					e.what());
			}
			waypoint_index++;
		}
		path->finalize();
	}

	void parseFixedRotation(const XmlElement& rotation, radar::Platform* platform, const params::RotationAngleUnit unit)
	{
		math::RotationPath* path = platform->getRotationPath();
		try
		{
			const RealType start_az_deg = get_child_real_type(rotation, "startazimuth");
			const RealType start_el_deg = get_child_real_type(rotation, "startelevation");
			const RealType rate_az_deg_s = get_child_real_type(rotation, "azimuthrate");
			const RealType rate_el_deg_s = get_child_real_type(rotation, "elevationrate");
			const std::string owner = std::format("platform '{}' fixedrotation", platform->getName());

			rotation_warning_utils::maybe_warn_about_rotation_value(
				start_az_deg, unit, rotation_warning_utils::ValueKind::Angle, "XML", owner, "startazimuth");
			rotation_warning_utils::maybe_warn_about_rotation_value(
				start_el_deg, unit, rotation_warning_utils::ValueKind::Angle, "XML", owner, "startelevation");
			rotation_warning_utils::maybe_warn_about_rotation_value(
				rate_az_deg_s, unit, rotation_warning_utils::ValueKind::Rate, "XML", owner, "azimuthrate");
			rotation_warning_utils::maybe_warn_about_rotation_value(
				rate_el_deg_s, unit, rotation_warning_utils::ValueKind::Rate, "XML", owner, "elevationrate");
			const math::RotationCoord start =
				rotation_angle_utils::external_rotation_to_internal(start_az_deg, start_el_deg, 0.0, unit);
			const math::RotationCoord rate =
				rotation_angle_utils::external_rotation_rate_to_internal(rate_az_deg_s, rate_el_deg_s, 0.0, unit);

			path->setConstantRate(start, rate);
			LOG(logging::Level::DEBUG, "Added fixed rotation to platform {}", platform->getName());
		}
		catch (XmlException& e)
		{
			LOG(logging::Level::FATAL, "Failed to set fixed rotation for platform {}. {}", platform->getName(),
				e.what());
			throw XmlException("Failed to set fixed rotation for platform " + platform->getName());
		}
	}

	/// Parses a transmitter after its operation mode has already been determined.
	radar::Transmitter* parseTransmitterWithMode(const XmlElement& transmitter, radar::Platform* platform,
												 ParserContext& ctx, const ReferenceLookup& refs,
												 const radar::OperationMode mode)
	{
		const std::string name = XmlElement::getSafeAttribute(transmitter, "name");
		const SimId id = assign_id_from_attribute("transmitter '" + name + "'", ObjectType::Transmitter);
		const XmlElement pulsed_mode_element = transmitter.childElement("pulsed_mode", 0);
		const bool is_pulsed = mode == radar::OperationMode::PULSED_MODE;

		auto transmitter_obj = std::make_unique<radar::Transmitter>(platform, name, mode, id);

		const SimId waveform_id =
			resolve_reference_id(transmitter, "waveform", "transmitter '" + name + "'", *refs.waveforms);
		fers_signal::RadarSignal* wave = ctx.world->findWaveform(waveform_id);
		if (wave == nullptr)
		{
			throw XmlException("Waveform ID '" + std::to_string(waveform_id) + "' not found for transmitter '" + name +
							   "'");
		}
		validate_fmcw_waveform(*wave, "Waveform '" + wave->getName() + "'");
		validate_waveform_mode_match(*wave, mode, "Transmitter '" + name + "'");
		transmitter_obj->setWave(wave);

		if (is_pulsed)
		{
			transmitter_obj->setPrf(get_child_real_type(pulsed_mode_element, "prf"));
		}

		const SimId antenna_id =
			resolve_reference_id(transmitter, "antenna", "transmitter '" + name + "'", *refs.antennas);
		const antenna::Antenna* ant = ctx.world->findAntenna(antenna_id);
		if (ant == nullptr)
		{
			throw XmlException("Antenna ID '" + std::to_string(antenna_id) + "' not found for transmitter '" + name +
							   "'");
		}
		transmitter_obj->setAntenna(ant);

		const SimId timing_id =
			resolve_reference_id(transmitter, "timing", "transmitter '" + name + "'", *refs.timings);
		transmitter_obj->setTiming(resolve_timing_instance(timing_id, ctx, "transmitter '" + name + "'"));

		RealType pri = is_pulsed ? (1.0 / transmitter_obj->getPrf()) : 0.0;
		auto schedule = parseSchedule(transmitter, name, is_pulsed, pri);
		if (wave->isFmcwFamily())
		{
			validate_fmcw_schedule(schedule, *wave, "Transmitter '" + name + "'");
		}
		if (!schedule.empty())
		{
			transmitter_obj->setSchedule(std::move(schedule));
		}

		ctx.world->add(std::move(transmitter_obj));
		return ctx.world->getTransmitters().back().get();
	}

	radar::Transmitter* parseTransmitter(const XmlElement& transmitter, radar::Platform* platform, ParserContext& ctx,
										 const ReferenceLookup& refs)
	{
		const std::string name = XmlElement::getSafeAttribute(transmitter, "name");
		const radar::OperationMode mode = parse_mode_elements(transmitter, "Transmitter '" + name + "'");
		return parseTransmitterWithMode(transmitter, platform, ctx, refs, mode);
	}

	/// Parses a receiver after its operation mode has already been determined.
	radar::Receiver* parseReceiverWithMode(const XmlElement& receiver, radar::Platform* platform, ParserContext& ctx,
										   const ReferenceLookup& refs, const radar::OperationMode mode)
	{
		const std::string name = XmlElement::getSafeAttribute(receiver, "name");
		const SimId id = assign_id_from_attribute("receiver '" + name + "'", ObjectType::Receiver);
		const XmlElement pulsed_mode_element = receiver.childElement("pulsed_mode", 0);
		const bool is_pulsed = mode == radar::OperationMode::PULSED_MODE;

		auto receiver_obj = std::make_unique<radar::Receiver>(platform, name, next_seed(*ctx.master_seeder), mode, id);

		const SimId ant_id = resolve_reference_id(receiver, "antenna", "receiver '" + name + "'", *refs.antennas);
		const antenna::Antenna* antenna = ctx.world->findAntenna(ant_id);
		if (antenna == nullptr)
		{
			throw XmlException("Antenna ID '" + std::to_string(ant_id) + "' not found for receiver '" + name + "'");
		}
		receiver_obj->setAntenna(antenna);

		if (!receiver.childElement("noise_temp", 0).isValid())
		{
			LOG(logging::Level::DEBUG, "Receiver '{}' does not specify noise temperature",
				receiver_obj->getName().c_str());
		}
		else
		{
			try
			{
				receiver_obj->setNoiseTemperature(get_child_real_type(receiver, "noise_temp"));
			}
			catch (const XmlException&)
			{
				LOG(logging::Level::WARNING, "Receiver '{}' has an empty noise temperature; using default.",
					receiver_obj->getName().c_str());
			}
		}

		if (is_pulsed)
		{
			const RealType window_length = get_child_real_type(pulsed_mode_element, "window_length");
			if (window_length <= 0)
			{
				throw XmlException("<window_length> must be positive for receiver '" + name + "'");
			}

			const RealType prf = get_child_real_type(pulsed_mode_element, "prf");
			if (prf <= 0)
			{
				throw XmlException("<prf> must be positive for receiver '" + name + "'");
			}

			const RealType window_skip = get_child_real_type(pulsed_mode_element, "window_skip");
			if (window_skip < 0)
			{
				throw XmlException("<window_skip> must not be negative for receiver '" + name + "'");
			}
			receiver_obj->setWindowProperties(window_length, prf, window_skip);
		}
		const SimId timing_id = resolve_reference_id(receiver, "timing", "receiver '" + name + "'", *refs.timings);
		receiver_obj->setTiming(resolve_timing_instance(timing_id, ctx, "receiver '" + name + "'"));

		if (get_attribute_bool(receiver, "nodirect", false))
		{
			receiver_obj->setFlag(radar::Receiver::RecvFlag::FLAG_NODIRECT);
			LOG(logging::Level::DEBUG, "Ignoring direct signals for receiver '{}'", receiver_obj->getName().c_str());
		}

		if (get_attribute_bool(receiver, "nopropagationloss", false))
		{
			receiver_obj->setFlag(radar::Receiver::RecvFlag::FLAG_NOPROPLOSS);
			LOG(logging::Level::DEBUG, "Ignoring propagation losses for receiver '{}'",
				receiver_obj->getName().c_str());
		}

		RealType pri = is_pulsed ? (1.0 / receiver_obj->getWindowPrf()) : 0.0;
		auto schedule = parseSchedule(receiver, name, is_pulsed, pri);
		if (!schedule.empty())
		{
			receiver_obj->setSchedule(std::move(schedule));
		}

		ctx.world->add(std::move(receiver_obj));
		return ctx.world->getReceivers().back().get();
	}

	radar::Receiver* parseReceiver(const XmlElement& receiver, radar::Platform* platform, ParserContext& ctx,
								   const ReferenceLookup& refs)
	{
		const std::string name = XmlElement::getSafeAttribute(receiver, "name");
		const radar::OperationMode mode = parse_mode_elements(receiver, "Receiver '" + name + "'");
		return parseReceiverWithMode(receiver, platform, ctx, refs, mode);
	}

	void parseMonostatic(const XmlElement& monostatic, radar::Platform* platform, ParserContext& ctx,
						 const ReferenceLookup& refs)
	{
		const std::string name = XmlElement::getSafeAttribute(monostatic, "name");
		const radar::OperationMode monostatic_mode = parse_mode_elements(monostatic, "Monostatic '" + name + "'");
		radar::Transmitter* trans = parseTransmitterWithMode(monostatic, platform, ctx, refs, monostatic_mode);
		radar::Receiver* recv = parseReceiverWithMode(monostatic, platform, ctx, refs, monostatic_mode);
		if (trans->getMode() != monostatic_mode || recv->getMode() != monostatic_mode)
		{
			throw XmlException("Monostatic '" + name + "' parsed inconsistent transmitter/receiver modes.");
		}
		if (trans->getSignal() != nullptr)
		{
			validate_waveform_mode_match(*trans->getSignal(), trans->getMode(),
										 "Monostatic '" + trans->getName() + "'");
		}
		trans->setAttached(recv);
		recv->setAttached(trans);
	}

	void parseTarget(const XmlElement& target, radar::Platform* platform, ParserContext& ctx)
	{
		const std::string name = XmlElement::getSafeAttribute(target, "name");
		const SimId id = assign_id_from_attribute("target '" + name + "'", ObjectType::Target);

		const XmlElement rcs_element = target.childElement("rcs", 0);
		if (!rcs_element.isValid())
		{
			throw XmlException("<rcs> element is required in <target>!");
		}

		const std::string rcs_type = XmlElement::getSafeAttribute(rcs_element, "type");
		std::unique_ptr<radar::Target> target_obj;
		const unsigned seed = next_seed(*ctx.master_seeder);

		if (rcs_type == "isotropic")
		{
			target_obj = radar::createIsoTarget(platform, name, get_child_real_type(rcs_element, "value"), seed, id);
		}
		else if (rcs_type == "file")
		{
			// Defer to dependency-injected file loader
			target_obj = ctx.loaders.loadFileTarget(platform, name,
													XmlElement::getSafeAttribute(rcs_element, "filename"), seed, id);
		}
		else
		{
			throw XmlException("Unsupported RCS type: " + rcs_type);
		}

		if (const XmlElement model = target.childElement("model", 0); model.isValid())
		{
			if (const std::string model_type = XmlElement::getSafeAttribute(model, "type"); model_type == "constant")
			{
				target_obj->setFluctuationModel(std::make_unique<radar::RcsConst>());
			}
			else if (model_type == "chisquare" || model_type == "gamma")
			{
				target_obj->setFluctuationModel(
					std::make_unique<radar::RcsChiSquare>(target_obj->getRngEngine(), get_child_real_type(model, "k")));
			}
			else
			{
				throw XmlException("Unsupported model type: " + model_type);
			}
		}

		LOG(logging::Level::DEBUG, "Added target {} with RCS type {} to platform {}", name, rcs_type,
			platform->getName());
		ctx.world->add(std::move(target_obj));
	}

	void parsePlatformElements(const XmlElement& platform, ParserContext& ctx, radar::Platform* plat,
							   const std::function<void(const XmlElement&, std::string_view)>& register_name,
							   const ReferenceLookup& refs)
	{
		auto parseChildrenWithRefs = [&](const std::string& elementName, auto parseFunc)
		{
			unsigned index = 0;
			while (true)
			{
				const XmlElement element = platform.childElement(elementName, index++);
				if (!element.isValid())
					break;
				register_name(element, elementName);
				parseFunc(element, plat, ctx, refs);
			}
		};

		auto parseChildrenWithoutRefs = [&](const std::string& elementName, auto parseFunc)
		{
			unsigned index = 0;
			while (true)
			{
				const XmlElement element = platform.childElement(elementName, index++);
				if (!element.isValid())
					break;
				register_name(element, elementName);
				parseFunc(element, plat, ctx);
			}
		};

		parseChildrenWithRefs("monostatic", parseMonostatic);
		parseChildrenWithRefs("transmitter", parseTransmitter);
		parseChildrenWithRefs("receiver", parseReceiver);
		parseChildrenWithoutRefs("target", parseTarget);
	}

	void parsePlatform(const XmlElement& platform, ParserContext& ctx,
					   const std::function<void(const XmlElement&, std::string_view)>& register_name,
					   const ReferenceLookup& refs)
	{
		std::string name = XmlElement::getSafeAttribute(platform, "name");
		const SimId id = assign_id_from_attribute("platform '" + name + "'", ObjectType::Platform);
		auto plat = std::make_unique<radar::Platform>(name, id);

		parsePlatformElements(platform, ctx, plat.get(), register_name, refs);

		if (const XmlElement motion_path = platform.childElement("motionpath", 0); motion_path.isValid())
		{
			parseMotionPath(motion_path, plat.get());
		}

		const XmlElement rot_path = platform.childElement("rotationpath", 0);
		if (const XmlElement fixed_rot = platform.childElement("fixedrotation", 0);
			rot_path.isValid() && fixed_rot.isValid())
		{
			LOG(logging::Level::ERROR,
				"Both <rotationpath> and <fixedrotation> are declared for platform {}. Only <rotationpath> will be "
				"used.",
				plat->getName());
			parseRotationPath(rot_path, plat.get(), ctx.parameters.rotation_angle_unit);
		}
		else if (rot_path.isValid())
		{
			parseRotationPath(rot_path, plat.get(), ctx.parameters.rotation_angle_unit);
		}
		else if (fixed_rot.isValid())
		{
			parseFixedRotation(fixed_rot, plat.get(), ctx.parameters.rotation_angle_unit);
		}

		ctx.world->add(std::move(plat));
	}

	void collectIncludeElements(const XmlDocument& doc, const fs::path& currentDir, std::vector<fs::path>& includePaths)
	{
		unsigned index = 0;
		while (true)
		{
			XmlElement include_element = doc.getRootElement().childElement("include", index++);
			if (!include_element.isValid())
				break;

			std::string include_filename = include_element.getText();
			if (include_filename.empty())
			{
				LOG(logging::Level::ERROR, "<include> element is missing the filename!");
				continue;
			}

			fs::path include_path = currentDir / include_filename;
			includePaths.push_back(include_path);

			XmlDocument included_doc;
			if (!included_doc.loadFile(include_path.string()))
			{
				LOG(logging::Level::ERROR, "Failed to load included XML file: {}", include_path.string());
				continue;
			}

			collectIncludeElements(included_doc, include_path.parent_path(), includePaths);
		}
	}

	bool addIncludeFilesToMainDocument(const XmlDocument& mainDoc, const fs::path& currentDir)
	{
		std::vector<fs::path> include_paths;
		collectIncludeElements(mainDoc, currentDir, include_paths);
		bool did_combine = false;

		for (const auto& include_path : include_paths)
		{
			XmlDocument included_doc;
			if (!included_doc.loadFile(include_path.string()))
			{
				throw XmlException("Failed to load included XML file: " + include_path.string());
			}

			mergeXmlDocuments(mainDoc, included_doc);
			did_combine = true;
		}

		removeIncludeElements(mainDoc);
		return did_combine;
	}

	void validateXml(const bool didCombine, const XmlDocument& mainDoc)
	{
		LOG(logging::Level::DEBUG, "Validating the{}XML file...", didCombine ? " combined " : " ");
		if (!mainDoc.validateWithDtd(fers_xml_dtd))
		{
			LOG(logging::Level::FATAL, "{} XML file failed DTD validation!", didCombine ? "Combined" : "Main");
			throw XmlException("XML file failed DTD validation!");
		}
		LOG(logging::Level::DEBUG, "{} XML file passed DTD validation.", didCombine ? "Combined" : "Main");

		if (!mainDoc.validateWithXsd(fers_xml_xsd))
		{
			LOG(logging::Level::FATAL, "{} XML file failed XSD validation!", didCombine ? "Combined" : "Main");
			throw XmlException("XML file failed XSD validation!");
		}
		LOG(logging::Level::DEBUG, "{} XML file passed XSD validation.", didCombine ? "Combined" : "Main");
	}

	void processParsedDocument(const XmlDocument& doc, ParserContext& ctx)
	{
		const XmlElement root = doc.getRootElement();
		if (root.name() != "simulation")
		{
			throw XmlException("Root element is not <simulation>!");
		}

		std::unordered_map<std::string, std::string> name_registry;
		name_registry.reserve(64); // TODO: reserve 64?
		const auto register_name = [&](const XmlElement& element, const std::string_view kind)
		{
			const std::string name = XmlElement::getSafeAttribute(element, "name");
			const auto [iter, inserted] = name_registry.emplace(name, std::string(kind));
			if (!inserted)
			{
				throw XmlException("Duplicate name '" + name + "' found for " + std::string(kind) +
								   "; previously used by " + iter->second + ".");
			}
		};

		try
		{
			ctx.parameters.simulation_name = XmlElement::getSafeAttribute(root, "name");
			if (!ctx.parameters.simulation_name.empty())
			{
				LOG(logging::Level::INFO, "Simulation name set to: {}", ctx.parameters.simulation_name);
			}
		}
		catch (const XmlException&)
		{
			LOG(logging::Level::WARNING, "No 'name' attribute found in <simulation> tag. KML name will default.");
		}

		parseParameters(root.childElement("parameters", 0), ctx.parameters);

		params::params = ctx.parameters;

		auto parseElements =
			[](const XmlElement& parent, const std::string& elementName, ParserContext& parser_ctx, auto parseFunction)
		{
			unsigned index = 0;
			while (true)
			{
				XmlElement element = parent.childElement(elementName, index++);
				if (!element.isValid())
					break;
				parseFunction(element, parser_ctx);
			}
		};

		parseElements(root, "waveform", ctx,
					  [&](const XmlElement& p, ParserContext& c)
					  {
						  register_name(p, "waveform");
						  parseWaveform(p, c);
					  });

		parseElements(root, "timing", ctx,
					  [&](const XmlElement& p, ParserContext& c)
					  {
						  register_name(p, "timing");
						  parseTiming(p, c);
					  });

		parseElements(root, "antenna", ctx,
					  [&](const XmlElement& p, ParserContext& c)
					  {
						  register_name(p, "antenna");
						  parseAntenna(p, c);
					  });

		std::unordered_map<std::string, SimId> waveform_refs;
		std::unordered_map<std::string, SimId> antenna_refs;
		std::unordered_map<std::string, SimId> timing_refs;
		waveform_refs.reserve(ctx.world->getWaveforms().size());
		antenna_refs.reserve(ctx.world->getAntennas().size());
		timing_refs.reserve(ctx.world->getTimings().size());

		for (const auto& [id, waveform] : ctx.world->getWaveforms())
			waveform_refs.emplace(waveform->getName(), id);
		for (const auto& [id, antenna] : ctx.world->getAntennas())
			antenna_refs.emplace(antenna->getName(), id);
		for (const auto& [id, timing] : ctx.world->getTimings())
			timing_refs.emplace(timing->getName(), id);

		const ReferenceLookup refs{&waveform_refs, &antenna_refs, &timing_refs};

		parseElements(root, "platform", ctx,
					  [&](const XmlElement& p, ParserContext& c)
					  {
						  register_name(p, "platform");
						  parsePlatform(p, c, register_name, refs);
					  });

		const RealType start_time = ctx.parameters.start;
		const RealType end_time = ctx.parameters.end;
		const RealType dt_sim = 1.0 / (ctx.parameters.rate * ctx.parameters.oversample_ratio);
		const auto num_samples = static_cast<size_t>(std::ceil((end_time - start_time) / dt_sim));

		for (const auto& receiver : ctx.world->getReceivers())
		{
			if (receiver->getMode() == radar::OperationMode::CW_MODE ||
				receiver->getMode() == radar::OperationMode::FMCW_MODE)
			{
				receiver->prepareStreamingData(num_samples);
			}
		}

		ctx.world->scheduleInitialEvents();

		LOG(logging::Level::DEBUG, "Initial Event Queue State:\n{}", ctx.world->dumpEventQueue());
	}

	AssetLoaders createDefaultAssetLoaders()
	{
		return {.loadWaveform = [](const std::string& name, const fs::path& pulse_path, RealType power,
								   RealType carrierFreq, SimId id)
				{ return serial::loadWaveformFromFile(name, pulse_path.string(), power, carrierFreq, id); },
				.loadXmlAntenna = [](const std::string& name, const std::string& filename, SimId id)
				{ return std::make_unique<antenna::XmlAntenna>(name, filename, id); },
				.loadH5Antenna = [](const std::string& name, const std::string& filename, SimId id)
				{ return std::make_unique<antenna::H5Antenna>(name, filename, id); },
				.loadFileTarget = [](radar::Platform* platform, const std::string& name, const std::string& filename,
									 unsigned seed, SimId id)
				{ return radar::createFileTarget(platform, name, filename, seed, id); }};
	}
}
