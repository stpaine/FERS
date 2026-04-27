// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file antenna_factory.cpp
 * @brief Implementation of the Antenna class and its derived classes.
 */

#include "antenna/antenna_factory.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <limits>
#include <optional>
#include <stdexcept>

#include "antenna_pattern_dtd.h"
#include "antenna_pattern_xsd.h"
#include "core/config.h"
#include "core/logging.h"
#include "core/portable_utils.h"
#include "math/geometry_ops.h"
#include "serial/libxml_wrapper.h"

using logging::Level;
using math::SVec3;
using math::Vec3;

namespace
{
	/// Unit used for XML antenna axis sample angles.
	enum class AxisUnit
	{
		Radians, ///< Axis samples are expressed in radians.
		Degrees, ///< Axis samples are expressed in degrees.
	};

	/// Gain format used for XML antenna axis sample values.
	enum class AxisGainFormat
	{
		Linear, ///< Gain samples are already linear.
		DBi, ///< Gain samples are expressed in dBi.
	};

	/// Parsed metadata attributes for one XML antenna axis.
	struct AxisMetadata
	{
		AxisUnit unit{AxisUnit::Radians}; ///< Parsed angle unit.
		AxisGainFormat format{AxisGainFormat::Linear}; ///< Parsed gain format.
		antenna::XmlAntenna::AxisSymmetry symmetry{antenna::XmlAntenna::AxisSymmetry::Mirrored}; ///< Parsed symmetry.
		bool unit_explicit{}; ///< True when unit was explicitly set.
		bool format_explicit{}; ///< True when format was explicitly set.
		bool symmetry_explicit{}; ///< True when symmetry was explicitly set.
	};

	/// Result from loading one XML antenna gain axis.
	struct AxisLoadResult
	{
		RealType max_gain{}; ///< Maximum linear gain found on the axis.
		antenna::XmlAntenna::AxisSymmetry symmetry{
			antenna::XmlAntenna::AxisSymmetry::Mirrored}; ///< Axis symmetry mode.
		std::size_t sample_count{}; ///< Number of parsed samples.
	};

	/// Returns a lowercase copy of a string.
	std::string toLowerCopy(std::string value)
	{
		std::transform(value.begin(), value.end(), value.begin(),
					   [](const unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
		return value;
	}

	/// Parses a finite real value and reports the provided context on failure.
	RealType parseRealValue(const std::string_view text, const std::string_view context)
	{
		const std::string value(text);
		size_t index = 0;
		double parsed = 0.0;

		try
		{
			parsed = std::stod(value, &index);
		}
		catch (const std::exception&)
		{
			throw std::runtime_error("Invalid numeric value for " + std::string(context) + ": '" + value + "'.");
		}

		while (index < value.size() && (std::isspace(static_cast<unsigned char>(value[index])) != 0))
		{
			++index;
		}

		if (index != value.size())
		{
			throw std::runtime_error("Invalid numeric value for " + std::string(context) + ": '" + value + "'.");
		}

		if (!std::isfinite(parsed))
		{
			throw std::runtime_error("Non-finite numeric value for " + std::string(context) + ".");
		}

		return static_cast<RealType>(parsed);
	}

	/// Returns the XML token for an axis unit.
	const char* axisUnitName(const AxisUnit unit) noexcept { return unit == AxisUnit::Degrees ? "deg" : "rad"; }

	/// Returns the XML token for an axis gain format.
	const char* axisFormatName(const AxisGainFormat format) noexcept
	{
		return format == AxisGainFormat::DBi ? "dBi" : "linear";
	}

	/// Returns the XML token for an axis symmetry mode.
	const char* axisSymmetryName(const antenna::XmlAntenna::AxisSymmetry symmetry) noexcept
	{
		return symmetry == antenna::XmlAntenna::AxisSymmetry::None ? "none" : "mirrored";
	}

	/// Parses optional metadata attributes from one XML antenna axis.
	AxisMetadata parseAxisMetadata(const XmlElement& axisXml)
	{
		AxisMetadata metadata;
		const std::string axis_name(axisXml.name());

		if (const auto unit_attr = XmlElement::getOptionalAttribute(axisXml, "unit"); unit_attr.has_value())
		{
			metadata.unit_explicit = true;
			const std::string value = toLowerCopy(*unit_attr);
			if (value == "rad")
			{
				metadata.unit = AxisUnit::Radians;
			}
			else if (value == "deg")
			{
				metadata.unit = AxisUnit::Degrees;
			}
			else
			{
				throw std::runtime_error("Unsupported unit '" + *unit_attr + "' on <" + axis_name + "> axis.");
			}
		}

		if (const auto format_attr = XmlElement::getOptionalAttribute(axisXml, "format"); format_attr.has_value())
		{
			metadata.format_explicit = true;
			const std::string value = toLowerCopy(*format_attr);
			if (value == "linear")
			{
				metadata.format = AxisGainFormat::Linear;
			}
			else if (value == "dbi")
			{
				metadata.format = AxisGainFormat::DBi;
			}
			else
			{
				throw std::runtime_error("Unsupported format '" + *format_attr + "' on <" + axis_name + "> axis.");
			}
		}

		if (const auto symmetry_attr = XmlElement::getOptionalAttribute(axisXml, "symmetry"); symmetry_attr.has_value())
		{
			metadata.symmetry_explicit = true;
			const std::string value = toLowerCopy(*symmetry_attr);
			if (value == "mirrored")
			{
				metadata.symmetry = antenna::XmlAntenna::AxisSymmetry::Mirrored;
			}
			else if (value == "none")
			{
				metadata.symmetry = antenna::XmlAntenna::AxisSymmetry::None;
			}
			else
			{
				throw std::runtime_error("Unsupported symmetry '" + *symmetry_attr + "' on <" + axis_name + "> axis.");
			}
		}

		return metadata;
	}

	/**
	 * @brief Compute the sinc function.
	 *
	 * @param theta The angle for which to compute the sinc function.
	 * @return The value of the sinc function at the given angle theta.
	 */
	RealType sinc(const RealType theta) noexcept
	{
		if (std::abs(theta) < EPSILON)
		{
			return 1.0;
		}
		return std::sin(theta) / theta;
	}

	/**
	 * @brief Compute the Bessel function of the first kind.
	 *
	 * @param x The value for which to compute the Bessel function.
	 * @return The value of the Bessel function of the first kind at the given value x.
	 */
	RealType j1C(const RealType x) noexcept { return x == 0 ? 1.0 : core::besselJ1(x) / x; }

	/**
	 * @brief Load antenna gain axis data from an XML element.
	 *
	 * @param set The interpolation set to store the gain axis data.
	 * @param axisXml The XML element containing the gain axis data.
	 */
	AxisLoadResult loadAntennaGainAxis(const interp::InterpSet* set, const XmlElement& axisXml)
	{
		if (!axisXml.isValid())
		{
			throw std::runtime_error("XML antenna pattern is missing a required gain axis.");
		}

		const AxisMetadata metadata = parseAxisMetadata(axisXml);
		const std::string axis_name(axisXml.name());
		XmlElement sample = axisXml.childElement("gainsample");
		if (!sample.isValid())
		{
			throw std::runtime_error("XML antenna <" + axis_name + "> axis must contain at least one <gainsample>.");
		}

		RealType min_angle = std::numeric_limits<RealType>::max();
		RealType max_angle = std::numeric_limits<RealType>::lowest();
		RealType max_gain = 0.0;
		std::size_t sample_count = 0;

		while (sample.isValid())
		{
			XmlElement angle_element = sample.childElement("angle", 0);

			if (XmlElement gain_element = sample.childElement("gain", 0);
				angle_element.isValid() && gain_element.isValid())
			{
				const RealType raw_angle = parseRealValue(angle_element.getText(), "<" + axis_name + "> sample angle");
				const RealType raw_gain = parseRealValue(gain_element.getText(), "<" + axis_name + "> sample gain");
				const RealType angle = metadata.unit == AxisUnit::Degrees ? raw_angle * (PI / 180.0) : raw_angle;
				const RealType gain =
					metadata.format == AxisGainFormat::DBi ? std::pow(10.0, raw_gain / 10.0) : raw_gain;

				if (!std::isfinite(angle))
				{
					throw std::runtime_error("Converted <" + axis_name + "> sample angle is non-finite.");
				}
				if (!std::isfinite(gain))
				{
					throw std::runtime_error("Converted <" + axis_name + "> sample gain is non-finite.");
				}
				if (gain < 0.0)
				{
					throw std::runtime_error("Converted <" + axis_name + "> sample gain must be non-negative.");
				}

				set->insertSample(angle, gain);
				min_angle = std::min(min_angle, angle);
				max_angle = std::max(max_angle, angle);
				max_gain = std::max(max_gain, gain);
				++sample_count;
			}
			else if (sample.name() == "gainsample")
			{
				throw std::runtime_error("Each <gainsample> in <" + axis_name + "> must contain <angle> and <gain>.");
			}

			sample = XmlElement(sample.getNode()->next);
		}

		AxisLoadResult result;
		result.max_gain = max_gain;
		result.symmetry = metadata.symmetry;
		result.sample_count = sample_count;

		if (metadata.symmetry_explicit)
		{
			if (result.symmetry == antenna::XmlAntenna::AxisSymmetry::Mirrored && min_angle < 0.0)
			{
				throw std::runtime_error("XML antenna <" + axis_name +
										 "> axis uses symmetry='mirrored' but defines negative sample angles.");
			}
			if (result.symmetry == antenna::XmlAntenna::AxisSymmetry::None && !(min_angle < 0.0 && max_angle > 0.0))
			{
				throw std::runtime_error(
					"XML antenna <" + axis_name +
					"> axis uses symmetry='none' but does not span both negative and positive angles.");
			}
		}
		else if (min_angle < 0.0)
		{
			if (!(min_angle < 0.0 && max_angle > 0.0))
			{
				throw std::runtime_error("XML antenna <" + axis_name +
										 "> axis contains negative sample angles but does not span both sides of zero "
										 "for direct signed-angle lookup.");
			}
			result.symmetry = antenna::XmlAntenna::AxisSymmetry::None;
		}

		const char* unit_source = metadata.unit_explicit ? "explicit" : "legacy default";
		const char* format_source = metadata.format_explicit ? "explicit" : "legacy default";
		const char* symmetry_source = metadata.symmetry_explicit
			? "explicit"
			: (result.symmetry == antenna::XmlAntenna::AxisSymmetry::None ? "auto-detected" : "legacy default");

		LOG(Level::INFO,
			"XML antenna axis '{}' using unit='{}' ({}) format='{}' ({}) symmetry='{}' ({}) with {} samples.",
			axis_name, axisUnitName(metadata.unit), unit_source, axisFormatName(metadata.format), format_source,
			axisSymmetryName(result.symmetry), symmetry_source, result.sample_count);

		return result;
	}
}

namespace antenna
{
	void Antenna::setEfficiencyFactor(const RealType loss) noexcept
	{
		if (loss > 1)
		{
			LOG(Level::INFO, "Using greater than unity antenna efficiency.");
		}
		_loss_factor = loss;
	}

	RealType Antenna::getAngle(const SVec3& angle, const SVec3& refangle) noexcept
	{
		SVec3 normangle(angle);
		normangle.length = 1;
		return std::acos(dotProduct(Vec3(normangle), Vec3(refangle)));
	}

	RealType Gaussian::getGain(const SVec3& angle, const SVec3& refangle, RealType /*wavelength*/) const noexcept
	{
		const SVec3 a = angle - refangle;
		return std::exp(-a.azimuth * a.azimuth * _azscale) * std::exp(-a.elevation * a.elevation * _elscale);
	}

	RealType Sinc::getGain(const SVec3& angle, const SVec3& refangle, RealType /*wavelength*/) const noexcept
	{
		const RealType theta = getAngle(angle, refangle);
		const RealType sinc_val = sinc(_beta * theta);
		const RealType gain_pattern = std::pow(std::abs(sinc_val), _gamma);
		return _alpha * gain_pattern * getEfficiencyFactor();
	}

	RealType SquareHorn::getGain(const SVec3& angle, const SVec3& refangle, const RealType wavelength) const noexcept
	{
		const RealType ge = 4 * PI * std::pow(_dimension, 2) / std::pow(wavelength, 2);
		const RealType x = PI * _dimension * std::sin(getAngle(angle, refangle)) / wavelength;
		return ge * std::pow(sinc(x), 2) * getEfficiencyFactor();
	}

	RealType Parabolic::getGain(const SVec3& angle, const SVec3& refangle, const RealType wavelength) const noexcept
	{
		const RealType ge = std::pow(PI * _diameter / wavelength, 2);
		const RealType x = PI * _diameter * std::sin(getAngle(angle, refangle)) / wavelength;
		return ge * std::pow(2 * j1C(x), 2) * getEfficiencyFactor();
	}

	std::optional<RealType> XmlAntenna::lookupAxisGain(const interp::InterpSet* set, const RealType angle,
													   const AxisSymmetry symmetry) const noexcept
	{
		if (symmetry == AxisSymmetry::Mirrored)
		{
			return set->getValueAt(std::abs(angle));
		}
		return set->getValueAt(angle);
	}

	RealType XmlAntenna::getGain(const SVec3& angle, const SVec3& refangle, RealType /*wavelength*/) const
	{
		const SVec3 delta_angle = angle - refangle;

		const std::optional<RealType> azi_value =
			lookupAxisGain(_azi_samples.get(), delta_angle.azimuth, _azi_symmetry);

		if (const std::optional<RealType> elev_value =
				lookupAxisGain(_elev_samples.get(), delta_angle.elevation, _elev_symmetry);
			azi_value && elev_value)
		{
			return *azi_value * *elev_value * _max_gain * getEfficiencyFactor();
		}

		LOG(Level::FATAL, "Could not get antenna gain value");
		throw std::runtime_error("Could not get antenna gain value");
	}

	void XmlAntenna::loadAntennaDescription(const std::string_view filename)
	{
		_filename = filename;
		XmlDocument doc;
		if (!doc.loadFile(std::string(filename)))
		{
			LOG(Level::FATAL, "Could not load antenna description {}", filename.data());
			throw std::runtime_error("Could not load antenna description");
		}
		(void)doc.validateWithDtd(antenna_pattern_dtd);
		(void)doc.validateWithXsd(antenna_pattern_xsd);

		const XmlElement root(doc.getRootElement());
		const AxisLoadResult elev_result = loadAntennaGainAxis(_elev_samples.get(), root.childElement("elevation", 0));
		const AxisLoadResult azi_result = loadAntennaGainAxis(_azi_samples.get(), root.childElement("azimuth", 0));

		_elev_symmetry = elev_result.symmetry;
		_azi_symmetry = azi_result.symmetry;

		_max_gain = std::max(azi_result.max_gain, elev_result.max_gain);
		if (!std::isfinite(_max_gain) || _max_gain <= 0.0)
		{
			throw std::runtime_error("XML antenna pattern peak linear gain must be greater than zero.");
		}
		_elev_samples->divide(_max_gain);
		_azi_samples->divide(_max_gain);
	}

	RealType H5Antenna::getGain(const SVec3& angle, const SVec3& refangle, RealType /*wavelength*/) const
	{
		constexpr RealType two_pi = 2.0 * PI;

		const SVec3& pattern_angle = angle - refangle;

		const double ex1 = (pattern_angle.azimuth + PI) / two_pi;
		const double ey1 = (pattern_angle.elevation + PI) / two_pi;

		const auto calc_grid_point = [](const double value, const std::size_t size)
		{
			const double grid_size = static_cast<double>(size - 1);
			const double x1 = std::floor(value * grid_size) / grid_size;
			const double x2 = std::min(x1 + 1.0 / static_cast<double>(size), 1.0);
			return std::pair{x1, x2};
		};

		const std::size_t size_azi = _pattern.size();
		const std::size_t size_elev = _pattern[0].size();

		LOG(logging::Level::TRACE, "Size of pattern: {} x {}", size_azi, size_elev);

		const auto [x1, x2] = calc_grid_point(ex1, size_azi);
		const auto [y1, y2] = calc_grid_point(ey1, size_elev);

		const double t = (ex1 - x1) / (x2 - x1);
		const double u = (ey1 - y1) / (y2 - y1);

		const auto calc_array_index = [](const double value, const std::size_t size)
		{ return std::min(static_cast<std::size_t>(std::floor(value * static_cast<double>(size))), size - 1); };

		const std::size_t arr_x = calc_array_index(x1, size_azi);
		const std::size_t arr_y = calc_array_index(y1, size_elev);

		const RealType interp = (1.0 - t) * (1.0 - u) * _pattern[arr_x][arr_y] +
			t * (1.0 - u) * _pattern[(arr_x + 1) % size_azi][arr_y] +
			t * u * _pattern[(arr_x + 1) % size_azi][(arr_y + 1) % size_elev] +
			(1.0 - t) * u * _pattern[arr_x][(arr_y + 1) % size_elev];

		return interp * getEfficiencyFactor();
	}
}
