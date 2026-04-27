// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file antenna_factory.h
 * @brief Header file defining various types of antennas and their gain patterns.
 */

#pragma once

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "core/config.h"
#include "core/logging.h"
#include "core/sim_id.h"
#include "interpolation/interpolation_set.h"
#include "math/geometry_ops.h"

namespace serial
{
	/// Reads a 2D antenna gain pattern from the named dataset.
	std::vector<std::vector<RealType>> readPattern(const std::string& name, const std::string& datasetName);
}

namespace antenna
{
	/**
	 * @class Antenna
	 * @brief Abstract base class representing an antenna.
	 */
	class Antenna
	{
	public:
		/**
		 * @brief Constructs an Antenna object with the given name.
		 *
		 * @param name The name of the antenna.
		 */
		explicit Antenna(std::string name, const SimId id = 0) noexcept :
			_loss_factor(1), _id(id == 0 ? SimIdGenerator::instance().generateId(ObjectType::Antenna) : id),
			_name(std::move(name))
		{
		}

		virtual ~Antenna() = default;

		Antenna(const Antenna&) = delete;

		Antenna& operator=(const Antenna&) = delete;

		Antenna(Antenna&&) = default;

		Antenna& operator=(Antenna&&) = default;

		/**
		 * @brief Computes the gain of the antenna based on the input angle and reference angle.
		 *
		 * @param angle The angle at which the gain is to be computed.
		 * @param refangle The reference angle.
		 * @param wavelength The wavelength of the signal.
		 * @return The gain of the antenna at the specified angle and wavelength.
		 */
		[[nodiscard]] virtual RealType getGain(const math::SVec3& angle, const math::SVec3& refangle,
											   RealType wavelength) const = 0;

		/**
		 * @brief Retrieves the efficiency factor of the antenna.
		 *
		 * @return The efficiency factor of the antenna.
		 */
		[[nodiscard]] RealType getEfficiencyFactor() const noexcept { return _loss_factor; }

		/**
		 * @brief Retrieves the name of the antenna.
		 *
		 * @return The name of the antenna.
		 */
		[[nodiscard]] std::string getName() const noexcept { return _name; }

		/**
		 * @brief Retrieves the unique ID of the antenna.
		 *
		 * @return The antenna SimId.
		 */
		[[nodiscard]] SimId getId() const noexcept { return _id; }

		/**
		 * @brief Computes the noise temperature of the antenna based on the angle.
		 *
		 * @param angle The angle at which the noise temperature is to be computed.
		 * @return The noise temperature of the antenna.
		 */
		// TODO: Implement noise temperature calculation
		[[nodiscard]] virtual RealType getNoiseTemperature(const math::SVec3& /*angle*/) const noexcept { return 0; }

		/**
		 * @brief Sets the efficiency factor of the antenna.
		 *
		 * @param loss The new efficiency factor.
		 */
		void setEfficiencyFactor(RealType loss) noexcept;

		/**
		 * @brief Sets the name of the antenna.
		 *
		 * @param name The new name of the antenna.
		 */
		void setName(std::string name) noexcept { _name = std::move(name); }

	protected:
		/**
		 * @brief Computes the angle between the input and reference angles.
		 *
		 * @param angle The input angle.
		 * @param refangle The reference angle.
		 * @return The computed angle.
		 */
		static RealType getAngle(const math::SVec3& angle, const math::SVec3& refangle) noexcept;

	private:
		RealType _loss_factor; ///< Efficiency factor of the antenna.
		SimId _id; ///< Unique ID for this antenna.
		std::string _name; ///< Name of the antenna.
	};

	/**
	 * @class Isotropic
	 * @brief Represents an isotropic antenna with uniform gain in all directions.
	 *
	 * This class models an ideal isotropic antenna, which has a directivity of 1 (0 dB).
	 */
	class Isotropic final : public Antenna
	{
	public:
		/**
		 * @brief Constructs an Isotropic antenna with the given name.
		 *
		 * @param name The name of the antenna.
		 */
		explicit Isotropic(const std::string_view name, const SimId id = 0) : Antenna(name.data(), id) {}

		~Isotropic() override = default;

		Isotropic(const Isotropic&) = delete;

		Isotropic& operator=(const Isotropic&) = delete;

		Isotropic(Isotropic&&) = delete;

		Isotropic& operator=(Isotropic&&) = delete;

		/**
		 * @brief Computes the gain of the isotropic antenna.
		 *
		 * @return The gain of the antenna, which is equal to the efficiency factor.
		 */
		[[nodiscard]] RealType getGain(const math::SVec3& /*angle*/, const math::SVec3& /*refangle*/,
									   RealType /*wavelength*/) const override
		{
			// David Young: Isotropic antennas have a directivity of 1 (or 0 dB),
			// therefore, the gain of the antenna is the efficiency factor
			return getEfficiencyFactor();
		}
	};

	/**
	 * @class Sinc
	 * @brief Represents a sinc function-based antenna gain pattern.
	 *
	 * This antenna has a gain pattern defined by a sinc function, with customizable parameters.
	 */
	class Sinc final : public Antenna
	{
	public:
		/**
		 * @brief Constructs a Sinc antenna with the given parameters.
		 *
		 * @param name The name of the antenna.
		 * @param alpha The alpha parameter.
		 * @param beta The beta parameter.
		 * @param gamma The gamma parameter.
		 */
		Sinc(const std::string_view name, const RealType alpha, const RealType beta, const RealType gamma,
			 const SimId id = 0) : Antenna(name.data(), id), _alpha(alpha), _beta(beta), _gamma(gamma)
		{
		}

		~Sinc() override = default;

		Sinc(const Sinc&) = delete;

		Sinc& operator=(const Sinc&) = delete;

		Sinc(Sinc&&) = delete;

		Sinc& operator=(Sinc&&) = delete;

		/** @brief Gets the alpha parameter of the sinc function. */
		[[nodiscard]] RealType getAlpha() const noexcept { return _alpha; }

		/** @brief Gets the beta parameter of the sinc function. */
		[[nodiscard]] RealType getBeta() const noexcept { return _beta; }

		/** @brief Gets the gamma parameter of the sinc function. */
		[[nodiscard]] RealType getGamma() const noexcept { return _gamma; }

		/**
		 * @brief Computes the gain of the sinc antenna based on the input parameters.
		 *
		 * @param angle The angle at which the gain is to be computed.
		 * @param refangle The reference angle.
		 * @param wavelength The wavelength of the signal.
		 * @return The computed gain of the antenna.
		 */
		[[nodiscard]] RealType getGain(const math::SVec3& angle, const math::SVec3& refangle,
									   RealType wavelength) const noexcept override;

		/**
		 * @brief Sets the alpha parameter of the sinc function.
		 *
		 * @param alpha The new alpha parameter.
		 */
		void setAlpha(RealType alpha) noexcept { _alpha = alpha; }

		/**
		 * @brief Sets the beta parameter of the sinc function.
		 *
		 * @param beta The new beta parameter.
		 */
		void setBeta(RealType beta) noexcept { _beta = beta; }

		/**
		 * @brief Sets the gamma parameter of the sinc function.
		 *
		 * @param gamma The new gamma parameter.
		 */
		void setGamma(RealType gamma) noexcept { _gamma = gamma; }

	private:
		RealType _alpha; ///< Parameter defining the shape of the gain pattern.
		RealType _beta; ///< Parameter defining the shape of the gain pattern.
		RealType _gamma; ///< Parameter defining the shape of the gain pattern.
	};

	/**
	 * @class Gaussian
	 * @brief Represents a Gaussian-shaped antenna gain pattern.
	 *
	 * This antenna has a gain pattern that follows a Gaussian distribution.
	 */
	class Gaussian final : public Antenna
	{
	public:
		/**
		 * @brief Constructs a Gaussian antenna with the given parameters.
		 *
		 * @param name The name of the antenna.
		 * @param azscale The azimuth scale factor.
		 * @param elscale The elevation scale factor.
		 */
		Gaussian(const std::string_view name, const RealType azscale, const RealType elscale, const SimId id = 0) :
			Antenna(name.data(), id), _azscale(azscale), _elscale(elscale)
		{
		}

		~Gaussian() override = default;

		Gaussian(const Gaussian&) = delete;

		Gaussian& operator=(const Gaussian&) = delete;

		Gaussian(Gaussian&&) = delete;

		Gaussian& operator=(Gaussian&&) = delete;

		/**
		 * @brief Computes the gain of the Gaussian antenna.
		 *
		 * @param angle The angle at which the gain is to be computed.
		 * @param refangle The reference angle.
		 * @param wavelength The wavelength of the signal.
		 * @return The computed gain of the antenna.
		 */
		[[nodiscard]] RealType getGain(const math::SVec3& angle, const math::SVec3& refangle,
									   RealType wavelength) const noexcept override;

		/** @brief Gets the azimuth scale factor. */
		[[nodiscard]] RealType getAzimuthScale() const noexcept { return _azscale; }

		/** @brief Gets the elevation scale factor. */
		[[nodiscard]] RealType getElevationScale() const noexcept { return _elscale; }

		/**
		 * @brief Sets the azimuth scale factor of the Gaussian function.
		 *
		 * @param azscale The new azimuth scale factor.
		 */
		void setAzimuthScale(RealType azscale) noexcept { _azscale = azscale; }

		/**
		 * @brief Sets the elevation scale factor of the Gaussian function.
		 *
		 * @param elscale The new elevation scale factor.
		 */
		void setElevationScale(RealType elscale) noexcept { _elscale = elscale; }

	private:
		RealType _azscale; ///< Azimuth scale factor.
		RealType _elscale; ///< Elevation scale factor.
	};

	/**
	 * @class SquareHorn
	 * @brief Represents a square horn antenna.
	 *
	 * This antenna models a square horn with a specific dimension.
	 */
	class SquareHorn final : public Antenna
	{
	public:
		/**
		 * @brief Constructs a SquareHorn antenna with the given dimension.
		 *
		 * @param name The name of the antenna.
		 * @param dimension The dimension of the square horn.
		 */
		SquareHorn(const std::string_view name, const RealType dimension, const SimId id = 0) :
			Antenna(name.data(), id), _dimension(dimension)
		{
		}

		~SquareHorn() override = default;

		SquareHorn(const SquareHorn&) = delete;

		SquareHorn& operator=(const SquareHorn&) = delete;

		SquareHorn(SquareHorn&&) = delete;

		SquareHorn& operator=(SquareHorn&&) = delete;

		/**
		 * @brief Computes the gain of the square horn antenna.
		 *
		 * @param angle The angle at which the gain is to be computed.
		 * @param refangle The reference angle.
		 * @param wavelength The wavelength of the signal.
		 * @return The computed gain of the antenna.
		 */
		[[nodiscard]] RealType getGain(const math::SVec3& angle, const math::SVec3& refangle,
									   RealType wavelength) const noexcept override;

		/** @brief Gets the dimension of the square horn. */
		[[nodiscard]] RealType getDimension() const noexcept { return _dimension; }

		/**
		 * @brief Sets the dimension of the square horn.
		 *
		 * @param dimension The new dimension of the square horn.
		 */
		void setDimension(RealType dimension) noexcept { _dimension = dimension; }

	private:
		RealType _dimension; ///< Dimension of the square horn.
	};

	/**
	 * @class Parabolic
	 * @brief Represents a parabolic reflector antenna.
	 *
	 * This antenna models a parabolic reflector with a specific diameter.
	 */
	class Parabolic final : public Antenna
	{
	public:
		/**
		 * @brief Constructs a Parabolic antenna with the given diameter.
		 *
		 * @param name The name of the antenna.
		 * @param diameter The diameter of the parabolic reflector.
		 */
		Parabolic(const std::string_view name, const RealType diameter, const SimId id = 0) :
			Antenna(name.data(), id), _diameter(diameter)
		{
		}

		~Parabolic() override = default;

		Parabolic(const Parabolic&) = delete;

		Parabolic& operator=(const Parabolic&) = delete;

		Parabolic(Parabolic&&) = delete;

		Parabolic& operator=(Parabolic&&) = delete;

		/**
		 * @brief Computes the gain of the parabolic antenna.
		 *
		 * @param angle The angle at which the gain is to be computed.
		 * @param refangle The reference angle.
		 * @param wavelength The wavelength of the signal.
		 * @return The computed gain of the antenna.
		 */
		[[nodiscard]] RealType getGain(const math::SVec3& angle, const math::SVec3& refangle,
									   RealType wavelength) const noexcept override;

		/** @brief Gets the diameter of the parabolic reflector. */
		[[nodiscard]] RealType getDiameter() const noexcept { return _diameter; }

		/**
		 * @brief Sets the diameter of the parabolic reflector.
		 *
		 * @param diameter The new diameter of the parabolic reflector.
		 */
		void setDiameter(RealType diameter) noexcept { _diameter = diameter; }

	private:
		RealType _diameter; ///< Diameter of the parabolic reflector.
	};

	/**
	 * @class XmlAntenna
	 * @brief Represents an antenna whose gain pattern is defined by an XML file.
	 *
	 * This class models an antenna where the gain pattern is read from an XML file.
	 */
	class XmlAntenna final : public Antenna
	{
	public:
		/// Symmetry mode for one-dimensional XML antenna gain axes.
		enum class AxisSymmetry
		{
			Mirrored, ///< Mirror positive-axis samples onto negative angles.
			None, ///< Use the axis samples exactly as provided.
		};

		/**
		 * @brief Constructs an XmlAntenna with the specified name and XML configuration file.
		 *
		 * The constructor loads the azimuth and elevation gain patterns from the provided XML file.
		 *
		 * @param name The name of the antenna.
		 * @param filename The path to the XML file containing the antenna's gain pattern data.
		 */
		XmlAntenna(const std::string_view name, const std::string_view filename, const SimId id = 0) :
			Antenna(name.data(), id), _azi_samples(std::make_unique<interp::InterpSet>()),
			_elev_samples(std::make_unique<interp::InterpSet>())
		{
			loadAntennaDescription(filename);
		}

		~XmlAntenna() override = default;

		XmlAntenna(const XmlAntenna&) = delete;

		XmlAntenna& operator=(const XmlAntenna&) = delete;

		XmlAntenna(XmlAntenna&&) = delete;

		XmlAntenna& operator=(XmlAntenna&&) = delete;

		/**
		 * @brief Computes the gain of the antenna based on the input angle and reference angle.
		 *
		 * @param angle The angle at which the gain is to be computed.
		 * @param refangle The reference angle.
		 * @param wavelength The wavelength of the signal (not used in this antenna type).
		 * @return The gain of the antenna at the specified angle.
		 * @throws std::runtime_error If gain values cannot be retrieved from the interpolation sets.
		 */
		[[nodiscard]] RealType getGain(const math::SVec3& angle, const math::SVec3& refangle,
									   RealType wavelength) const override;

		/** @brief Gets the filename of the antenna description. */
		[[nodiscard]] const std::string& getFilename() const noexcept { return _filename; }

		/** @brief Gets the maximum gain of the antenna. */
		[[nodiscard]] RealType getMaxGain() const noexcept { return _max_gain; }

		/** @brief Gets the interpolation set for azimuth gain samples. */
		[[nodiscard]] const interp::InterpSet* getAzimuthSamples() const noexcept { return _azi_samples.get(); }

		/** @brief Gets the interpolation set for elevation gain samples. */
		[[nodiscard]] const interp::InterpSet* getElevationSamples() const noexcept { return _elev_samples.get(); }

	private:
		/// Looks up a gain value from an XML antenna axis interpolation set.
		[[nodiscard]] std::optional<RealType> lookupAxisGain(const interp::InterpSet* set, RealType angle,
															 AxisSymmetry symmetry) const noexcept;

		/**
		 * @brief Loads the antenna gain pattern from the specified XML file.
		 *
		 * @param filename The path to the XML file containing the antenna's gain pattern data.
		 * @throws std::runtime_error If the XML file cannot be loaded or parsed.
		 */
		void loadAntennaDescription(std::string_view filename);

		std::string _filename; ///< The original filename for the antenna description.
		RealType _max_gain{}; ///< The maximum gain of the antenna.
		std::unique_ptr<interp::InterpSet> _azi_samples; ///< Interpolation set for azimuth gain samples.
		std::unique_ptr<interp::InterpSet> _elev_samples; ///< Interpolation set for elevation gain samples.
		AxisSymmetry _azi_symmetry{AxisSymmetry::Mirrored}; ///< Lookup mode for azimuth gain samples.
		AxisSymmetry _elev_symmetry{AxisSymmetry::Mirrored}; ///< Lookup mode for elevation gain samples.
	};

	/**
	 * @class H5Antenna
	 * @brief Represents an antenna whose gain pattern is loaded from a HDF5 file.
	 *
	 * This class models an antenna with a gain pattern defined in an HDF5 file. The gain pattern is stored in a
	 * `Pattern` object, which is used to compute the antenna's gain based on the input angle and reference angle.
	 */
	class H5Antenna final : public Antenna
	{
	public:
		/**
		 * @brief Constructs a H5Antenna with the specified name and gain pattern file.
		 *
		 * @param name The name of the antenna.
		 * @param filename The path to the file containing the antenna's gain pattern.
		 */
		H5Antenna(const std::string_view name, const std::string& filename, const SimId id = 0) :
			Antenna(name.data(), id), _pattern(serial::readPattern(filename, "antenna")), _filename(filename)
		{
		}

		~H5Antenna() override = default;

		H5Antenna(const H5Antenna&) = delete;

		H5Antenna& operator=(const H5Antenna&) = delete;

		H5Antenna(H5Antenna&&) = delete;

		H5Antenna& operator=(H5Antenna&&) = delete;

		/**
		 * @brief Computes the gain of the antenna based on the input angle and reference angle.
		 *
		 * @param angle The angle at which the gain is to be computed.
		 * @param refangle The reference angle.
		 * @return The gain of the antenna at the specified angle.
		 */
		[[nodiscard]] RealType getGain(const math::SVec3& angle, const math::SVec3& refangle,
									   RealType /*wavelength*/) const override;

		/** @brief Gets the filename of the antenna description. */
		[[nodiscard]] const std::string& getFilename() const noexcept { return _filename; }

		/** @brief Gets the gain pattern object. */
		[[nodiscard]] const std::vector<std::vector<RealType>>& getPattern() const noexcept { return _pattern; }

	private:
		std::vector<std::vector<RealType>> _pattern; ///< The 2D pattern data.
		std::string _filename; ///< The original filename for the antenna description.
	};
}
