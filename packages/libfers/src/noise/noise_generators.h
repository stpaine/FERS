// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file noise_generators.h
 * @brief Header file for noise generator classes.
 */

#pragma once

#include <functional>
#include <memory>
#include <random>
#include <vector>

#include "core/config.h"
#include "falpha_branch.h"

namespace noise
{
	/**
	 * @class NoiseGenerator
	 * @brief Abstract base class for noise generators.
	 */
	class NoiseGenerator
	{
	public:
		NoiseGenerator() = default;

		virtual ~NoiseGenerator() = default;

		NoiseGenerator(const NoiseGenerator&) = delete;

		NoiseGenerator& operator=(const NoiseGenerator&) = delete;

		NoiseGenerator(NoiseGenerator&&) = delete;

		NoiseGenerator& operator=(NoiseGenerator&&) = delete;

		/**
		 * @brief Pure virtual method to generate a noise sample.
		 *
		 * @return A noise sample of type RealType.
		 */
		virtual RealType getSample() = 0;
	};

	/**
	 * @class WgnGenerator
	 * @brief Generates white Gaussian noise.
	 */
	class WgnGenerator final : public NoiseGenerator
	{
	public:
		/**
		 * @brief Constructor to initialize the WGN generator with a given standard deviation.
		 *
		 * @param rngEngine The random number engine to use for generation.
		 * @param stddev The standard deviation of the generated Gaussian noise. Default is 1.0.
		 */
		explicit WgnGenerator(std::mt19937& rngEngine, const RealType stddev = 1.0) noexcept :
			_rng_engine(rngEngine), _dist(0.0, stddev), _stddev(stddev)
		{
		}

		/**
		 * @brief Generates a sample of white Gaussian noise.
		 *
		 * @return A noise sample of type RealType.
		 */
		RealType getSample() noexcept override { return _dist(_rng_engine.get()); }

	private:
		std::reference_wrapper<std::mt19937> _rng_engine; ///< Reference to the RNG engine.
		std::normal_distribution<> _dist; ///< Normal distribution for generating Gaussian noise.
		RealType _stddev; ///< Standard deviation of the generated noise.
	};

	/**
	 * @class GammaGenerator
	 * @brief Generates Gamma-distributed noise.
	 */
	class GammaGenerator final : public NoiseGenerator
	{
	public:
		/**
		 * @brief Constructor to initialize the Gamma generator with a shape parameter.
		 *
		 * @param rngEngine The random number engine to use for generation.
		 * @param k The shape parameter of the Gamma distribution.
		 */
		explicit GammaGenerator(std::mt19937& rngEngine, const RealType k) noexcept :
			_rng_engine(rngEngine), _dist(k, 1.0)
		{
		}

		/**
		 * @brief Generates a sample of Gamma noise.
		 *
		 * @return A noise sample of type RealType.
		 */
		RealType getSample() noexcept override { return _dist(_rng_engine.get()); }

	private:
		std::reference_wrapper<std::mt19937> _rng_engine; ///< Reference to the RNG engine.
		std::gamma_distribution<> _dist; ///< Gamma distribution for generating noise.
	};

	/**
	 * @class MultirateGenerator
	 * @brief Generates multirate noise using a hierarchical tree structure.
	 */
	class MultirateGenerator final : public NoiseGenerator
	{
	public:
		/**
		 * @brief Constructor to initialize the multirate generator.
		 *
		 * @param rngEngine The random number engine to use for generation.
		 * @param alpha The scaling parameter that controls the noise properties.
		 * @param branches The number of branches in the tree structure.
		 */
		MultirateGenerator(std::mt19937& rngEngine, RealType alpha, unsigned branches);

		/**
		 * @brief Generates a multirate noise sample.
		 *
		 * @return A noise sample of type RealType.
		 */
		RealType getSample() override { return _topbranch->getSample() * _scale; }

		/**
		 * @brief Skips a number of samples in the noise sequence.
		 *
		 * @param samples The number of samples to skip.
		 */
		void skipSamples(long long samples) noexcept;

		/**
		 * @brief Resets the noise generator state.
		 */
		void reset() noexcept;

	private:
		RealType _scale; ///< Scaling factor for the noise.
		std::unique_ptr<FAlphaBranch> _topbranch; ///< Pointer to the top branch in the tree structure.
		std::reference_wrapper<std::mt19937> _rng_engine; ///< Reference to the RNG engine.

		/**
		 * @brief Helper method to create the hierarchical tree structure.
		 *
		 * @param fAlpha Fractional alpha value.
		 * @param fInt Integer part of the scaling factor.
		 * @param branches The number of branches in the tree.
		 */
		void createTree(RealType fAlpha, int fInt, unsigned branches);
	};

	/**
	 * @class ClockModelGenerator
	 * @brief Generates noise using a clock model with multiple rates.
	 *
	 * This is useful for simulating clock jitter or other similar phenomena.
	 */
	class ClockModelGenerator final : public NoiseGenerator
	{
	public:
		/**
		 * @brief Constructor to initialize the clock model generator.
		 *
		 * @param rngEngine The random number engine to use for generation.
		 * @param alpha Vector of scaling parameters for the noise.
		 * @param inWeights Vector of weights for each rate process.
		 * @param frequency The base frequency of the clock model.
		 * @param phaseOffset The phase offset of the generated noise.
		 * @param freqOffset The frequency offset of the generated noise.
		 * @param branches The number of branches in each rate process.
		 */
		ClockModelGenerator(std::mt19937& rngEngine, const std::vector<RealType>& alpha,
							const std::vector<RealType>& inWeights, RealType frequency, RealType phaseOffset,
							RealType freqOffset, int branches) noexcept;

		/**
		 * @brief Generates a clock model noise sample.
		 *
		 * @return A noise sample of type RealType.
		 */
		RealType getSample() override;

		/**
		 * @brief Skips a number of samples in the noise sequence.
		 *
		 * @param samples The number of samples to skip.
		 */
		void skipSamples(long long samples);

		/**
		 * @brief Resets the noise generator state.
		 */
		void reset();

		/**
		 * @brief Checks if the noise generator is enabled.
		 *
		 * @return True if the generator is enabled, false otherwise.
		 */
		[[nodiscard]] bool enabled() const;

	private:
		std::reference_wrapper<std::mt19937> _rng_engine; ///< Reference to the RNG engine.
		std::vector<std::unique_ptr<MultirateGenerator>> _generators; ///< Multirate noise generators.
		std::vector<RealType> _weights; ///< Weights for each rate process.
		RealType _phase_offset; ///< Phase offset for the noise.
		RealType _freq_offset; ///< Frequency offset for the noise.
		RealType _frequency; ///< Base frequency of the clock model.
		unsigned long _count = 0; ///< Sample counter.
	};
}
