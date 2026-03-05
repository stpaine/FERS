// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file timing.cpp
 * @brief Implementation of timing sources.
 */

#include "timing.h"

#include <stdexcept>

#include "core/logging.h"
#include "prototype_timing.h"

using logging::Level;

namespace timing
{
	Timing::Timing(std::string name, const unsigned seed, const SimId id) noexcept :
		_name(std::move(name)), _id(id == 0 ? SimIdGenerator::instance().generateId(ObjectType::Timing) : id),
		_rng(seed), _seed(seed)
	{
	}

	// NOLINTNEXTLINE(readability-make-member-function-const)
	void Timing::skipSamples(const long long samples) noexcept
	{
		if (_enabled && _model)
		{
			_model->skipSamples(samples);
		}
	}

	void Timing::initializeModel(const PrototypeTiming* timing) noexcept
	{
		if (_model)
		{
			LOG(Level::WARNING, "Timing source '{}' already initialized. Skipping re-initialization.", _name);
			return;
		}

		_prototype = timing;
		_frequency = timing->getFrequency();

		std::normal_distribution normal_dist{0.0, 1.0};

		_freq_offset = timing->getFreqOffset().value_or(0);
		if (const std::optional<RealType> random_freq_stdev = timing->getRandomFreqOffsetStdev(); random_freq_stdev)
		{
			LOG(Level::INFO, "Timing source '{}': applying random frequency offset with stdev {} Hz.", _name,
				random_freq_stdev.value());
			_freq_offset += normal_dist(_rng) * random_freq_stdev.value();
		}

		_phase_offset = timing->getPhaseOffset().value_or(0);
		if (const std::optional<RealType> random_phase_stdev = timing->getRandomPhaseOffsetStdev(); random_phase_stdev)
		{
			LOG(Level::INFO, "Timing source '{}': applying random phase offset with stdev {} radians.", _name,
				random_phase_stdev.value());
			_phase_offset += normal_dist(_rng) * random_phase_stdev.value();
		}

		timing->copyAlphas(_alphas, _weights);

		_model = std::make_unique<noise::ClockModelGenerator>(_rng, _alphas, _weights, _frequency, _phase_offset,
															  _freq_offset, 15);

		if (timing->getFrequency() == 0.0f)
		{
			LOG(Level::INFO, "Timing source frequency not set, results could be incorrect.");
		}

		_sync_on_pulse = timing->getSyncOnPulse();
		_enabled = true;
	}

	std::unique_ptr<Timing> Timing::clone() const
	{
		if (!_prototype)
		{
			LOG(Level::FATAL, "Cannot clone a Timing object that has not been initialized from a prototype.");
			throw std::logic_error("Cannot clone a Timing object that has not been initialized from a prototype.");
		}
		auto new_timing = std::make_unique<Timing>(_name, _seed, _id);
		new_timing->initializeModel(_prototype);
		return new_timing;
	}
}
