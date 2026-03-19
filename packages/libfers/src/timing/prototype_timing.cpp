// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file prototype_timing.cpp
 * @brief Implementation file for the PrototypeTiming class.
 */

#include "prototype_timing.h"

#include "core/logging.h"

using logging::Level;

namespace timing
{
	void PrototypeTiming::setAlpha(const RealType alpha, const RealType weight) noexcept
	{
		_alphas.emplace_back(alpha);
		_weights.emplace_back(weight);
	}

	void PrototypeTiming::copyAlphas(std::vector<RealType>& alphas, std::vector<RealType>& weights) const noexcept
	{
		alphas = _alphas;
		weights = _weights;
	}
}
