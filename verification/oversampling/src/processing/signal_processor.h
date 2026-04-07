// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// Upstream reference:
//   packages/libfers/src/processing/signal_processor.h

#pragma once

#include <memory>
#include <random>
#include <span>
#include <vector>

#include "core/config.h"

namespace serial
{
	class Response;
}

namespace processing
{
	void renderWindow(std::vector<ComplexType>& window, RealType length, RealType start, RealType fracDelay,
					  std::span<const std::unique_ptr<serial::Response>> responses);

	void applyThermalNoise(std::span<ComplexType> window, RealType noiseTemperature, std::mt19937& rngEngine);

	RealType quantizeAndScaleWindow(std::span<ComplexType> window);
}

