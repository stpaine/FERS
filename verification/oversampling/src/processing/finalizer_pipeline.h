// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// Upstream reference:
//   packages/libfers/src/processing/finalizer_pipeline.h
//
// Local simplification:
//   only the oversampling-specific downsampling + quantization step is required here.

#pragma once

#include <vector>

#include "core/config.h"

namespace processing::pipeline
{
	RealType applyDownsamplingAndQuantization(std::vector<ComplexType>& buffer);
}

