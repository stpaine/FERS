// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// Upstream reference:
//   packages/libfers/src/processing/finalizer_pipeline.cpp

#include "finalizer_pipeline.h"

#include "core/parameters.h"
#include "processing/signal_processor.h"
#include "signal/dsp_filters.h"

namespace processing::pipeline
{
	RealType applyDownsamplingAndQuantization(std::vector<ComplexType>& buffer)
	{
		if (params::oversampleRatio() > 1)
		{
			buffer = std::move(fers_signal::downsample(buffer));
		}
		return quantizeAndScaleWindow(buffer);
	}
}

