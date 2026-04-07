// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// Upstream reference:
//   packages/libfers/src/interpolation/interpolation_point.h

#pragma once

#include "core/config.h"

namespace interp
{
	struct InterpPoint
	{
		RealType power{};
		RealType time{};
		RealType delay{};
		RealType phase{};
	};
}

