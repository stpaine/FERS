// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// Upstream reference:
//   packages/libfers/src/serial/response.cpp

#include "response.h"

#include "signal/radar_signal.h"

namespace serial
{
	void Response::addInterpPoint(const interp::InterpPoint& point) { _points.push_back(point); }

	std::vector<ComplexType> Response::renderBinary(RealType& rate, unsigned& size, const RealType fracWinDelay) const
	{
		rate = _wave->getRate();
		return _wave->render(_points, size, fracWinDelay);
	}
}

