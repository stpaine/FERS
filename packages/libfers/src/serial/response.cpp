// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file response.cpp
 * @brief Implementation of the Response class
 */

#include "response.h"

#include "core/sim_id.h"
#include "libxml_wrapper.h"
#include "radar/radar_obj.h"
#include "radar/transmitter.h"
#include "signal/radar_signal.h"

using interp::InterpPoint;

namespace serial
{
	SimId Response::getTransmitterId() const noexcept { return _transmitter->getId(); }

	void Response::addInterpPoint(const InterpPoint& point) { _points.push_back(point); }

	std::vector<ComplexType> Response::renderBinary(RealType& rate, unsigned& size, const RealType fracWinDelay) const
	{
		rate = _wave->getRate();
		return _wave->render(_points, size, fracWinDelay);
	}

	std::vector<ComplexType> Response::renderSlice(const RealType outputRate, const RealType outputStartTime,
												   const std::size_t sampleCount, const RealType fracWinDelay) const
	{
		return _wave->renderSlice(_points, outputStartTime, outputRate, sampleCount, fracWinDelay);
	}

	RealType Response::sampleRate() const noexcept { return _wave->getRate(); }

	unsigned Response::sampleCount() const noexcept { return _wave->getSampleCount(); }
}
