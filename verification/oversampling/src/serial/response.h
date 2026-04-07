// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// Upstream reference:
//   packages/libfers/src/serial/response.h
//
// Local simplification:
//   transmitter ownership is removed because the harness only needs render timing.

#pragma once

#include <vector>

#include "core/config.h"
#include "interpolation/interpolation_point.h"

namespace fers_signal
{
	class RadarSignal;
}

namespace serial
{
	class Response
	{
	public:
		explicit Response(const fers_signal::RadarSignal* wave) noexcept : _wave(wave) {}

		[[nodiscard]] RealType startTime() const noexcept { return _points.empty() ? 0.0 : _points.front().time; }
		[[nodiscard]] RealType endTime() const noexcept { return _points.empty() ? 0.0 : _points.back().time; }
		[[nodiscard]] RealType getLength() const noexcept { return endTime() - startTime(); }

		void addInterpPoint(const interp::InterpPoint& point);
		std::vector<ComplexType> renderBinary(RealType& rate, unsigned& size, RealType fracWinDelay) const;
		[[nodiscard]] const std::vector<interp::InterpPoint>& points() const noexcept { return _points; }

	private:
		const fers_signal::RadarSignal* _wave;
		std::vector<interp::InterpPoint> _points;
	};
}

