// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// Upstream reference:
//   packages/libfers/src/signal/radar_signal.h
//
// Local simplification:
//   SimId support is removed because the isolation harness only audits signal behavior.

#pragma once

#include <memory>
#include <optional>
#include <span>
#include <string>
#include <tuple>
#include <vector>

#include "core/config.h"

namespace interp
{
	struct InterpPoint;
}

namespace fers_signal
{
	class Signal
	{
	public:
		virtual ~Signal() = default;
		Signal() = default;

		Signal(const Signal&) = delete;
		Signal& operator=(const Signal&) = delete;
		Signal(Signal&&) = default;
		Signal& operator=(Signal&&) = default;

		void clear() noexcept;
		void load(std::span<const ComplexType> inData, unsigned samples, RealType sampleRate);
		[[nodiscard]] RealType getRate() const noexcept { return _rate; }

		virtual std::vector<ComplexType> render(const std::vector<interp::InterpPoint>& points, unsigned& size,
												double fracWinDelay) const;

	private:
		std::vector<ComplexType> _data;
		unsigned _size{0};
		RealType _rate{0};

		[[nodiscard]] constexpr std::tuple<double, double, double, int>
		calculateWeightsAndDelays(std::vector<interp::InterpPoint>::const_iterator iter,
								  std::vector<interp::InterpPoint>::const_iterator next, double sampleTime, double idelay,
								  double fracWinDelay) const noexcept;

		ComplexType performConvolution(int i, const double* filt, int filtLength, double amplitude,
									   int iSampleUnwrap) const noexcept;
	};

	class RadarSignal
	{
	public:
		RadarSignal(std::string name, RealType power, RealType carrierfreq, RealType length, std::unique_ptr<Signal> signal);

		~RadarSignal() = default;
		RadarSignal(const RadarSignal&) = delete;
		RadarSignal& operator=(const RadarSignal&) = delete;
		RadarSignal(RadarSignal&&) noexcept = delete;
		RadarSignal& operator=(RadarSignal&&) noexcept = delete;

		void setFilename(const std::string& filename) noexcept { _filename = filename; }
		[[nodiscard]] const std::optional<std::string>& getFilename() const noexcept { return _filename; }
		[[nodiscard]] RealType getPower() const noexcept { return _power; }
		[[nodiscard]] RealType getCarrier() const noexcept { return _carrierfreq; }
		[[nodiscard]] const std::string& getName() const noexcept { return _name; }
		[[nodiscard]] RealType getRate() const noexcept { return _signal->getRate(); }
		[[nodiscard]] RealType getLength() const noexcept { return _length; }
		[[nodiscard]] const Signal* getSignal() const noexcept { return _signal.get(); }

		std::vector<ComplexType> render(const std::vector<interp::InterpPoint>& points, unsigned& size,
										RealType fracWinDelay) const;

	private:
		std::string _name;
		RealType _power;
		RealType _carrierfreq;
		RealType _length;
		std::unique_ptr<Signal> _signal;
		std::optional<std::string> _filename;
	};

	class CwSignal final : public Signal
	{
	public:
		CwSignal() = default;
		~CwSignal() override = default;

		std::vector<ComplexType> render(const std::vector<interp::InterpPoint>& points, unsigned& size,
										RealType fracWinDelay) const override;
	};
}

