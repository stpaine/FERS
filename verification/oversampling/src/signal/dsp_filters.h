// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2007-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// Upstream reference:
//   packages/libfers/src/signal/dsp_filters.h

#pragma once

#include <memory>
#include <span>
#include <vector>

#include "core/config.h"

namespace fers_signal
{
	void upsample(std::span<const ComplexType> in, unsigned size, std::span<ComplexType> out);
	std::vector<ComplexType> downsample(std::span<const ComplexType> in);

	class DspFilter
	{
	public:
		DspFilter() = default;
		virtual ~DspFilter() = default;

		virtual RealType filter(RealType sample) = 0;
		virtual void filter(std::span<RealType> samples) = 0;

		DspFilter(const DspFilter&) = delete;
		DspFilter& operator=(const DspFilter&) = delete;
		DspFilter(DspFilter&&) noexcept = default;
		DspFilter& operator=(DspFilter&&) noexcept = default;
	};

	class IirFilter final : public DspFilter
	{
	public:
		IirFilter(const RealType* denCoeffs, const RealType* numCoeffs, unsigned order) noexcept;
		~IirFilter() override = default;

		RealType filter(RealType sample) noexcept override;
		void filter(std::span<RealType> samples) noexcept override;

	private:
		std::vector<RealType> _a;
		std::vector<RealType> _b;
		std::vector<RealType> _w;
		unsigned _order{};
	};

	class FirFilter final : public DspFilter
	{
	public:
		explicit FirFilter(std::span<const RealType> coeffs) :
			_filter(coeffs.begin(), coeffs.end()), _w(coeffs.size()), _order(coeffs.size())
		{
		}

		~FirFilter() override = default;

		RealType filter(RealType) override { return 0; }
		void filter(std::span<RealType>) noexcept override {}
		void filter(std::vector<ComplexType>& samples) const;

	private:
		std::vector<RealType> _filter;
		std::vector<RealType> _w;
		unsigned _order{};
	};

	class DecadeUpsampler
	{
	public:
		DecadeUpsampler();
		~DecadeUpsampler() = default;

		void upsample(RealType sample, std::span<RealType> out) const;

		DecadeUpsampler(const DecadeUpsampler&) = delete;
		DecadeUpsampler& operator=(const DecadeUpsampler&) = delete;
		DecadeUpsampler(DecadeUpsampler&&) noexcept = default;
		DecadeUpsampler& operator=(DecadeUpsampler&&) noexcept = default;

	private:
		std::unique_ptr<IirFilter> _filter;
	};
}

