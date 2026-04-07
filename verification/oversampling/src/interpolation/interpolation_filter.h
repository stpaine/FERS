// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// Upstream reference:
//   packages/libfers/src/interpolation/interpolation_filter.h

#pragma once

#include <expected>
#include <span>
#include <vector>

#include "core/config.h"

namespace interp
{
	class InterpFilter
	{
	public:
		InterpFilter(const InterpFilter&) = delete;
		InterpFilter(InterpFilter&&) = delete;
		InterpFilter& operator=(const InterpFilter&) = delete;
		InterpFilter& operator=(InterpFilter&&) = delete;
		~InterpFilter() = default;

		static constexpr RealType sinc(const RealType x) noexcept
		{
			return x == 0.0 ? 1.0 : std::sin(x * PI) / (x * PI);
		}

		[[nodiscard]] std::expected<RealType, std::string> kaiserWinCompute(RealType x) const noexcept;
		[[nodiscard]] std::expected<RealType, std::string> interpFilter(RealType x) const noexcept;
		[[nodiscard]] std::span<const RealType> getFilter(RealType delay) const;
		static InterpFilter& getInstance() noexcept;

	private:
		InterpFilter();

		RealType _alpha;
		RealType _beta = 5;
		RealType _bessel_beta;
		int _length;
		int _table_filters;
		std::vector<RealType> _filter_table;
	};
}

