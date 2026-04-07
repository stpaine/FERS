// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2024-present FERS Contributors (see AUTHORS.md).
//
// Upstream reference:
//   packages/libfers/src/core/config.h

#pragma once

#include <complex>
#include <limits>
#include <numbers>

using RealType = double;
using ComplexType = std::complex<RealType>;

constexpr RealType PI = std::numbers::pi_v<RealType>;
constexpr RealType EPSILON = std::numeric_limits<RealType>::epsilon();

