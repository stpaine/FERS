// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2025-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file rendering_job.h
 * @brief Defines the data packet for asynchronous receiver finalization.
 */

#pragma once

#include <memory>
#include <vector>

#include "config.h"
#include "radar/transmitter.h"
#include "serial/response.h"

namespace core
{
	/**
	 * @struct RenderingJob
	 * @brief A data packet containing all information needed to process one receive window.
	 *
	 * This packet is created by the main simulation loop when a pulsed receiver's
	 * window ends. It is then passed to a dedicated finalizer thread for
	 * processing, decoupling the physics simulation from the expensive rendering
	 * and I/O tasks.
	 */
	struct RenderingJob
	{
		/// The ideal, jitter-free start time of the receive window.
		RealType ideal_start_time;

		/// The duration of the receive window in seconds.
		RealType duration;

		/// A list of all Response objects that overlap with this window.
		std::vector<std::unique_ptr<serial::Response>> responses;

		/// A list of all CW transmitters that were active during this window.
		std::vector<radar::Transmitter*> active_cw_sources;
	};
}
