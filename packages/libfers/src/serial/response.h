// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file response.h
 * @brief Classes for managing radar signal responses.
 */

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "core/config.h"
#include "core/sim_id.h"
#include "interpolation/interpolation_point.h"

class XmlElement;

namespace fers_signal
{
	class RadarSignal;
}

namespace radar
{
	class Transmitter;
}

namespace serial
{
	/**
	 * @class Response
	 * @brief Manages radar signal responses from a transmitter.
	 */
	class Response
	{
	public:
		/**
		 * @brief Constructor for the `Response` class.
		 *
		 * @param wave Pointer to the radar signal object.
		 * @param transmitter Pointer to the transmitter object.
		 */
		Response(const fers_signal::RadarSignal* wave, const radar::Transmitter* transmitter) noexcept :
			_transmitter(transmitter), _wave(wave)
		{
		}

		~Response() = default;
		Response(const Response&) = delete;
		Response& operator=(const Response&) = delete;
		Response(Response&&) = delete;
		Response& operator=(Response&&) = delete;

		/**
		 * @brief Retrieves the start time of the response.
		 *
		 * @return Start time as a `RealType`. Returns 0.0 if no points are present.
		 */
		[[nodiscard]] RealType startTime() const noexcept { return _points.empty() ? 0.0 : _points.front().time; }

		/**
		 * @brief Retrieves the end time of the response.
		 *
		 * @return End time as a `RealType`. Returns 0.0 if no points are present.
		 */
		[[nodiscard]] RealType endTime() const noexcept { return _points.empty() ? 0.0 : _points.back().time; }

		/**
		 * @brief Adds an interpolation point to the response.
		 *
		 * @param point The interpolation point to be added.
		 * @throws std::logic_error If the new point has a time earlier than the last point.
		 */
		void addInterpPoint(const interp::InterpPoint& point);

		/**
		 * @brief Renders the response in binary format.
		 *
		 * @param rate Output parameter for the signal rate.
		 * @param size Output parameter for the size of the binary data.
		 * @param fracWinDelay Delay factor applied during windowing.
		 * @return A vector of `ComplexType` representing the binary data.
		 */
		std::vector<ComplexType> renderBinary(RealType& rate, unsigned& size, RealType fracWinDelay) const;

		/// Renders a bounded absolute-time response slice on the requested output grid.
		[[nodiscard]] std::vector<ComplexType> renderSlice(RealType outputRate, RealType outputStartTime,
														   std::size_t sampleCount, RealType fracWinDelay) const;

		/// Returns the waveform native sample rate.
		[[nodiscard]] RealType sampleRate() const noexcept;

		/// Returns the waveform native sample count.
		[[nodiscard]] unsigned sampleCount() const noexcept;

		/**
		 * @brief Retrieves the length of the response.
		 *
		 * @return The length of the response as a `RealType`.
		 */
		[[nodiscard]] RealType getLength() const noexcept { return endTime() - startTime(); }

		/**
		 * @brief Retrieves the ID of the associated transmitter.
		 *
		 * @return The transmitter SimId.
		 */
		[[nodiscard]] SimId getTransmitterId() const noexcept;

	private:
		const radar::Transmitter* _transmitter; ///< Pointer to the transmitter object.
		const fers_signal::RadarSignal* _wave; ///< Pointer to the radar signal object.
		std::vector<interp::InterpPoint> _points; ///< Vector of interpolation points.
	};
}
