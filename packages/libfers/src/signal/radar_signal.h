// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2006-2008 Marc Brooker and Michael Inggs
// Copyright (c) 2008-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file radar_signal.h
 * @brief Classes for handling radar waveforms and signals.
 */

#pragma once

#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include "core/config.h"
#include "core/sim_id.h"

namespace interp
{
	struct InterpPoint;
}

namespace fers_signal
{
	/// Sweep direction for a linear FMCW chirp.
	enum class FmcwChirpDirection
	{
		Up, ///< Instantaneous baseband frequency increases over the chirp.
		Down ///< Instantaneous baseband frequency decreases over the chirp.
	};

	/// Converts a chirp direction to the schema token.
	[[nodiscard]] std::string_view fmcwChirpDirectionToken(FmcwChirpDirection direction) noexcept;

	/// Parses a schema chirp direction token.
	[[nodiscard]] FmcwChirpDirection parseFmcwChirpDirection(std::string_view direction);

	/**
	 * @class Signal
	 * @brief Class for handling radar waveform signal data.
	 */
	class Signal
	{
	public:
		virtual ~Signal() = default;

		Signal() = default;

		Signal(const Signal&) = delete;

		Signal& operator=(const Signal&) = delete;

		Signal(Signal&&) = default;

		Signal& operator=(Signal&&) = default;

		/**
		 * @brief Clears the internal signal data.
		 */
		void clear() noexcept;

		/**
		 * @brief Loads complex radar waveform data.
		 *
		 * @param inData The input span of complex signal data.
		 * @param samples The number of samples in the input data.
		 * @param sampleRate The sample rate of the input data.
		 */
		void load(std::span<const ComplexType> inData, unsigned samples, RealType sampleRate);

		/**
		 * @brief Gets the sample rate of the signal.
		 *
		 * @return The sample rate of the signal.
		 */
		[[nodiscard]] RealType getRate() const noexcept { return _rate; }

		/**
		 * @brief Renders the signal data based on interpolation points.
		 *
		 * @param points A vector of interpolation points used to render the signal.
		 * @param size Reference to store the size of the rendered data.
		 * @param fracWinDelay Fractional window delay to apply during rendering.
		 * @return A vector of rendered complex signal data.
		 */
		virtual std::vector<ComplexType> render(const std::vector<interp::InterpPoint>& points, unsigned& size,
												double fracWinDelay) const;

	private:
		std::vector<ComplexType> _data; ///< The complex signal data.
		unsigned _size{0}; ///< The size of the signal data.
		RealType _rate{0}; ///< The sample rate of the signal.

		/**
		 * @brief Calculates weights and delays for rendering.
		 *
		 * @param iter Iterator pointing to the current interpolation point.
		 * @param next Iterator pointing to the next interpolation point.
		 * @param sampleTime Current sample time.
		 * @param idelay Integer delay value.
		 * @param fracWinDelay Fractional window delay to apply.
		 * @return A tuple containing amplitude, phase, fractional delay, and sample unwrap index.
		 */
		[[nodiscard]] std::tuple<double, double, double, int>
		calculateWeightsAndDelays(std::vector<interp::InterpPoint>::const_iterator iter,
								  std::vector<interp::InterpPoint>::const_iterator next, double sampleTime,
								  double idelay, double fracWinDelay) const noexcept;

		/**
		 * @brief Performs convolution with a filter.
		 *
		 * @param i Index of the current sample.
		 * @param filt Pointer to the filter coefficients.
		 * @param filtLength Length of the filter.
		 * @param amplitude Amplitude scaling factor.
		 * @param iSampleUnwrap Unwrapped sample index for the convolution.
		 * @return The result of the convolution for the given sample.
		 */
		ComplexType performConvolution(int i, const double* filt, int filtLength, double amplitude,
									   int iSampleUnwrap) const noexcept;
	};

	/**
	 * @class RadarSignal
	 * @brief Class representing a radar signal with associated properties.
	 */
	class RadarSignal
	{
	public:
		/**
		 * @brief Constructs a RadarSignal object.
		 *
		 * @param name The name of the radar signal.
		 * @param power The power of the radar signal.
		 * @param carrierfreq The carrier frequency of the radar signal.
		 * @param length The length of the radar signal.
		 * @param signal A unique pointer to the `Signal` object containing the waveform data.
		 * @throws std::runtime_error if the signal is null.
		 */
		RadarSignal(std::string name, RealType power, RealType carrierfreq, RealType length,
					std::unique_ptr<Signal> signal, const SimId id = 0);

		~RadarSignal() = default;

		RadarSignal(const RadarSignal&) noexcept = delete;

		RadarSignal& operator=(const RadarSignal&) noexcept = delete;

		RadarSignal(RadarSignal&&) noexcept = delete;

		RadarSignal& operator=(RadarSignal&&) noexcept = delete;

		/**
		 * @brief Sets the filename associated with this signal.
		 * @param filename The source filename.
		 */
		void setFilename(const std::string& filename) noexcept { _filename = filename; }

		/**
		 * @brief Gets the filename associated with this signal.
		 * @return The source filename, if one was set.
		 */
		[[nodiscard]] const std::optional<std::string>& getFilename() const noexcept { return _filename; }

		/**
		 * @brief Gets the power of the radar signal.
		 *
		 * @return The power of the radar signal.
		 */
		[[nodiscard]] RealType getPower() const noexcept { return _power; }

		/**
		 * @brief Gets the carrier frequency of the radar signal.
		 *
		 * @return The carrier frequency of the radar signal.
		 */
		[[nodiscard]] RealType getCarrier() const noexcept { return _carrierfreq; }

		/**
		 * @brief Gets the name of the radar signal.
		 *
		 * @return The name of the radar signal.
		 */
		[[nodiscard]] const std::string& getName() const noexcept { return _name; }

		/**
		 * @brief Gets the unique ID of the radar signal.
		 *
		 * @return The radar signal SimId.
		 */
		[[nodiscard]] SimId getId() const noexcept { return _id; }

		/**
		 * @brief Gets the sample rate of the radar signal.
		 *
		 * @return The sample rate of the radar signal.
		 */
		[[nodiscard]] RealType getRate() const noexcept { return _signal->getRate(); }

		/**
		 * @brief Gets the length of the radar signal.
		 *
		 * @return The length of the radar signal.
		 */
		[[nodiscard]] RealType getLength() const noexcept { return _length; }

		/**
		 * @brief Gets the underlying signal object.
		 * @return A const pointer to the Signal object.
		 */
		[[nodiscard]] const Signal* getSignal() const noexcept { return _signal.get(); }

		/// Returns true when this signal is a continuous-wave signal.
		[[nodiscard]] bool isCw() const noexcept;

		/// Returns true when this signal is an FMCW linear chirp signal.
		[[nodiscard]] bool isFmcwChirp() const noexcept;

		/// Gets the FMCW chirp implementation, if this signal owns one.
		[[nodiscard]] const class FmcwChirpSignal* getFmcwChirpSignal() const noexcept;

		/**
		 * @brief Renders the radar signal.
		 *
		 * @param points A vector of interpolation points.
		 * @param size Reference to store the size of the rendered data.
		 * @param fracWinDelay Fractional window delay to apply during rendering.
		 * @return A vector of rendered complex radar signal data.
		 */
		std::vector<ComplexType> render(const std::vector<interp::InterpPoint>& points, unsigned& size,
										RealType fracWinDelay) const;

	private:
		std::string _name; ///< The name of the radar signal.
		SimId _id; ///< Unique ID for this radar signal.
		RealType _power; ///< The power of the radar signal.
		RealType _carrierfreq; ///< The carrier frequency of the radar signal.
		RealType _length; ///< The length of the radar signal.
		std::unique_ptr<Signal> _signal; ///< The `Signal` object containing the radar signal data.
		std::optional<std::string> _filename; ///< The original filename for file-based signals.
	};

	/// Continuous-wave signal implementation.
	class CwSignal final : public Signal
	{
	public:
		CwSignal() = default;

		~CwSignal() override = default;

		CwSignal(const CwSignal&) noexcept = delete;

		CwSignal& operator=(const CwSignal&) noexcept = delete;

		CwSignal(CwSignal&&) noexcept = delete;

		CwSignal& operator=(CwSignal&&) noexcept = delete;

		/**
		 * @brief Renders the signal data. For CW signals, this is a no-op.
		 * @return An empty vector of complex signal data.
		 */
		std::vector<ComplexType> render(const std::vector<interp::InterpPoint>& points, unsigned& size,
										RealType fracWinDelay) const override;
	};

	/// FMCW linear chirp signal implementation.
	class FmcwChirpSignal final : public Signal
	{
	public:
		/// Constructs an FMCW chirp signal with timing and sweep parameters.
		FmcwChirpSignal(RealType chirp_bandwidth, RealType chirp_duration, RealType chirp_period,
						RealType start_frequency_offset = 0.0, std::optional<std::size_t> chirp_count = std::nullopt,
						FmcwChirpDirection direction = FmcwChirpDirection::Up);

		~FmcwChirpSignal() override = default;

		FmcwChirpSignal(const FmcwChirpSignal&) noexcept = delete;

		FmcwChirpSignal& operator=(const FmcwChirpSignal&) noexcept = delete;

		FmcwChirpSignal(FmcwChirpSignal&&) noexcept = delete;

		FmcwChirpSignal& operator=(FmcwChirpSignal&&) noexcept = delete;

		/// Gets the chirp bandwidth in hertz.
		[[nodiscard]] RealType getChirpBandwidth() const noexcept { return _chirp_bandwidth; }

		/// Gets the chirp duration in seconds.
		[[nodiscard]] RealType getChirpDuration() const noexcept { return _chirp_duration; }

		/// Gets the chirp period in seconds.
		[[nodiscard]] RealType getChirpPeriod() const noexcept { return _chirp_period; }

		/// Gets the start frequency offset relative to carrier in hertz.
		[[nodiscard]] RealType getStartFrequencyOffset() const noexcept { return _start_frequency_offset; }

		/// Gets the optional finite chirp count.
		[[nodiscard]] const std::optional<std::size_t>& getChirpCount() const noexcept { return _chirp_count; }

		/// Gets the chirp rate in hertz per second.
		[[nodiscard]] RealType getChirpRate() const noexcept { return _chirp_rate; }

		/// Gets the signed chirp rate in hertz per second.
		[[nodiscard]] RealType getSignedChirpRate() const noexcept
		{
			return isDownChirp() ? -_chirp_rate : _chirp_rate;
		}

		/// Gets the FMCW sweep direction.
		[[nodiscard]] FmcwChirpDirection getDirection() const noexcept { return _direction; }

		/// Returns true when this chirp sweeps downward.
		[[nodiscard]] bool isDownChirp() const noexcept { return _direction == FmcwChirpDirection::Down; }

		/// Returns the active chirp index for a time since the segment start.
		[[nodiscard]] std::optional<std::size_t> activeChirpIndexAt(RealType time_since_segment_start) const noexcept;

		/// Returns true when the signal is inside an active chirp at the specified time.
		[[nodiscard]] bool isActiveAt(RealType time_since_segment_start) const noexcept
		{
			return activeChirpIndexAt(time_since_segment_start).has_value();
		}

		/// Computes baseband phase for a time inside a chirp.
		[[nodiscard]] RealType basebandPhaseForChirpTime(RealType chirp_time) const noexcept;

		/// Computes instantaneous baseband phase at a time since segment start.
		[[nodiscard]] std::optional<RealType>
		instantaneousBasebandPhase(RealType time_since_segment_start) const noexcept;

		/// Renders an FMCW waveform from interpolation points.
		std::vector<ComplexType> render(const std::vector<interp::InterpPoint>& points, unsigned& size,
										RealType fracWinDelay) const override;

	private:
		RealType _chirp_bandwidth{}; ///< Chirp bandwidth in hertz.
		RealType _chirp_duration{}; ///< Active chirp duration in seconds.
		RealType _chirp_period{}; ///< Chirp repetition period in seconds.
		RealType _start_frequency_offset{}; ///< Start frequency offset relative to carrier in hertz.
		std::optional<std::size_t> _chirp_count; ///< Optional finite chirp count.
		RealType _chirp_rate{}; ///< Frequency sweep rate in hertz per second.
		FmcwChirpDirection _direction{FmcwChirpDirection::Up}; ///< Frequency sweep direction.
	};
}
