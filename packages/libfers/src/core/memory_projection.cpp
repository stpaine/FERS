// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "memory_projection.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <format>
#include <limits>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>
#include <unordered_map>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

#include "core/logging.h"
#include "core/parameters.h"
#include "core/sim_id.h"
#include "core/world.h"
#include "radar/receiver.h"
#include "radar/transmitter.h"
#include "timing/timing.h"

namespace core
{
	namespace
	{
		using radar::OperationMode;

		constexpr std::uint64_t max_uint64 = std::numeric_limits<std::uint64_t>::max(); ///< Saturating byte-count cap.

		/**
		 * @brief Checks whether a receiver stores a full-duration streaming IQ buffer.
		 * @param receiver The receiver to inspect.
		 * @return True for CW and FMCW receivers, false otherwise.
		 */
		[[nodiscard]] bool isStreamingReceiver(const radar::Receiver& receiver) noexcept
		{
			return receiver.getMode() == OperationMode::CW_MODE || receiver.getMode() == OperationMode::FMCW_MODE;
		}

		/**
		 * @brief Converts a finite floating-point count to an integer by rounding up.
		 * @param value The floating-point count to convert.
		 * @param overflowed Set to true when the input cannot be represented as `uint64_t`.
		 * @return The rounded-up count, clamped to `uint64_t` max on overflow.
		 */
		[[nodiscard]] std::uint64_t ceilToUint64(const RealType value, bool& overflowed) noexcept
		{
			if (!std::isfinite(value))
			{
				overflowed = true;
				return max_uint64;
			}
			if (value <= 0.0)
			{
				return 0;
			}
			if (value >= static_cast<RealType>(max_uint64))
			{
				overflowed = true;
				return max_uint64;
			}
			const RealType nearest = std::round(value);
			const RealType tolerance = 1.0e-12 * std::max<RealType>(1.0, std::abs(nearest));
			if (std::abs(value - nearest) <= tolerance)
			{
				return static_cast<std::uint64_t>(nearest);
			}
			return static_cast<std::uint64_t>(std::ceil(value));
		}

		/**
		 * @brief Adds two byte projections while preserving overflow state.
		 * @param lhs The left-hand byte projection.
		 * @param rhs The right-hand byte projection.
		 * @return The summed projection, clamped to `uint64_t` max on overflow.
		 */
		[[nodiscard]] ByteProjection addBytes(const ByteProjection lhs, const ByteProjection rhs) noexcept
		{
			ByteProjection result{};
			result.overflowed = lhs.overflowed || rhs.overflowed;
			if (max_uint64 - lhs.bytes < rhs.bytes)
			{
				result.bytes = max_uint64;
				result.overflowed = true;
				return result;
			}
			result.bytes = lhs.bytes + rhs.bytes;
			return result;
		}

		/**
		 * @brief Multiplies two byte-count factors while preserving prior overflow state.
		 * @param lhs The left-hand factor.
		 * @param rhs The right-hand factor.
		 * @param input_overflowed True if an upstream calculation has already overflowed.
		 * @return The product projection, clamped to `uint64_t` max on overflow.
		 */
		[[nodiscard]] ByteProjection multiplyBytes(const std::uint64_t lhs, const std::uint64_t rhs,
												   const bool input_overflowed = false) noexcept
		{
			ByteProjection result{};
			result.overflowed = input_overflowed;
			if (lhs != 0 && rhs > max_uint64 / lhs)
			{
				result.bytes = max_uint64;
				result.overflowed = true;
				return result;
			}
			result.bytes = lhs * rhs;
			return result;
		}

		/**
		 * @brief Subtracts a byte projection from an RSS value without underflowing.
		 * @param lhs The baseline byte count.
		 * @param rhs The byte projection to subtract.
		 * @return The clamped difference with the input overflow flag preserved.
		 */
		[[nodiscard]] ByteProjection subtractClamped(const std::uint64_t lhs, const ByteProjection rhs) noexcept
		{
			return {.bytes = lhs > rhs.bytes ? lhs - rhs.bytes : 0, .overflowed = rhs.overflowed};
		}

		/**
		 * @brief Converts an oversampled sample count to the rendered output sample count.
		 * @param oversampled_samples The sample count at the internal simulation rate.
		 * @return The sample count after applying the configured oversample ratio.
		 */
		[[nodiscard]] std::uint64_t downsampledSampleCount(const std::uint64_t oversampled_samples) noexcept
		{
			const unsigned ratio = params::oversampleRatio();
			if (ratio <= 1)
			{
				return oversampled_samples;
			}
			return oversampled_samples / ratio;
		}

		/**
		 * @brief Counts samples required for a duration at a given sample rate.
		 * @param duration_seconds Duration of the interval in seconds.
		 * @param sample_rate_hz Sample rate used for the interval.
		 * @param overflowed Set to true if the count cannot be represented as `uint64_t`.
		 * @return The rounded-up sample count for the interval.
		 */
		[[nodiscard]] std::uint64_t countSamplesForDuration(const RealType duration_seconds,
															const RealType sample_rate_hz, bool& overflowed) noexcept
		{
			return ceilToUint64(duration_seconds * sample_rate_hz, overflowed);
		}

		/**
		 * @brief Counts evenly spaced start times within an inclusive time range.
		 * @param first_start First candidate start time.
		 * @param last_start Last allowed start time.
		 * @param step_seconds Spacing between successive starts.
		 * @param overflowed Set to true if the count cannot be represented as `uint64_t`.
		 * @return The number of starts in range, or zero for an invalid range.
		 */
		[[nodiscard]] std::uint64_t countStartsInRange(const RealType first_start, const RealType last_start,
													   const RealType step_seconds, bool& overflowed) noexcept
		{
			if (first_start > last_start || step_seconds <= 0.0 || !std::isfinite(step_seconds))
			{
				return 0;
			}
			return ceilToUint64(std::floor((last_start - first_start) / step_seconds) + 1.0, overflowed);
		}

		/**
		 * @brief Projects the number of receive windows emitted by a pulsed receiver.
		 * @param receiver The pulsed receiver to inspect.
		 * @param overflowed Set to true if any window count arithmetic overflows.
		 * @return The projected number of receive windows during the simulation.
		 */
		[[nodiscard]] std::uint64_t countPulsedWindows(const radar::Receiver& receiver, bool& overflowed)
		{
			const RealType prf = receiver.getWindowPrf();
			if (prf <= 0.0 || !std::isfinite(prf))
			{
				return 0;
			}

			const RealType sim_end = params::endTime();
			const RealType step_seconds = 1.0 / prf;
			const auto& schedule = receiver.getSchedule();
			std::uint64_t total = 0;

			if (schedule.empty())
			{
				const RealType first_start = receiver.getWindowStart(0);
				if (first_start >= sim_end)
				{
					return 0;
				}
				const auto count = countStartsInRange(first_start, sim_end, step_seconds, overflowed);
				return count;
			}

			RealType next_requested = receiver.getWindowStart(0);
			bool counted_any_window = false;
			for (const auto& period : schedule)
			{
				const RealType period_end = std::min(period.end, sim_end);
				if (period_end < period.start || next_requested > period_end)
				{
					continue;
				}

				const RealType first_start = std::max(next_requested, period.start);
				if (!counted_any_window && first_start >= sim_end)
				{
					break;
				}

				const auto count = countStartsInRange(first_start, period_end, step_seconds, overflowed);
				if (count == 0)
				{
					continue;
				}

				const auto added = addBytes({.bytes = total}, {.bytes = count});
				total = added.bytes;
				overflowed = overflowed || added.overflowed;
				counted_any_window = true;
				next_requested = first_start + static_cast<RealType>(count) * step_seconds;
			}

			return total;
		}

		/**
		 * @brief Reads the process resident set size from the current platform when supported.
		 * @return The current resident set size in bytes, or `std::nullopt` when unavailable.
		 */
		[[nodiscard]] std::optional<std::uint64_t> currentResidentSetBytes() noexcept
		{
#if defined(__linux__)
			long page_size = sysconf(_SC_PAGESIZE);
			if (page_size <= 0)
			{
				return std::nullopt;
			}

			FILE* statm = std::fopen("/proc/self/statm", "r");
			if (statm == nullptr)
			{
				return std::nullopt;
			}

			unsigned long resident_pages = 0;
			const int scanned = std::fscanf(statm, "%*s %lu", &resident_pages);
			std::fclose(statm);
			if (scanned != 1)
			{
				return std::nullopt;
			}

			const auto pages = static_cast<std::uint64_t>(resident_pages);
			const auto bytes = multiplyBytes(pages, static_cast<std::uint64_t>(page_size));
			if (bytes.overflowed)
			{
				return max_uint64;
			}
			return bytes.bytes;
#elif defined(__APPLE__) && defined(__MACH__)
			mach_task_basic_info_data_t info{};
			mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
			if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info), &count) !=
				KERN_SUCCESS)
			{
				return std::nullopt;
			}
			return static_cast<std::uint64_t>(info.resident_size);
#else
			return std::nullopt;
#endif
		}

		/**
		 * @brief Converts a byte projection to the JSON shape used by memory reports.
		 * @param projection The byte projection to serialize.
		 * @return A JSON object containing raw bytes, formatted bytes, and overflow state.
		 */
		[[nodiscard]] nlohmann::json byteProjectionToJson(const ByteProjection projection)
		{
			return {{"bytes", projection.bytes},
					{"human", formatByteSize(projection.bytes)},
					{"overflowed", projection.overflowed}};
		}

		/**
		 * @brief Converts an optional byte projection to JSON, preserving unknown values.
		 * @param projection The optional byte projection to serialize.
		 * @return A JSON object with nullable bytes when the projection is unavailable.
		 */
		[[nodiscard]] nlohmann::json optionalByteProjectionToJson(const std::optional<ByteProjection>& projection)
		{
			if (!projection.has_value())
			{
				return {{"bytes", nullptr}, {"human", "unknown"}, {"overflowed", false}};
			}
			return byteProjectionToJson(*projection);
		}
	}

	std::vector<std::shared_ptr<timing::Timing>> collectCwPhaseNoiseTimings(const World& world)
	{
		std::unordered_map<SimId, std::shared_ptr<timing::Timing>> unique_timings;

		for (const auto& transmitter_ptr : world.getTransmitters())
		{
			if (!transmitter_ptr->isStreamingMode())
			{
				continue;
			}
			unique_timings.try_emplace(transmitter_ptr->getTiming()->getId(), transmitter_ptr->getTiming());
		}

		for (const auto& receiver_ptr : world.getReceivers())
		{
			if (!isStreamingReceiver(*receiver_ptr))
			{
				continue;
			}
			unique_timings.try_emplace(receiver_ptr->getTiming()->getId(), receiver_ptr->getTiming());
		}

		std::vector<std::shared_ptr<timing::Timing>> timings;
		timings.reserve(unique_timings.size());
		for (const auto& entry : unique_timings)
		{
			timings.push_back(entry.second);
		}
		return timings;
	}

	SimulationMemoryProjection projectSimulationMemory(const World& world)
	{
		SimulationMemoryProjection projection{};
		projection.duration_seconds = std::max<RealType>(0.0, params::endTime() - params::startTime());
		projection.oversample_ratio = params::oversampleRatio();
		projection.simulation_sample_rate_hz = params::rate() * static_cast<RealType>(projection.oversample_ratio);

		bool sample_count_overflowed = false;
		projection.streaming_sample_count = countSamplesForDuration(
			projection.duration_seconds, projection.simulation_sample_rate_hz, sample_count_overflowed);

		const RealType phase_noise_lookup_start = world.earliestPhaseNoiseLookupStart();
		bool phase_noise_count_overflowed = false;
		projection.phase_noise_sample_count =
			countSamplesForDuration(std::max<RealType>(0.0, params::endTime() - phase_noise_lookup_start),
									projection.simulation_sample_rate_hz, phase_noise_count_overflowed);
		if (projection.phase_noise_sample_count != max_uint64)
		{
			++projection.phase_noise_sample_count;
		}
		if (phase_noise_count_overflowed)
		{
			projection.phase_noise_sample_count = max_uint64;
		}

		const auto timings = collectCwPhaseNoiseTimings(world);
		projection.phase_noise_timing_count = static_cast<std::uint64_t>(timings.size());
		for (const auto& timing : timings)
		{
			if (timing && timing->isEnabled())
			{
				++projection.enabled_phase_noise_timing_count;
			}
		}

		projection.phase_noise_lookup =
			multiplyBytes(projection.phase_noise_sample_count,
						  projection.enabled_phase_noise_timing_count * static_cast<std::uint64_t>(sizeof(RealType)),
						  sample_count_overflowed);

		for (const auto& receiver_ptr : world.getReceivers())
		{
			const auto& receiver = *receiver_ptr;
			if (isStreamingReceiver(receiver))
			{
				++projection.streaming_receiver_count;
				const auto capacity = static_cast<std::uint64_t>(receiver.getStreamingData().capacity());
				projection.allocated_streaming_iq_buffers =
					addBytes(projection.allocated_streaming_iq_buffers,
							 multiplyBytes(capacity, static_cast<std::uint64_t>(sizeof(ComplexType))));
				const auto rendered_samples = receiver.isDechirpEnabled()
					? projection.streaming_sample_count
					: downsampledSampleCount(projection.streaming_sample_count);
				projection.rendered_hdf5_sample_count =
					addBytes({.bytes = projection.rendered_hdf5_sample_count},
							 {.bytes = rendered_samples, .overflowed = sample_count_overflowed})
						.bytes;
				continue;
			}

			if (receiver.getMode() == OperationMode::PULSED_MODE)
			{
				++projection.pulsed_receiver_count;
				bool window_count_overflowed = false;
				const auto windows = countPulsedWindows(receiver, window_count_overflowed);
				projection.pulsed_window_count = addBytes({.bytes = projection.pulsed_window_count},
														  {.bytes = windows, .overflowed = window_count_overflowed})
													 .bytes;

				bool pulsed_sample_count_overflowed = false;
				const auto raw_window_samples = countSamplesForDuration(
					receiver.getWindowLength(), projection.simulation_sample_rate_hz, pulsed_sample_count_overflowed);
				const auto rendered_window_samples = downsampledSampleCount(raw_window_samples);
				const auto rendered_samples = multiplyBytes(windows, rendered_window_samples,
															window_count_overflowed || pulsed_sample_count_overflowed);
				projection.rendered_hdf5_sample_count =
					addBytes({.bytes = projection.rendered_hdf5_sample_count}, rendered_samples).bytes;
			}
		}

		projection.streaming_iq_buffers =
			multiplyBytes(projection.streaming_sample_count,
						  projection.streaming_receiver_count * static_cast<std::uint64_t>(sizeof(ComplexType)),
						  sample_count_overflowed);

		projection.rendered_hdf5_payload =
			multiplyBytes(projection.rendered_hdf5_sample_count, 2ULL * static_cast<std::uint64_t>(sizeof(RealType)));

		projection.current_resident_set = currentResidentSetBytes();
		if (projection.current_resident_set.has_value())
		{
			projection.resident_baseline =
				subtractClamped(*projection.current_resident_set, projection.allocated_streaming_iq_buffers);
			projection.projected_total_footprint =
				addBytes(addBytes(addBytes(projection.phase_noise_lookup, projection.streaming_iq_buffers),
								  projection.rendered_hdf5_payload),
						 *projection.resident_baseline);
		}

		return projection;
	}

	std::string formatByteSize(const std::uint64_t bytes)
	{
		constexpr std::array units = {"B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB"};
		auto value = static_cast<long double>(bytes);
		std::size_t unit_index = 0;
		while (value >= 1024.0L && unit_index + 1 < units.size())
		{
			value /= 1024.0L;
			++unit_index;
		}
		if (unit_index == 0)
		{
			return std::format("{} {}", bytes, units[unit_index]);
		}
		return std::format("{:.2f} {}", static_cast<double>(value), units[unit_index]);
	}

	std::string memoryProjectionToJsonString(const SimulationMemoryProjection& projection)
	{
		const nlohmann::json result = {
			{"duration_seconds", projection.duration_seconds},
			{"simulation_sample_rate_hz", projection.simulation_sample_rate_hz},
			{"oversample_ratio", projection.oversample_ratio},
			{"sample_counts",
			 {{"streaming_samples", projection.streaming_sample_count},
			  {"phase_noise_samples_per_enabled_timing", projection.phase_noise_sample_count},
			  {"rendered_hdf5_samples", projection.rendered_hdf5_sample_count},
			  {"pulsed_windows", projection.pulsed_window_count}}},
			{"object_counts",
			 {{"phase_noise_timing_sources", projection.phase_noise_timing_count},
			  {"enabled_phase_noise_timing_sources", projection.enabled_phase_noise_timing_count},
			  {"streaming_receivers", projection.streaming_receiver_count},
			  {"pulsed_receivers", projection.pulsed_receiver_count}}},
			{"phase_noise_lookups", byteProjectionToJson(projection.phase_noise_lookup)},
			{"streaming_iq_buffers", byteProjectionToJson(projection.streaming_iq_buffers)},
			{"allocated_streaming_iq_buffers", byteProjectionToJson(projection.allocated_streaming_iq_buffers)},
			{"rendered_hdf5_dataset_payload", byteProjectionToJson(projection.rendered_hdf5_payload)},
			{"current_resident_set",
			 projection.current_resident_set.has_value()
				 ? nlohmann::json{{"bytes", *projection.current_resident_set},
								  {"human", formatByteSize(*projection.current_resident_set)}}
				 : nlohmann::json{{"bytes", nullptr}, {"human", "unknown"}}},
			{"resident_baseline", optionalByteProjectionToJson(projection.resident_baseline)},
			{"projected_total_footprint", optionalByteProjectionToJson(projection.projected_total_footprint)}};
		return result.dump(2);
	}

	void logSimulationMemoryProjection(const World& world)
	{
		const auto projection = projectSimulationMemory(world);
		const std::string resident_baseline =
			projection.resident_baseline.has_value() ? formatByteSize(projection.resident_baseline->bytes) : "unknown";
		const std::string total = projection.projected_total_footprint.has_value()
			? formatByteSize(projection.projected_total_footprint->bytes)
			: "unknown";

		LOG(logging::Level::DEBUG,
			"Projected simulation footprint: phase_noise_lookup_memory={} ({} enabled timing sources x {} samples), "
			"streaming_iq_buffer_memory={} ({} receivers x {} samples), rendered_hdf5_dataset_payload={} "
			"({} output samples), resident_baseline={} (current RSS minus allocated streaming IQ buffers), "
			"projected_total_footprint={}.",
			formatByteSize(projection.phase_noise_lookup.bytes), projection.enabled_phase_noise_timing_count,
			projection.phase_noise_sample_count, formatByteSize(projection.streaming_iq_buffers.bytes),
			projection.streaming_receiver_count, projection.streaming_sample_count,
			formatByteSize(projection.rendered_hdf5_payload.bytes), projection.rendered_hdf5_sample_count,
			resident_baseline, total);

		constexpr std::uint64_t one_gib = 1024ULL * 1024ULL * 1024ULL;
		for (const auto& receiver_ptr : world.getReceivers())
		{
			const auto& receiver = *receiver_ptr;
			if (!receiver.isDechirpEnabled())
			{
				continue;
			}
			const auto sample_count = receiver.getStreamingData().empty()
				? projection.streaming_sample_count
				: static_cast<std::uint64_t>(receiver.getStreamingData().size());
			const auto projected_payload =
				multiplyBytes(sample_count, 2ULL * static_cast<std::uint64_t>(sizeof(RealType)));
			if (projected_payload.bytes > one_gib || projected_payload.overflowed)
			{
				const auto gib = static_cast<long double>(projected_payload.bytes) / static_cast<long double>(one_gib);
				// TODO: UI workflows should have disk+memory usage stats on the SimulationView page (shows before user
				// runs the sim)
				LOG(logging::Level::WARNING,
					"Receiver {} is outputting RF-rate IF data. Projected file size is {:.2f} GiB. This is expected "
					"for V1 native dechirping, but ensure you have sufficient disk space. Future versions will support "
					"IF-rate decimation.",
					receiver.getName(), static_cast<double>(gib));
			}
		}
	}
}
