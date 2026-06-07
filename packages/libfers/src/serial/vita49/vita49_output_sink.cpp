// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "serial/vita49/vita49_output_sink.h"

#include <algorithm>
#include <chrono>
#include <optional>
#include <stdexcept>

#include "core/parameters.h"
#include "serial/vita49/udp_sender.h"

namespace serial::vita49
{
	namespace
	{
		constexpr auto kStatsEmitInterval = std::chrono::milliseconds(250);
		constexpr auto kPacketTraceEmitInterval = std::chrono::milliseconds(250);
		constexpr std::size_t kPacketTraceBatchSize = 64;

		[[nodiscard]] std::uint64_t defaultEpochNanoseconds()
		{
			const auto now = std::chrono::system_clock::now().time_since_epoch();
			return static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(now).count());
		}

		[[nodiscard]] core::Vita49Timestamp toCoreTimestamp(const Timestamp& timestamp) noexcept
		{
			return core::Vita49Timestamp{.integer_seconds = timestamp.integer_seconds,
										 .fractional_picoseconds = timestamp.fractional_picoseconds};
		}

		[[nodiscard]] std::optional<core::Vita49Timestamp>
		tryCoreTimestampFromEpoch(const std::uint64_t epoch_unix_nanoseconds,
								  const RealType sample_time_seconds) noexcept
		{
			try
			{
				return toCoreTimestamp(timestampFromEpoch(epoch_unix_nanoseconds, sample_time_seconds));
			}
			catch (...)
			{
				return std::nullopt;
			}
		}
	}

	Vita49OutputSink::Vita49OutputSink(std::unique_ptr<DatagramSender> sender,
									   core::ReceiverOutputTelemetryCallback telemetry_callback) :
		_telemetry_callback(std::move(telemetry_callback)), _provided_sender(std::move(sender))
	{
	}

	Vita49OutputSink::~Vita49OutputSink()
	{
		if (!_finalized)
		{
			try
			{
				(void)finalize();
			}
			catch (...)
			{
			}
		}
	}

	void Vita49OutputSink::initializeRun(const core::OutputConfig& config, std::string simulation_name)
	{
		std::scoped_lock const lock(_mutex);
		if (config.mode != core::OutputMode::Vita49Udp)
		{
			throw std::invalid_argument("Vita49OutputSink requires VITA output mode");
		}
		if (config.vita49.host.empty() || config.vita49.port == 0)
		{
			throw std::invalid_argument("VITA output requires destination host and port");
		}
		if (config.vita49.queue_depth == 0)
		{
			throw std::invalid_argument("VITA output queue depth must be positive");
		}

		_config = config;
		_simulation_name = std::move(simulation_name);
		const auto epoch_ns = config.vita49.epoch_unix_nanoseconds.value_or(defaultEpochNanoseconds());
		_packetizer =
			std::make_unique<Vita49Packetizer>(epoch_ns, config.vita49.adc_fullscale, config.vita49.max_udp_payload);
		_sender = std::make_unique<PacedSender>(
			_provided_sender ? std::move(_provided_sender) : std::make_unique<UdpSender>(), config.vita49.queue_depth);
		_sender->open(config.vita49.host, config.vita49.port);
		_sender->start(params::startTime());
		_last_stats_emit = std::chrono::steady_clock::time_point::min();
		_last_packet_trace_emit = std::chrono::steady_clock::now();
		_pending_packet_traces.clear();
		_trace_sequence = 0;
		_initialized = true;
		_finalized = false;
		emitTelemetry({}, true);
	}

	std::uint32_t Vita49OutputSink::registerStream(const core::ReceiverStreamDescriptor& stream)
	{
		std::scoped_lock const lock(_mutex);
		const auto stream_id = _registry.registerStream(stream);
		if (!_streams.contains(stream_id))
		{
			StreamState state;
			state.descriptor = stream;
			state.stats.receiver_id = stream.receiver_id;
			state.stats.receiver_name = stream.receiver_name;
			state.stats.stream_id = stream_id;
			state.stats.mode = stream.mode.empty() ? "unknown" : stream.mode;
			state.stats.sample_rate = stream.sample_rate;
			state.stats.reference_frequency = stream.reference_frequency;
			_streams.emplace(stream_id, std::move(state));
			emitTelemetry({}, true);
		}
		return stream_id;
	}

	void Vita49OutputSink::openStream(const std::uint32_t stream_id, const RealType first_sample_time)
	{
		std::scoped_lock const lock(_mutex);
		emitContext(stream_id, first_sample_time, true, false);
		stateFor(stream_id).opened = true;
		emitTelemetry({}, true);
	}

	void Vita49OutputSink::submitBlock(const core::ReceiverSampleBlock& block)
	{
		std::scoped_lock const lock(_mutex);
		if (!_initialized || !_packetizer)
		{
			throw std::logic_error("VITA output sink has not been initialized");
		}
		emitTelemetry(consumeSenderDropsLocked(), false);
		const auto stream_id = registerStream(block.stream);
		auto& state = stateFor(stream_id);
		if (!state.opened)
		{
			openStream(stream_id, block.first_sample_time);
		}

		const RealType sample_rate = block.sample_rate > 0.0 ? block.sample_rate : state.stats.sample_rate;
		auto result = _packetizer->packetize(block, stream_id, state.packet_counts, state.sample_loss_pending);
		state.sample_loss_pending = false;
		for (auto& packet : result.packets)
		{
			const auto packet_sample_count = packet.sample_count;
			const auto packet_first_sample_time = packet.first_sample_time;
			const auto packet_end_sample_time = sample_rate > 0.0
				? packet_first_sample_time + static_cast<RealType>(packet_sample_count) / sample_rate
				: packet_first_sample_time;
			const auto packet_over_range = packet.over_range;
			const auto packet_timestamp = toCoreTimestamp(packet.timestamp);
			const auto packet_end_timestamp =
				tryCoreTimestampFromEpoch(_packetizer->epochUnixNanoseconds(), packet_end_sample_time);
			if (!enqueuePacket(std::move(packet)))
			{
				continue;
			}
			state.stats.samples_emitted += packet_sample_count;
			++state.stats.packets_emitted;
			if (!state.stats.first_sample_time.has_value())
			{
				state.stats.first_sample_time = packet_first_sample_time;
				state.stats.first_timestamp = packet_timestamp;
			}
			state.stats.end_sample_time = packet_end_sample_time;
			state.stats.end_timestamp = packet_end_timestamp;
			if (packet_over_range)
			{
				state.over_range_pending = true;
			}
		}
		state.stats.over_range_count += result.over_range_count;
	}

	void Vita49OutputSink::emitContextHeartbeat(const RealType simulation_time)
	{
		std::scoped_lock const lock(_mutex);
		emitTelemetry(consumeSenderDropsLocked(), false);
		for (auto& [stream_id, state] : _streams)
		{
			if (!state.closed && simulation_time - state.last_context_time >= 1.0)
			{
				emitContext(stream_id, simulation_time, false, false);
			}
		}
	}

	void Vita49OutputSink::closeStream(const std::uint32_t stream_id)
	{
		std::scoped_lock const lock(_mutex);
		emitTelemetry(consumeSenderDropsLocked(), false);
		auto& state = stateFor(stream_id);
		if (!state.closed)
		{
			emitContext(stream_id, state.last_context_time, false, true);
			state.closed = true;
			emitTelemetry({}, true);
		}
	}

	core::OutputStats Vita49OutputSink::finalize()
	{
		std::scoped_lock const lock(_mutex);
		if (_finalized)
		{
			auto stats = snapshotStatsLocked();
			emitTelemetry({}, true);
			return stats;
		}

		if (_sender)
		{
			_sender->flush();
			emitTelemetry(consumeSenderDropsLocked(), false);
		}

		for (auto& [stream_id, state] : _streams)
		{
			if (!state.closed)
			{
				emitContext(stream_id, state.last_context_time, false, true);
				state.closed = true;
			}
		}

		if (_sender)
		{
			_sender->stop();
			emitTelemetry(consumeSenderDropsLocked(), false);
		}

		core::OutputStats stats{.mode = core::OutputMode::Vita49Udp,
								.epoch_unix_nanoseconds = _packetizer
									? std::optional<std::uint64_t>(_packetizer->epochUnixNanoseconds())
									: std::nullopt,
								.streams = {}};
		for (auto& [stream_id, state] : _streams)
		{
			if (_sender)
			{
				state.stats.late_packet_count = _sender->latePacketCount(stream_id);
			}
			stats.streams.push_back(state.stats);
		}
		_finalized = true;
		emitTelemetry({}, true);
		return stats;
	}

	core::OutputStats Vita49OutputSink::snapshotStats() const
	{
		std::scoped_lock const lock(_mutex);
		return snapshotStatsLocked();
	}

	Vita49OutputSink::StreamState& Vita49OutputSink::stateFor(const std::uint32_t stream_id)
	{
		const auto found = _streams.find(stream_id);
		if (found == _streams.end())
		{
			throw std::out_of_range("Unknown VITA stream ID");
		}
		return found->second;
	}

	const Vita49OutputSink::StreamState& Vita49OutputSink::stateFor(const std::uint32_t stream_id) const
	{
		const auto found = _streams.find(stream_id);
		if (found == _streams.end())
		{
			throw std::out_of_range("Unknown VITA stream ID");
		}
		return found->second;
	}

	core::OutputStats Vita49OutputSink::snapshotStatsLocked() const
	{
		core::OutputStats stats{.mode = core::OutputMode::Vita49Udp,
								.epoch_unix_nanoseconds = _packetizer
									? std::optional<std::uint64_t>(_packetizer->epochUnixNanoseconds())
									: std::nullopt,
								.streams = {}};
		for (const auto& [stream_id, state] : _streams)
		{
			auto stream_stats = state.stats;
			if (_sender)
			{
				stream_stats.late_packet_count = _sender->latePacketCount(stream_id);
			}
			stats.streams.push_back(std::move(stream_stats));
		}
		return stats;
	}

	std::vector<core::ReceiverOutputPacketTrace> Vita49OutputSink::consumeSenderDropsLocked()
	{
		std::vector<core::ReceiverOutputPacketTrace> traces;
		if (!_sender)
		{
			return traces;
		}

		for (const auto& dropped : _sender->consumeDroppedDatagrams())
		{
			applyDropped(dropped);
			if (_config.vita49.packet_trace_enabled)
			{
				traces.push_back(makeDropTrace(dropped));
			}
		}
		return traces;
	}

	bool Vita49OutputSink::enqueuePacket(SerializedPacket&& packet)
	{
		if (!_sender)
		{
			throw std::logic_error("VITA paced sender is unavailable");
		}

		std::vector<core::ReceiverOutputPacketTrace> traces;
		std::optional<core::ReceiverOutputPacketTrace> sent_trace;
		if (_config.vita49.packet_trace_enabled)
		{
			sent_trace = makeTrace(packet, packet.context_packet ? "context" : "data");
		}
		const auto result = _sender->enqueue(std::move(packet));
		if (_config.vita49.packet_trace_enabled && result.dropped)
		{
			traces.push_back(makeDropTrace(*result.dropped));
		}
		if (result.dropped)
		{
			applyDropped(*result.dropped);
		}
		if (result.enqueued && sent_trace.has_value())
		{
			traces.push_back(std::move(*sent_trace));
		}
		emitTelemetry(std::move(traces), false);
		return result.enqueued;
	}

	void Vita49OutputSink::emitTelemetry(std::vector<core::ReceiverOutputPacketTrace> packets, const bool force_stats)
	{
		if (!_telemetry_callback)
		{
			return;
		}

		std::optional<core::OutputStats> stats;
		const auto now = std::chrono::steady_clock::now();
		if (force_stats || _last_stats_emit == std::chrono::steady_clock::time_point::min() ||
			now - _last_stats_emit >= kStatsEmitInterval)
		{
			stats = snapshotStatsLocked();
			_last_stats_emit = now;
		}

		for (auto& packet : packets)
		{
			packet.sequence = ++_trace_sequence;
			_pending_packet_traces.push_back(std::move(packet));
		}

		std::vector<core::ReceiverOutputPacketTrace> packet_batch;
		const bool trace_interval_elapsed = _last_packet_trace_emit == std::chrono::steady_clock::time_point::min() ||
			now - _last_packet_trace_emit >= kPacketTraceEmitInterval;
		if (!_pending_packet_traces.empty() &&
			(force_stats || _pending_packet_traces.size() >= kPacketTraceBatchSize || trace_interval_elapsed))
		{
			packet_batch = std::move(_pending_packet_traces);
			_pending_packet_traces.clear();
			_last_packet_trace_emit = now;
		}

		if (!stats.has_value() && packet_batch.empty())
		{
			return;
		}
		_telemetry_callback(stats, packet_batch);
	}

	core::ReceiverOutputPacketTrace Vita49OutputSink::makeTrace(const SerializedPacket& packet, std::string event) const
	{
		return core::ReceiverOutputPacketTrace{.sequence = 0,
											   .event = std::move(event),
											   .stream_id = packet.stream_id,
											   .byte_count = packet.bytes.size(),
											   .sample_count = packet.sample_count,
											   .first_sample_time = packet.first_sample_time,
											   .timestamp = toCoreTimestamp(packet.timestamp),
											   .data_packet = packet.data_packet,
											   .context_packet = packet.context_packet,
											   .dropped = false,
											   .over_range = packet.over_range,
											   .sample_loss = packet.sample_loss};
	}

	core::ReceiverOutputPacketTrace Vita49OutputSink::makeDropTrace(const DroppedDatagram& dropped) const
	{
		return core::ReceiverOutputPacketTrace{.sequence = 0,
											   .event = "drop",
											   .stream_id = dropped.stream_id,
											   .byte_count = 0,
											   .sample_count = dropped.sample_count,
											   .first_sample_time = 0.0,
											   .timestamp = std::nullopt,
											   .data_packet = dropped.data_packet,
											   .context_packet = dropped.context_packet,
											   .dropped = true,
											   .over_range = false,
											   .sample_loss = true};
	}

	void Vita49OutputSink::emitContext(const std::uint32_t stream_id, const RealType simulation_time,
									   const bool stream_open, const bool stream_close)
	{
		if (!_packetizer)
		{
			throw std::logic_error("VITA packetizer is unavailable");
		}
		auto& state = stateFor(stream_id);
		const RealType context_time = simulation_time <= -1.0e200 ? 0.0 : simulation_time;
		const auto timestamp = timestampFromEpoch(_packetizer->epochUnixNanoseconds(), context_time);
		const ContextBuildRequest request{.stream = state.descriptor,
										  .stream_id = stream_id,
										  .simulation_name = _simulation_name,
										  .adc_fullscale = _packetizer->adcFullscale(),
										  .timestamp = timestamp,
										  .packet_count = state.packet_counts.next(),
										  .valid_data = true,
										  .calibrated_time = true,
										  .reference_lock = true,
										  .over_range = state.over_range_pending,
										  .sample_loss = state.sample_loss_pending,
										  .stream_open = stream_open,
										  .stream_close = stream_close};
		const auto context = Vita49ContextBuilder::build(request);
		auto packet = _packetizer->makeContextPacket(context);
		packet.first_sample_time = context_time;
		if (enqueuePacket(std::move(packet)))
		{
			++state.stats.context_packets;
		}
		state.last_context_time = context_time;
		state.sample_loss_pending = false;
		state.over_range_pending = false;
	}

	void Vita49OutputSink::applyDropped(const DroppedDatagram& dropped)
	{
		if (dropped.stream_id == 0 || !_streams.contains(dropped.stream_id))
		{
			return;
		}
		auto& state = stateFor(dropped.stream_id);
		++state.stats.packets_dropped;
		state.stats.samples_dropped += dropped.sample_count;
		if (dropped.data_packet || (!dropped.context_packet && dropped.sample_count > 0))
		{
			state.stats.packets_emitted -= std::min<std::uint64_t>(state.stats.packets_emitted, 1u);
			state.stats.samples_emitted -= std::min(state.stats.samples_emitted, dropped.sample_count);
		}
		else if (dropped.context_packet)
		{
			state.stats.context_packets -= std::min<std::uint64_t>(state.stats.context_packets, 1u);
		}
		state.sample_loss_pending = true;
	}

	std::unique_ptr<core::ReceiverOutputSink>
	makeVita49OutputSink(core::ReceiverOutputTelemetryCallback telemetry_callback)
	{
		return std::make_unique<Vita49OutputSink>(nullptr, std::move(telemetry_callback));
	}
}
