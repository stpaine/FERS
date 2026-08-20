// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#pragma once

#include <chrono>
#include <memory>
#include <mutex>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "core/receiver_output.h"
#include "serial/vita49/paced_sender.h"
#include "serial/vita49/stream_registry.h"
#include "serial/vita49/vita49_context_builder.h"
#include "serial/vita49/vita49_packetizer.h"

namespace serial::vita49
{
	class Vita49OutputSink final : public core::ReceiverOutputSink
	{
	public:
		explicit Vita49OutputSink(std::unique_ptr<DatagramSender> sender = nullptr,
								  core::ReceiverOutputTelemetryCallback telemetry_callback = nullptr);
		~Vita49OutputSink() override;
		Vita49OutputSink(const Vita49OutputSink&) = delete;
		Vita49OutputSink& operator=(const Vita49OutputSink&) = delete;
		Vita49OutputSink(Vita49OutputSink&&) = delete;
		Vita49OutputSink& operator=(Vita49OutputSink&&) = delete;

		void initializeRun(const core::OutputConfig& config, std::string simulation_name) override;
		std::uint32_t registerStream(const core::ReceiverStreamDescriptor& stream) override;
		void openStream(std::uint32_t stream_id, RealType first_sample_time) override;
		void submitBlock(const core::ReceiverSampleBlock& block) override;
		void submitBlocks(std::span<const core::ReceiverSampleBlock> blocks) override;
		void emitContextHeartbeat(RealType simulation_time) override;
		void closeStream(std::uint32_t stream_id) override;
		core::OutputStats finalize() override;
		[[nodiscard]] core::OutputStats snapshotStats() const override;

	private:
		struct StreamState
		{
			core::ReceiverStreamDescriptor descriptor;
			core::ReceiverStreamStats stats;
			PacketCountSequencer packet_counts;
			bool opened = false;
			bool closed = false;
			bool sample_loss_pending = false;
			bool over_range_pending = false;
			RealType last_context_time = -1.0e300;
		};

		struct PendingContext
		{
			std::uint32_t stream_id = 0;
			RealType simulation_time = 0.0;
			bool stream_open = false;
			bool stream_close = false;
		};

		[[nodiscard]] StreamState& stateFor(std::uint32_t stream_id);
		[[nodiscard]] const StreamState& stateFor(std::uint32_t stream_id) const;
		[[nodiscard]] core::OutputStats snapshotStatsLocked() const;
		void ensurePacketizer();
		void startPacing();
		void appendPendingContexts(std::vector<SerializedPacket>& packets,
								   std::vector<std::uint32_t>& context_stream_ids);
		[[nodiscard]] bool enqueuePacket(SerializedPacket&& packet);
		[[nodiscard]] bool enqueuePackets(std::vector<SerializedPacket> packets);
		void emitTelemetry(std::vector<core::ReceiverOutputPacketTrace> packets = {}, bool force_stats = false);
		[[nodiscard]] std::vector<core::ReceiverOutputPacketTrace> consumeSenderDropsLocked();
		[[nodiscard]] core::ReceiverOutputPacketTrace makeTrace(const SerializedPacket& packet,
																std::string event) const;
		[[nodiscard]] core::ReceiverOutputPacketTrace makeDropTrace(const DroppedDatagram& dropped) const;
		[[nodiscard]] SerializedPacket buildContextPacket(std::uint32_t stream_id, RealType simulation_time,
														  bool stream_open, bool stream_close);
		void emitContext(std::uint32_t stream_id, RealType simulation_time, bool stream_open, bool stream_close);
		void applyDropped(const DroppedDatagram& dropped);

		core::OutputConfig _config;
		std::string _simulation_name;
		core::ReceiverOutputTelemetryCallback _telemetry_callback;
		std::unique_ptr<DatagramSender> _provided_sender;
		StreamRegistry _registry;
		std::unique_ptr<Vita49Packetizer> _packetizer;
		std::unique_ptr<PacedSender> _sender;
		std::unordered_map<std::uint32_t, StreamState> _streams;
		std::vector<PendingContext> _pending_contexts;
		mutable std::recursive_mutex _mutex;
		std::chrono::steady_clock::time_point _last_stats_emit = std::chrono::steady_clock::time_point::min();
		std::chrono::steady_clock::time_point _last_packet_trace_emit = std::chrono::steady_clock::time_point::min();
		std::vector<core::ReceiverOutputPacketTrace> _pending_packet_traces;
		std::uint64_t _trace_sequence = 0;
		bool _initialized = false;
		bool _finalized = false;
		bool _pacing_started = false;
	};

	[[nodiscard]] std::unique_ptr<core::ReceiverOutputSink>
	makeVita49OutputSink(core::ReceiverOutputTelemetryCallback telemetry_callback = nullptr);
}
