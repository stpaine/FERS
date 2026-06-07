// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "serial/vita49/paced_sender.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <stdexcept>

#if defined(__i386__) || defined(__x86_64__)
#include <immintrin.h>
#endif

namespace serial::vita49
{
	namespace
	{
		constexpr auto kCoarsePacingSleep = std::chrono::milliseconds(1);
		constexpr auto kFinePacingSpin = std::chrono::microseconds(200);

		void cpuPause() noexcept
		{
#if defined(__i386__) || defined(__x86_64__)
			_mm_pause();
#elif defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
			__asm__ __volatile__("yield" ::: "memory");
#else
			std::atomic_signal_fence(std::memory_order_seq_cst);
#endif
		}
	}

	PacedSender::PacedSender(std::unique_ptr<DatagramSender> sender, const std::size_t queue_depth) :
		_sender(std::move(sender)), _queue_depth(queue_depth)
	{
		if (!_sender)
		{
			throw std::invalid_argument("PacedSender requires a datagram sender");
		}
		if (queue_depth == 0)
		{
			throw std::invalid_argument("PacedSender queue depth must be positive");
		}
	}

	PacedSender::~PacedSender() { stop(); }

	void PacedSender::open(const std::string& host, const std::uint16_t port) { _sender->open(host, port); }

	void PacedSender::start(const RealType simulation_epoch_time)
	{
		std::scoped_lock const lock(_mutex);
		if (_started)
		{
			return;
		}
		_simulation_epoch_time = simulation_epoch_time;
		_steady_epoch = std::chrono::steady_clock::now();
		_stopping = false;
		_started = true;
		_thread = std::thread([this] { run(); });
	}

	EnqueueResult PacedSender::enqueue(SerializedPacket packet)
	{
		std::unique_lock lock(_mutex);
		if (!_started)
		{
			throw std::logic_error("PacedSender must be started before enqueue");
		}
		if (_stopping)
		{
			return EnqueueResult{
				.enqueued = false,
				.dropped = std::nullopt,
			};
		}

		_cv.wait(lock, [this] { return _stopping || queuedOrSendingCount() < _queue_depth; });
		if (_stopping)
		{
			return EnqueueResult{
				.enqueued = false,
				.dropped = std::nullopt,
			};
		}

		_queue.push_back(std::move(packet));
		_cv.notify_one();
		return EnqueueResult{.enqueued = true, .dropped = std::nullopt};
	}

	void PacedSender::flush()
	{
		std::unique_lock lock(_mutex);
		_cv.wait(lock, [this] { return _queue.empty() && !_send_in_progress; });
	}

	void PacedSender::stop()
	{
		{
			std::scoped_lock const lock(_mutex);
			if (!_started && !_thread.joinable())
			{
				_sender->close();
				return;
			}
			_stopping = true;
			_cv.notify_all();
		}

		if (_thread.joinable())
		{
			_thread.join();
		}

		{
			std::scoped_lock const lock(_mutex);
			_started = false;
			_stopping = false;
			_send_in_progress = false;
			_cv.notify_all();
		}
		_sender->close();
	}

	std::uint64_t PacedSender::latePacketCount(const std::uint32_t stream_id) const
	{
		std::scoped_lock const lock(_mutex);
		const auto found = _late_packets.find(stream_id);
		return found == _late_packets.end() ? 0 : found->second;
	}

	std::uint64_t PacedSender::sentPacketCount(const std::uint32_t stream_id) const
	{
		std::scoped_lock const lock(_mutex);
		const auto found = _sent_packets.find(stream_id);
		return found == _sent_packets.end() ? 0 : found->second;
	}

	std::uint64_t PacedSender::sendFailureCount(const std::uint32_t stream_id) const
	{
		std::scoped_lock const lock(_mutex);
		const auto found = _send_failures.find(stream_id);
		return found == _send_failures.end() ? 0 : found->second;
	}

	std::uint64_t PacedSender::droppedDataPacketCount(const std::uint32_t stream_id) const
	{
		std::scoped_lock const lock(_mutex);
		const auto found = _dropped_data_packets.find(stream_id);
		return found == _dropped_data_packets.end() ? 0 : found->second;
	}

	std::uint64_t PacedSender::droppedContextPacketCount(const std::uint32_t stream_id) const
	{
		std::scoped_lock const lock(_mutex);
		const auto found = _dropped_context_packets.find(stream_id);
		return found == _dropped_context_packets.end() ? 0 : found->second;
	}

	std::uint64_t PacedSender::droppedSampleCount(const std::uint32_t stream_id) const
	{
		std::scoped_lock const lock(_mutex);
		const auto found = _dropped_samples.find(stream_id);
		return found == _dropped_samples.end() ? 0 : found->second;
	}

	std::vector<DroppedDatagram> PacedSender::consumeDroppedDatagrams()
	{
		std::scoped_lock const lock(_mutex);
		auto result = std::move(_pending_dropped_datagrams);
		_pending_dropped_datagrams.clear();
		return result;
	}

	void PacedSender::run()
	{
		std::unique_lock lock(_mutex);
		while (true)
		{
			if (_queue.empty())
			{
				if (_stopping)
				{
					return;
				}
				_cv.wait(lock, [this] { return _stopping || !_queue.empty(); });
				continue;
			}

			const auto due = dueTime(_queue.front());
			waitUntilDue(lock, due);

			auto packet = std::move(_queue.front());
			_queue.pop_front();
			_send_in_progress = true;
			const auto now = std::chrono::steady_clock::now();
			lock.unlock();
			sendOneUnlocked(std::move(packet), now);
			lock.lock();
			_send_in_progress = false;
			_cv.notify_all();
		}
	}

	void PacedSender::waitUntilDue(std::unique_lock<std::mutex>& lock, const std::chrono::steady_clock::time_point due)
	{
		while (true)
		{
			const auto now = std::chrono::steady_clock::now();
			if (now >= due)
			{
				return;
			}

			const auto remaining = due - now;
			if (remaining > kCoarsePacingSleep)
			{
				const auto coarse_sleep =
					std::chrono::duration_cast<std::chrono::steady_clock::duration>(kCoarsePacingSleep);
				const auto fine_spin = std::chrono::duration_cast<std::chrono::steady_clock::duration>(kFinePacingSpin);
				const auto wait_duration = std::min(remaining - fine_spin, coarse_sleep);
				_cv.wait_for(lock, wait_duration);
				continue;
			}

			lock.unlock();
			while (std::chrono::steady_clock::now() < due)
			{
				cpuPause();
			}
			lock.lock();
		}
	}

	void PacedSender::sendOneUnlocked(SerializedPacket packet, const std::chrono::steady_clock::time_point now)
	{
		const auto due = dueTime(packet);
		try
		{
			_sender->send(packet.bytes);
		}
		catch (...)
		{
			std::scoped_lock const lock(_mutex);
			++_send_failures[packet.stream_id];
			recordDroppedUnlocked(packet);
			_pending_dropped_datagrams.push_back(makeDroppedDatagram(packet));
			return;
		}
		std::scoped_lock const lock(_mutex);
		++_sent_packets[packet.stream_id];
		if (now > due + std::chrono::milliseconds(1))
		{
			++_late_packets[packet.stream_id];
		}
	}

	void PacedSender::recordDroppedUnlocked(const SerializedPacket& packet)
	{
		if (packet.data_packet || (!packet.context_packet && packet.sample_count > 0))
		{
			++_dropped_data_packets[packet.stream_id];
			_dropped_samples[packet.stream_id] += packet.sample_count;
			return;
		}
		if (packet.context_packet)
		{
			++_dropped_context_packets[packet.stream_id];
		}
	}

	DroppedDatagram PacedSender::makeDroppedDatagram(const SerializedPacket& packet) const noexcept
	{
		return DroppedDatagram{.stream_id = packet.stream_id,
							   .sample_count = packet.sample_count,
							   .data_packet = packet.data_packet,
							   .context_packet = packet.context_packet};
	}

	std::size_t PacedSender::queuedOrSendingCount() const noexcept
	{
		return _queue.size() + (_send_in_progress ? 1u : 0u);
	}

	std::chrono::steady_clock::time_point PacedSender::dueTime(const SerializedPacket& packet) const
	{
		const auto seconds = packet.first_sample_time - _simulation_epoch_time;
		const auto nanos = static_cast<std::int64_t>(seconds * 1'000'000'000.0);
		return _steady_epoch + std::chrono::nanoseconds(nanos);
	}

}
