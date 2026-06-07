// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>

namespace serial::vita49
{
	class DatagramSender
	{
	public:
		virtual ~DatagramSender() = default;
		virtual void open(const std::string& host, std::uint16_t port) = 0;
		virtual void send(std::span<const std::uint8_t> bytes) = 0;
		virtual void close() noexcept = 0;
	};

	class UdpSender final : public DatagramSender
	{
	public:
		UdpSender() = default;
		~UdpSender() override;

		UdpSender(const UdpSender&) = delete;
		UdpSender& operator=(const UdpSender&) = delete;
		UdpSender(UdpSender&& other) noexcept;
		UdpSender& operator=(UdpSender&& other) noexcept;

		void open(const std::string& host, std::uint16_t port) override;
		void send(std::span<const std::uint8_t> bytes) override;
		void close() noexcept override;

	private:
		int _socket = -1;
		std::unique_ptr<std::byte[]> _address;
		std::size_t _address_size = 0;
	};
}
