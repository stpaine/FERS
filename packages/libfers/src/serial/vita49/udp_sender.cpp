// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2026-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

#include "serial/vita49/udp_sender.h"

#include <bit>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>

#ifdef _WIN32
#include <winsock2.h>
#include <ws2tcpip.h>
#else
#include <netdb.h>
#include <sys/socket.h>
#include <unistd.h>
#endif

namespace serial::vita49
{
	namespace
	{
		constexpr int kRequestedUdpSendBufferBytes = 4 * 1024 * 1024;

		void tuneSendBuffer(const int socket_fd) noexcept
		{
			const int send_buffer_bytes = kRequestedUdpSendBufferBytes;
#ifdef _WIN32
			const auto* option_value = std::bit_cast<const char*>(&send_buffer_bytes);
#else
			const void* option_value = &send_buffer_bytes;
#endif
			(void)::setsockopt(socket_fd, SOL_SOCKET, SO_SNDBUF, option_value, sizeof(send_buffer_bytes));
		}
	}

	UdpSender::~UdpSender() { close(); }

	UdpSender::UdpSender(UdpSender&& other) noexcept :
		_socket(other._socket), _address(std::move(other._address)), _address_size(other._address_size)
	{
		other._socket = -1;
		other._address_size = 0;
	}

	UdpSender& UdpSender::operator=(UdpSender&& other) noexcept
	{
		if (this != &other)
		{
			close();
			_socket = other._socket;
			_address = std::move(other._address);
			_address_size = other._address_size;
			other._socket = -1;
			other._address_size = 0;
		}
		return *this;
	}

	void UdpSender::open(const std::string& host, const std::uint16_t port)
	{
		close();
		if (host.empty() || port == 0)
		{
			throw std::invalid_argument("VITA UDP destination requires non-empty host and non-zero port");
		}

#ifdef _WIN32
		WSADATA wsa_data;
		if (WSAStartup(MAKEWORD(2, 2), &wsa_data) != 0)
		{
			throw std::runtime_error("WSAStartup failed for VITA UDP sender");
		}
#endif

		addrinfo hints{};
		hints.ai_family = AF_UNSPEC;
		hints.ai_socktype = SOCK_DGRAM;
		hints.ai_protocol = IPPROTO_UDP;

		addrinfo* result = nullptr;
		const auto port_string = std::to_string(port);
		const int gai = getaddrinfo(host.c_str(), port_string.c_str(), &hints, &result);
		if (gai != 0)
		{
			throw std::runtime_error("VITA UDP destination resolution failed: " + host + ":" + port_string);
		}
		std::unique_ptr<addrinfo, decltype(&freeaddrinfo)> const result_guard(result, freeaddrinfo);

		for (auto* candidate = result; candidate != nullptr; candidate = candidate->ai_next)
		{
#ifdef _WIN32
			const int fd =
				static_cast<int>(::socket(candidate->ai_family, candidate->ai_socktype, candidate->ai_protocol));
#else
			const int fd = ::socket(candidate->ai_family, candidate->ai_socktype, candidate->ai_protocol);
#endif
			if (fd < 0)
			{
				continue;
			}

			tuneSendBuffer(fd);
			_socket = fd;
			_address_size = candidate->ai_addrlen;
			auto storage = std::make_unique<std::byte[]>(_address_size);
			std::memcpy(storage.get(), candidate->ai_addr, _address_size);
			_address = std::move(storage);
			return;
		}

		throw std::runtime_error("VITA UDP socket creation failed");
	}

	void UdpSender::send(const std::span<const std::uint8_t> bytes)
	{
		if (_socket < 0 || _address == nullptr)
		{
			throw std::runtime_error("VITA UDP sender is not open");
		}
		if (bytes.empty())
		{
			return;
		}

#ifdef _WIN32
		const auto* payload = std::bit_cast<const char*>(bytes.data());
		const auto sent = ::sendto(_socket, payload, static_cast<int>(bytes.size()), 0,
								   static_cast<const sockaddr*>(static_cast<const void*>(_address.get())),
								   static_cast<int>(_address_size));
#else
		const auto sent = ::sendto(_socket, bytes.data(), bytes.size(), 0,
								   static_cast<const sockaddr*>(static_cast<const void*>(_address.get())),
								   static_cast<socklen_t>(_address_size));
#endif
		if (sent < 0 || static_cast<std::size_t>(sent) != bytes.size())
		{
			throw std::runtime_error("VITA UDP send failed");
		}
	}

	void UdpSender::close() noexcept
	{
		if (_socket >= 0)
		{
#ifdef _WIN32
			::closesocket(_socket);
			WSACleanup();
#else
			::close(_socket);
#endif
			_socket = -1;
		}
		_address.reset();
		_address_size = 0;
	}
}
