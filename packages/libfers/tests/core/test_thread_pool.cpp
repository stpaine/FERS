#include <atomic>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <future>
#include <stdexcept>
#include <thread>

#include "core/thread_pool.h"

using namespace std::chrono_literals;

TEST_CASE("ThreadPool executes tasks and returns results", "[core][thread_pool]")
{
	pool::ThreadPool pool(2);

	std::atomic<int> sum{0};
	auto future_value = pool.enqueue(
		[&sum]
		{
			sum.fetch_add(2);
			return 7;
		});
	auto future_void = pool.enqueue([&sum] { sum.fetch_add(3); });

	REQUIRE(future_value.get() == 7);
	future_void.get();
	pool.wait();

	REQUIRE(sum.load() == 5);
}

TEST_CASE("ThreadPool wait tracks pending tasks", "[core][thread_pool]")
{
	pool::ThreadPool pool(2);

	std::promise<void> gate;
	std::shared_future<void> gate_future = gate.get_future().share();
	std::atomic<int> started{0};

	auto task = [&started, gate_future]
	{
		started.fetch_add(1);
		gate_future.wait();
	};

	auto future_a = pool.enqueue(task);
	auto future_b = pool.enqueue(task);

	for (int i = 0; i < 100 && started.load() < 2; ++i)
	{
		std::this_thread::sleep_for(1ms);
	}

	REQUIRE(started.load() == 2);
	REQUIRE(pool.getAvailableThreads() == 0);

	gate.set_value();
	future_a.get();
	future_b.get();
	pool.wait();

	REQUIRE(pool.getAvailableThreads() == 2);
}

TEST_CASE("ThreadPool propagates exceptions through futures", "[core][thread_pool]")
{
	pool::ThreadPool pool(1);

	auto future = pool.enqueue([] { throw std::runtime_error("boom"); });

	REQUIRE_THROWS_AS(future.get(), std::runtime_error);
	pool.wait();
}

TEST_CASE("ThreadPool handles zero-thread pools", "[core][thread_pool]")
{
	pool::ThreadPool pool(0);

	REQUIRE(pool.getAvailableThreads() == 0);
	pool.wait();
}

TEST_CASE("ThreadPool exception logging not directly testable", "[core][thread_pool]")
{
	// TODO: The worker try/catch path is not reliably triggered because packaged_task stores exceptions.
	REQUIRE(true);
}
