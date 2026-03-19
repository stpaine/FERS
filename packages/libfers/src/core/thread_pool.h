// SPDX-License-Identifier: GPL-2.0-only
//
// Copyright (c) 2024-present FERS Contributors (see AUTHORS.md).
//
// See the GNU GPLv2 LICENSE file in the FERS project root for more information.

/**
 * @file thread_pool.h
 * @brief A simple thread pool implementation.
 */

#pragma once

#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

namespace pool
{
	/**
	 * @class ThreadPool
	 * @brief A simple thread pool implementation.
	 */
	class ThreadPool
	{
	public:
		/**
		 * @brief Constructs a ThreadPool with a specified number of threads.
		 * @param numThreads The number of threads in the pool.
		 */
		explicit ThreadPool(unsigned numThreads);

		/**
		 * @brief Destroys the ThreadPool, joining all threads.
		 */
		~ThreadPool();

		ThreadPool(const ThreadPool&) = delete;
		ThreadPool& operator=(const ThreadPool&) = delete;
		ThreadPool(ThreadPool&&) = delete;
		ThreadPool& operator=(ThreadPool&&) = delete;

		/**
		 * @brief Enqueues a task to be executed by the thread pool.
		 * @tparam F The type of the function to be executed.
		 * @tparam Args The types of the arguments to the function.
		 * @param f The function to be executed.
		 * @param args The arguments to the function.
		 * @return A future that will hold the result of the function.
		 * @throws std::runtime_error if the thread pool is stopped.
		 */
		template <class F, class... Args>
		std::future<std::invoke_result_t<F, Args...>> enqueue(F&& f, Args&&... args)
		{
			using ReturnType = std::invoke_result_t<F, Args...>;

			auto task = std::make_shared<std::packaged_task<ReturnType()>>(
				std::bind(std::forward<F>(f), std::forward<Args>(args)...));

			std::future<ReturnType> res = task->get_future();
			{
				std::unique_lock lock(_queue_mutex);
				if (_stop)
				{
					throw std::runtime_error("enqueue on stopped ThreadPool");
				}

				++_pending_tasks;
				_tasks.emplace([task] { (*task)(); });
				_condition.notify_one();
			}
			return res;
		}

		/**
		 * @brief Waits for all tasks in the thread pool to finish.
		 */
		void wait()
		{
			std::unique_lock lock(_queue_mutex);
			_done_condition.wait(lock, [this] { return _pending_tasks == 0; });
		}

		/**
		 * @brief Returns the number of threads available for executing tasks.
		 * @return The number of available threads.
		 */
		[[nodiscard]] unsigned getAvailableThreads();

	private:
		using Task = std::function<void()>;

		std::vector<std::thread> _workers; ///< Vector of worker threads.
		std::queue<Task> _tasks; ///< Queue of tasks to be executed.
		std::mutex _queue_mutex; ///< Mutex for synchronizing access to the task queue.
		std::condition_variable _condition; ///< Condition variable for task notification.
		std::condition_variable _done_condition; ///< Condition variable for task completion notification.
		std::atomic<bool> _stop = false; ///< Flag indicating whether the thread pool is stopped.
		std::atomic<unsigned> _pending_tasks = 0; ///< Count of pending tasks.
	};
}
