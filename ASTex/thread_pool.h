/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/



/*
 * IMPORTANT : The ThreadPool code (thread_pool.h and thread_pool.cpp) is
 * based on "A Simple c++11 threadpool implementation" found on github
 * (https://github.com/progschj/ThreadPool, latest commit : 9a42ec1 )
 * (c) 2012 Jakob Progsch, Václav Zeman
 * It has been modified to fit to our purposes.
 * A copy of its license is provided in the following lines.
 */

/****************************************************************************
*Copyright (c) 2012 Jakob Progsch, Václav Zeman                             *
*                                                                           *
*This software is provided 'as-is', without any express or implied          *
*warranty. In no event will the authors be held liable for any damages      *
*arising from the use of this software.                                     *
*                                                                           *
*Permission is granted to anyone to use this software for any purpose,      *
*including commercial applications, and to alter it and redistribute it     *
*freely, subject to the following restrictions:                             *
*                                                                           *
*1. The origin of this software must not be misrepresented; you must not    *
*claim that you wrote the original software. If you use this software       *
*in a product, an acknowledgment in the product documentation would be      *
*appreciated but is not required.                                           *
*                                                                           *
*2. Altered source versions must be plainly marked as such, and must not be *
*misrepresented as being the original software.                             *
*                                                                           *
*3. This notice may not be removed or altered from any source               *
*distribution.                                                              *
****************************************************************************/

#ifndef _UTILS_THREADPOOL_H_
#define _UTILS_THREADPOOL_H_

#include <ASTex/dll.h>
#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <cstdint>
#include <iostream>

#if defined(_MSC_VER)
#define ASTEX_TLS __declspec( thread )
#else
#define ASTEX_TLS __thread
#endif


//#ifndef _UTILS_THREADPOOL_CPP_
extern ASTEX_TLS uint32_t thread_index_;
//#endif

class ASTEX_API ThreadPool
{
	// no copy
	ThreadPool(const ThreadPool&) = delete;
	// no move
	ThreadPool(ThreadPool&&) = delete;
	// no affectation
	ThreadPool operator=(const ThreadPool&) = delete;
	// no move affectation
	ThreadPool operator=(ThreadPool&&) = delete;
public:

	using Future = std::future<void>;

	
	ThreadPool(uint32_t nb_ww);
	
	~ThreadPool();

#if defined(_MSC_VER) && _MSC_VER < 1900
	using PackagedTask = std::shared_ptr<std::packaged_task<void()>>; // avoiding a MSVC 2013 Bug
#else
	using PackagedTask = std::packaged_task<void()>;
#endif

	template <class F, class... Args>
	Future enqueue(const F& f, Args&&... args);


	/**
	 * @brief get the number of currently working thread for parallel algos
	 */
	inline uint32_t nb_workers() const
	{
		return nb_working_workers_;
	}

	/**
	* @brief get the number of threads that could be used for parallel algos
	*/
	inline uint32_t max_nb_workers() const
	{
		return uint32_t(workers_.size());
	}

	/**
	 * @brief set nb working threads for parallel algos ( no param = full power)
	 * @param nb [0,nb_max_workers()] (for 0 parallel algo are replaced by normal version)
	 */
	void set_nb_workers(uint32_t nb = 0xffffffff);

private:
#pragma warning(push)
#pragma warning(disable:4251)
	// need to keep track of threads so we can join them
	std::vector<std::thread> workers_;
	// the task queue
	std::queue<PackagedTask> tasks_;

	// synchronization for task queuing
	std::mutex queue_mutex_;
	std::condition_variable condition_;
	// is the thread
	bool stop_;

	// limit usage to the n-th first workers
	uint32_t nb_working_workers_;
	std::mutex running_mutex_;
	std::condition_variable condition_running_;
#pragma warning(pop)
};




// add new work item to the pool


template <class F, class... Args>
ThreadPool::Future ThreadPool::enqueue(const F& f, Args&&... args)
{
#if defined(_MSC_VER) && _MSC_VER < 1900
	PackagedTask task = std::make_shared<std::packaged_task<void()>>(std::bind(f, std::forward<Args>(args)...));
	std::future<void> res = task->get_future();
#else

	PackagedTask task([&, f]() -> void
	{
		return f(std::forward<Args>(args)...);
	});
	std::future<void> res = task.get_future();

#endif

	{
		std::unique_lock<std::mutex> lock(queue_mutex_);
		// don't allow enqueueing after stopping the pool
		if (stop_)
		{
			std::cerr << "Enqueue on stopped ThreadPool."<< std::endl;
			return res;
		}
		// Push work back on the queue
		tasks_.push(std::move(task));
	}
	// Notify a thread that there is new work to perform
	condition_.notify_one();
	return res;
}






/**
 * @brief thread index [0..nb_workers] for use in code of lambdas
 */
ASTEX_API uint32_t current_thread_index();

ASTEX_API ThreadPool* internal_thread_pool();

ASTEX_API ThreadPool* external_thread_pool();

/**
 * launch an external thread
 */
template <class F, class... Args>
ThreadPool::Future launch_thread(const F& f, Args&&... args)
{
	return external_thread_pool()->enqueue(f,args...);
}



/**
* @brief get the number of threads that are launched by foreach_xxx
*/
inline uint16_t nb_launched_threads()
{
	return uint16_t(std::thread::hardware_concurrency());
}



#endif 
