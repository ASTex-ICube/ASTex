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



#define _UTILS_THREADPOOL_CPP_

#include "thread_pool.h"

ASTEX_TLS uint32_t thread_index_;

ASTEX_API ThreadPool* internal_thread_pool()
{
	// thread safe accoring to http://stackoverflow.com/questions/8102125/is-local-static-variable-initialization-thread-safe-in-c11
	static ThreadPool pool_i(std::thread::hardware_concurrency());
	return &pool_i;
}

ASTEX_API ThreadPool* external_thread_pool()
{
	// thread safe accoring to http://stackoverflow.com/questions/8102125/is-local-static-variable-initialization-thread-safe-in-c11
	static ThreadPool pool_e(32);
	return &pool_e;
}


ASTEX_API uint32_t current_thread_index()
{
	return thread_index_;
}




ThreadPool::~ThreadPool()
{
	nb_working_workers_ = uint32_t(workers_.size());
	condition_running_.notify_all();
	{
		std::unique_lock<std::mutex> lock(queue_mutex_);
		stop_ = true;
	}
#if !(defined(CGOGN_WIN_VER) && (CGOGN_WIN_VER <= 61))
	condition_running_.notify_all();
	condition_.notify_all();
#endif
	for(std::thread &worker: workers_)
		worker.join();
}


ThreadPool::ThreadPool(uint32_t nb_ww)
	:  stop_(false)
{
	this->nb_working_workers_ = nb_ww;

	for(uint32_t i = 0u; i< nb_ww; ++i)
	{
		workers_.emplace_back(
		[this, i] () -> void
		{
			thread_index_ = i;
			for(;;)
			{
				while (i >= this->nb_working_workers_)
				{
					std::unique_lock<std::mutex> lock(this->running_mutex_);
					this->condition_running_.wait(lock);
				}

				std::unique_lock<std::mutex> lock(this->queue_mutex_);
				this->condition_.wait(
					lock,
					[this] { return this->stop_ || !this->tasks_.empty(); }
				);

				if (this->stop_ && this->tasks_.empty())
				{
					return;
				}

				if (i < this->nb_working_workers_)
				{
					auto task = std::move(this->tasks_.front());
					this->tasks_.pop();
					lock.unlock();
#if defined(_MSC_VER) && _MSC_VER < 1900
					(*task)();
#else
					task();
#endif
				}
				else
				{
					lock.unlock();
					condition_.notify_one();
				}
			}
		});
	}
}

void ThreadPool::set_nb_workers(uint32_t nb )
{
	if (nb == 0xffffffff)
		nb_working_workers_ = uint32_t(workers_.size());
	else
		nb_working_workers_ = std::min(uint32_t(workers_.size()), nb);

	condition_running_.notify_all();
}
