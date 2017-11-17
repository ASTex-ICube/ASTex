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

#include<vector>

#ifndef _ASTEX_STORE_H_
#define _ASTEX_STORE_H_


namespace ASTex
{


template <typename T>
class Store
{
public:

	using Self = Store<T>;
	using value_type = T;

protected:

	static const uint32_t CHUNK_SIZE = 4096;
	// vector of block pointers
	std::vector<T*> table_data_;
	uint32_t size_;

public:

	inline Store() :
		size_(0)
	{
		table_data_.reserve(1024u);
	}

	~Store()
	{
		for(auto chunk : table_data_)
			delete[] chunk;
	}

	uint32_t element_size() const
	{
		return sizeof(std::declval<T>());
	}

	/**
	 * @brief get the capacity of the array
	 * @return number of allocated lines
	 */
	uint32_t capacity() const
	{
		return uint32(table_data_.size())*CHUNK_SIZE;
	}

	void clear()
	{
		for(auto chunk : table_data_)
			delete[] chunk;
		table_data_.clear();
		table_data_.shrink_to_fit();
		table_data_.reserve(1024u);
		size_ = 0;
	}

	void shrink_to_fit()
	{
		for (uint32_i= (size/CHUNK_SIZE)+1 ; i<table_data_.size(); ++i)
				delete[] table_data_[i];
	}


	inline void push_back(const T& v)
	{
		uint32_t j = size_%CHUNK_SIZE;
		if (j == 0)
			table_data_.push_back(new T[CHUNK_SIZE]);
		(table_data_.back())[j] = v;
		++size_;
	}

	inline Self& operator << (const T& v)
	{
		push_back(v);
		return *this;
	}

	inline void push_back(const T&& v)
	{
		uint32_t j = size_%CHUNK_SIZE;
		if (j == 0)
			table_data_.push_back(new T[CHUNK_SIZE]);
		(table_data_.back())[j] = std::move(v);
		++size_;
	}

	inline Self& operator << (const T&& v)
	{
		push_back(v);
		return *this;
	}


	template< class... Args >
	void emplace_back(Args&&... args )
	{
		uint32_t j = size_%CHUNK_SIZE;
		if (j == 0)
			table_data_.push_back(new T[CHUNK_SIZE]);
		(table_data_.back())[j] = T(args...);
		++size_;
	}


	inline void pop_back()
	{
		assert(size_>0);
		--size_;
	}

	inline T& back()
	{
		assert(size_>0);
		uint32_t sz = size_-1;
		return table_data_[sz / CHUNK_SIZE][sz % CHUNK_SIZE];
	}

	inline T& operator[](uint32_t i)
	{
		assert(i / CHUNK_SIZE < table_data_.size());
		return table_data_[i / CHUNK_SIZE][i % CHUNK_SIZE];
	}

	inline const T& operator[](uint32_t i) const
	{
		assert(i / CHUNK_SIZE < table_data_.size());
		return table_data_[i / CHUNK_SIZE][i % CHUNK_SIZE];
	}

	/**
	 * @brief remove warning do not preserve order of element
	 * @param i index
	 */
	inline void remove(uint32_t i)
	{
		assert(i / CHUNK_SIZE < table_data_.size());
		table_data_[i / CHUNK_SIZE][i % CHUNK_SIZE] = back();
		pop_back();

	}

	inline void set_all_values(const T& v)
	{
		for (T* chunk : table_data_)
		{
			for(uint32_t i = 0; i < CHUNK_SIZE; ++i)
				*chunk++ = v;
		}
	}


	class iterator
	{
		const std::vector<T*>& td_;
		uint32_t sz_;
		uint32_t i_;

	public:

		inline iterator(const Self* S, uint32_t i) :
			td_(S->table_data_),sz_(S->size_),i_(i)
		{}

		inline iterator(const iterator& it) :
			td_(it.td_),
			sz_(it.sz_),
			i_(it.i_)
		{}

		inline iterator& operator=(const iterator& it)
		{
			td_ = it.td_;
			sz_ = it.sz_;
			i_ = it.i_;
			return *this;
		}

		inline iterator& operator++()
		{
			++i_;
			return *this;
		}

		inline T& operator*()
		{
			return td_[i_ / CHUNK_SIZE][i_ % CHUNK_SIZE];
		}

		inline bool operator!=(iterator it) const
		{
			return i_ != it.i_;
		}

	};

	inline iterator begin()
	{
		return iterator(this,0);
	}

	inline iterator end()
	{
		return iterator(this, size_);
	}

	class const_iterator
	{
		const std::vector<T*>& td_;
		uint32_t sz_;
		uint32_t i_;

	public:

		inline const_iterator(const Self* S, uint32_t i) :
			td_(S->table_data_),sz_(S->size_),i_(i)
		{}

		inline const_iterator(const const_iterator& it) :
			td_(it.td_),
			sz_(it.sz_),
			i_(it.i_)
		{}

		inline const_iterator& operator=(const iterator& it)
		{
			td_ = it.td_;
			sz_ = it.sz_;
			i_ = it.i_;
			return *this;
		}

		inline const_iterator& operator++()
		{
			++i_;
			return *this;
		}

		inline const T& operator*()
		{
			return td_[i_ / CHUNK_SIZE][i_ % CHUNK_SIZE];
		}

		inline bool operator!=(const_iterator it) const
		{
			return i_ != it.i_;
		}

	};

	inline const_iterator begin() const
	{
		return const_iterator(this,0);
	}

	inline const_iterator end() const
	{
		return const_iterator(this, size_);
	}

};












} //ASTex

#endif
