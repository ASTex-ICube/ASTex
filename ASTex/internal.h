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



#ifndef __ASTEX_INTERNAL__
#define __ASTEX_INTERNAL__

#include <type_traits>
#include <algorithm>
#include <string>
#include <cassert>
#include <cstdint>

namespace ASTex
{

namespace IO
{

inline std::string remove_path(const std::string& full_path_name)
{
	std::size_t x = full_path_name.rfind('/');
	if (x==std::string::npos)
	{
#ifdef WIN32
		x = full_path_name.rfind('\\');
		if (x == std::string::npos)
#endif
			return full_path_name;
	}
	return full_path_name.substr(x+1, std::string::npos);
}

inline std::string remove_ext(const std::string& file_name)
{
	std::size_t x = file_name.rfind('.');
	if (x==std::string::npos)
	{
		return file_name;
	}
	return file_name.substr(0,x);
}

}




//template <typename T>
//struct function_traits : public function_traits<decltype(&T::operator())>
//{};

//template <typename ClassType, typename ReturnType, typename... Args>
//struct function_traits<ReturnType(ClassType::*)(Args...) const>
//// we specialize for pointers to member function
//{
//	static const size_t arity = sizeof...(Args);
//	// arity is the number of arguments.

//	using result_type = ReturnType;

//	template <size_t i>
//	struct arg
//	{
//		static_assert(i < sizeof...(Args), "Trying to access to an argument whose index is higher than the function arity.");
//		using type = typename std::tuple_element<i, std::tuple<Args...>>::type;
//		// the i-th argument is equivalent to the i-th tuple element of a tuple
//		// composed of those arguments.
//	};
//};


class AnyType
{};

class NoneType
{};


// specialization for lambda functions
template <typename T>
struct function_traits : public function_traits<decltype(&T::operator())>
{};

// General case
template <typename ReturnType, typename... Args>
struct function_traits<ReturnType(Args...)>
{
	static const size_t arity = sizeof...(Args);

	using result_type = ReturnType;

	template <size_t i,  class Enable = void>
	struct arg
	{
		using type = NoneType;
	};

	template <size_t i>
	struct arg<i, typename std::enable_if< (i < sizeof...(Args))>::type>
	{
		using type = typename std::tuple_element<i, std::tuple<Args...>>::type;
	};
};


// specialization for function pointers
template <typename ReturnType, typename... Args>
struct function_traits<ReturnType(*)(Args...)> : public function_traits<ReturnType(Args...)> {};

// specialization for function references
template <typename ReturnType, typename... Args>
struct function_traits<ReturnType(&)(Args...)> : public function_traits<ReturnType(Args...)> {};

// specialization for member function pointers
template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits<ReturnType(ClassType::*)(Args...)>: public function_traits<ReturnType(Args...)> {};

// specialization for const member function pointers
template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits<ReturnType(ClassType::*)(Args...) const> : public function_traits<ReturnType(Args...)>{};







template<typename F>
using func_return_type = typename function_traits<F>::result_type;

template<typename F>
using return_bool = std::is_same<func_return_type<F>, bool>;



template <int N, typename FUNC, typename T,typename... Args>
struct check_param;

template <int N, typename FUNC, typename T>
struct check_param<N,FUNC,T>
{
	using TFR = typename function_traits<FUNC>::result_type;
	static const bool value = std::is_same<T,TFR>::value || std::is_same<AnyType,T>::value;
};

template <int N, typename FUNC, typename T, typename... Args>
struct check_param
{
	using TFP = typename function_traits<FUNC>::template arg<N>::type;
	static const bool local_value = std::is_same<T,TFP>::value || std::is_same<AnyType,T>::value;
	static const bool value = local_value && check_param<N+1,FUNC, Args...>::value;
};


template <typename FUNC, typename T, typename... Args>
struct check_signature
{
	static const bool nb_value = function_traits<FUNC>::arity == sizeof...(Args);
	static const bool value = nb_value && check_param<0,FUNC,T,Args...>::value;
};




/// P,i,j,t -> bool

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,T1,int32_t,int32_t,uint16_t,bool>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, T1& p1, int p2, int p3, uint16_t p4)
{
	return func(p1, p2, p3, p4);
}

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,const T1&, int32_t,int32_t,uint16_t,bool>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, const T1& p1, int p2, int p3, uint16_t p4)
{
	return func(p1, p2, p3, p4);
}


/// i,j,t -> bool

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,int32_t,int32_t,uint16_t,bool>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, T1 /*p1*/, int p2, int p3, uint16_t p4)
{
	return func(p2, p3, p4);
}



/// P,i,j -> bool

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,T1,int32_t,int32_t,bool>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, T1& p1, int p2, int p3, uint16_t)
{
	return func(p1,p2, p3);
}

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,const T1&, int32_t,int32_t,bool>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, const T1& p1, int p2, int p3, uint16_t)
{
	return func(p1,p2, p3);
}


/// P,t -> bool

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,T1,uint16_t,bool>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, T1& p1, int, int, uint16_t p4)
{
	return func(p1, p4);
}

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,const T1&, uint16_t,bool>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, const T1& p1, int, int, uint16_t p4)
{
	return func(p1, p4);
}


/// i,j -> bool
template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,int32_t,int32_t,bool>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, T1, int p2, int p3, uint16_t)
{
	return func(p2, p3);
}


/// P -> bool
template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,T1&,bool>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, T1& p1, int, int, uint16_t)
{
	return func(p1);
}

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,const T1&, bool>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, const T1& p1, int, int, uint16_t)
{
	return func(p1);
}


/// P,i,j,t -> void

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,T1&,int32_t,int32_t,uint16_t,void>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, T1& p1, int p2, int p3, uint16_t p4)
{
	func(p1, p2, p3, p4);
	return true;
}

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC, const T1&, int32_t,int32_t,uint16_t,void>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, const T1& p1, int p2, int p3, uint16_t p4)
{
	func(p1, p2, p3, p4);
	return true;
}


/// i,j,t -> void

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,int32_t,int32_t,uint16_t,void>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, T1 /*p1*/, int p2, int p3, uint16_t p4)
{
	func(p2, p3, p4);
	return true;
}


/// P,i,j -> void

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,T1&,int32_t,int32_t,void>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, T1& p1, int p2, int p3, uint16_t)
{
	func(p1,p2, p3);
	return true;
}

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,const T1&, int32_t,int32_t,void>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, const T1& p1, int p2, int p3, uint16_t)
{
	func(p1,p2, p3);
	return true;
}


/// P,t -> void

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,T1&,uint16_t,void>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, T1& p1, int, int, uint16_t p4)
{
	func(p1, p4);
	return true;
}

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,const T1&, uint16_t,void>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, const T1& p1, int, int, uint16_t p4)
{
	func(p1, p4);
	return true;
}

/// i,j -> void
template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,int32_t,int32_t,void>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, T1, int p2, int p3, uint16_t)
{
	func(p2, p3);
	return true;
}


/// P -> void

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC,T1&,void>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, T1& p1, int, int, uint16_t)
{
	func(p1);
	return true;
}

template < typename T1, typename FUNC,
typename std::enable_if<check_signature<FUNC, const T1&, void>::value>::type* = nullptr>
inline bool four_param_binder(const FUNC& func, const T1& p1, int, int, uint16_t)
{
	func(p1);
	return true;
}



inline void assert_msg(bool b, const std::string& msg)
{
	if (!b)
		std::cerr << msg << std::endl;
	assert(b);
}


} //ASTex

#endif

