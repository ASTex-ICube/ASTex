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



// hvError.h: interface for the runtime Error class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ERROR_H__09002DDD_2472_43B9_B7B6_FCB6FF1B6B0D__INCLUDED_)
#define AFX_ERROR_H__09002DDD_2472_43B9_B7B6_FCB6FF1B6B0D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <iostream>
#include <cstdlib>

inline void hvFatal(const char *str)
{
    std::cout<<"Internal error:"<<str<<"\n";
//    char buffer[16];
//	std::cin.getline(buffer,16);
    std::abort();
//    std::exit();
//    std::terminate();
}

#endif // !defined(AFX_ARRAY1_H__09002DDD_2472_43B9_B7B6_FCB6FF1B6B0D__INCLUDED_)
