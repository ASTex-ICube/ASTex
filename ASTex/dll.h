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



#ifndef __ASTEX_DLL_H_
#define __ASTEX_DLL_H_

/**
* \brief Linkage declaration for ASTex symbols.
*/
#ifdef WIN32 
#ifdef BUILD_STATIC
#define ASTEX_API
#endif
#ifndef ASTEX_API
#if defined ASTEX_DLL_EXPORT
#define ASTEX_API __declspec(dllexport)
#else
#define ASTEX_API __declspec(dllimport)
#endif
#endif
#else
#define ASTEX_API
#endif

#endif // __ASTEX_DLL_H_
