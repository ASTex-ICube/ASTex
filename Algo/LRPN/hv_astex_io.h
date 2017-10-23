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



#ifndef _HV_ASTEX_IO_H_
#define _HV_ASTEX_IO_H_

#include "hv/hvPictRGB.h"
#include "hv/hvPict.h"
#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>

namespace ASTex
{
template<typename T>
void load_to_hv(const std::string& filename, hview::hvPictRGB<T>& pict)
{
	ImageRGB<T> img;
	img.load(filename);
	pict.reset(img.width(),img.height(), hview::hvColRGB<T>(0));
	auto src = img.getDataPtr();
	auto dst = pict.data();
	memcpy(reinterpret_cast<uint8_t*>(dst),reinterpret_cast<uint8_t*>(src),img.width()*img.height()*3*sizeof(T));
}

template<typename T>
void save_from_hv(const std::string& filename, hview::hvPictRGB<T>& pict)
{
	ImageRGB<T> img = create_from_buffer<ImageRGB<T>>(pict.sizeX(), pict.sizeY(), reinterpret_cast<char*>(pict.data()), false);
	img.save(filename);
}

template<typename T>
void load_to_hv(const std::string& filename, hview::hvPict<T>& pict)
{
	ImageGray<T> img;
	img.load(filename);
	pict.reset(img.width(),img.height());
	auto src = img.getDataPtr();
	auto dst = pict.data();
	memcpy(reinterpret_cast<uint8_t*>(dst),reinterpret_cast<uint8_t*>(src),img.width()*img.height()*3*sizeof(T));
}

template<typename T>
void save_from_hv(const std::string& filename, hview::hvPict<T>& pict)
{
	ImageGray<T> img = create_from_buffer<ImageGray<T>>(pict.sizeX(), pict.sizeX(), reinterpret_cast<char*>(pict.data()), false);
	img.save(filename);
}

}

#endif
