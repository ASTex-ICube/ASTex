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



#ifndef __ASTEX_IMAGE_EASY_IO__
#define __ASTEX_IMAGE_EASY_IO__

#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>
#include <ASTex/image_rgba.h>

namespace ASTex
{

template <typename PIX>
inline auto NICE(const PIX& p) -> typename std::enable_if<pixel_traits<PIX>::dim ==4, typename Eigen::Matrix<double,1,pixel_traits<PIX>::dim>>::type
{
	return eigenPixel<double>(p).transpose();
}

template <typename PIX>
inline auto NICE(const PIX& p) -> typename std::enable_if<pixel_traits<PIX>::dim == 1, PIX>::type
{
  return double(p);
}


namespace IO
{

void ASTEX_API save01_in_u8(const ImageGrayd& img, const std::string& filename);

void ASTEX_API save01_in_u8(const ImageGrayf& img, const std::string& filename);

void ASTEX_API save01_in_u8(const ImageRGBd& img, const std::string& filename);

void ASTEX_API save01_in_u8(const ImageRGBf& img, const std::string& filename);

void ASTEX_API save01_in_u8(const ImageRGBAd& img, const std::string& filename);

void ASTEX_API save01_in_u8(const ImageRGBAf& img, const std::string& filename);



bool ASTEX_API loadu8_in_01(ImageGrayd& img, const std::string& filename);

bool ASTEX_API loadu8_in_01(ImageGrayf& img, const std::string& filename);

bool ASTEX_API loadu8_in_01(ImageRGBd& img, const std::string& filename);

bool ASTEX_API loadu8_in_01(ImageRGBf& img, const std::string& filename);

bool ASTEX_API loadu8_in_01(ImageRGBAd& img, const std::string& filename);

bool ASTEX_API loadu8_in_01(ImageRGBAf& img, const std::string& filename);

}
}

#endif
