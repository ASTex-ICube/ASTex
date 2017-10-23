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



#ifndef __ASTEX_CONVERSION__
#define __ASTEX_CONVERSION__

#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>
#include <ASTex/image_rgba.h>

namespace ASTex
{

/**
 * @brief convert a vector components pixel into average scalar pixel
 * @tparam TSCR input image type
 * @tparam TDST output image type
 * @param v_in input pixel
 * @param v_out output pixel
 */
template <typename TSRC, typename TDST>
void averageVec2Gray( const typename TSRC::PixelType& v_in, typename TDST::PixelType& v_out)
{
	double tempo = v_in[0];

	for(uint32_t i=1; i<TSRC::NB_CHANNELS; ++i)
	{
		tempo += v_in[i];
	}

	v_out= typename TDST::PixelType(tempo/TSRC::NB_CHANNELS);
}

/**
 * @brief convert a scalar pixel into &vector components pixel by copying values
 * @tparam TSCR input image type
 * @tparam TDST output image type
 * @param v_in input pixel
 * @param v_out output pixel
 */
template <typename TSRC, typename TDST>
void copyGray2Vec( const typename TSRC::PixelType& v_in, typename TDST::PixelType& v_out)
{
	for(uint32_t i=0; i<TDST::NB_CHANNELS; ++i)
	{
		v_out[i] = v_in;
	}
}

























/**
 * @brief image conversion function
 * @tparam TSCR input pointer image type (deduced)
 * @tparam TDST output pointer image type (deduced)
 * @param src input image pointer
 * @param dst output image pointer
 * @param convFunc pixel conversion function void f(const T1& in, T2& out)
 */
template <typename TSRC, typename TDST, typename CONV>
void convertImage(typename TSRC::Pointer src, typename TDST::Pointer dst, CONV convFunc)
{
	dst->SetRegions(src->GetLargestPossibleRegion());
	dst->Allocate();

	uint32_t W = src->width();
	uint32_t H = src->height();

	for(uint32_t j=0; j<H  ; ++j)
	{
		for(uint32_t i=0; i< W ; ++i)
		{
			convFunc(src->pixel(i,j), dst->pixel(i,j));
		}
	}
}




















}

#endif



