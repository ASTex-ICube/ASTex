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






template<typename T>
/**
 * @brief interlace 3 scalar images into rgb image
 * @param im_r
 * @param im_g
 * @param im_b
 * @return
 */
ImageRGB<T> interlace(const ImageGray<T>& im_r, const ImageGray<T>& im_g, const ImageGray<T>& im_b)
{
	using GRAY = typename ImageGray<T>::PixelType;
	using RGB = typename ImageRGB<T>::PixelType;
	using OUT_IT = typename ImageRGB<T>::Iterator;
	using IN_IT  = typename ImageGray<T>::ConstIterator;

	ImageRGB<T> res(im_r.width(),im_r.height());

	IN_IT ir = im_r.beginConstIterator();
	IN_IT ig = im_g.beginConstIterator();
	IN_IT ib = im_b.beginConstIterator();
	IN_IT ia = im_a.beginConstIterator();
	for (OUT_IT irgb = res.beginIterator(); !irgb.IsAtEnd(); ++irgb)
	{
		irgb.Value() = RGBA(ir.Value(),ig.Value(),ib.Value());
		++ir; ++ig; ++ib;
	}
}

template<typename T>
ImageRGBA<T> interlace(const ImageGray<T>& im_r, const ImageGray<T>& im_g, const ImageGray<T>& im_b, const ImageGray<T>& im_a)
{
	using GRAY = typename ImageGray<T>::PixelType;
	using RGBA = typename ImageRGBA<T>::PixelType;
	using OUT_IT = typename ImageRGBA<T>::Iterator;
	using IN_IT  = typename ImageGray<T>::ConstIterator;

	ImageRGBA<T> res(im_r.width(),im_r.height());

	IN_IT ir = im_r.beginConstIterator();
	IN_IT ig = im_g.beginConstIterator();
	IN_IT ib = im_b.beginConstIterator();
	IN_IT ia = im_a.beginConstIterator();
	for (OUT_IT irgba = res.beginIterator(); !irgba.IsAtEnd(); ++irgba)
	{
		irgba.Value() = RGBA(ir.Value(),ig.Value(),ib.Value(),ia.Value());
		++ir; ++ig; ++ib; ++ia;
	}
}


template<typename T>
const deinterlace(const ImageRGB<T> im_rgb, ImageGray<T>& im_r, ImageGray<T>& im_g, ImageGray<T>& im_b)
{
	using GRAY = typename ImageGray<T>::PixelType;
	using RGB = typename ImageRGB<T>::PixelType;
	using IN_IT = typename ImageRGB<T>::ConstIterator;
	using OUT_IT  = typename ImageGray<T>::Iterator;

	OUT_IT ir = im_r.beginIterator();
	OUT_IT ig = im_g.beginIterator();
	OUT_IT ib = im_b.beginIterator();
	OUT_IT ia = im_a.beginIterator();
	for (IN_IT irgb = im_rgb.beginIterator(); !irgb.IsAtEnd(); ++irgb)
	{
		const RGB& P = irgb.Value();
		ir.Value() = P.GetRed();
		ig.Value() = P.GetGreen();
		ib.Value() = P.GetBlue();
		++ir; ++ig; ++ib;
	}
}

template<typename T>
const deinterlace(const ImageRGBA<T> im_rgba, ImageGray<T>& im_r, ImageGray<T>& im_g, ImageGray<T>& im_b, ImageGray<T>& im_a)
{
	using GRAY = typename ImageGray<T>::PixelType;
	using RGBA = typename ImageRGBA<T>::PixelType;
	using IN_IT = typename ImageRGBA<T>::ConstIterator;
	using OUT_IT  = typename ImageGray<T>::Iterator;

	OUT_IT ir = im_r.beginIterator();
	OUT_IT ig = im_g.beginIterator();
	OUT_IT ib = im_b.beginIterator();
	OUT_IT ia = im_a.beginIterator();
	for (IN_IT irgba = im_rgba.beginIterator(); !irgba.IsAtEnd(); ++irgba)
	{
		const RGBA& P = irgba.Value();
		ir.Value() = P.GetRed();
		ig.Value() = P.GetGreen();
		ib.Value() = P.GetBlue();
		ia.Value() = P.GetAlpha();
		++ir; ++ig; ++ib; ++ia;
	}
}

































}

#endif



