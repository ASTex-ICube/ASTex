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



#ifndef __ASTEX_IMAGE_RGBA__
#define __ASTEX_IMAGE_RGBA__

#include "itkImage.h"
#include <ASTex/image_common.h>

namespace ASTex
{

/**
 * @brief The Pixel RGBA class (just to add nice constructor)
 * @tparam CHANNEL_TYPE (uchar/char/ushort/short/.../float/double)
 */
template<typename CHANNEL_TYPE>
class RGBA: public itk::RGBAPixel<CHANNEL_TYPE>
{
	using Inherit =  itk::RGBAPixel<CHANNEL_TYPE>;
	using Self =  RGBA<CHANNEL_TYPE>;
public:
	RGBA() {}

	RGBA(const Inherit& itkrgba):
		Inherit(itkrgba)
	{}

	RGBA(CHANNEL_TYPE r, CHANNEL_TYPE g, CHANNEL_TYPE b, CHANNEL_TYPE a)
	{
		Inherit::Set(r,g,b,a);
	}

	explicit RGBA(CHANNEL_TYPE v)
	{
		Inherit::Set(v,v,v,v);
	}

	template<typename PIX>
	explicit RGBA(const PIX& p, typename std::enable_if<std::is_base_of<itk::FixedArray<CHANNEL_TYPE,4>,PIX>::value>::type* = nullptr)
	{
		Inherit::Set(CHANNEL_TYPE(p[0]),CHANNEL_TYPE(p[1]),CHANNEL_TYPE(p[2]),CHANNEL_TYPE(p[3]));
	}

	/**
	 * @brief convert a reference to itk::RGBAPixel into a reference to RGBA
	 */
	static Self& convert(Inherit& itkrgba)
	{
		return *(reinterpret_cast<Self*>(&itkrgba));
	}

	/**
	 * @brief convert a const reference to itk::RGBAPixel into a const reference to RGBA
	 */
	static const Self& convert(const Inherit& itkrgba)
	{
		return *(reinterpret_cast<const Self*>(&itkrgba));
	}

	inline Self& operator += (const Self& p)
	{
		Inherit::operator += (p);
		return *this;
	}

	inline Self& operator -= (const Self& p)
	{
		Inherit::operator -= (p);
		return *this;
	}

	inline Self& operator *= (CHANNEL_TYPE s)
	{
		Inherit::operator *= (s);
		return *this;
	}

	inline Self& operator /= (CHANNEL_TYPE s)
	{
		Inherit::operator[](0) /= s;
		Inherit::operator[](1) /= s;
		Inherit::operator[](2) /= s;
		Inherit::operator[](3) /= s;
		return *this;
	}


	inline Self operator + (const Self& p) const
	{
		return Inherit::operator + (p);
	}

	inline Self operator - (const Self& p) const
	{
		return Inherit::operator - (p);
	}


	inline Self operator * (CHANNEL_TYPE s) const
	{
		Self res = Self(Inherit::GetRed()*s,Inherit::GetGreen()*s,Inherit::GetBlue()*s,Inherit::GetAlpha()*s);
		return res;
	}

	inline Self operator / (CHANNEL_TYPE s) const
	{
		Self res = Self(Inherit::GetRed()/s,Inherit::GetGreen()/s,Inherit::GetBlue()/s,Inherit::GetAlpha()/s);
		return res;
	}

};

template<typename CHANNEL_TYPE>
inline RGBA<CHANNEL_TYPE> operator * (CHANNEL_TYPE s, const RGBA<CHANNEL_TYPE>& rgba)
{
	return rgba*s;
}


/**
 * @brief The RGBA image class
 * @tparam CHANNEL_TYPE (uchar/char/ushort/short/.../float/double)
 */
template<typename CHANNEL_TYPE>
class ImageRGBABase : public ImageBase
{
protected:
	typename itk::Image< itk::RGBAPixel<CHANNEL_TYPE> >::Pointer itk_img_;


public:

	typedef itk::Image< itk::RGBAPixel<CHANNEL_TYPE> >	   ItkImg;
	typedef itk::ImageRegionIteratorWithIndex<ItkImg>      IteratorIndexed;
	typedef itk::ImageRegionIterator<ItkImg>               Iterator;
	typedef itk::ImageRegionConstIteratorWithIndex<ItkImg> ConstIteratorIndexed;
	typedef itk::ImageRegionConstIterator<ItkImg>          ConstIterator;

	static const uint32_t NB_CHANNELS = 4;
	typedef itk::RGBAPixel<CHANNEL_TYPE> PixelType;
	typedef itk::RGBAPixel<double> DoublePixelType;
	typedef RGBA<CHANNEL_TYPE> ASTexPixelType;
	typedef CHANNEL_TYPE DataType;

	ImageRGBABase():
		itk_img_(NULL)
	{}

	ImageRGBABase(typename itk::Image< itk::RGBAPixel<CHANNEL_TYPE> >::Pointer itk_im):
		itk_img_(itk_im)
	{}

protected:

	inline PixelType* getPixelsPtr()
	{
		return this->itk_img_->GetPixelContainer()->GetBufferPointer();
	}

	inline DataType* getDataPtr()
	{
		return getPixelsPtr()->GetDataPointer();
	}

	inline const PixelType* getPixelsPtr() const
	{
		return this->itk_img_->GetPixelContainer()->GetBufferPointer();
	}

	inline const DataType* getDataPtr() const
	{
		return getPixelsPtr()->GetDataPointer();
	}

};




using RGBAu8  = RGBA< uint8_t >;
using RGBA8   = RGBA< int8_t >;
using RGBAu16 = RGBA< uint16_t >;
using RGBA16  = RGBA< int16_t >;
using RGBAu32 = RGBA< uint32_t >;
using RGBA32  = RGBA< int32_t >;
using RGBAu64 = RGBA< uint64_t >;
using RGBA64  = RGBA< int64_t >;
using RGBAf   = RGBA< float >;
using RGBAd   = RGBA< double >;



template <class T> using ImageRGBA = ImageCommon< ImageRGBABase< T >, false >;
using ImageRGBAu8  = ImageRGBA< uint8_t >;
using ImageRGBA8   = ImageRGBA< int8_t >;
using ImageRGBAu16 = ImageRGBA< uint16_t >;
using ImageRGBA16  = ImageRGBA< int16_t >;
using ImageRGBAu32 = ImageRGBA< uint32_t >;
using ImageRGBA32  = ImageRGBA< int32_t >;
using ImageRGBAu64 = ImageRGBA< uint64_t >;
using ImageRGBA64  = ImageRGBA< int64_t >;
using ImageRGBAf   = ImageRGBA< float >;
using ImageRGBAd   = ImageRGBA< double >;


// ConstImages are Images that only support const operations, they are not modifiable.
// But they could be construct from Itk::Image::ConstPointer -> usefull in Filter writing.

template <class T> using ConstImageRGBA = ImageCommon< ImageRGBABase< T >, true >;
using ConstImageRGBAu8  = ConstImageRGBA< uint8_t >;
using ConstImageRGBA8   = ConstImageRGBA< int8_t >;
using ConstImageRGBAu16 = ConstImageRGBA< uint16_t >;
using ConstImageRGBA16  = ConstImageRGBA< int16_t >;
using ConstImageRGBAu32 = ConstImageRGBA< uint32_t >;
using ConstImageRGBA32  = ConstImageRGBA< int32_t >;
using ConstImageRGBAu64 = ConstImageRGBA< uint64_t >;
using ConstImageRGBA64  = ConstImageRGBA< int64_t >;
using ConstImageRGBAf   = ConstImageRGBA< float >;
using ConstImageRGBAd   = ConstImageRGBA< double >;





}


#endif



