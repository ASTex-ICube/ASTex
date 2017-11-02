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



#ifndef __ASTEX_IMAGE_RGB__
#define __ASTEX_IMAGE_RGB__

#include <ASTex/image_common.h>


namespace ASTex
{

// undefine RGB macro (exist when compile on Windows!!)
#ifdef RGB
#undef RGB
#endif


template<typename CHANNEL_TYPE>
inline itk::RGBPixel<CHANNEL_TYPE> itkRGBPixel(CHANNEL_TYPE r, CHANNEL_TYPE g, CHANNEL_TYPE b)
{
	return itk::RGBPixel<CHANNEL_TYPE>(std::array<CHANNEL_TYPE,3>({{r,g,b}}).data());
}




/**
 * @brief The Pixem RGB class (just to add nice constructor)
 * @tparam CHANNEL_TYPE (uchar/char/ushort/short/.../float/double)
 */
template<typename CHANNEL_TYPE>
class RGB: public itk::RGBPixel<CHANNEL_TYPE>
{
	using Inherit =  itk::RGBPixel<CHANNEL_TYPE>;
	using Self =  RGB<CHANNEL_TYPE>;
public:
	RGB() {}

	RGB(const Inherit& itkrgb):
		Inherit(itkrgb)
	{}

	RGB(CHANNEL_TYPE r, CHANNEL_TYPE g, CHANNEL_TYPE b)
	{
		Inherit::Set(r,g,b);
	}

	explicit RGB(CHANNEL_TYPE v)
	{
		Inherit::Set(v,v,v);
	}

	template<typename PIX>
	explicit RGB(const PIX& p, typename std::enable_if<std::is_base_of<itk::FixedArray<CHANNEL_TYPE,3>,PIX>::value>::type* = nullptr)
	{
		Inherit::Set(CHANNEL_TYPE(p[0]),CHANNEL_TYPE(p[1]),CHANNEL_TYPE(p[2]));
	}


	/**
	 * @brief convert a reference to itk::RGBPixel into a reference to RGB
	 */
	static Self& convert(Inherit& itkrgb)
	{
		return *(reinterpret_cast<Self*>(&itkrgb));
	}

	/**
	 * @brief convert a const reference to itk::RGBPixel into a const reference to RGB
	 */
	static const Self& convert(const Inherit& itkrgb)
	{
		return *(reinterpret_cast<const Self*>(&itkrgb));
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
		Self res = Self(Inherit::GetRed()*s,Inherit::GetGreen()*s,Inherit::GetBlue()*s);
		return res;
	}

	inline Self operator / (CHANNEL_TYPE s) const
	{
		Self res = Self(Inherit::GetRed()/s,Inherit::GetGreen()/s,Inherit::GetBlue()/s);
		return res;
	}
};

template<typename CHANNEL_TYPE>
inline RGB<CHANNEL_TYPE> operator * (CHANNEL_TYPE s, const RGB<CHANNEL_TYPE>& rgb)
{
	return rgb*s;
}

template<typename T>
itk::RGBPixel<T> blend(const itk::RGBPixel<T>& c1, const itk::RGBPixel<T>& c2, double alpha)
{
	double beta = 1.0-alpha;
	double r = double(c1.GetRed())*alpha+double(c2.GetRed())*beta;
	double g = double(c1.GetGreen())*alpha+double(c2.GetGreen())*beta;
	double b = double(c1.GetBlue())*alpha+double(c2.GetBlue())*beta;
	return  itkRGBPixel((T)r, (T)g, (T)b);
}


/**
 * @brief The RGB Image class
 * @tparam CHANNEL_TYPE (uchar/char/ushort/short/.../float/double)
 */
template<typename CHANNEL_TYPE>
class ImageRGBBase : public ImageBase
{
protected:
	typename itk::Image< itk::RGBPixel<CHANNEL_TYPE> >::Pointer itk_img_;


public:

	typedef itk::Image< itk::RGBPixel<CHANNEL_TYPE> >	   ItkImg;
	typedef itk::ImageRegionIteratorWithIndex<ItkImg>      IteratorIndexed;
	typedef itk::ImageRegionIterator<ItkImg>               Iterator;
	typedef itk::ImageRegionConstIteratorWithIndex<ItkImg> ConstIteratorIndexed;
	typedef itk::ImageRegionConstIterator<ItkImg>          ConstIterator;


	static const uint32_t NB_CHANNELS = 3;
	typedef itk::RGBPixel<CHANNEL_TYPE> PixelType;
	typedef RGB<double> DoublePixelType;
	typedef RGB<CHANNEL_TYPE> ASTexPixelType;
	typedef CHANNEL_TYPE DataType;


	ImageRGBBase():
		itk_img_(NULL)
	{}

	ImageRGBBase(typename itk::Image< itk::RGBPixel<CHANNEL_TYPE> >::Pointer itk_im):
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


using RGBu8  = RGB< uint8_t >;
using RGB8   = RGB< int8_t >;
using RGBu16 = RGB< uint16_t >;
using RGB16  = RGB< int16_t >;
using RGBu32 = RGB< uint32_t >;
using RGB32  = RGB< int32_t >;
using RGBu64 = RGB< uint64_t >;
using RGB64  = RGB< int64_t >;
using RGBf   = RGB< float >;
using RGBd   = RGB< double >;



template <class T> using ImageRGB = ImageCommon< ImageRGBBase< T >, false >;
using ImageRGBu8  = ImageRGB< uint8_t >;
using ImageRGB8   = ImageRGB< int8_t >;
using ImageRGBu16 = ImageRGB< uint16_t >;
using ImageRGB16  = ImageRGB< int16_t >;
using ImageRGBu32 = ImageRGB< uint32_t >;
using ImageRGB32  = ImageRGB< int32_t >;
using ImageRGBu64 = ImageRGB< uint64_t >;
using ImageRGB64  = ImageRGB< int64_t >;
using ImageRGBf   = ImageRGB< float >;
using ImageRGBd   = ImageRGB< double >;


// ConstImages are Images that only support const operations, they are not modifiable.
// But they could be construct from Itk::Image::ConstPointer -> usefull in Filter writing.

template <class T> using ConstImageRGB = ImageCommon< ImageRGBBase< T >, true >;
using ConstImageRGBu8  = ConstImageRGB< uint8_t >;
using ConstImageRGB8   = ConstImageRGB< int8_t >;
using ConstImageRGBu16 = ConstImageRGB< uint16_t >;
using ConstImageRGB16  = ConstImageRGB< int16_t >;
using ConstImageRGBu32 = ConstImageRGB< uint32_t >;
using ConstImageRGB32  = ConstImageRGB< int32_t >;
using ConstImageRGBu64 = ConstImageRGB< uint64_t >;
using ConstImageRGB64  = ConstImageRGB< int64_t >;
using ConstImageRGBf   = ConstImageRGB< float >;
using ConstImageRGBd   = ConstImageRGB< double >;





}

#endif



