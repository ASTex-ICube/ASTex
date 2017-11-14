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
#include <array>


namespace ASTex
{

// undefine RGB macro (exist when compile on Windows!!)
#ifdef RGB
#undef RGB
#endif



/**
 * @brief The RGB Image class
 * @tparam CHANNEL_TYPE (uchar/char/ushort/short/.../float/double)
 */
template<typename CHANNEL_TYPE>
class ImageRGBBase : public ImageBase
{
public:
	using PixelType =            itk::RGBPixel<CHANNEL_TYPE>;
	using ItkImg =               itk::Image< PixelType >;
	using IteratorIndexed =      itk::ImageRegionIteratorWithIndex<ItkImg>;
	using Iterator =             itk::ImageRegionIterator<ItkImg>;
	using ConstIteratorIndexed = itk::ImageRegionConstIteratorWithIndex<ItkImg> ;
	using ConstIterator =        itk::ImageRegionConstIterator<ItkImg>;
	using DoublePixelEigen =     Eigen::Vector3d;
	using DataType =             CHANNEL_TYPE;

	static const uint32_t NB_CHANNELS = 3;

protected:
	typename ItkImg::Pointer itk_img_;

public:
	ImageRGBBase():
		itk_img_(NULL)
	{}

	ImageRGBBase(typename itk::Image< itk::RGBPixel<CHANNEL_TYPE> >::Pointer itk_im):
		itk_img_(itk_im)
	{}

	/**
	 * @brief itkPixel
	 * @return a pixel of value (r/g/b)
	 */
	inline static PixelType itkPixel(CHANNEL_TYPE r, CHANNEL_TYPE g, CHANNEL_TYPE b)
	{
		return PixelType(std::array<CHANNEL_TYPE, 3>({ { r,g,b } }).data());
	}

	/**
	 * @brief itkPixel
	 * @return a pixel of value (r/r/r)
	 */
	inline static PixelType itkPixel(CHANNEL_TYPE r)
	{
		return PixelType(std::array<CHANNEL_TYPE, 3>({ { r,r,r } }).data());
	}

	/**
	 * @brief itkPixelNorm
	 * @param r normalized value [0,1]
	 * @param g normalized value [0,1]
	 * @param b normalized value [0,1]
	 * @return a pixel RGB
	 */
	template <bool B = true>
	inline static auto itkPixelNorm(double r, double g, double b) -> typename std::enable_if<B && std::is_arithmetic<CHANNEL_TYPE>::value, PixelType>::type
	{
		if (std::is_floating_point<CHANNEL_TYPE>::value)
			return itkPixel(r,g,b);

		if (std::is_unsigned<CHANNEL_TYPE>::value)
			return itkPixel(r * std::numeric_limits<CHANNEL_TYPE>::max(), g * std::numeric_limits<CHANNEL_TYPE>::max(), b * std::numeric_limits<CHANNEL_TYPE>::max());

		return itkPixel(r*(double(std::numeric_limits<CHANNEL_TYPE>::max()) - double(std::numeric_limits<CHANNEL_TYPE>::lowest())) + std::numeric_limits<CHANNEL_TYPE>::lowest(),
			g*(double(std::numeric_limits<CHANNEL_TYPE>::max()) - double(std::numeric_limits<CHANNEL_TYPE>::lowest())) + std::numeric_limits<CHANNEL_TYPE>::lowest(), 
			b*(double(std::numeric_limits<CHANNEL_TYPE>::max()) - double(std::numeric_limits<CHANNEL_TYPE>::lowest())) + std::numeric_limits<CHANNEL_TYPE>::lowest());
	}

	/**
	 * @brief eigenPixel
	 * @param p a itk::Pixel
	 * @return an Eigen::Vector3d with same value
	 */
	inline static DoublePixelEigen eigenPixel(const PixelType& p)
	{
		if (std::is_same<CHANNEL_TYPE,double>::value)
			return *(reinterpret_cast<const DoublePixelEigen*>(&p));
		return DoublePixelEigen(p[0],p[1],p[2]);
	}

	/**
	 * @brief itkPixel
	 * @param p an Eigen::Vector3d
	 * @return an itk::Pixel with same values of p
	 */
	inline static PixelType itkPixel(const DoublePixelEigen& p)
	{
		if (std::is_same<CHANNEL_TYPE,double>::value)
			return *(reinterpret_cast<const PixelType*>(&p));
		return itkPixel(CHANNEL_TYPE(p[0]),CHANNEL_TYPE(p[1]),CHANNEL_TYPE(p[2]));
	}

	inline static DoublePixelEigen eigenPixel(double v)
	{
		return DoublePixelEigen(v,v,v);
	}


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


template<typename CHANNEL_TYPE>
inline itk::RGBPixel<CHANNEL_TYPE> itkRGBPixel(CHANNEL_TYPE r, CHANNEL_TYPE g, CHANNEL_TYPE b)
{
	return ImageRGBBase<CHANNEL_TYPE>::itkPixel(r,g,b);
}

template<typename CHANNEL_TYPE>
inline itk::RGBPixel<CHANNEL_TYPE> itkRGBPixel(CHANNEL_TYPE r)
{
	return ImageRGBBase<CHANNEL_TYPE>::itkPixel(r);;
}


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



