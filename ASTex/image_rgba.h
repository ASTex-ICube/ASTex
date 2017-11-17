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
 * @brief The RGBA image class
 * @tparam CHANNEL_TYPE (uchar/char/ushort/short/.../float/double)
 */
template<typename CHANNEL_TYPE>
class ImageRGBABase : public ImageBase
{
public:
	using PixelType =            itk::RGBAPixel<CHANNEL_TYPE>;
	using ItkImg =               itk::Image< PixelType >;
	using IteratorIndexed =      itk::ImageRegionIteratorWithIndex<ItkImg>;
	using Iterator =             itk::ImageRegionIterator<ItkImg>;
	using ConstIteratorIndexed = itk::ImageRegionConstIteratorWithIndex<ItkImg> ;
	using ConstIterator =        itk::ImageRegionConstIterator<ItkImg>;
	using DoublePixelEigen =     Eigen::Vector4d;
	using DataType =             CHANNEL_TYPE;

	static const uint32_t NB_CHANNELS = 4;

protected:
	typename ItkImg::Pointer itk_img_;

public:
	ImageRGBABase():
		itk_img_(NULL)
	{}

	ImageRGBABase(typename itk::Image< itk::RGBAPixel<CHANNEL_TYPE> >::Pointer itk_im):
		itk_img_(itk_im)
	{}

	/**
	 * @brief itkPixel
	 * @return a pixel of value (r/g/b/a)
	 */
	inline static PixelType itkPixel(CHANNEL_TYPE r, CHANNEL_TYPE g, CHANNEL_TYPE b, CHANNEL_TYPE a)
	{
		return PixelType(std::array<CHANNEL_TYPE, 4>({ { r,g,b,a } }).data());
	}

	/**
	 * @brief itkPixel
	 * @return a pixel of value (r/r/r/r)
	 */
	inline static PixelType itkPixel(CHANNEL_TYPE r)
	{
		return PixelType(std::array<CHANNEL_TYPE, 4>({ { r,r,r,r } }).data());
	}

	/**
	 * @brief itkPixelNorm
	 * @param r normalized value [0,1]
	 * @param g normalized value [0,1]
	 * @param b normalized value [0,1]
	 * @param a normalized value [0,1]
	 * @return a pixel RGBA
	 */
	template <bool B = true>
	inline static auto itkPixelNorm(double r, double g, double b, double a) -> typename std::enable_if<B && std::is_arithmetic<CHANNEL_TYPE>::value, PixelType>::type
	{
		if (std::is_floating_point<CHANNEL_TYPE>::value)
			return itkPixel(CHANNEL_TYPE(r), CHANNEL_TYPE(g), CHANNEL_TYPE(b), CHANNEL_TYPE(a));

		if (std::is_unsigned<CHANNEL_TYPE>::value)
			return itkPixel(r * std::numeric_limits<CHANNEL_TYPE>::max(), g * std::numeric_limits<CHANNEL_TYPE>::max(), b * std::numeric_limits<CHANNEL_TYPE>::max(), a * std::numeric_limits<CHANNEL_TYPE>::max());

		return itkPixel(r*(double(std::numeric_limits<CHANNEL_TYPE>::max()) - double(std::numeric_limits<CHANNEL_TYPE>::lowest())) + std::numeric_limits<CHANNEL_TYPE>::lowest(),
			g*(double(std::numeric_limits<CHANNEL_TYPE>::max()) - double(std::numeric_limits<CHANNEL_TYPE>::lowest())) + std::numeric_limits<CHANNEL_TYPE>::lowest(),
			b*(double(std::numeric_limits<CHANNEL_TYPE>::max()) - double(std::numeric_limits<CHANNEL_TYPE>::lowest())) + std::numeric_limits<CHANNEL_TYPE>::lowest(),
			a*(double(std::numeric_limits<CHANNEL_TYPE>::max()) - double(std::numeric_limits<CHANNEL_TYPE>::lowest())) + std::numeric_limits<CHANNEL_TYPE>::lowest());
	}

	/**
	 * @brief eigenPixel
	 * @param p a itk::Pixel
	 * @return an Eigen::Vector4d with same value
	 */
	inline static DoublePixelEigen eigenPixel(const PixelType& p)
	{
		if (std::is_same<CHANNEL_TYPE,double>::value)
			return *(reinterpret_cast<const DoublePixelEigen*>(&p));
		return DoublePixelEigen(p[0],p[1],p[2],p[3]);
	}

	/**
	 * @brief itkPixel
	 * @param p an Eigen::Vector4d
	 * @return an itk::Pixel with same values of p
	 */
	inline static PixelType itkPixel(const DoublePixelEigen& p)
	{
		if (std::is_same<CHANNEL_TYPE,double>::value)
			return *(reinterpret_cast<const PixelType*>(&p));
		return itkPixel(CHANNEL_TYPE(p[0]),CHANNEL_TYPE(p[1]),CHANNEL_TYPE(p[2]),CHANNEL_TYPE(p[3]));
	}

	inline static DoublePixelEigen eigenPixel(double v)
	{
		return DoublePixelEigen(v,v,v,v);
	}

	inline static DoublePixelEigen normalized(const PixelType& p)
	{
		return DoublePixelEigen(ASTex::normalized(p[0]), ASTex::normalized(p[1]), ASTex::normalized(p[2]), ASTex::normalized(p[3]));
	}

	inline static PixelType unnormalized(const DoublePixelEigen& p)
	{
		return itkPixel(ASTex::unnormalized<DataType>(p[0]), ASTex::unnormalized<DataType>(p[1]), ASTex::unnormalized<DataType>(p[2]), ASTex::unnormalized<DataType>(p[2]));
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
inline itk::RGBAPixel<CHANNEL_TYPE> itkRGBAPixel(CHANNEL_TYPE r, CHANNEL_TYPE g, CHANNEL_TYPE b, CHANNEL_TYPE a)
{
	return ImageRGBABase<CHANNEL_TYPE>::itkPixel(r,g,b,a);
}

template<typename CHANNEL_TYPE>
inline itk::RGBAPixel<CHANNEL_TYPE> itkRGBAPixel(CHANNEL_TYPE r)
{
	return ImageRGBABase<CHANNEL_TYPE>::itkPixel(r);;
}


template<typename CHANNEL_TYPE>
inline itk::RGBAPixel<CHANNEL_TYPE> itkRGBAPixel(const Eigen::Vector4d& v)
{
	return ImageRGBABase<CHANNEL_TYPE>::itkPixel(v);
}


template<typename CHANNEL_TYPE>
inline Eigen::Vector4d eigenRGBAPixel(const itk::RGBAPixel<CHANNEL_TYPE>& p)
{
	return ImageRGBABase<CHANNEL_TYPE>::eigenPixel(p);
}



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



