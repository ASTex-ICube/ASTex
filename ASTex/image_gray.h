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



#ifndef __ASTEX_IMAGE_GRAY__
#define __ASTEX_IMAGE_GRAY__

#include "itkImage.h"
#include <ASTex/image_common.h>

namespace ASTex
{

template<typename CHANNEL_TYPE>
class ImageGrayBase : public ImageBase
{
public:
	using PixelType =            CHANNEL_TYPE;
	using ItkImg =               itk::Image< PixelType >;
	using IteratorIndexed =      itk::ImageRegionIteratorWithIndex<ItkImg>;
	using Iterator =             itk::ImageRegionIterator<ItkImg>;
	using ConstIteratorIndexed = itk::ImageRegionConstIteratorWithIndex<ItkImg> ;
	using ConstIterator =        itk::ImageRegionConstIterator<ItkImg>;
	using DoublePixelEigen =     double;
	using LongPixelEigen =       int64_t;
	using DataType =             CHANNEL_TYPE;
	template<typename S>
	using EigenVector =          S;

	static const uint32_t NB_CHANNELS = 1;

protected:
	typename ItkImg::Pointer itk_img_;

public:
	ImageGrayBase():
		itk_img_(NULL)
	{}

	ImageGrayBase(typename itk::Image< CHANNEL_TYPE >::Pointer itk_im):
		itk_img_(itk_im)
	{}

	/**
	 * @brief itkPixel
	 * @return a pixel of value (v)
	 */
	inline static PixelType itkPixel(CHANNEL_TYPE v)
	{
		return PixelType(v);
	}

	/**
	 * @brief itkPixelNorm
	 * @param v normalized value [0,1]
	 * @return a pixel V
	 */
	template <bool B=true>
	inline static auto itkPixelNorm(double v)->typename std::enable_if<B && std::is_arithmetic<CHANNEL_TYPE>::value,PixelType>::type
	{
		assert((v>=0.0)&&(v<=1.0));
		if (std::is_floating_point<CHANNEL_TYPE>::value)
			return PixelType(v);

		if (std::is_unsigned<CHANNEL_TYPE>::value)
			return PixelType(v * std::numeric_limits<CHANNEL_TYPE>::max());

		return PixelType(v*(double(std::numeric_limits<CHANNEL_TYPE>::max()) - double(std::numeric_limits<CHANNEL_TYPE>::lowest())) + std::numeric_limits<CHANNEL_TYPE>::lowest());
	}

	/**
	 * @brief eigenDoublePixel
	 * @param p a itk::Pixel
	 * @return a double with same value
	 */
	inline static DoublePixelEigen eigenDoublePixel(const PixelType& p)
	{
		return DoublePixelEigen(p);
	}

	inline static LongPixelEigen eigenLongPixel(const PixelType& p)
	{
		return LongPixelEigen(p);
	}

	template<typename S>
	inline static S eigenPixel(const PixelType& p)
	{
		if (std::is_same<PixelType,S>::value)
			return *(reinterpret_cast<const S*>(&p));
		return S(p);
	}



	inline static DoublePixelEigen normalized(const PixelType& p)
	{
		return normalized(p);
	}

	inline static PixelType unnormalized(const DoublePixelEigen& p)
	{
		return ASTex::unnormalized<DataType>(p);
	}

protected:

	inline PixelType* getPixelsPtr()
	{
		return this->itk_img_->GetPixelContainer()->GetBufferPointer();
	}

	inline DataType* getDataPtr()
	{
		return getPixelsPtr();
	}

	inline const PixelType* getPixelsPtr() const
	{
		return this->itk_img_->GetPixelContainer()->GetBufferPointer();
	}

	inline const DataType* getDataPtr() const
	{
		return getPixelsPtr();
	}
};

template<typename CHANNEL_TYPE>
inline CHANNEL_TYPE itkPixel(const double& v)
{
	return CHANNEL_TYPE(v);
}


template<typename S, typename CHANNEL_TYPE>
inline typename ImageGrayBase<CHANNEL_TYPE>::template EigenVector<S> eigenPixel(const CHANNEL_TYPE& p)
{
	return ImageGrayBase<CHANNEL_TYPE>::template eigenPixel<S>(p);
}



//template<typename CHANNEL_TYPE>
//inline double eigenDoublePixel(const CHANNEL_TYPE& p)
//{
//	return doubl(p);
//}

//template<typename CHANNEL_TYPE>
//inline auto eigenLongPixel(const CHANNEL_TYPE& p) -> typename std::enable_if<!std::is_floating_point<CHANNEL_TYPE>::value,int64_t>::type
//{
//	return int64_t(p);
//}



template <class T> using ImageGray = ImageCommon< ImageGrayBase< T >, false >;
using ImageGrayu8  = ImageGray< uint8_t >;
using ImageGray8   = ImageGray< int8_t >;
using ImageGrayu16 = ImageGray< uint16_t >;
using ImageGray16  = ImageGray< int16_t >;
using ImageGrayu32 = ImageGray< uint32_t >;
using ImageGray32  = ImageGray< int32_t >;
using ImageGrayu64 = ImageGray< uint64_t >;
using ImageGray64  = ImageGray< int64_t >;
using ImageGrayf   = ImageGray< float >;
using ImageGrayd   = ImageGray< double >;
using ImageGraycf  = ImageGray< std::complex<float> >;
using ImageGraycd  = ImageGray< std::complex<double> >;


// ConstImages are Images that only support const operations, they are not modifiable.
// But they could be construct from Itk::Image::ConstPointer -> usefull in Filter writing.

template <class T> using ConstImageGray = ImageCommon< ImageGrayBase< T >, true >;
using ConstImageGrayu8  = ConstImageGray< uint8_t >;
using ConstImageGray8   = ConstImageGray< int8_t >;
using ConstImageGrayu16 = ConstImageGray< uint16_t >;
using ConstImageGray16  = ConstImageGray< int16_t >;
using ConstImageGrayu32 = ConstImageGray< uint32_t >;
using ConstImageGray32  = ConstImageGray< int32_t >;
using ConstImageGrayu64 = ConstImageGray< uint64_t >;
using ConstImageGray64  = ConstImageGray< int64_t >;
using ConstImageGrayf   = ConstImageGray< float >;
using ConstImageGrayd   = ConstImageGray< double >;
using ConstImageGraycf  = ConstImageGray< std::complex<float> >;
using ConstImageGraycd  = ConstImageGray< std::complex<double> >;


class ImagePixelCompactIndex : public ImageGrayu32
{
public:
	inline static uint32_t pos_to_ui32(const Index& p) { return p[0] | p[1]<<16;}
	inline static  Index ui32_to_pos(uint32_t idx) { return gen_index(idx & 0xffff,idx >> 16);}
	inline Index iget(int x, int y) { return ui32_to_pos(this->pixelAbsolute(x,y));	}
	inline void iset(int x, int y, int a, int b) {	this->pixelAbsolute(x,y) = a | b<<16;}
	inline void iset(int x, int y, const Index& p)	{ this->pixelAbsolute(x,y) = pos_to_ui32(p);}
};

}

#endif



