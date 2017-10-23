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



