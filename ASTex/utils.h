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



#ifndef __ASTEX_UTILS__
#define __ASTEX_UTILS__

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <cstring>

#ifdef WIN32
#include<windows.h>
#endif

namespace ASTex
{


//#ifdef WIN32
//inline void create_directory(const std::string& path)
//{
//	CreateDirectory( path.c_str(), NULL);
//}
//#else
//inline void create_directory(const std::string& path)
//{
//	system(std::string("mkdir -p "+path).c_str());
//}
//#endif

#ifdef WIN32
inline bool create_directory(const std::string& path)
{
	return CreateDirectory( path.c_str(), NULL) == 0;
}
#else
inline bool create_directory(const std::string& path)
{
	if (system(std::string("test -d "+path).c_str()) == 0)
		return false;
	system(std::string("mkdir -p "+path).c_str());
	return true;
}
#endif


template <typename SCALAR, typename SCALAR2>
SCALAR clamp_scalar(const SCALAR n, const SCALAR2 lower, const SCALAR2 upper)
{
  return std::max(SCALAR(lower), std::min(n, SCALAR(upper)));
}

////Random in [A,B]
template <typename SCALAR,typename SCALAR2>
inline SCALAR clampedRandom (const SCALAR sa, const SCALAR2 sb)
{
	double a(sa);
	double b(sb);
	double tmp = double(std::rand())/RAND_MAX; // [0,1]
	double r = a+tmp*(b-a);
	return SCALAR(r);
}



template <typename SCALAR>
inline bool compare_scalar(const SCALAR i, const SCALAR j, const SCALAR EPSILON = 0.0001)
{
	return std::abs(j-i) < EPSILON;
}


template <typename IMG>
void crop_image(const IMG& input, IMG& output, int i_min, int i_max, int j_min, int j_max)
{
	assert(0<=i_min && i_min < i_max && i_max < input.width() && 0<=j_min && j_min < j_max && j_max < input.height());

	for(int x = i_min; x <= i_max; x++)
	{
		for(int y = j_min; y <= j_max; y++)
		{
			output.pixelAbsolute(x-i_min,y-j_min)=input.pixelAbsolute(x,y);
		}
	}
}

template <typename IMG>
int fulfill_crop_image(const IMG& crop, IMG& output, int i_min, int j_min)
{
	assert(0 <= i_min && i_min <= output.width()-crop.width() && 0 <= j_min && j_min <= output.height() - crop.height());

	for(int x = 0; x < crop.width(); x++)
	{
		for(int y = 0; y <crop.height(); y++)
		{
			output.pixelAbsolute( i_min+x , j_min+y ) = crop.pixelAbsolute(x,y);
		}
	}
	return 0;
}




template <typename ISRC, typename IDST>
void distance_invert_expo (const ISRC& src, IDST& dst, double dev)
{

	dst.for_all_pixels([&](typename IDST::PixelType& p, int x, int y)
	{
		p = std::exp(-src.pixelAbsolute(x,y)/dev);
	});
}

template <typename ISRC, typename IDST>
void image_log (const ISRC& src, IDST& dst, double ratio)
{
	dst.for_all_pixels([&](typename IDST::PixelType& p, int x, int y)
	{
		p = std::sqrt((src.pixelAbsolute(x,y)*ratio ));
	});
}

template <typename ISRC, typename IDST>
void image_expo (const ISRC& src, IDST& dst, double dev)
{
	dst.for_all_pixels([&](typename IDST::PixelType& p, int x, int y)
	{
		p = 1.0 - std::exp(-src.pixelAbsolute(x,y)/dev);
	});
}

template <typename IMG, typename MASK>
typename IMG::ASTexPixelType compute_mean(const IMG& img, const MASK& mask)
{
	using Pix = typename IMG::ASTexPixelType;

	Pix mean(0);
	long int nbpix = 0;

	img.for_all_pixels([&] (const Pix& p)
	{
		mean += p;
		nbpix++;
	}, mask);

	mean /= nbpix;
	return mean;
}

template <typename IMG>
typename IMG::ASTexPixelType compute_mean(const IMG& img)
{
	return compute_mean(img,[](int,int){return true;});
}





template <typename IMG, typename MASK, typename SUP>
auto compute_max(const IMG& img, const MASK& mask, const SUP& sup) ->
		typename IMG::ASTexPixelType
{
	using Pix = typename IMG::ASTexPixelType;

	Pix max_val = img.pixelAbsolute(0,0);
	img.for_all_pixels([&](const Pix& p)
	{
		if (sup(p,max_val))
			 max_val = p;
	},mask);

	return max_val;
}


template <typename IMG, typename MASK>
auto compute_max(const IMG& img, const MASK& mask) ->
		typename std::enable_if<std::is_arithmetic<typename IMG::ASTexPixelType>::value, typename IMG::ASTexPixelType>::type
{
	using Pix = typename IMG::ASTexPixelType;
	return compute_max(img, mask, [](const Pix& p, const Pix& q){return p>q;});
}

template <typename IMG>
typename IMG::ASTexPixelType compute_max(const IMG& img)
{
	using Pix = typename IMG::ASTexPixelType;
	return compute_max(img,
					[](int,int){return true;},
					[](const Pix& p, const Pix& q){return p>q;});
}





template <typename IMG, typename MASK, typename SUP>
auto compute_min(const IMG& img, const MASK& mask, const SUP& inf) ->
		typename IMG::ASTexPixelType
{
	using Pix = typename IMG::ASTexPixelType;

	Pix min_val = img.pixelAbsolute(0,0);
	img.for_all_pixels([&](const Pix& p)
	{
		if (inf(p,min_val))
			 min_val = p;
	},mask);

	return min_val;
}


template <typename IMG, typename MASK>
auto compute_min(const IMG& img, const MASK& mask) ->
		typename std::enable_if<std::is_arithmetic<typename IMG::ASTexPixelType>::value, typename IMG::ASTexPixelType>::type
{
	using Pix = typename IMG::ASTexPixelType;
	return compute_min(img, mask, [](const Pix& p, const Pix& q){return p<q;});
}

template <typename IMG>
typename IMG::ASTexPixelType compute_min(const IMG& img)
{
	using Pix = typename IMG::ASTexPixelType;
	return compute_min(img,
					[](int,int){return true;},
					[](const Pix& p, const Pix& q){return p<q;});
}

template <typename I>
double mse(const I& i1, const I& i2, int x1, int y1, int x2, int y2, int neighborhoodRadius, int periodicity=false)
{
    double error=0.0;
    static typename I::PixelType ms_zero;
    typename I::PixelType diff=ms_zero;
    size_t arraySize = sizeof(typename I::PixelType) / sizeof(typename I::DataType);
    typename I::DataType    *a_pixi1 = new typename I::DataType[arraySize],
                            *a_pixi2 = new typename I::DataType[arraySize];

    unsigned hit=0;
    for(int dy=-neighborhoodRadius; dy<=neighborhoodRadius; ++dy)
        for(int dx=-neighborhoodRadius; dx<=neighborhoodRadius; ++dx)
        {
            int xx1 = periodicity ? (x1+dx)%i1.width() : x1+dx;
            int xx2 = periodicity ? (x2+dx)%i2.width() : x2+dx;
            int yy1 = periodicity ? (y1+dy)%i1.height() : y1+dy;
            int yy2 = periodicity ? (y2+dy)%i2.height() : y2+dy;

            if(xx1 >= 0 && xx1<i1.width() &&
               xx2 >= 0 && xx2<i2.width() &&
               yy1 >= 0 && yy1<i1.height() &&
               yy2 >= 0 && yy1<i2.height())
            {
                typename I::PixelType pixi1=i1.pixelAbsolute(xx1, yy1), pixi2=i2.pixelAbsolute(xx2, yy2);
                assert(sizeof(typename I::DataType) <= 64
                       && "mse: using Image types of data types with sizes higher than 64 bits is not allowed");
                std::memcpy(a_pixi1, &pixi1, sizeof(typename I::PixelType));
                std::memcpy(a_pixi2, &pixi2, sizeof(typename I::PixelType));
                if(!std::is_floating_point<typename I::DataType>::value)
                {
                    for(unsigned i=0; i<arraySize; ++i)
                    {
                        error += (int64_t(a_pixi1[i]) - a_pixi2[i]) * (int64_t(a_pixi1[i]) - a_pixi2[i]);
                    }
                }
                else
                {
                    for(unsigned i=0; i<arraySize; ++i)
                    {
                        error += (double(a_pixi1[i]) - a_pixi2[i]) * (double(a_pixi1[i]) - a_pixi2[i]);
                    }
                }
                ++hit;
            }
        }
    assert(hit!=0 && "mse: error with parameters");
    delete[](a_pixi1);
    delete[](a_pixi2);
    return error/arraySize/hit;
}

}

#endif
