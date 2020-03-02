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

#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>
#include <ASTex/image_rgba.h>

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
void nullify_image(IMG& image)
{
	static typename IMG::PixelType pix_zero;
	image.for_all_pixels([&] (typename IMG::PixelType &pix)
	{
		pix = pix_zero;
	});
}


template <typename IMG>
void crop_image(const IMG& input, IMG& output, int i_min, int i_max, int j_min, int j_max)
{
	static typename IMG::PixelType pix_zero;
	assert(0<=i_min && i_min < i_max && i_max < input.width() && 0<=j_min && j_min < j_max && j_max < input.height());
	output.initItk(i_max-i_min+1, j_max-j_min+1, true);
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
typename IMG::PixelType compute_mean(const IMG& img, const MASK& mask)
{
	using Pix = typename IMG::DoublePixelEigen;

	Pix mean(0);
	long int nbpix = 0;

	img.for_all_pixels([&] (const Pix& p)
	{
		mean += eigenPixel<double>(p);
		nbpix++;
	}, mask);

	mean /= nbpix;
	return IMG::itkPixel(mean);
}

template <typename IMG>
typename IMG::PixelType compute_mean(const IMG& img)
{
	return compute_mean(img,[](int,int){return true;});
}





template <typename IMG, typename MASK, typename SUP>
auto compute_max(const IMG& img, const MASK& mask, const SUP& sup) ->
		typename IMG::PixelType
{
	using Pix = typename IMG::PixelType;

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
		typename std::enable_if<std::is_arithmetic<typename IMG::PixelType>::value, typename IMG::PixelType>::type
{
	using Pix = typename IMG::PixelType;
	return compute_max(img, mask, [](const Pix& p, const Pix& q){return p>q;});
}

template <typename IMG>
typename IMG::PixelType compute_max(const IMG& img)
{
	using Pix = typename IMG::PixelType;
	return compute_max(img,
					[](int,int){return true;},
					[](const Pix& p, const Pix& q){return p>q;});
}





template <typename IMG, typename MASK, typename SUP>
auto compute_min(const IMG& img, const MASK& mask, const SUP& inf) ->
		typename IMG::PixelType
{
	using Pix = typename IMG::PixelType;

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
	using Pix = typename IMG::PixelType;
	return compute_min(img, mask, [](const Pix& p, const Pix& q){return p<q;});
}

template <typename IMG>
typename IMG::PixelType compute_min(const IMG& img)
{
	using Pix = typename IMG::PixelType;
	return compute_min(img,
					[](int,int){return true;},
					[](const Pix& p, const Pix& q){return p<q;});
}

template <typename I>
double mse(const I& i1, const I& i2, int x1, int y1, int x2, int y2, int neighborhoodRadius=0, bool periodicity=false)
{
	double error=0.0;
	static typename I::PixelType ms_zero; //A way to always get a zero-valued pixel
	size_t arraySize = sizeof(typename I::PixelType) / sizeof(typename I::DataType); //Finding out the nb of channels
	typename I::DataType    *a_pixi1 = new typename I::DataType[arraySize],
							*a_pixi2 = new typename I::DataType[arraySize]; //Creating a dummy pixel you can use [] on

	unsigned hit=0;
	for(int dy=-neighborhoodRadius; dy<=neighborhoodRadius; ++dy)
		for(int dx=-neighborhoodRadius; dx<=neighborhoodRadius; ++dx)
		{
			int xx1 = periodicity ? (x1+dx+4*i1.width())%i1.width() : x1+dx;
			int xx2 = periodicity ? (x2+dx+4*i2.width())%i2.width() : x2+dx;
			int yy1 = periodicity ? (y1+dy+4*i1.height())%i1.height() : y1+dy;
			int yy2 = periodicity ? (y2+dy+4*i2.height())%i2.height() : y2+dy;

			if(xx1 >= 0 && xx1<i1.width() &&
			   xx2 >= 0 && xx2<i2.width() &&
			   yy1 >= 0 && yy1<i1.height() &&
			   yy2 >= 0 && yy2<i2.height())
			{
				typename I::PixelType pixi1=i1.pixelAbsolute(xx1, yy1), pixi2=i2.pixelAbsolute(xx2, yy2);
				assert(sizeof(typename I::DataType) <= 64
					   && "mse: using Image types of data types with sizes higher than 64 bits is not allowed");
				std::memcpy(a_pixi1, &pixi1, sizeof(typename I::PixelType)); //Filling the dummy pixel
				std::memcpy(a_pixi2, &pixi2, sizeof(typename I::PixelType));
				if(!std::is_floating_point<typename I::DataType>::value)
				{
					for(unsigned i=0; i<arraySize; ++i)
					{//computing the error per channel
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
	if(hit==0)
		return std::numeric_limits<double>::infinity();
	delete[](a_pixi1);
	delete[](a_pixi2);
	return error/arraySize/hit; //It returns a double but you may want an array,
								//so either extract the channels from your image or make another mse
}

template<typename I>
void printImage(const I &image, std::ostream &out_stream, bool show_position = false)
{
	image.for_all_pixels([&] (const typename I::PixelType &pix, int x, int y)
	{
		if(show_position)
			out_stream << "(" << x << ", " << y << "): " << pix << std::endl;
		else
			out_stream << pix << std::endl;
	});
	return;
}

template<typename I>
typename I::PixelType bilinear_interpolation(const I& image,
											 double x, double y, bool periodicity) //Peut-être rendre plus lisible ?
{
	size_t pixelSize = sizeof(typename I::PixelType)/sizeof(typename I::DataType);
	typename I::DataType *pixelArray = new typename I::DataType[pixelSize];
	auto turnPixelIntoArithmeticArray = [&] (const typename I::PixelType &pix) -> double*
	{
		double *arithmeticArray = new double[pixelSize];
		std::memcpy(pixelArray, &pix, pixelSize*sizeof(typename I::DataType));
		for(unsigned i=0; i<pixelSize; ++i)
		{
			arithmeticArray[i] = double(pixelArray[i]);
		}
		return arithmeticArray;
	};
	auto turnArithmeticArrayIntoPixel = [&] (double *arithmeticArray) -> typename I::PixelType
	{
		typename I::PixelType pix;
		for(unsigned i=0; i<pixelSize; ++i)
		{
			pixelArray[i] = typename I::DataType(arithmeticArray[i]);
		}
		std::memcpy(&pix, pixelArray, pixelSize*sizeof(typename I::DataType));
		return pix;
	};
	typename I::PixelType outputValue;
	if(!periodicity)
	{
		if(x == image.width()-1)
		{
			if(y == image.height()-1)
				outputValue=image.pixelAbsolute(x, y);
			else
			{
				double *q0 = turnPixelIntoArithmeticArray(image.pixelAbsolute(int(x), int(y)));
				double *q1 = turnPixelIntoArithmeticArray(image.pixelAbsolute(int(x), int(y+1)));
				double *qInterp = new double[pixelSize];
				double y_y1 = y - int(y);
				for(unsigned i=0; i<pixelSize; ++i)
				{
					qInterp[i] = (1.0-y_y1) * q0[i] + y_y1*q1[i];
				}
				outputValue=turnArithmeticArrayIntoPixel(qInterp);
				delete[] qInterp;
				delete[] q0;
				delete[] q1;
			}
		}
		else if(y == image.height()-1)
		{
			double *q0 = turnPixelIntoArithmeticArray(image.pixelAbsolute(int(x), int(y)));
			double *q1 = turnPixelIntoArithmeticArray(image.pixelAbsolute(int(x+1), int(y)));
			double *qInterp = new double[pixelSize];
			double x_x1 = x - int(x);
			for(unsigned i=0; i<pixelSize; ++i)
			{
				qInterp[i] = (1.0-x_x1) * q0[i] + x_x1*q1[i];
			}
			outputValue=turnArithmeticArrayIntoPixel(qInterp);
			delete[] qInterp;
			delete[] q0;
			delete[] q1;
		}
		else
		{
			double *q00 = turnPixelIntoArithmeticArray(image.pixelAbsolute(int(x), int(y)));
			double *q10 = turnPixelIntoArithmeticArray(image.pixelAbsolute(int(x+1), int(y)));
			double *q01 = turnPixelIntoArithmeticArray(image.pixelAbsolute(int(x), int(y+1)));
			double *q11 = turnPixelIntoArithmeticArray(image.pixelAbsolute(int(x+1), int(y+1)));
			double *qInterp = new double[pixelSize];
			double x_x1 = x - int(x);
			double y_y1 = y - int(y);
			for(unsigned i=0; i<pixelSize; ++i)
			{
				qInterp[i] = (1.0-x_x1) * (1.0-y_y1) * q00[i]
							+ x_x1	    * (1.0-y_y1) * q10[i]
							+(1.0-x_x1) * y_y1		 * q01[i]
							+x_x1		* y_y1		 * q11[i];
			}
			outputValue=turnArithmeticArrayIntoPixel(qInterp);
			delete[] qInterp;
			delete[] q00;
			delete[] q10;
			delete[] q01;
			delete[] q11;
		}
	}
	else
	{
		double *q00 = turnPixelIntoArithmeticArray(image.pixelAbsolute(int(x)%image.width(), int(y)%image.height()));
		double *q10 = turnPixelIntoArithmeticArray(image.pixelAbsolute(int(x+1)%image.width(), int(y)%image.height()));
		double *q01 = turnPixelIntoArithmeticArray(image.pixelAbsolute(int(x)%image.width(), int(y+1)%image.height()));
		double *q11 = turnPixelIntoArithmeticArray(image.pixelAbsolute(int(x+1)%image.width(), int(y+1)%image.height()));
		double *qInterp = new double[pixelSize];
		double x_x1 = x - int(x);
		double y_y1 = y - int(y);
		for(unsigned i=0; i<pixelSize; ++i)
		{
			qInterp[i] = (1.0-x_x1) * (1.0-y_y1) * q00[i]
						+ x_x1	    * (1.0-y_y1) * q10[i]
						+(1.0-x_x1) * y_y1		 * q01[i]
						+x_x1		* y_y1		 * q11[i];
		}
		outputValue=turnArithmeticArrayIntoPixel(qInterp);
		delete[] qInterp;
		delete[] q00;
		delete[] q10;
		delete[] q01;
		delete[] q11;
	}
	delete[] pixelArray;
	return outputValue;
}

template<typename I, typename MASK_TYPE, typename std::enable_if<std::is_floating_point<MASK_TYPE>::value>::type> //I need masks able to
typename I::PixelType normalize(I &image, const ImageGray<MASK_TYPE> &mask)
{
	using PixelType = typename I::PixelType;
	double pixelCount = 0;
	PixelType norm = I::zero();
	image.for_all_pixels([&] (PixelType &pix, int x, int y)
	{
		double maskPix = double(mask.pixelAbsolute(x, y));
		if(maskPix>0)
		{
			pixelCount += maskPix;
			norm += pix * maskPix;
		}
	});
	image.for_all_pixels([&] (PixelType &pix, int x, int y)
	{
		double maskPix = double(mask.pixelAbsolute(x, y));
		if(maskPix>0)
		{
			pixelCount += maskPix;
			norm += pix * maskPix;
		}
	});
}

}

#endif
