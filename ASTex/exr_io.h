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



#include <iostream>
#include <iomanip>

#include <ASTex/image_rgba.h>
#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>



#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-register"
#pragma clang diagnostic ignored "-Wunused-parameter"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#elif defined(_MSC_VER)
#pragma warning( push)
#pragma warning ( disable : 4101 )
#endif


#include <ImfHeader.h>
#include <ImfFrameBuffer.h>
#include <ImfOutputFile.h>
#include <ImfInputFile.h>
#include <ImfChannelList.h>
#include <ImfStringAttribute.h>



#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#pragma warning( pop )
#endif



//#include "dll.h"

namespace ASTex
{
namespace IO
{
namespace EXR
{


/**
 * @brief save image in exr Y float 32bit
 * @param img image to save
 * @param out_name filename
 */
void ASTEX_API save(const ImageGrayf& img, const std::string& out_name);

/**
 * @brief save image in exr RGB float 32bit
 * @param img image to save
 * @param out_name filename
 */
void ASTEX_API save(const ImageRGBf& img, const std::string& out_name);

/**
 * @brief save image in exr RGBA float 32bit
 * @param img image to save
 * @param out_name filename
 */
void ASTEX_API save(const ImageRGBAf& img, const std::string& out_name);

/**
 * @brief save image in exr Y float 32bit
 * @param img image to save
 * @param out_name filename
 */
void ASTEX_API save(const ImageGrayd& img, const std::string& out_name);

/**
 * @brief save image in exr RGB float 32bit
 * @param img image to save
 * @param out_name filename
 */
void ASTEX_API save(const ImageRGBd& img, const std::string& out_name);

/**
 * @brief save image in exr RGBA float 32bit
 * @param img image to save
 * @param out_name filename
 */
void ASTEX_API save(const ImageRGBAd& img, const std::string& out_name);


/**
 * @brief load exr Y float 32bit
 * @param img output image loaded (only float/double components)
 * @param in_name file name
 * @return true of load ok
 */
template <typename T>
auto load(ImageGray<T>& img, const std::string in_name)-> typename std::enable_if<std::is_floating_point<T>::value, bool>::type
{
	Imf::InputFile file (in_name.c_str());
	const Imf::Channel *chY = file.header().channels().findChannel ("Y");
	if (chY == NULL)
	{
		std::cerr << "Error no Y canal" << std::endl;
		return false;
	}

	Imath::Box2i dw = file.header().dataWindow();
	int width = dw.max.x - dw.min.x + 1;
	int height = dw.max.y - dw.min.y + 1;

	float* pixels = new float[width*height];

	Imf::FrameBuffer frameBuffer;
	frameBuffer.insert ("Y", Imf::Slice (Imf::FLOAT,(char *) (pixels - dw.min.x - dw.min.y * width),
		sizeof (float) , sizeof (float) * width,	1, 1, FLT_MAX));

	file.setFrameBuffer (frameBuffer);
	file.readPixels (dw.min.y, dw.max.y);

	img.initItk(width,height);

	float* ptr = pixels;
	img.for_all_pixels([&] (typename ImageGray<T>::PixelType& P )
	{
		P = *ptr++;
	});

	delete[] pixels;
	return true;
}

/**
 * @brief load exr RGB float 32bit
 * @param img output image loaded (only float/double components)
 * @param in_name file name
 * @return true of load ok
 */
template < typename T>
auto load(ImageRGB<T>& img, const std::string in_name) -> typename std::enable_if<std::is_floating_point<T>::value, bool>::type
{
	Imf::InputFile file (in_name.c_str());
	const Imf::Channel *chR = file.header().channels().findChannel ("R");
	const Imf::Channel *chG = file.header().channels().findChannel ("G");
	const Imf::Channel *chB = file.header().channels().findChannel ("B");
	if ((chR == NULL)||(chG == NULL)||(chB == NULL))
	{
		std::cerr << "Missing RBG canals" << std::endl;
		return false;
	}

	Imath::Box2i dw = file.header().dataWindow();
	int width = dw.max.x - dw.min.x + 1;
	int height = dw.max.y - dw.min.y + 1;

	float* pixels = new float[3*width*height];

	Imf::FrameBuffer frameBuffer;
	const std::size_t PSZ = 3*sizeof(float);
	frameBuffer.insert ("R", Imf::Slice (Imf::FLOAT,(char *) (pixels - dw.min.x - dw.min.y * width),
		PSZ, PSZ*width,	1, 1, FLT_MAX));
	frameBuffer.insert ("G", Imf::Slice (Imf::FLOAT,(char *) (pixels + 1 - dw.min.x - dw.min.y * width),
		PSZ, PSZ*width,	1, 1, FLT_MAX));
	frameBuffer.insert ("B", Imf::Slice (Imf::FLOAT,(char *) (pixels + 2 - dw.min.x - dw.min.y * width),
		PSZ, PSZ*width,	1, 1, FLT_MAX));

	file.setFrameBuffer (frameBuffer);
	file.readPixels (dw.min.y, dw.max.y);

	img.initItk(width,height);

	float* ptr = pixels;
	img.for_all_pixels([&] (typename ImageRGB<T>::PixelType& P )
	{
		P[0] = *ptr++;
		P[1] = *ptr++;
		P[2] = *ptr++;
	});

	delete[] pixels;
	return true;
}

/**
 * @brief load exr RGBA float 32bit
 * @param img output image loaded (only float/double components)
 * @param in_name file name
 * @return true of load ok
 */
template < typename T>
auto load(ImageRGBA<T>& img, const std::string in_name) -> typename std::enable_if<std::is_floating_point<T>::value, bool>::type
{
	Imf::InputFile file (in_name.c_str());
	const Imf::Channel *chR = file.header().channels().findChannel ("R");
	const Imf::Channel *chG = file.header().channels().findChannel ("G");
	const Imf::Channel *chB = file.header().channels().findChannel ("B");
	const Imf::Channel *chA = file.header().channels().findChannel ("A");
	if ((chR == NULL)||(chG == NULL)||(chB == NULL)||(chA == NULL))
	{
		std::cerr << "Missing RBGA canals" << std::endl;
		return false;
	}

	Imath::Box2i dw = file.header().dataWindow();
	int width = dw.max.x - dw.min.x + 1;
	int height = dw.max.y - dw.min.y + 1;

	float* pixels = new float[4*width*height];

	Imf::FrameBuffer frameBuffer;
	const std::size_t PSZ = 4*sizeof(float);
	frameBuffer.insert ("R", Imf::Slice (Imf::FLOAT,(char *) (pixels - dw.min.x - dw.min.y * width),
		PSZ, PSZ*width,	1, 1, FLT_MAX));
	frameBuffer.insert ("G", Imf::Slice (Imf::FLOAT,(char *) (pixels + 1 - dw.min.x - dw.min.y * width),
		PSZ, PSZ*width,	1, 1, FLT_MAX));
	frameBuffer.insert ("B", Imf::Slice (Imf::FLOAT,(char *) (pixels + 2 - dw.min.x - dw.min.y * width),
		PSZ, PSZ*width,	1, 1, FLT_MAX));
	frameBuffer.insert ("A", Imf::Slice (Imf::FLOAT,(char *) (pixels + 3 - dw.min.x - dw.min.y * width),
		PSZ, PSZ*width,	1, 1, FLT_MAX));

	file.setFrameBuffer (frameBuffer);
	file.readPixels (dw.min.y, dw.max.y);

	img.initItk(width,height);

	float* ptr = pixels;
	img.for_all_pixels([&] (typename ImageRGBA<T>::PixelType& P )
	{
		P[0] = *ptr++;
		P[1] = *ptr++;
		P[2] = *ptr++;
		P[3] = *ptr++;
	});

	delete[] pixels;
	return true;
}

}
}
}



