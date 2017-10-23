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

#include "exr_io.h"

namespace ASTex
{
namespace IO
{
namespace EXR
{

void ASTEX_API save_data_yf(int width, int height, const float* pixels, const std::string& out_name)
{
	Imf::Header header (width, height);
	header.insert ("comments", Imf::StringAttribute ("Gray ASTex image saved (in Y channel)"));
	header.channels().insert ("Y", Imf::Channel(Imf::FLOAT));

	Imf::OutputFile file (out_name.c_str(), header);
	Imf::FrameBuffer frameBuffer;

	const std::size_t PSZ = sizeof(float);
	frameBuffer.insert ("Y", Imf::Slice(Imf::FLOAT, (char *) pixels, PSZ, PSZ * width));
	file.setFrameBuffer (frameBuffer);
	file.writePixels (height);
}

void ASTEX_API save_data_rgbf(int width, int height, const float* pixels, const std::string& out_name)
{
	Imf::Header header (width, height);
	header.insert ("comments", Imf::StringAttribute ("RGB ASTex image saved"));
	header.channels().insert ("R", Imf::Channel(Imf::FLOAT));
	header.channels().insert ("G", Imf::Channel(Imf::FLOAT));
	header.channels().insert ("B", Imf::Channel(Imf::FLOAT));

	Imf::OutputFile file (out_name.c_str(), header);
	Imf::FrameBuffer frameBuffer;

	const std::size_t PSZ = 3*sizeof(float);
	frameBuffer.insert ("R", Imf::Slice(Imf::FLOAT, (char *) pixels, PSZ, PSZ * width));
	frameBuffer.insert ("G", Imf::Slice(Imf::FLOAT, (char *) (pixels+1), PSZ, PSZ * width));
	frameBuffer.insert ("B", Imf::Slice(Imf::FLOAT, (char *) (pixels+2), PSZ, PSZ * width));
	file.setFrameBuffer (frameBuffer);
	file.writePixels (height);
}

void ASTEX_API save_data_rgbaf(int width, int height, const float* pixels, const std::string& out_name)
{
	Imf::Header header (width, height);
	header.insert ("comments", Imf::StringAttribute ("RGBA ASTex image saved"));
	header.channels().insert ("R", Imf::Channel(Imf::FLOAT));
	header.channels().insert ("G", Imf::Channel(Imf::FLOAT));
	header.channels().insert ("B", Imf::Channel(Imf::FLOAT));
	header.channels().insert ("A", Imf::Channel(Imf::FLOAT));

	Imf::OutputFile file (out_name.c_str(), header);
	Imf::FrameBuffer frameBuffer;

	const std::size_t PSZ = 4*sizeof(float);
	frameBuffer.insert ("R", Imf::Slice(Imf::FLOAT, (char *) pixels, PSZ, PSZ * width));
	frameBuffer.insert ("G", Imf::Slice(Imf::FLOAT, (char *) (pixels+1), PSZ, PSZ * width));
	frameBuffer.insert ("B", Imf::Slice(Imf::FLOAT, (char *) (pixels+2), PSZ, PSZ * width));
	frameBuffer.insert ("A", Imf::Slice(Imf::FLOAT, (char *) (pixels+3), PSZ, PSZ * width));
	file.setFrameBuffer (frameBuffer);
	file.writePixels (height);
}


void ASTEX_API save(const ImageGrayf& img, const std::string& out_name)
{
	int width = img.width();
	int height = img.height();

	const float* ptr = &(img.pixelAbsolute(0,0));
	save_data_yf(width,height,ptr,out_name);
}


void ASTEX_API save(const ImageRGBf& img, const std::string& out_name)
{
	int width = img.width();
	int height = img.height();
	const float* ptr = &(img.pixelAbsolute(0,0)[0]);
	save_data_rgbf(width,height,ptr,out_name);

}


void ASTEX_API save(const ImageRGBAf& img, const std::string& out_name)
{
	int width = img.width();
	int height = img.height();
	const float* ptr = &(img.pixelAbsolute(0,0)[0]);
	save_data_rgbaf(width,height,ptr,out_name);
}



void ASTEX_API save(const ImageGrayd& img, const std::string& out_name)
{
	int width = img.width();
	int height = img.height();

	float* pixels = new float[width*height];
	float* ptr= pixels;
	img.for_all_pixels([&] (const ImageGrayd::PixelType& P )
	{
		*ptr++ = float(P);
	});

	save_data_yf(width,height,pixels,out_name);
	delete[] pixels;
}


void ASTEX_API save(const ImageRGBd& img, const std::string& out_name)
{
	int width = img.width();
	int height = img.height();

	float* pixels = new float[3*width*height];
	float* ptr = pixels;
	img.for_all_pixels([&] (const ImageRGBd::PixelType& P )
	{
		*ptr++ = float(P[0]);
		*ptr++ = float(P[1]);
		*ptr++ = float(P[2]);
	});

	save_data_rgbf(width,height,pixels,out_name);
	delete[] pixels;
}


void ASTEX_API save(const ImageRGBAd& img, const std::string& out_name)
{
	int width = img.width();
	int height = img.height();

	float* pixels = new float[4*width*height];
	float* ptr = pixels;
	img.for_all_pixels([&] (const ImageRGBAd::PixelType& P )
	{
		*ptr++ = float(P[0]);
		*ptr++ = float(P[1]);
		*ptr++ = float(P[2]);
		*ptr++ = float(P[3]);
	});

	save_data_rgbaf(width,height,pixels,out_name);
	delete[] pixels;
}


}
}
}



