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
#include <ASTex/image_rgb.h>
#include <ASTex/image_rgba.h>


#include <ASTex/special_io.h>
#include <ASTex/easy_io.h>

//#include <ASTex/image_merging.h>

using namespace ASTex;


int main()
{
	ImageRGBu8 image;

	bool ok = image.load(TEMPO_PATH+"simpleRGB.png");
	if (!ok)
		return 1;


	image.pixelAbsolute(0,0) = itkRGBPixel(255,127,52);
	image.pixelAbsolute(1,0) = itkRGBPixel(20,104,51);

	//ImageRGBu8::DoublePixelEigen dp1 = image.pixelEigenAbsolute(0,0);
	//ImageRGBu8::DoublePixelEigen dp2 = image.pixelEigenAbsolute(1,0);
//	ImageRGBu8::DoublePixelEigen dp = (image.pixelEigenAbsolute(0, 0) + image.pixelEigenAbsolute(1, 0))/2;

	image.pixelEigenAbsoluteWrite(0, 1) = (image.pixelEigenAbsolute(0, 0) + image.pixelEigenAbsolute(1, 0)) / 2;


	std::cout << image.pixelAbsolute(0,1) << std::endl;

	return 0;



//	ImageRGBu8 imx;

//	auto hm = Assembler1D::into(imx);
//	hm <<image << 3 << image << gen_region(50,50,150,150) << Assembler1D::HorizontalFlush;
//	imx.save(TEMPO_PATH+"h2simple.png");

//	Assembler1D::into(imx) <<image << 3 << image << gen_region(50,50,200,200) << 2 << image << gen_region(0,0,150,150) << Assembler1D::VerticalFlush;
//	imx.save(TEMPO_PATH+"v3simple.png");

//	Assembler2D::into(imx) <<image << 1 << image << Assembler2D::EndLine(2)<<image << 1 << image << Assembler2D::FinalFlush;
//	imx.save(TEMPO_PATH+"4simple.png");

	int W = image.width()/4;
	int H = image.height()/4;

	image.setCenter(image.width()/2, image.height()/2);

	for( int j= 0; j<H  ; ++j)
	{
		for( int i= 0; i< W ; ++i)
		{
			if (j%2)
				image.pixelRelative(i,j) = ImageRGBu8::itkPixel(255,128,0);
			else
				image.pixelRelative(i,j)[2]=0;
		}
	}


	for( int j= 0; j<H  ; ++j)
	{
		for( int i= 0; i< W ; ++i)
		{

			ImageRGBu8::PixelType& P = image.pixelAbsolute(i,j);
			P[0] = 0;
			P[1] = 30;

			auto P2 = image.pixelAbsolute(i,j);
			P2[2] = 200;

		}
	}

	image.save(TEMPO_PATH+"out.png");

	ImageRGBAd im(512,512);
	im.for_all_pixels([](ImageRGBAd::PixelType& P, int /*x*/, int y)
	{
		P = ImageRGBAd::itkPixel(1.0,0.3,0,(511-y%255)/255.0);
	});
	IO::save01_in_u8(im,TEMPO_PATH + "out_rgba.png");



  return EXIT_SUCCESS;
}

