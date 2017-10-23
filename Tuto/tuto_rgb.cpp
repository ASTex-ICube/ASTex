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

using namespace ASTex;


int main()
{
	ImageRGBu8 image;

	bool ok = image.load(TEMPO_PATH+"simpleRGB.png");
	if (!ok)
		return 1;


	int W = image.width()/4;
	int H = image.height()/4;

	image.setCenter(image.width()/2, image.height()/2);

	for( int j= 0; j<H  ; ++j)
	{
		for( int i= 0; i< W ; ++i)
		{
			if (j%2)
				image.pixelRelative(i,j) = RGBu8(255,128,0);
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

			RGBu8 Q;
			Q = P2*(uint8_t)(3) + P2*(uint8_t)(2);

		}
	}

	image.save(TEMPO_PATH+"out.png");

	ImageRGBAd im(512,512);
	im.for_all_pixels([](ImageRGBAd::PixelType& P, int /*x*/, int y)
	{
		P = ImageRGBAd::ASTexPixelType(1.0,0.3,0,(511-y%255)/255.0);
	});
	IO::save01_in_u8(im,TEMPO_PATH + "out_rgba.png");



  return EXIT_SUCCESS;
}

