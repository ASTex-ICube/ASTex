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
#include <chrono>

#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>
#include <ASTex/store.h>


using namespace ASTex;

struct Pix
{
	uint16_t x;
	uint16_t y;
	uint8_t g;
	inline Pix() {}
	inline Pix(uint16_t i, uint16_t j, uint8_t c):x(i), y(j), g(c) {}
};



int main()
{
	ImageRGBu8 image(600,600,true);

	Store<Index> si;

	si << gen_index(1,1);
	si << gen_index(2,2);
	si << gen_index(3,1);
	si << gen_index(4,2);

	for(int i=6 ; i < 295; i+=3)
		for(int j=1; j<599; ++j)
			si << gen_index(j,i);

	si.remove(5);
	si.remove(7);
	si.remove(9);

	for(const auto& p: si)
		image.pixelAbsolute(p)= ImageRGBu8::itkPixelNorm(1,0,0);


	Store<Pix> sp;

	for(int i=300 ; i < 595; i+=3)
		for(int j=1; j<599; ++j)
			sp.emplace_back(j,i,((i*600+j)%255));


	for(const auto& p: sp)
		image.pixelAbsolute(p.x,p.y)= ImageRGBu8::itkPixel(p.g,p.g,p.g);


	image.save(TEMPO_PATH+"store_out.png");

	return EXIT_SUCCESS;
}

