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

#include <ASTex/image_indexed.h>

using namespace ASTex;



int main()
{
	ImageIndexedu8 img;
	bool ok = img.loadIndexedPNG(TEMPO_PATH+"small4.png");
	if (!ok)
		return 1;

	std::vector<ImageIndexedu8::PaletteColorType>& pal = img.palette();

	for (size_t i= 0; i<pal.size(); ++i)
	{
		std::cout << pal[i] << std::endl;
	}



	uint32_t HW = std::min(img.height(),img.width());

	for(uint32_t j=0; j<HW; ++j)
	{
		img.pixelAbsolute(j,j)=3;
	}

	for (ImageIndexedu8::IteratorIndexed it = img.beginIteratorIndexed(); !it.IsAtEnd(); ++it)
	{
		std::cout << "Pixel["<<it.GetIndex()[0]<< ","<< it.GetIndex()[1]<<"] >>> index = "<<int(it.Get()) << " >>> color = "<<img.color(it.Get())<< std::endl;
	}

	img.save(TEMPO_PATH+"indexed2.png");

	img.saveIndexedPNG(TEMPO_PATH+"ind4.png");

	ImageRGBu8 img2 = img.createRGB();
	img2.save(TEMPO_PATH+"indexedRGB.png");

	return EXIT_SUCCESS;

}

