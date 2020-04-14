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
#include <ASTex/image_gray.h>
#include <ASTex/histogram.h>

using namespace ASTex;


int main()
{
	ImageGrayu8 image;

	bool ok = image.load(TEMPO_PATH+"simpleGray.png");
	if (!ok)
		return 1;

	std::cout << " Taille =" << image.width() << "/" << image.height() << std::endl;

	uint32_t W = image.width()/4;
	uint32_t H = image.height()/4;

	for(uint32_t j=0; j<H  ; ++j)
	{
		for(uint32_t i=0; i< W ; ++i)
		{
			image.pixelAbsolute(i,j) = 128;
		}
	}

	HistogramGrayu8 h(image);

	double mean = h.mean();
	double sigma = h.variance();

	std::cout << "Mean "<< mean << std::endl;
	std::cout << "Variance " << sigma << std::endl;

	image.save(TEMPO_PATH+"simpleG2.png");

  return EXIT_SUCCESS;
}



