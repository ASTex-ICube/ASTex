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



#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>
#include <ASTex/easy_io.h>

#include "wang_tiles.h"

using namespace ASTex;


int main(int argc, char** argv)
{
	std::string fn = TEMPO_PATH+"quilting_input8.png";
	int tw = 100;
	int gen_sz = 10;
	int nbcol = 2;

	if (argc>4)
	{
		fn = std::string(argv[1]);
		tw = atoi(argv[2]);
		gen_sz = atoi(argv[3]);
		nbcol = atoi(argv[4]);
	}
	else
	{
		std::cout << argv[0]<< "input tile_width generated_width nb_sides ...ssssss using default"<< std::endl;
	}

	ImageRGBu8 im;
	im.load(fn);

	auto start_chrono = std::chrono::system_clock::now();

	auto wang = WangTilesGenerator<ImageRGBu8>::create(im,tw,nbcol);

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "wang tile timing: " << elapsed_seconds.count() << " s." << std::endl;


	ImageRGBu8 ti = wang.all_tiles();
	ti.save(TEMPO_PATH+"wang_tiles.png");

	ImageRGBu8 gen = wang.compose(gen_sz, gen_sz);
	gen.save(TEMPO_PATH+"wang_generated.png");

	return EXIT_SUCCESS;
}
