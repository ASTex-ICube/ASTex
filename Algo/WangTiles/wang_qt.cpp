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

#include "imageviewer.h"

using namespace ASTex;

void app_mouse_clicked(int /*button*/, int /*x*/, int /*y*/, int /*id*/) {}
void app_key_pressed(int /*code*/, char /*key*/, int /*id*/) {}

int main(int argc, char** argv)
{
	QApplication app(argc, argv);
	std::string fn = TEMPO_PATH+"quilting_input8.png";
	int tw = 100;
	int gen_sz = 1000;

	if (argc>=4)
	{
		fn = std::string(argv[1]);
		tw = atoi(argv[2]);
		gen_sz = atoi(argv[3]);
	}
	else
	{
		std::cout << argv[0]<< " tile_width generated_width   using default"<< std::endl;
	}

	ImageRGBu8 im;
	im.load(fn);

	auto start_chrono = std::chrono::system_clock::now();

	auto wang = WangTilesGenerator<ImageRGBu8,3>::create(im,tw);

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "wang tile timing: " << elapsed_seconds.count() << " s." << std::endl;


	ImageRGBu8 ti = wang.all_tiles();
	auto v1 = image_viewer(ti);

	ImageRGBu8 gen = wang.compose(gen_sz/tw, gen_sz/tw);
	auto v2 = image_viewer(gen,"gen", &app);


	return app.exec();
}
