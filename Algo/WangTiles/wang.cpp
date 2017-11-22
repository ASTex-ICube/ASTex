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
//	std::string fn = "C:/Users/thery/Desktop/blue_rust.png";
//	std::string fn = "/Users/thery/Desktop/blue_rust.png";
	std::string fn = "/tmp/blue_rust.png";
//	std::string fn = "/tmp/ASTex_data/quilting_input8.png";

	QApplication app(argc, argv);

	ImageRGBu8 im;
	im.load(fn);

	auto start_chrono = std::chrono::system_clock::now();

//	WangTilesGenerator<ImageRGBu8,4> wta(im,120);
//	auto wang = wta.create();
	auto wang = WangTilesGenerator<ImageRGBu8,2>::create(im,200);

	std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start_chrono;
	std::cout << "wang : " << elapsed_seconds.count() << " s." << std::endl;


	ImageRGBu8 ti = wang.all_tiles();
	auto v1 = image_viewer(ti,"merge", &app);

	ImageRGBu8 gen = wang.compose(15, 10);
	auto v2 = image_viewer(gen,"gen", &app);

//	ImageViewer imgv("tiled", &app);
//	imgv.update(gen);
//	imgv.show();

	return app.exec();
}
