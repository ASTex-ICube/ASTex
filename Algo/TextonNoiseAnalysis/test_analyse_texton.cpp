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




#include "analyse_texton.h"
#include <fstream>
#include <ASTex/utils.h>
#include <ASTex/easy_io.h>


int main(int argc, char** argv)
{

	if( argc < 3 )
	{
            std::cerr << "Usage: " << std::endl;
            std::cerr << argv[0] << " <Source_file ; Outdir >" << std::endl;
            return EXIT_FAILURE;
	}

    // Get the arguments
    std::string source_dir = argv[1];
    std::string out_dir_base = argv[2];

    // Get the name of the file without extention en creation the folder for the res
    std::string name_file = ASTex::IO::remove_ext(ASTex::IO::remove_path(source_dir));
    std::string out_dir = out_dir_base+name_file+"/";

    ASTex::create_directory(out_dir);

    // LOAD INPUT
    ASTex::ImageRGBd input_color;
    ASTex::IO::loadu8_in_01(input_color,source_dir);

	std::vector<ASTex::ImageRGBd> images_alpha_visu;
	std::vector<ASTex::ImageRGBd> images_alpha_visu_bounded;

	ASTex::analyse_texton(input_color, images_alpha_visu, images_alpha_visu_bounded);

	for( const ASTex::ImageRGBd& im : images_alpha_visu)
		ASTex::IO::save01_in_u8(im, out_dir+"Texton_plus_mean_"+std::to_string(im.width())+".png");
	for( const ASTex::ImageRGBd& im : images_alpha_visu_bounded)
		ASTex::IO::save01_in_u8(im, out_dir+"Texton_plus_mean_bounded_"+std::to_string(im.width())+".png");

	return EXIT_SUCCESS;
 }
