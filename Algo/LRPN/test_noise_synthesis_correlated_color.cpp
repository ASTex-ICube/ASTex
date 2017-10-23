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
#include <ASTex/internal.h>
#include "noise_synthesis_correlated_color.h"

int main( int argc, char ** argv )
{
	if( argc < 1 )
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " Source " << std::endl;
		return EXIT_FAILURE;
	}
	std::string filename_source = argv[1];
	std::string base_dir = ASTex::IO::remove_ext(filename_source) + "_noise_filter/";

	int nb_clusters = 2;

	ASTex::Synthesis_corel(filename_source, base_dir, nb_clusters);

	return EXIT_SUCCESS;
}

