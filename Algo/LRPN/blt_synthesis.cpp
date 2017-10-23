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





#include<ASTex/utils.h>
#include<ASTex/internal.h>
#include "noise_filter.h"
#include "noise_kmean.h"
#include "noise_synthesis_correlated_color.h"

extern void almost_wang_tiles( const std::string& filename_source, const std::string& base_dir, const int nb_clusters);
extern void biscalenoisepatchexg(const std::string& filename_source, const std::string& base_dir, int NCLUSTERS, int CONTENTS, int NSCALES, int NROTANGLES, float TSCALE, int SZ_MULT);


const int filt_size = 16;
const int size_fft = 8;
const float sig_freq = 0.03f;
const int nb_step_filt = 5;


int main( int argc, char ** argv )
{
	if (argc < 2)
	{
		std::cerr << ASTex::IO::remove_path(std::string(argv[0])) << " texture_input [size_mult=2] [nb_cluster=2]" << std::endl;
#ifdef WIN32
		std::getchar();
#endif
		return 1;
	}

	std::string filename_source = argv[1];
	std::string base_dir = ASTex::IO::remove_ext(filename_source) + "_Generated/";

	int nb_clusters = 2;
	int sz_mult = 2;

	if (argc == 3)
	{
		sz_mult = atoi(argv[2]);
		
	}
	if( argc == 4 )
	{
		sz_mult = atoi(argv[2]);
		nb_clusters = atoi(argv[3]);
	}

	if (!ASTex::create_directory(base_dir))
		std::cout << "warning directory " << base_dir << " already exists!" << std::endl;

	std::cout << "NOISE FILTERING" << std::endl;
	ASTex::RPnoise(filename_source, base_dir, filt_size, size_fft, sig_freq, nb_step_filt);

	std::cout << "NOISE KMEAN" << std::endl;
	ASTex::KMnoise(filename_source, base_dir, nb_clusters, size_fft);

	std::cout << "NOISE SYNTHESIS" << std::endl;
	ASTex::Synthesis_corel(filename_source, base_dir, nb_clusters);

	std::cout << "WANG TILE" << std::endl;
	almost_wang_tiles(filename_source, base_dir, nb_clusters);

	std::cout << "BISCALE PATCH EXCHANGE" << std::endl;
	biscalenoisepatchexg(filename_source, base_dir, nb_clusters, 15, 3, 5, 3,sz_mult);

	return EXIT_SUCCESS;
}
