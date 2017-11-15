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



#include <string>
#include <cmath>
#include <ASTex/special_io.h>
#include <ASTex/fourier.h>
#include <ASTex/local_spectrum.h>
#include <ASTex/utils.h>
#include <ASTex/distances_maps.h>


namespace ASTex
{


template<typename SEEDS>
int kmean(const typename SEEDS::InputType& input, std::vector<ImageGrayd>& output_dists, std::vector<ImageGrayd>& output_bin ,
		  const int nb_clusters, const int welchsize, const int welchstep, const int fft_size)
{
	srand(time(NULL));
	const int im_w = input.width();
	const int im_h = input.height();

	for (int c=0; c<nb_clusters; ++c)
	{
		output_dists[c].initItk(im_w,im_h,true);
		output_bin[c].initItk(im_w,im_h,true);
//		masks_spectrum_visu[c].initItk(fft_size,fft_size,true);
//		masks_spectrum_save[c].initItk(fft_size,fft_size,true);

	}

	SEEDS seeds(input, nb_clusters, welchsize, welchstep, fft_size);

	ImageGray16 clusters(im_w, im_h, true);

	std::vector< ImageGrayd > seeds_distance_maps; 	// allocate distance maps
	seeds_distance_maps.resize(nb_clusters);

	std::vector< std::vector<double> > distance_inter_seed; 	//
	distance_inter_seed.resize( nb_clusters , std::vector<double>( nb_clusters ) );

	std::vector< double > distance_intra_classe; 	//
	distance_intra_classe.resize(nb_clusters);

	std::vector< double > average_dist_intra_class;
	average_dist_intra_class.resize(nb_clusters);

//	ImageSpectrald sp_tmp(fft_size,fft_size,true);

	for (int c = 0; c < nb_clusters; ++c)
	{
		seeds_distance_maps[c].initItk(im_w,im_h,true);
	}

	// iterate

	int iter = 0;
	const int max_iter = 20;
	int nb_changes = 1;

	while (iter < max_iter && nb_changes > 0)
	{
		++iter;
		nb_changes = 0;

		// compute distance maps
		seeds.compute(seeds_distance_maps);

		// compute clusters
		for (int y=0; y<im_h; ++y)
		{
			for (int x=0; x<im_w; ++x)
			{
				int c_min = 0;
				double dist_min = seeds_distance_maps[0].pixelAbsolute(x,y);
				for (int c = 1; c < nb_clusters; ++c)
				{
					double dist = seeds_distance_maps[c].pixelAbsolute(x,y);
					if (dist < dist_min)
					{
						c_min = c;
						dist_min = dist;
					}
				}
				if (clusters.pixelAbsolute(x,y) != c_min)
				{
					nb_changes++;
				}
				clusters.pixelAbsolute(x,y) = c_min;
			}
		}

		// update seeds
		seeds.average(clusters);
	}
	// KMEAN TERMINE

//	double average= 0;

//	// PRINT DISTANCE ET SPECTRE
//	for (int c = 0; c < nb_clusters; ++c)
//	{
//		average += get_average(seeds_distance_maps[c]);
//	}
//	average/=nb_clusters;

	for (int c = 0; c < nb_clusters; ++c)
	{
		output_dists[c].initItk(im_w,im_h,true);
	}

	for (int y=0; y<im_h; ++y)
	{
		for (int x=0; x<im_w; ++x)
		{

			double sum_dist=0;
			for (int c = 0; c < nb_clusters; ++c)
				sum_dist += seeds_distance_maps[c].pixelAbsolute(x,y);

			for (int c = 0; c < nb_clusters; ++c)
				output_dists[c].pixelAbsolute(x,y) = seeds_distance_maps[c].pixelAbsolute(x,y) / sum_dist;

		}
	}




	for (int y=0; y<im_h; ++y)
	{
		for (int x=0; x<im_w; ++x)
		{
			 output_bin[clusters.pixelAbsolute(x,y)].pixelAbsolute(x,y)=1;
		}
	}

	return 0;
}



/**
 * @brief KMnoise
 * @param filename_source
 * @param base_dir
 * @param nb_clusters
 * @param fft_size
 */
void KMnoise(const std::string& filename_source, const std::string& base_dir, const int nb_clusters, const int fft_size);

}
