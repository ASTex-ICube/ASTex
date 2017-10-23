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



#include <cmath>
#include <ASTex/special_io.h>
#include <ASTex/fourier.h>
#include <ASTex/local_spectrum.h>
#include <ASTex/utils.h>
#include <ASTex/distances_maps.h>

namespace ASTex
{


void distance_invert_expo (const ImageGrayd& s, ImageGrayd& d, double dev)
{
	for_all_pixels(s,d,[dev](double P, double& Q)
	{
		Q = std::exp(-P/dev);
	});
}

void Image_log (const ImageGrayd& s, ImageGrayd& d, double ratio)
{
	for_all_pixels(s,d,[ratio](double P, double& Q)
	{
		Q = std::sqrt(P*ratio);
	});
}

void Image_expo (const ImageGrayd& s, ImageGrayd& d, double dev)
{
	for_all_pixels(s,d,[dev](double P, double& Q)
	{
		Q = 1.0 - std::exp(-P/dev);
	});
}

double get_average(const ImageGrayd& s)
{
	double sum=0;

	s.for_all_pixels([&sum] (double p)
	{
		sum += p;
	});

	return sum/double(s.width()*s.height());

}

double get_max(const ImageGrayd& s)
{
	double max=s.pixelAbsolute(0,0);

	s.for_all_pixels([&max] (double p)
	{
		if (max < p) max = p;
	});
	return max;
}



int kmean(const ImageGrayd& input, std::vector<ImageGrayd>& output_dists, std::vector<ImageGrayd>& output_bin ,
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

	Fourier::LocalSpectrum lsp (input,fft_size,'p');
	lsp.welch_post_loading(welchsize, welchstep);

	ImageGray16 clusters(im_w, im_h, true);

	std::vector< ImageSpectrald > seeds; // seeds = spectra
	seeds.resize(nb_clusters);

	std::vector<int> clusters_size; // nb elements in each cluster
	clusters_size.resize(nb_clusters);

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
		seeds[c].initItk(fft_size,fft_size,true); // ?true?
		clusters_size[c] = -1;
		seeds_distance_maps[c].initItk(im_w,im_h,true);
	}

	// random seed init
	for (int c = 0; c < nb_clusters; ++c)
	{
		int x = rand()% im_w;
		int y = rand()% im_h;
		seeds[c].copy_pixels(lsp.spectrum(x,y));
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
		for (int c = 0; c < nb_clusters; ++c)
			lsp.distance_map_to_spectum_uniform_weights(seeds_distance_maps[c],seeds[c]);

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
		for (int c = 0; c < nb_clusters; ++c)
		{
			seeds[c].zero();
			clusters_size[c] = 0;
		}

		for (int y=0; y<im_h; ++y)
		{
			for (int x=0; x<im_w; ++x)
			{
				const int c = clusters.pixelAbsolute(x,y);
				seeds[c] += lsp.spectrum(x,y);
				clusters_size[c] ++;
			}
		}

		for (int c = 0; c < nb_clusters; ++c)
		{
			if (clusters_size[c] != 0)
			{
				seeds[c] /= double(clusters_size[c]);
			}
		}
	}
	// KMEAN TERMINE
	
	double average= 0;
	double max_spec = 0;

	// PRINT DISTANCE ET SPECTRE
	for (int c = 0; c < nb_clusters; ++c)
	{
		average += get_average(seeds_distance_maps[c]);
		max_spec = std::max(max_spec,get_max(seeds[c]));
	}
	average/=nb_clusters;
	
	for (int c = 0; c < nb_clusters; ++c)
	{
		output_dists[c].initItk(im_w,im_h,true);
//		masks_spectrum_visu[c].initItk(fft_size,fft_size,true);
//		masks_spectrum_save[c].initItk(fft_size,fft_size,true);

////		distance_invert_expo(seeds_distance_maps[c],output_dists[c],average);
//		Image_log(seeds[c],masks_spectrum_save[c],1/double(fft_size));
//		Image_expo(seeds[c],masks_spectrum_visu[c],1/max_spec);
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



void KMnoise(const std::string& filename_source, const std::string& base_dir, const int nb_clusters, const int fft_size)
{
	std::string name_file = IO::remove_ext(IO::remove_path(filename_source));

	// LOAD INPUT
	ImageGrayd input;
	IO::load_RGB_2_lightness(input, filename_source);

	// COMPUTE MASKS
	std::vector<ImageGrayd> masks_dist(nb_clusters);
	std::vector<ImageGrayd> masks_bin(nb_clusters);

	// TODO check param                              V                 V
	kmean(input,masks_dist,masks_bin,nb_clusters,(fft_size*2)-1,4,fft_size );

	for (int c=0; c<nb_clusters; ++c)
	{
		IO::save(masks_dist[c],base_dir+name_file + "_mask_dist"+std::to_string(c)+".png");
		IO::save(masks_bin[c],base_dir+name_file+"_mask_bin"+std::to_string(c)+".png");

	}
}

}
