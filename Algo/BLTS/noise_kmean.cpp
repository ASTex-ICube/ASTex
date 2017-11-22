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


#include "noise_kmean.h"

namespace ASTex
{


class SpectralSeeds
{
public:
	using InputType = ImageGrayd;
	using Type = ImageSpectrald;


	SpectralSeeds(const InputType& input, const int nb_clusters, const int welchsize, const int welchstep, const int fft_size):
		nbc(nb_clusters)
	{
		lsp = new Fourier::LocalSpectrum(input,fft_size,'p');
		lsp->welch_post_loading(welchsize, welchstep);
		seeds.resize(nbc);
		clusters_size.resize(nbc);

		im_w = input.width();
		im_h = input.height();

		for (int c = 0; c < nbc; ++c)
		{
			seeds[c].initItk(fft_size,fft_size,true); // ?true?
			int x = rand()% im_w;
			int y = rand()% im_h;
			seeds[c].copy_pixels(lsp->spectrum(x,y));
			clusters_size[c] = -1;
		}

	}

	~SpectralSeeds()
	{
		delete lsp;
	}

	void compute(std::vector< ImageGrayd >& seeds_distance_maps)
	{
		for (int c = 0; c < nbc; ++c)
			lsp->distance_map_to_spectum_uniform_weights(seeds_distance_maps[c],seeds[c]);

	}

	void average(ImageGray16& clusters)
	{
		// update seeds
		for (int c = 0; c < nbc; ++c)
		{
			seeds[c].zero();
			clusters_size[c] = 0;
		}

		for (int y=0; y<im_h; ++y)
		{
			for (int x=0; x<im_w; ++x)
			{
				const int c = clusters.pixelAbsolute(x,y);
				seeds[c] += lsp->spectrum(x,y);
				clusters_size[c] ++;
			}
		}

		for (int c = 0; c < nbc; ++c)
		{
			if (clusters_size[c] != 0)
			{
				seeds[c] /= double(clusters_size[c]);
			}
		}
	}

protected:
	std::vector< Type > seeds; // seeds = spectra

	std::vector<int> clusters_size; // nb elements in each cluster

	int nbc;

	Fourier::LocalSpectrum* lsp;

	int im_w;
	int im_h;

};




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
	kmean<SpectralSeeds>(input,masks_dist,masks_bin,nb_clusters,(fft_size*2)-1,4,fft_size );

	for (int c=0; c<nb_clusters; ++c)
	{
		IO::save(masks_dist[c],base_dir+name_file + "_mask_dist"+std::to_string(c)+".png");
		IO::save(masks_bin[c],base_dir+name_file+"_mask_bin"+std::to_string(c)+".png");

	}
}

}
