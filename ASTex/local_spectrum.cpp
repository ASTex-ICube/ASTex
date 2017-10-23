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




#include <fstream>



#include <algorithm>
#include "local_spectrum.h"

//#include "output_geoffrey.h"

//includes pour la convert d'espace de couleur
#include <ASTex/colorspace_filters.h>
#include "itkCastImageFilter.h"

#include <ASTex/utils.h>

//Include pour le crop d'image ( pour Fourrier sur Crop de l'image)
// #include "itkCropImageFilter.h"
// #include <itkImage.h>

#include <ASTex/fourier.h>

#include <cmath>

namespace ASTex
{

namespace Fourier
{






double LocalSpectrum::g(double x, double y, double sigma)
{
	return 1.0 / (2.0 * M_PI * sigma * sigma) *
			std::exp(-(x * x + y * y) / (2 * sigma * sigma));
}


double LocalSpectrum::g_non_norm(double x, double y, double sigma)
{
	return   std::exp(-(x * x + y * y) / (2 * sigma * sigma));
}

double LocalSpectrum::hamming(double x, double y, double window_size)
{
	double coef = (2*M_PI) / (window_size-1);
	double rx = std::abs(x)<window_size/2 ? (0.54+0.46*cos(coef*x)): 0;
	double ry = std::abs(y)<window_size/2 ? (0.54+0.46*cos(coef*y)): 0;
	return rx*ry ;
}






LocalSpectrum::LocalSpectrum(const ImageGrayd& input, int windowsize, char windowtype, std::string name):
	m_windowsize(windowsize), m_windowtype(windowtype), m_name(name){
	m_input_width = input.width();
	m_input_height = input.height();

	// compute the window mask
	std::vector<std::vector< double>> W;
	W.resize( m_windowsize , std::vector< double>( m_windowsize ) );

	m_left_bnd = -(m_windowsize-1) / 2;
	m_right_bnd = (m_windowsize) / 2;

	for(int y = m_left_bnd; y <= m_right_bnd; y++)
	{
		for(int x = m_left_bnd; x <= m_right_bnd ; x++)
		{
			switch (windowtype)
			{
			case 'g' : 	    	//GAUSSIAN
				W[x - m_left_bnd][y - m_left_bnd] = LocalSpectrum::g_non_norm(x, y, m_windowsize/3);
				break;
			case 'h' : //HAMMING
				W[x - m_left_bnd][y - m_left_bnd] = LocalSpectrum::hamming(x,y, m_windowsize);
				break;
			case 'p' : //PORTE
				W[x - m_left_bnd][y - m_left_bnd] =  1;
//				W[x - m_left_bnd][y - m_left_bnd] =  (1/double(windowsize));
				break;
			default : //PORTE
				W[x - m_left_bnd][y - m_left_bnd] = 1;
				break;
			}
		}
	}

	// compute local spectra

	m_spectrum.resize(m_input_width*m_input_height);
	for(ImageSpectrald& sp: m_spectrum)
		sp.initItk(m_windowsize,true);
			

	ImageGrayd temp_crop;
	temp_crop.initItk(m_windowsize,m_windowsize,true);

	input.for_all_pixels([&] (int i, int j)
	{
		const int ic = std::max(std::min(i,m_input_width- m_right_bnd-1),-m_left_bnd); // the crop is not centered on i when too close to image boundaries
		const int jc = std::max(std::min(j,m_input_height- m_right_bnd-1),-m_left_bnd); // the crop is not centered on j when too close to image boundaries

		ImageSpectrald modulus;
		ImageSpectrald phase;

		crop_image(input,temp_crop,ic+m_left_bnd,ic+m_right_bnd,jc+m_left_bnd,jc+m_right_bnd);

		for (int l = 0; l < m_windowsize	; ++l)
			for (int k = 0; k < m_windowsize; ++k)
			{
				temp_crop.pixelAbsolute(k,l) *= W[k][l];
			}


		//FFT
		fftForwardModulusAndPhase(temp_crop, modulus, phase);

		for (int l = 0; l < m_windowsize; ++l)
			for (int k = 0; k < m_windowsize; ++k)
				spectrum(i,j).pixelAbsolute(k,l)= modulus.pixelAbsolute(k,l);

	});
}

LocalSpectrum::LocalSpectrum(const ImageGrayd& input, int windowsize, int welch_windowsize, int welch_step): m_windowsize(windowsize), m_windowtype('w')
{
	m_input_width = input.width();
	m_input_height = input.height();

	m_left_bnd = -(m_windowsize-1) / 2;
	m_right_bnd = (m_windowsize) / 2;

	m_spectrum.resize(m_input_width*m_input_height);
	for(ImageSpectrald& sp: m_spectrum)
		sp.initItk(m_windowsize,true);


	ImageGrayd temp_crop;
	temp_crop.initItk(welch_windowsize,welch_windowsize,true);

	int welch_left_bnd = -(welch_windowsize-1) / 2;
	int welch_right_bnd = (welch_windowsize) / 2;

	input.for_all_pixels([&] (int i, int j)
	{
		const int ic = std::max(std::min(i,m_input_width - welch_right_bnd-1),-welch_left_bnd); // the crop is not centered on i when too close to image boundaries
		const int jc = std::max(std::min(j,m_input_height - welch_right_bnd-1),-welch_left_bnd); // the crop is not centered on j when too close to image boundaries
		ImageSpectrald modulus;
		crop_image(input,temp_crop,ic+welch_left_bnd,ic+welch_right_bnd,jc+welch_left_bnd,jc+welch_right_bnd);
		welch(temp_crop, modulus, welch_step);

		modulus.for_all_pixels([&] (double M, int k, int l)
		{
			spectrum(i,j).pixelAbsolute(k,l)= M;
		});

	});
}


void LocalSpectrum::welch_post_loading(int welch_windowsize, int welch_step)
{
	
	int welch_left_bnd = -(welch_windowsize-1) / 2;
	int welch_right_bnd = (welch_windowsize) / 2;
	if (m_windowsize>1)
	{
		welch_left_bnd += m_windowsize/2 - 1;
		welch_right_bnd -= m_windowsize/2;
	}

	std::vector<ImageSpectrald> spectrum_res(m_input_width*m_input_height);
	for(ImageSpectrald& sp: spectrum_res)
		sp.initItk(m_windowsize,true);


	//Tous les pixels de l'image
	for (int j=0; j< m_input_height; ++j)
	{
		for (int i=0; i< m_input_width; ++i)
		{
			ImageSpectrald& sp_res = spectrum_res[i + j*m_input_width];

			//Toutes les fréquences
			sp_res.for_all_pixels([&] (double& P, int fx, int fy)
			{
				P=0; // ?
				int card = 0;
				double sum = 0;
				//Tous les pixels du voisinage
				for (int k = welch_left_bnd; k <= welch_right_bnd; k+=welch_step)
				{
					for (int l = welch_left_bnd; l <= welch_right_bnd;  l+=welch_step)
					{
						const int isum = std::max(0,std::min(m_input_width-1,i+k));
						const int jsum = std::max(0,std::min(m_input_height-1,j+l));
						double val = spectrum(isum,jsum).pixelAbsolute(fx,fy);
						sum+= val*val;
						card ++;
					}
				}
				P = std::sqrt(sum/card);
			});
		}
	}

	m_spectrum.swap(spectrum_res);
}



void LocalSpectrum::mosaic (ImageSpectrald& output)
{
	for (int i = -m_left_bnd; i < m_input_width- m_right_bnd; i+=m_windowsize)
	{
		for (int j = -m_left_bnd; j < m_input_height-m_right_bnd; j+=m_windowsize)
		{
			fulfill_crop_image(spectrum(i,j), output, i+m_left_bnd,j+m_left_bnd);
		}
	}
}


void LocalSpectrum::distance_map_to_pixel_linear_weights (ImageGrayd& output, int i_ref, int j_ref)
{
	for (int j = 0; j < m_input_height; ++j)
	{
		for (int i = 0; i < m_input_width; ++i)
		{
			double dist = distance_spectrum_to_spectrum_linear_weights(spectrum(i_ref,j_ref), spectrum(i,j));
			output.pixelAbsolute(i,j) = dist;
		}
	}
}

void LocalSpectrum::distance_map_to_pixel_uniform_weights (ImageGrayd& output, int i_ref, int j_ref)
{
	for (int j = 0; j < m_input_height; ++j)
	{
		for (int i = 0; i < m_input_width; ++i)
		{
			double dist = distance_spectrum_to_spectrum_uniform_weights(spectrum(i_ref,j_ref), spectrum(i,j));
			output.pixelAbsolute(i,j) = dist;
		}
	}
}

void LocalSpectrum::distance_map_to_spectum_linear_weights(ImageGrayd& output,const ImageSpectrald& sp_ref)
{
	for (int j = 0; j < m_input_height; ++j)
	{
		for (int i = 0; i < m_input_width; ++i)
		{
			double dist = distance_spectrum_to_spectrum_linear_weights(sp_ref,spectrum(i,j));
			output.pixelAbsolute(i,j) = dist;
		}
	}
}



void LocalSpectrum::distance_map_to_spectum_uniform_weights(ImageGrayd& output,const ImageSpectrald& sp_ref)
{
	for (int j = 0; j < m_input_height; ++j)
	{
		for (int i = 0; i < m_input_width; ++i)
		{
			double dist = distance_spectrum_to_spectrum_uniform_weights(sp_ref,spectrum(i,j));
			output.pixelAbsolute(i,j) = dist;
		}
	}
}



void LocalSpectrum::distance_mosaic (ImageGrayd& output, int step)
{
	int left_step = -(step-1) / 2;
	int right_step = (step) / 2;

	for (int j = -m_left_bnd; j < m_input_height-m_right_bnd; j+=m_windowsize)
	{
		for (int i = -m_left_bnd; i < m_input_width- m_right_bnd; i+=m_windowsize)
		{
			for (int l= left_step; l<= right_step; ++l)
			{
				for(int k= left_step; k<= right_step; ++k)
				{
					if (i+k>= 0 && i+k< m_input_width && j+l>= 0 && j+l< m_input_height )
					{
						double dist = distance_spectrum_to_spectrum_linear_weights(spectrum(i,j),spectrum(i+k,j+l));
						output.pixelAbsolute(i+k,j+l)= dist;
					}
				}
			}
		}
	}
}

void LocalSpectrum::compute_local_guidance_map (int filtre_size)
{
	int left_step = -(filtre_size-1) / 2;
	int right_step = (filtre_size) / 2;

	m_local_guidance.resize( m_input_width * m_input_height );

	for (ImageGrayd& lgm : m_local_guidance)
		lgm.initItk(filtre_size,filtre_size,true);

	for (int j = 0; j < m_input_height; ++j)
	{
		for (int i = 0; i < m_input_width; ++i)
		{
			for (int l= left_step; l<= right_step; ++l)
			{
				for(int k= left_step; k<= right_step; ++k)
				{
					if (i+k>= 0 && i+k< m_input_width && j+l>= 0 && j+l< m_input_height )
					{
						double dist = distance_squared_spectrum_to_spectrum_uniform_weights(spectrum(i,j),spectrum(i+k,j+l));
						local_guidance(i,j).pixelAbsolute(k-left_step, l-left_step)= dist;
					}

				}
			}
		}
	}
}

char LocalSpectrum::get_windowtype (){ return m_windowtype;}

char LocalSpectrum::get_windowsize (){ return m_windowsize;}

std::string LocalSpectrum::get_name (){ return m_name;}


} //namespace Fourier
} //namespace ASTex
